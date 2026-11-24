import os
import re
import sys
from pathlib import Path
import MDAnalysis
import freesasa
import MDAnalysis as mda
import pandas as pd
from MDAnalysis.analysis.rms import RMSD, RMSF
from MDAnalysis.analysis import align
import interfacea as ia
from openmm import app
import numpy as np
import DynaPIN.handling as hp
import DynaPIN.pdb_fixer as ptm
import random
import json
import subprocess
import time
import tempfile
from sys import platform
from shutil import make_archive
import shutil
import warnings

warnings.filterwarnings("ignore")
freesasa.setVerbosity(freesasa.silent)
handler = hp.tables_errors()


class dynapin:
    def __init__(self, trajectory_file, stride=1, split_models=True, chains=None, output_dir=None, topology_file=None 
                 ):
        """A class to perform Quality Control, Residue Based and Interaction Based analyses on MD simulations. Results of the analyses are printed as .csv files under a folder named 'tables'. MD outputs with any exception are transformed to .pdb file and the analyses are run through this .pdb file. Number of frames that will be fetched from initial input file to be analysed can be set by stride value.
        
        Keyword arguments:
        trajectory_file -- Path to the input file
        stride -- Number of frames to be skipped. By default, 1 which means no frames to be skipped.
        split_models -- Boolean. By default, False. If True, all the models will be splitted into a folder named 'models'.
        Return: None
        """
        if output_dir is None:
            base_for_name = os.path.splitext(os.path.basename(trajectory_file))[0]
            output_dir = base_for_name + '_db' + str(random.randrange(100,999))
        self.output_dir = output_dir
        self.job_path = os.path.abspath(self.output_dir)
        if not os.path.exists(self.job_path):
            os.makedirs(self.job_path, exist_ok=True)
        self.target_path = os.path.join(self.job_path, 'tables')
        if not os.path.exists(self.target_path):
            os.mkdir(self.target_path)

        def _resolve_path(p):
            if p is None:
                return None
            if os.path.isabs(p) and os.path.exists(p):
                return p
            if os.path.exists(p):
                return os.path.abspath(p)
            cand = os.path.join(self.job_path, p)
            if os.path.exists(cand):
                return os.path.abspath(cand)
            return p  

        trajectory_file = _resolve_path(trajectory_file)
        topology_file   = _resolve_path(topology_file) if topology_file else None

        handler.test_inp_path(trajectory_file)
        if topology_file:
            handler.test_inp_path(topology_file)
        handler.test_stride(stride)
        handler.test_split_models(split_models)

        #params
        self.trajectory_file_ = os.path.abspath(trajectory_file)
        self.split_models_ = split_models
        self.chains_ = chains
        self.qc_flag = False
        self.rb_flag = False
        self.ib_flag = False
        self.topology_file = topology_file
        self.stride = stride
        self.rmsd_data= None
        self.get_all_hph = None
        self.foldx_path = None
        self.dssp_flag = False

        def _normalize_chains(chains):
            if isinstance(chains, str):
                chains = [c.strip() for c in chains.replace(";", ",").split(",") if c.strip()]
            if isinstance(chains, (list, tuple)):
                chains = [c for c in chains if str(c).strip() != ""]
            return chains or None

        self.chains_ = _normalize_chains(self.chains_)

        if output_dir is None:
            file_name = trajectory_file.split('.')[0]
            output_dir = file_name + '_db' + str(random.randrange(100,999))

        self.output_dir = output_dir
        self.job_path = os.path.abspath(self.output_dir)


        if not os.path.exists(self.job_path): #create job folder
            os.mkdir(self.job_path)

        self.target_path = os.path.join(self.job_path, 'tables') 

        if not  os.path.exists(self.target_path): 
            os.mkdir(self.target_path)

        #check for preprocess ,define pdb_file
        file_name = os.path.splitext(os.path.basename(self.trajectory_file_))[0]
        file_ext = os.path.splitext(self.trajectory_file_)[1].lstrip('.')


        if file_ext in ('dcd', 'xtc', 'trr'):
            name = file_name.split("\\")[-1].split(".")[0]
            out_file = os.path.join(self.job_path, f"{name}.pdb")
            if not self.topology_file:
                if file_ext == 'dcd':
                    topo_prompt = 'Please provide topology file (.psf or .pdb) for the trajectory:\n'
                else:
                    topo_prompt = 'Please provide topology file (.tpr or .pdb) for the trajectory:\n'
                self.topology_file = _resolve_path(input(topo_prompt).strip())

            topo_ext = os.path.splitext(self.topology_file)[1].lower()
            allowed_top = ('.psf', '.pdb') if file_ext == 'dcd' else ('.gro', '.pdb', '.tpr')
            if topo_ext not in allowed_top:
                allowed_str = ' or '.join(ext.lstrip('.') for ext in allowed_top)
                raise ValueError(f"{file_ext.upper()} trajectories require a topology file in {allowed_str.upper()} format, got '{topo_ext or 'unknown'}'.")

            u = mda.Universe(self.topology_file, self.trajectory_file_)
            c = np.unique(u.atoms.segids)
            if all(str(x).strip() == "" for x in c):
                try:
                    c = np.unique(u.atoms.chainIDs)
                except AttributeError:
                    c = c
            if self.chains_ is None and len(c) > 2:
                sc = input("Select two chains to analyze (comma-separated, e.g. A,B): ")
                self.chains_ = _normalize_chains(sc)
            s = self.chains_
            
            self.pdb_file = os.path.join(self.job_path, self._preprocess_dcd(self.trajectory_file_, out_file, stride, self.topology_file, sc=s))

        elif file_ext == 'pdb':
            u = mda.Universe(self.trajectory_file_)
                
            segs = np.unique(u.atoms.segids)
            try:
                chainids = np.unique(u.atoms.chainIDs)
            except AttributeError:
                chainids = []
            chains_available = [c for c in segs if str(c).strip() != ""] or [c for c in chainids if str(c).strip() != ""]

            if self.chains_ is None and len(chains_available) > 2:
                sc = input("Select two chains to analyze (comma-separated, e.g. A,B): ")
                self.chains_ = _normalize_chains(sc)

            if self.chains_:
                sel = [f"(segid {c.upper()} or chainid {c.upper()})" for c in self.chains_]
                a = u.select_atoms(' or '.join(sel))
            else:
                a = u.atoms

            chain_len = len(chains_available) if chains_available else len(np.unique(u.atoms.segids))
            if self.stride != 1:
                name = file_name.split("\\")[-1]
                out_file = os.path.join(self.job_path, f"{name}_stride{stride}.pdb")
                
                with mda.Writer(out_file, u.atoms.n_atoms) as W:
                    for ts in u.trajectory[::stride]:
                        W.write(a)

                self.pdb_file = out_file
            else:
                if self.chains_:
                    name = file_name.split("\\")[-1]
                    out_file = os.path.join(self.job_path, f"{name}_ch_{','.join(self.chains_)}.pdb")
                    
                    with mda.Writer(out_file, u.atoms.n_atoms) as W:
                        for ts in u.trajectory:
                            W.write(a)

                    self.pdb_file = out_file
                else:
                    self.pdb_file = trajectory_file
        
        if self.chains_:
            handler.test_chain_sel(self.chains_, self.pdb_file)
            self.pdb_file = self._sel_chain(self.pdb_file, self.chains_, self.job_path)

        self.pdb_file = self._standardize_residue_names(self.pdb_file)

        # prepare models after chain selection/standardization so downstream analyses use the final structure
        if self.split_models_:
            self._prepare_models(self.pdb_file, self.job_path)

    def _standardize_residue_names(self, pdb_path):
        """Replaces common non-standard Histidine names with 'HIS' in a PDB file."""
        replacements = {
            "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
            "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
            "HISA": "HIS", "HISB": "HIS", "HISH": "HIS",
            "HISD": "HIS"
        }
        
        base_name = os.path.splitext(os.path.basename(pdb_path))[0]
        output_path = os.path.join(self.job_path, f"{base_name}_standardized.pdb")

        needs_replacement = False
        try:
            with open(pdb_path, 'r') as infile:
                content_check = infile.read()
                for old_name in replacements:
                    if f' {old_name} ' in content_check:
                        needs_replacement = True
                        break
        except FileNotFoundError:
            warnings.warn(f"PDB file not found at {pdb_path}. Skipping standardization.")
            return pdb_path

        if not needs_replacement:
            return pdb_path

        with open(pdb_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                if line.startswith(('ATOM', 'HETATM')):
                    res_name = line[17:20].strip()
                    if res_name in replacements:
                        new_line = line[:17] + 'HIS'.ljust(3) + line[20:]
                        outfile.write(new_line)
                    else:
                        outfile.write(line)
                else:
                    outfile.write(line)
        
        return output_path

    def _sel_chain(self, pdb, chains, job_path):
        u = mda.Universe(pdb)
        base = os.path.splitext(os.path.basename(pdb))[0]
        new_name = os.path.join(job_path, f"{base}_{','.join(chains)}.pdb")

        sel = []
        for el in chains:
            sel.append(f"(segid {el.upper()} or chainid {el.upper()})")

        t = u.select_atoms(' or '.join(sel))
        t.write(new_name, frames='all')

        return new_name

    def _get_params_(self):
        params = {
            'output_dir':self.output_dir,
            'trajectory_file': self.trajectory_file_,
            'topology_file': self.topology_file,
            'stride': self.stride,
            'split_models':self.split_models_,
            'chains': self.chains_,
            'QualityControl': {
                'Run': self.qc_flag,
                'rmsd_data':self.rmsd_data
                },
            'ResidueBased':{
                'Run': self.rb_flag,
                'FoldX_path': self.foldx_path,
            },
            'InteractionBased': { 
                'Run': self.ib_flag,
                'get_all_hph': self.get_all_hph
            }

        }

        json_path = os.path.join(self.job_path, 'table_params.json')
        with open(json_path, 'w+') as ofh:
            json.dump(params, ofh, indent=4)

        if self.split_models_:
            models_path = os.path.join(self.job_path, "models")
            archive_base = os.path.join(self.job_path, "models")
            make_archive(archive_base, "zip", root_dir=self.job_path, base_dir="models")
            shutil.rmtree(models_path)


    @staticmethod
    def _executable(name, comes_from):
        """Where a console script is, even when PATH does not know.

        DynaPIN is often started by a path rather than from an activated environment - a
        venv's bin directory, a notebook kernel, a scheduler - and a process spawned from
        there cannot see the console scripts installed beside the interpreter. So look
        next to the running interpreter first, then fall back to PATH.
        """
        beside = Path(sys.executable).parent / name
        if beside.is_file():
            return str(beside)
        found = shutil.which(name)
        if found:
            return found
        raise RuntimeError(
            f"{name} was not found. It ships with {comes_from}: "
            f"`pip install {comes_from}` into the same environment as DynaPIN.")

    @staticmethod
    def _pdb_splitmodel():
        return dynapin._executable("pdb_splitmodel", "pdb-tools")

    @staticmethod
    def _split_models(file, job_path):

        """ Creates folder 'models' and splits the trajectory into the models folder.
        
        Keyword arguments:
        file -- trajectory file
        Return: None
        """
        file_path = os.path.abspath(file)
        models_path = os.path.join(job_path, 'models')
        os.makedirs(models_path, exist_ok=True)

        # run it in models/ through cwd rather than chdir, which would leave the process
        # in the wrong directory if anything below raised
        result = subprocess.run([dynapin._pdb_splitmodel(), file_path], cwd=models_path,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True)
        produced = [f for f in os.listdir(models_path) if f.lower().endswith(".pdb")]
        if result.returncode != 0 or not produced:
            # this used to pass unnoticed: os.system's return code was ignored, so
            # models/ stayed empty and every CAPRI metric then quietly vanished from the
            # quality-control table behind a "expected for monomer systems" note
            raise RuntimeError(
                f"Splitting {file_path} into per-frame models produced nothing "
                f"(pdb_splitmodel exited {result.returncode}). "
                f"{result.stdout.strip()[:300]}")

    def _prepare_models(self, pdb_path, job_path):
        """Populate models folder from the final pdb (after chain selection/standardization).

        If the structure has multiple frames/models, split them.
        If it is a single-state structure, just copy it as model_00001.pdb.
        """
        models_path = os.path.join(job_path, 'models')
        if os.path.isdir(models_path):
            for f in os.listdir(models_path):
                fp = os.path.join(models_path, f)
                try:
                    if os.path.isfile(fp) or os.path.islink(fp):
                        os.remove(fp)
                    else:
                        shutil.rmtree(fp)
                except OSError:
                    continue

        u = mda.Universe(pdb_path)
        n_frames = len(u.trajectory)

        if n_frames > 1:
            self._split_models(pdb_path, job_path)
        else:
            os.makedirs(models_path, exist_ok=True)
            dest = os.path.join(models_path, "model_00001.pdb")
            shutil.copyfile(pdb_path, dest)

    @staticmethod
    def _preprocess_dcd(trajectory_file, output_file, stride, topology_file, sc=None, chains=None):


        """ Transforms input trajectory file into the .pdb file for given stride value.
        
        Keyword arguments:
        trajectory_file -- Input trajectory file
        output_file -- Output .pdb file
        stride -- Number of frames to be skipped.
        Return: None
        """
        u = mda.Universe(topology_file, trajectory_file)
        name = os.path.splitext(os.path.basename(output_file))[0]

        if chains is None:
            chains = []
        frames = list(range(0, len(u.trajectory), stride)) 

        with mda.Writer(output_file, multiframe=True) as W:
            for t in frames:
                u.trajectory[t]
                if sc:
                    sel = []
                    for i in sc:
                        sel.append(f"chainID {i.upper()}")

                    a = u.select_atoms(' or '.join(sel))
                else:
                    a = u.atoms
                
                W.write(a)

        ptm.main2(output_file, chains)

        return name + "_fixed.pdb"
        

    def run_quality_control(self, rmsd_data={'ref_struc':None, 'ref_frame':0}):
        """ The function to run Quality Control analyses in one line. Defines Quality Control class and runs the functions in it.
        
        Return: None
        """
        handler.test_rmsd_dict(rmsd_data)

        self.rmsd_data = rmsd_data
        self.qc_flag = True

        c = self.QualityControl(self.pdb_file, self.target_path, rmsd_data=rmsd_data, job_path=self.job_path,
                                stride=self.stride)
        c.quality_tbls_overres()
        
        dockq_df = pd.DataFrame()
        if len(c.segs) > 1:
            dockq_df = c.run_dockq()
        
        c.quality_tbls_overtime(dockq_df)


    def run_res_based(self, foldx_path=None):
        """ The function to run Residue Based analyses in one line. Defines Residue Based class and runs the functions in it.
        
        Return: None
        """
        self.foldx_path = foldx_path
        self.rb_flag = True
        self.dssp_flag = True
        c = self.ResidueBased(self.pdb_file, self.target_path, 
                              stride=self.stride, job_path=self.job_path, foldx_path=foldx_path, run_dssp=True)
        c.res_based_tbl()
        c.interface_table()


    def run_inter_based(self, get_all_hph=False):
        """The function to run Interaction Based analyses in one line. Defines Interaction Based class and runs the functions in it.
        
        Keyword arguments:
        get_all_hph -- True or False. Default is False. If True, the function calculates all the possible hydrophobic interactions.
        Return: None
        """
        self.ib_flag = True
        self.get_all_hph=get_all_hph
        c = self.InteractionBased(self.pdb_file, self.target_path, get_all_hph=get_all_hph, 
                                  stride=self.stride)
        c.int_based_tbl()

    
    class QualityControl:
        def __init__(self, pdb_file, target_path, job_path, rmsd_data, 
                     stride):
            """ A class to perform Quality Control analyses (RMSD, RG, and RMSF) of the trajectory. Initially runs private methods and defines RMSD, RG, and RMSF results of the trajectory.

            Return: None
            
            """
            self.stride = stride
            self.job_path = job_path
            self.u = mda.Universe(pdb_file, in_memory=True)

            handler.test_rmsd_refstruc(rmsd_data['ref_struc'])
            handler.test_rmsd_refframe(rmsd_data['ref_frame'], len(self.u.trajectory))
        
            self.target_path = target_path
            self.segs = self.u.segments
            self.protein = self.u.select_atoms('protein')

            self.rmsd_ref_struc = rmsd_data['ref_struc']
            self.rmsd_ref_frame = rmsd_data['ref_frame']
            self.rmsd_dict = self._calc_rmsd()
            self.rg_dict = self._calc_rg()
            self.rmsf_dict = self._calc_rmsf()

            self.header_overtime = ["Frame"]
            self.header_overres = ["Molecule", "Residue Number", "RMSF"]

        def _calc_rmsd(self):
            """ Calculates RMSD of each chain and the complex during the simulation. """
            n_frames = len(self.u.trajectory)
            if n_frames == 0:
                warnings.warn("No frames found in trajectory; RMSD values are not calculated.")
                return {}
            
            analysis = {}
            # Define the general selection for backbone (Protein + DNA/RNA)
            backbone_sel = "(backbone or nucleicbackbone)"

            # Load the reference universe ONCE outside the loop to save memory/time
            if self.rmsd_ref_struc:
                ref_universe = mda.Universe(self.rmsd_ref_struc)

                for chain in self.segs:
                    chainid = chain.segid
                    
                    current_chain_sel = f'(segid {chainid} or chainid {chainid}) and {backbone_sel}'
                    backbone_traj = self.u.select_atoms(current_chain_sel)
                    backbone_ref = ref_universe.select_atoms(current_chain_sel)

                    if backbone_traj.n_atoms == 0:
                        warnings.warn(f"No backbone atoms found for chain '{chainid}'. Skipping RMSD.")
                        continue

                    if backbone_traj.n_atoms != backbone_ref.n_atoms:
                        warnings.warn(f"Atom mismatch for Chain {chainid}: Traj={backbone_traj.n_atoms}, Ref={backbone_ref.n_atoms}. Skipping.")
                        continue
                    R1 = RMSD(backbone_traj, reference=backbone_ref, ref_frame=self.rmsd_ref_frame)
                    R1.run()
                    analysis[f"Chain {chainid}"] = [f"{x:.03f}" for x in R1.results.rmsd.T[2]]

                complex_traj = self.u.select_atoms(backbone_sel)

                complex_ref = ref_universe.select_atoms(backbone_sel)

                if complex_traj.n_atoms == 0:
                    warnings.warn("No backbone atoms found for the complex. Complex RMSD is not calculated.")
                elif complex_traj.n_atoms != complex_ref.n_atoms:
                    warnings.warn(f"Atom mismatch for Complex: Traj={complex_traj.n_atoms}, Ref={complex_ref.n_atoms}. Complex RMSD skipped.")
                else:
                    R1 = RMSD(complex_traj, reference=complex_ref, ref_frame=self.rmsd_ref_frame)
                    R1.run()
                    analysis["Complex"] = [f"{x:.03f}" for x in R1.results.rmsd.T[2]]
            else:
                 for chain in self.segs:
                     chainid = chain.segid
                     backbone = self.u.select_atoms(f'(segid {chainid} or chainid {chainid}) and {backbone_sel}')
                     if backbone.n_atoms == 0:
                         warnings.warn(f"No backbone atoms found for chain '{chainid}'. Skipping RMSD for this chain.")
                         continue
                     R1 = RMSD(backbone, reference=self.rmsd_ref_struc, ref_frame=self.rmsd_ref_frame)
                     R1.run()
                     analysis[f"Chain {chainid}"] = [f"{x:.03f}" for x in R1.results.rmsd.T[2]]
    
                 backbone = self.u.select_atoms(backbone_sel)
                 if backbone.n_atoms == 0:
                     warnings.warn("No backbone atoms found for the complex. Complex RMSD is not calculated.")
                 else:
                     R1 = RMSD(backbone, reference=self.rmsd_ref_struc, ref_frame=self.rmsd_ref_frame)
                     R1.run()
                     analysis["Complex"] = [f"{x:.03f}" for x in R1.results.rmsd.T[2]]                

            return analysis

        def _calc_rg(self):
            """Calculates RG of each chain and the complex during the simulation.
            
            Return: Dict: A dictionary that contains the chains and complexes as keys and their corresponding RMSD values as values.
            """
            n_frames = len(self.u.trajectory)
            if n_frames == 0:
                warnings.warn("No frames found in trajectory; RG values are not calculated.")
                return {}
            
            analysis = {}
            proteins = {}
            macro_sel = "(protein or nucleic)"
            for chain in self.segs:
                chainid = chain.segid
                ag = self.u.select_atoms(f"(segid {chainid} or chainid {chainid}) and {macro_sel}")
                if ag.n_atoms == 0:
                    print(f"[warn] No atoms found for chain '{chainid}'. "
                          f"Available segids: {list(set(self.u.atoms.segids))}")
                    continue
                proteins[chainid] = ag
                analysis[f"Chain {chainid}"] = []

            complex_ag = self.u.select_atoms(macro_sel)
            if complex_ag.n_atoms > 0:
                proteins["Complex"] = complex_ag
                analysis["Complex"] = []
            else:
                warnings.warn("No protein or nucleic atoms found for complex RG calculation; skipping Complex RG.")

            for _ in self.u.trajectory:
                for key, protein in proteins.items():
                    if key == "Complex":
                        analysis["Complex"].append(f"{protein.radius_of_gyration():.03f}")
                    else:
                        analysis[f"Chain {key}"].append(f"{protein.radius_of_gyration():.03f}")
            return analysis

        def _calc_rmsf(self):
            """ Calculates the RMSF of each residue for all simulation time.

            Return: Dict: A dictionary that contains the chain ids as keys and a list that has residue number in the index 0 and the corresponding RMSF values in the index 1 as values.
            """
            n_frames = len(self.u.trajectory)
            if n_frames == 0:
                warnings.warn("No frames found in trajectory; RMSF values are not calculated.")
                return {}
            
            analysis = {}
            for chain in self.segs:
                chainid = chain.segid
                analysis[chainid] = []

                backbone_sel = "(backbone or nucleicbackbone)"
                backbone_atoms = self.u.select_atoms(f'(segid {chainid} or chainid {chainid}) and {backbone_sel}')
                if backbone_atoms.n_atoms == 0:
                    print(f"[warn] No backbone atoms found for chain '{chainid}'. Skipping in RMSF.")
                    continue

                # allow single-frame inputs by returning zeros instead of skipping
                if n_frames == 1:
                    rmsf_vals = np.zeros(backbone_atoms.n_atoms)
                else:
                    rmsf = RMSF(backbone_atoms)
                    rmsf.run()
                    rmsf_vals = rmsf.results.rmsf

                resnums_np = np.array(backbone_atoms.resnums)
                unique_res = np.unique(resnums_np)
                per_res_rmsf = []
                for res in unique_res:
                    mask = resnums_np == res
                    per_res_rmsf.append(f"{np.mean(rmsf_vals[mask]):.03f}")

                analysis[chainid].append(unique_res)
                analysis[chainid].append(per_res_rmsf)

            return analysis

        
        def quality_tbls_overtime(self, dockq_df=None):
            """
            Build the over-time QC table (RMSD + RG) and, if provided, merge DockQ metrics.
            Writes a single 'qc_trajectory_metrics.csv' under self.target_path.
            """
            import os
            import pandas as pd
            import numpy as np

            key="Frame"

            if not hasattr(self, "rmsd_dict") or not hasattr(self, "rg_dict"):
                raise RuntimeError("rmsd_dict/rg_dict are not available; run QualityControl steps first.")

            def _series_len(vals):
                try:
                    return max(len(v) for v in vals.values() if v is not None)
                except ValueError:
                    return 0

            n_rows = max(_series_len(getattr(self, "rmsd_dict", {})), _series_len(getattr(self, "rg_dict", {})))
            
            qc_df = pd.DataFrame()

            if n_rows > 0:
                key_vals = np.arange(n_rows, dtype=int) * self.stride

                cols = [key]
                data_cols = {}
                all_names = sorted(set(self.rmsd_dict.keys()) | set(self.rg_dict.keys()))
                for name in all_names:
                    rmsd_vals = self.rmsd_dict.get(name)
                    rg_vals   = self.rg_dict.get(name)
                    if (rmsd_vals is None or len(rmsd_vals) == 0) and (rg_vals is None or len(rg_vals) == 0):
                        continue
                    cols.extend([f"{name} RMSD", f"{name} RG"])
                    def _pad(arr):
                        if arr is None:
                            return np.full(n_rows, np.nan)
                        arr_np = np.asarray(arr, dtype=float)
                        if len(arr_np) < n_rows:
                            pad = np.full(n_rows - len(arr_np), np.nan)
                            arr_np = np.concatenate([arr_np, pad])
                        return arr_np
                    data_cols[f"{name} RMSD"] = _pad(rmsd_vals)
                    data_cols[f"{name} RG"]   = _pad(rg_vals)

                qc_df = pd.DataFrame({key: key_vals, **data_cols}, columns=cols)

            if dockq_df is not None and not dockq_df.empty:
                if qc_df.empty:
                    qc_df = dockq_df
                else:
                    dq = dockq_df.copy()

                    sort_key = None
                    for cand in ("Frame", "frame"):
                        if cand in dq.columns:
                            sort_key = cand
                            break
                    if sort_key is not None:
                        dq = dq.sort_values(by=sort_key, kind="stable", na_position="last")

                    extra = [c for c in dq.columns if c != key and c not in qc_df.columns]
                    if extra:
                        qc_df = pd.merge(qc_df, dq[[key] + extra].drop_duplicates(subset=[key]), on=key, how="outer")

            out_path = os.path.join(self.target_path, "qc_trajectory_metrics.csv")
            os.makedirs(self.target_path, exist_ok=True)
            qc_df.to_csv(out_path, index=False)

        def run_dockq(self):
            """
            Run DockQ per frame and return a pandas DataFrame with metrics.
            No separate CSV is written by this function.
            """

            models_dir = os.path.join(self.job_path, "models")

            rx = re.compile(r'(\d+)\.pdb$', re.IGNORECASE)
            def _frame(fname):
                m = rx.search(fname)
                return int(m.group(1)) if m else None

            pdbs = [(f, _frame(f)) for f in os.listdir(models_dir) if f.lower().endswith(".pdb")]
            pdbs = [(f, n) for f, n in pdbs if n is not None]
            pdbs.sort(key=lambda x: x[1])
            if not pdbs:
                return pd.DataFrame()

            native_file, native_frame = pdbs[0]             
            native = os.path.join(models_dir, native_file)

            rows = [] 

            proc = subprocess.Popen(
                [dynapin._executable("DockQ", "DockQ"), native, native, "--short"],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            out, _ = proc.communicate()
            out = out or ""

            line = next((ln for ln in out.splitlines() if ln.startswith("DockQ")), None)
            if line:
                parts = line.split()
                def _f(x):
                    try: return float(x)
                    except: return None
                rows.append({
                    "Frame": native_frame - 1,          # 0 for first frame
                    "DockQ": _f(parts[1]),
                    "iRMSD": _f(parts[3]),
                    "lRMSD": _f(parts[5]),
                    "fnat":  _f(parts[7]),
                    "fnonnat": _f(parts[9]),
                    **({"clashes": _f(parts[13])} if len(parts) >= 14 and parts[12] == "clashes" else {}),
                    "DockQ_mapping": (parts[15])
                })


            for fname, frame in pdbs[1:]:
                mobile = os.path.join(models_dir, fname)
                proc = subprocess.Popen(
                    [dynapin._executable("DockQ", "DockQ"), mobile, native, "--short"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                out, _ = proc.communicate()
                out = out or ""

                line = next((ln for ln in out.splitlines() if ln.startswith("DockQ")), None)
                if not line:
                    continue
                parts = line.split()
                def _f(x):
                    try: return float(x)
                    except: return None
                rows.append({
                    "Frame": frame - 1,
                    "DockQ": _f(parts[1]),
                    "iRMSD": _f(parts[3]),
                    "lRMSD": _f(parts[5]),
                    "fnat":  _f(parts[7]),
                    "fnonnat": _f(parts[9]),
                    **({"clashes": _f(parts[13])} if len(parts) >= 14 and parts[12] == "clashes" else {}),
                    "DockQ_mapping": (parts[15])
                })


            df = pd.DataFrame(rows)
            if df.empty:
                return df
            df = df.sort_values(by=["Frame"], kind="stable", na_position="last")
            df = df.reset_index(drop=True)
            df["Frame"] = df.index * self.stride

            return df
                

        def quality_tbls_overres(self):
            """Reads RMSF results and writes them to the csv file with respect to chain id and residue number.

            Return: None
            """

            file = open(os.path.join(self.target_path, "qc_residue_rmsf.csv"), "w")

            print(",".join(self.header_overres), end="\n", file=file)

            for key, value in self.rmsf_dict.items():
                if not value or len(value) < 2:
                    warnings.warn(f"No RMSF data for {key}; skipping residue-level QC output for this chain.")
                    continue
                for index, val in enumerate(value[0]):
                    rms = value[1][index]

                    row = f"{key},{val},{rms}\n"

                    file.write(row)

            file.close()

    class ResidueBased:
        def __init__(self, pdb_path, target_path, 
                     stride, job_path, foldx_path, run_dssp):
            """ A class to perform Residue Based analyses such as core-rim, and  biophysical type classifications; Van der Waals, electrostatic, desolvation, and hydrogen bond energies either between same-chain or different chain residues. Also prints the interface class that residues had with the highest percentage for all simulation time. Initially calculates rASA values and residue energies.
            
            Keyword arguments:
            pdb_path -- Input .pdb file
            target_path -- Path to 'tables' folder
            Return: None
            """
            self.stride = stride
            self.foldx_path=foldx_path
            self.run_dssp = run_dssp               

            self.target_path = target_path
            self.header = ["Frame", "Chain", "Residue", "Residue Number","rASAc", "rASAm", "delta rASA", "SASA",
                           "Interface Label", "Residue Biophysical Type"]

            if self.foldx_path:
                self.header.extend(["Backbone Hbond Energy","Sidechain Hbond Energy", "Van der Waals Energy", "Electrostatic Energy",
                                    "Total Residue Energy"])

            if self.run_dssp:
                if platform == 'win32':
                    handler.check_wsl()
                
                handler.check_dssp(platform)
                self.header.append('Secondary Structure')

                self._run_dssp(job_path=job_path)
                self.dssp_data = self._run_dssp(job_path=job_path)
            
            self.header.append("\n")

            self.rasam_array = freesasa.structureArray(pdb_path,
                                                       {'separate-chains': True,
                                                        'separate-models': True})  
            self.rasac_array = freesasa.structureArray(pdb_path,
                                                       {'separate-chains': False,
                                                        'separate-models': True})  

            if self.foldx_path:
                self.energies = self._res_en(pdb_path, job_path, foldx_path)
            else:
                self.energies = {}
                warnings.warn("--foldx_path not provided. Skipping residue energy calculations.")


        @staticmethod
        def _calc_interface_class(delta_ras, rasc, rasm):
            """ Residues are classified according to the their rASA (relative accessible solvent area) values for both monomer (rASAm) and complex (rASAc) conformations by Levy et al. For more information about core-rim classification, please visit https://doi.org/10.1016/j.jmb.2010.09.028.

            This function classifies the residue according to the given rASA values. The classification labels are:
                0=Interior, 1=Surface, 2=Support, 3=Rim, 4=Core

            Keyword arguments:
            delta_ras -- rASAm - rASAc
            rasc -- rASA value for complex conformation
            rasm -- rASA value for monomer conformation
            Return: int: Interface Class label.
            """
            label = 0  # label will be used to extract static/dynamic interfaces
            if delta_ras == 0:
                if rasc < 0.25:
                    label = 0
                else:
                    label = 1

            if delta_ras > 0:
                if rasm < 0.25:
                    label = 2

                elif rasc > 0.25:
                    label = 3

                elif rasm > 0.25 > rasc:
                    label = 4

            return label

        @staticmethod
        def _calc_biophys_type(res_type):
            """ Classifies residues as +/-ly charged, hydrophobic or polar. Output labels are:
                hydrophobic: 0
                +ly charged: 1
                -ly charged: 2
                polar: 3

            Keyword arguments:
            res_type -- 3 letter residue code.
            Return: int: Biophysical Type label.
            """
            res_dict = {"GLY": 0, "ALA": 0, "PRO": 0, "VAL": 0, "LEU": 0, "ILE": 0, "MET": 0, "TRP": 0, "PHE": 0,
                        "SER": 3,
                        "THR": 3, "TYR": 3, "ASN": 3, "GLN": 3, "CYS": 3, "LYS": 1, "ARG": 1, "HIS": 1, "ASP": 2,
                        "GLU": 2}

            return res_dict[res_type]

        @staticmethod
        def _res_en(trajectory_file, job_path, foldx_path):
            """ Calculates the residue energies (Van der Waals, electrostatic, desolvation, and hydrogen bond for both same chain and different chain interactions) by running EvoEF1 for each frame.
            
            Keyword arguments:
            trajectory_file -- Input .pdb file
            Return: dict: Dictionary-in-dictionary with the order frame_number-residue-energy_values.
            """
            

            class ResEnergy:
                def __init__(self, name, resnum, chain, total, bb_hbond, sc_hbond, vdw, elec):
                    """A class to collect energy values of a residue. This class is defined for each residue. The energy values that ends with the letter 's' describes same-chain interaction, while letter 'd' stands for different chain interations.
                    
                    Keyword arguments:
                    vdwatt -- Van der Waals attraction energy
                    vdwrep -- Van der Waals repulsion energy
                    elec -- Electrostatic energy
                    HB -- description
                    desH -- Hydrophobic desolvation energy
                    desP -- Polar desolvation energy

                    Return: None
                    """
                    
                    self.res_name = name
                    self.resnum= resnum
                    self.chain = chain
                    self.total = total
                    self.bb_hbond = bb_hbond
                    self.sc_hbond = sc_hbond
                    self.vdw = vdw
                    self.elec = elec

                    self.hbond = self.bb_hbond + self.sc_hbond
            
            result = dict()

            # run FoldX in a controlled location
            inp_path = os.path.abspath(trajectory_file)
            current = os.getcwd()
            foldx_exe_path = foldx_path
            models_path = os.path.join(job_path, 'models')

            with tempfile.TemporaryDirectory() as output_path:
                try:
                    if not os.path.exists(models_path):
                        dynapin._split_models(inp_path, job_path)
                    os.chdir(models_path)

                    for i in os.listdir(models_path):
                        if i.split('.')[-1] == 'pdb' and "foldx" not in i:
                            subprocess.run([f"{foldx_exe_path}", '--command=SequenceDetail', f'--output-dir={output_path}', f'--pdb={i}'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                finally:
                    os.chdir(current)

                for file in os.listdir(output_path):
                    try:
                        frame = int(file.split('.')[0].split('_')[-2]) - 1
                    except ValueError:
                        frame = int(file.split('.')[0].split('_')[-1]) - 1
                    if frame not in result.keys():
                        result[frame] = dict()
                    with open(os.path.join(output_path, file), 'r') as fh:
                        for r in fh:
                            r = r.rstrip("\n")
                            splitted = r.split("\t")

                            try:
                                resname = splitted[1]
                                chain = splitted[2]
                                res_num = splitted[3]

                                total = float(splitted[8])
                                bb_hbond = float(splitted[9])
                                sc_hbond = float(splitted[10])
                                vdw = float(splitted[11])
                                elec = float(splitted[12])

                                res = ResEnergy(resname, res_num, chain, total, bb_hbond, sc_hbond, vdw, elec)

                                result[frame][f"{res.chain}{res.resnum}{res.res_name}"] = res

                            except IndexError:
                                continue

            return result
        
        @staticmethod
        def _run_dssp(job_path):
            """Run mkdssp for every model in parallel-friendly batches."""

            def extract_frame(filename):
                return int(filename.split('_')[-1].split('.')[0])

            models_path = os.path.join(job_path, 'models')
            if not os.path.isdir(models_path):
                raise FileNotFoundError("models directory not found; run split_models first.")

            pdb_files = [f for f in os.listdir(models_path) if f.lower().endswith('.pdb')]
            pdb_files.sort(key=extract_frame)

            records = []
            with tempfile.TemporaryDirectory(dir=job_path) as scratch:
                for file in pdb_files:
                    file_path = os.path.join(models_path, file)
                    file_name = os.path.splitext(file)[0]
                    frame = extract_frame(file)

                    tmp_input = None
                    with open(file_path, 'r', errors='ignore') as fh:
                        first_line = fh.readline()
                        if not first_line.startswith('CRYST1'):
                            fh.seek(0)
                            tmp_input = os.path.join(scratch, f"{file_name}_cryst1.pdb")
                            with open(tmp_input, 'w', newline='') as wh:
                                wh.write('CRYST1\n')
                                shutil.copyfileobj(fh, wh)
                    in_path = tmp_input or file_path

                    out_path = os.path.join(scratch, f"{file_name}.dssp")

                    cmd = ['mkdssp', '-i', in_path, '-o', out_path]
                    if platform == 'win32':
                        in_ws = in_path.replace('C:', '//mnt//c').replace('\\', '//')
                        out_ws = out_path.replace('C:', '//mnt//c').replace('\\', '//')
                        cmd = ['wsl', 'mkdssp', '-i', in_ws, '-o', out_ws]

                    subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

                    with open(out_path, 'r', errors='ignore') as fh:
                        parse_section = False
                        for line in fh:
                            if parse_section:
                                resnum = line[6:11].strip()
                                try:
                                    int(resnum)
                                except ValueError:
                                    continue
                                chain = line[11].strip()
                                resname = line[13].strip()
                                sec = line[16] if line[16] != " " else 'Null'
                                records.append((frame, resnum, chain, resname, sec))
                                continue
                            if line.startswith('  #'):
                                parse_section = True

            return pd.DataFrame.from_records(records, columns=['frame', 'resnum', 'chain', 'resname', 'secstruc'])
   
        def res_based_tbl(self):
            """Merges the residue based analyses results and writes to .csv file.
            
            
            Return: None
            """
            
            file = open(os.path.join(self.target_path, "res_trajectory_props.csv"), "w")  # open the output file

            file.write(",".join(self.header))

            rasac_array = [freesasa.calc(i).residueAreas() for  i in self.rasac_array]

            chain_len = len(rasac_array[0])
            rasm_array = list()
            for i in range(0, len(self.rasam_array), chain_len):
                group = self.rasam_array[i:i + chain_len]
                merged_dict = {}
                for obj in group:
                    dict_from_obj = freesasa.calc(obj).residueAreas()
                    merged_dict.update(dict_from_obj) 
                rasm_array.append(merged_dict)

            for frame, frame_obj in enumerate(rasac_array):  
                for chain, chain_mol in frame_obj.items():  
                    for resnum, val in chain_mol.items():
                        rasc = val.relativeTotal
                        rasm = rasm_array[frame][chain][resnum].relativeTotal
                        rasm_total = rasm_array[frame][chain][resnum].total
                        delta_ras = rasm - rasc
                        nbsa = (rasm * rasm_total) / 100
                        label = self._calc_interface_class(delta_ras, rasc, rasm)
                        try:
                            biophy_class = self._calc_biophys_type(val.residueType)
                        except:
                            biophy_class = None
                        
                        t = frame * self.stride
                        
                        row_list = [f"{round(t,3)}", chain, val.residueType, str(resnum), f'{rasc:.03f}', f'{rasm:.03f}', f'{delta_ras:.03f}', f'{nbsa:.03f}', str(label), str(biophy_class)]
                        
                        if self.foldx_path:
                            name = f"{chain.upper()}{val.residueNumber}{val.residueType}"
                            try:
                                obj = self.energies[frame][name]
                                row_list.extend([f'{obj.bb_hbond:.03f}',f'{obj.sc_hbond:.03f}',f'{obj.vdw:.03f}',f'{obj.elec:.03f}',f'{obj.total:.03f}'])
                            except KeyError:
                                continue

                        if self.run_dssp:
                            sec_structure = self.dssp_data.loc[(self.dssp_data['frame']==int(frame)+1)&(self.dssp_data['resnum']==str(resnum))&(self.dssp_data['chain']==chain), 'secstruc']
                            mlist = sec_structure.tolist()
                            try:
                                row_list.append(mlist[0])
                            except IndexError:
                                row_list.append("Null")
                        
                        file.write(",".join(row_list) + "\n")

            file.close()


        def interface_table(self):
            """Reads residue based csv file and writes the highest precentage of interface label of whole simulation for each residue. 
            
            Return: None
            """
            
            df = pd.read_csv(os.path.join(self.target_path, "res_trajectory_props.csv"),
                             usecols=["Chain", "Residue", "Residue Number", "Interface Label"])
            df["Residue Name"] = [a + str(b) for a, b in zip(df["Residue"], df["Residue Number"])]
            df["Interface Label"] = pd.to_numeric(df["Interface Label"], errors="coerce")
            df = df.dropna(subset=["Interface Label"])
            records = []
            label_ids = [0, 1, 2, 3, 4]
            label_names = {0: "Interior", 1: "Surface", 2: "Support", 3: "Rim", 4: "Core"}

            grouped = df.groupby(["Residue Name", "Chain"])
            for (res_name, chain), data in grouped:
                total = len(data)
                if total == 0:
                    continue
                percentages = {label: 0.0 for label in label_ids}
                counts = data["Interface Label"].value_counts()
                for label, count in counts.items():
                    int_label = int(label)
                    percentages.setdefault(int_label, 0.0)
                    percentages[int_label] = (count / total) * 100
                interface_score = sum(percentages.get(lbl, 0.0) for lbl in [2, 3, 4])
                row = {"Chain": chain, "Residue": res_name}
                for lbl in label_ids:
                    column_name = label_names[lbl]
                    row[column_name] = format(percentages.get(lbl, 0.0), ".2f")
                row["Interface_score"] = format(interface_score, ".2f")
                records.append(row)
            df2 = pd.DataFrame(records)
            if not df2.empty:
                df2.sort_values(by=["Chain", "Residue"], inplace=True)
            df2.to_csv(os.path.join(self.target_path, "res_interface_stats.csv"), index=False, mode="w+")

    class InteractionBased:
        def __init__(self, pdb_path, target_path, 
                      stride, get_all_hph=False):
            """ A class to calculate Hydrogen, Electrostatic, and Hydrophobic interactions with atom informations that contribute to the interaction, and writes interaction information for all frames to a .csv file.
            
            Keyword arguments:
            pdb_path -- Input .pdb file.
            target_path -- Path to the 'tables' fodler.
            get_all_hph -- Boolean. If True, prints all possible hydrophobic interactions. By default, prints only the closest interaction. 
            Return: None
            """

            self.stride = stride
            self.hph_status = get_all_hph
            self.pdb_path = pdb_path
            self._struct = app.PDBFile(pdb_path)
            self.target_path = target_path

        def _calc_interactions(self):
            """Calculates the interactions by Interfacea package.
            
            Return: list: A list of lists such as [interaction_table, frame_number].
            """
            positions = getattr(self._struct, "_positions", None)
            if positions is None:
                positions = []
            if len(positions) == 0 and getattr(self._struct, "positions", None) is not None:
                positions = [self._struct.positions]

            return_list = []
            for i, pos in enumerate(positions):
                self._struct.positions = pos
                self._struct.topology.createDisulfideBonds(self._struct.positions)
                mol = ia.Structure("pp.pdb", self._struct)
                
                analyzer = ia.InteractionAnalyzer(mol)
                analyzer.get_hbonds(strict=True)
                analyzer.get_hydrophobic()
                analyzer.get_ionic()
                table = analyzer.itable._table

                table.drop_duplicates(inplace=True)

                cha = table["chain_a"].tolist()
                chb = table["chain_b"].tolist()
                resa = table["resname_a"].tolist()
                resb = table["resname_b"].tolist()
                ida = table["resid_a"].tolist()
                idb = table["resid_b"].tolist()

                pairwises = [f"{z}{x}.{c}-{v}{b}.{n}" for z, x, c, v, b, n in zip(resa, ida, cha, resb, idb, chb)]
                resnamesa = [f"{a}{b}" for a, b in zip(resa, ida)]
                resnamesb = [f"{a}{b}" for a, b in zip(resb, idb)]

                table["pairwise"] = pairwises
                table.drop(["resname_a", "resname_b", "resid_a", "resid_b"], axis=1, inplace=True)

                table.insert(3, "residue_a", resnamesa)
                table.insert(4, "residue_b", resnamesb)

                return_list.append((analyzer.itable._table, i))  
            return return_list

        def int_based_tbl(self):
            """Reads list from _calc_interactions function and writes to .csv file with respect to the frame number.
            
            Return: None
            """
            
            handle = open(os.path.join(self.target_path, "int_pairwise_trajectory.csv"), "w+", newline="")
            return_list = self._calc_interactions()
            if len(return_list) == 0:
                warnings.warn("No coordinates available to compute interactions; interaction table will remain empty.")
                handle.close()
                return
            t = 'Frame,'
            handle.write(t + ",".join(return_list[0][0].columns.tolist()) + "\n")

            for tup in return_list:
                df = tup[0]
                frame = tup[1]

                frame = frame * self.stride

                df.insert(0, t.rstrip(','), round(frame,3))
                df.to_csv(handle, index=False, mode="w+", header=False)
            handle.close()
