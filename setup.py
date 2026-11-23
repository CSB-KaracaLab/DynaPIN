from setuptools import setup, find_packages

setup(
    name="DynaPIN",
    version="1.0.0",
    packages=find_packages(include=["DynaPIN", "DynaPIN.*"]),
    python_requires=">=3.10",
    # What DynaPIN imports, minus interfacea: that one is not on PyPI, and the copy used
    # here is vendored in interfacea_src, which environment.yml installs with pip -e.
    # Requiring it by name only made `pip install .` fail on something unresolvable.
    install_requires=[
        "MDAnalysis",
        "freesasa",
        "matplotlib",
        "numpy",
        "openmm",
        "pandas",
        "seaborn",
    ],
    entry_points={
        "console_scripts": [
            "dynapin=DynaPIN.dynapin:main",
        ],
    },
    include_package_data=True,
)
