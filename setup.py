from setuptools import setup, find_packages

setup(
    name="DynaPIN",
    version="1.0.0",
    packages=find_packages(include=["DynaPIN", "DynaPIN.*"]),
    python_requires=">=3.10",
    # interfacea is not on PyPI: the copy DynaPIN uses is vendored in interfacea_src and
    # installed from there by environment.yml ("pip: -e ./interfacea_src"). Naming it here
    # only made `pip install .` fail on a package that cannot be resolved.
    install_requires=[],
    entry_points={
        "console_scripts": [
            "dynapin=DynaPIN.dynapin:main",
        ],
    },
    include_package_data=True,
)
