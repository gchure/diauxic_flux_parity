from setuptools import setup, find_packages
import pathlib

__version__ = "0.0.1"

# The directory containing this file
HERE = pathlib.Path(__file__).parent

setup(
    name="diaux_flux_parity",
    version=__version__,
    description="Python utilities for simulating growth of self-replicating communities",
    license="MIT",
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
    ],
    author="Griffin Chure",
    author_email="griffinchure@gmail.com",
    packages=find_packages(
        exclude=('docs', 'doc', 'sandbox', 'dev', 'diaux.egg-info')),
    include_package_data=True,
    install_requires=[
        "matplotlib>=3.7.0",
        "numpy>=1.24.3",
        "pandas>=1.5.3",
        "scipy>=1.10.0",
        "seaborn>=0.13.2",
        "tqdm>=4.64.1"]
)