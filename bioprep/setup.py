from setuptools import setup, find_packages

setup(
    name="bioprep",
    version="0.1.0",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "biopython>=1.80",
        "pdbfixer>=1.9",
        "openmm>=8.0",
        "flask",
        "werkzeug",
        "numpy",
        "scipy",
        "scikit-learn"
    ],
    entry_points={
        "console_scripts": [
            "bioprep=bioprep.cli:main",
        ]
    },
    author="BioPrep Author",
    description="A CLI tool for preparing protein structures for molecular docking.",
)
