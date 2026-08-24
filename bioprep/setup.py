from setuptools import setup, find_packages

setup(
    name="bioprep",
    version="0.1.0",
    packages=find_packages(),
    # The web app lives inside the package so an installed copy can serve it.
    # find_packages() only collects .py files, so the templates and static
    # assets have to be listed explicitly or a pip install would produce a
    # Flask app with no templates to render.
    package_data={
        "bioprep": [
            "templates/*.html",
            "static/css/*.css",
            "static/js/*.js",
        ],
    },
    include_package_data=True,
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
            "bioprep-web=bioprep.webapp:main",
        ]
    },
    author="BioPrep Author",
    description="Prepare protein structures for molecular docking.",
)
