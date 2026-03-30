from setuptools import setup, find_packages

# Read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# Read requirements.txt for install_requires
with open("requirements.txt") as f:
    install_requires = f.read().splitlines()

setup(
    name="wasp-proteins-annotation",
    version="1.1.1",
    description="WASP: A computational pipeline for protein functional annotation using structural information.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Giorgia Del Missier",
    author_email="delmissiergiorgia@gmail.com",
    url="https://github.com/gioodm/WASP",
    license="MIT",
    packages=find_packages(include=["wasp", "wasp.*"]),
    include_package_data=True,
    install_requires=install_requires,
    entry_points={
        "console_scripts": [
            "wasp-run=wasp.main.run:main",
            "wasp-gem=wasp.gem.run_gem:main",
            "find-orphans=wasp.utils.gem_preprocessing.find_orphans:main",
            "run-rxn2code=wasp.utils.gem_preprocessing.rxn2code:main",
            "afdb-download=wasp.utils.fasta2pdb.afdb_downloader:main",
        ]
    },
    python_requires=">=3.9",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
