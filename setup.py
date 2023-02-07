from setuptools import setup, find_namespace_packages
import pdbe_covariations

setup(
    name="pdbe_covariations",
    version=pdbe_covariations.__version__,
    description="Package to calculate covariation pairs",
    project_urls={
        "Source code": "https://gitlab.ebi.ac.uk/pdbe/release/bioint",
        "Documentation": "https://gitlab.ebi.ac.uk/pdbe/release/bioint",
    },
    author="Protein Data Bank in Europe",
    author_email="pdbehelp@ebi.ac.uk",
    license="Apache License 2.0.",
    keywords="PDB FASTA UniProt hhsuite covariations UniClust gremlin sequence ",
    packages=find_namespace_packages(),
    zip_safe=False,
    include_package_data=True,
    install_requires=[
        "numpy",
        "requests",
        "pandas",
	"biopython",
        "gemmi"
    ],
    extras_require={
        "tests": ["pytest", "pytest-cov"],
    },
    entry_points={
        "console_scripts": [
            "covariations=pdbe_covariations.covariations:main",
        ]
    },
)
