"""
Setup Model
"""

from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / "README.MD").read_text(encoding="utf-8")

setup(
    name="MethyladenosineFinder",
    version="1.0.0",
    description="Find methyladenosine sits in PacBio",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zhuweix/MethyladenosineFinder",
    author="Zhuwei Xu",
    author_email="wszwei@gmail.com",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords="methyladenosine, methylated adenine, PacBio, computational biology",
    package_dir={"": "src",},
    package_data={"src": ["asset/*.data", "asset/*.csv", "asset/*.png"]},
    packages=find_packages(where="src"),
    python_requires=">=3.7, <4",
    install_requires=["numpy",
                      'pysam',
                      'biopython'
                      ],
    include_package_data=True,

    entry_points={
        'console_scripts': [
            'm6aworkflow = workflow.workflow:main',
            'm6apreprocess = workflow.preprocess:main',
            'm6abamtobed = workflow.bamtobed:main',
            'm6aworkflowlocal = workflow.workflow_local:main'
        ], },
    project_urls={
        "Bug Reports": "https://github.com/zhuweix/MethyladenosineFinder/issues",
        "Source": "https://github.com/zhuweix/MethyladenosineFinder",
        "Lab": "https://www.nichd.nih.gov/research/atNICHD/Investigators/clark",
    },
)