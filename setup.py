"""Package setup."""

import setuptools

version = "1.0.0"

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="FIREcaller",
    version=version,
    author="Jakub Lipiński",
    author_email="jakub@cellular-genomics.com",
    description="Python library for detecting Frequently Interacting REgions (FIREs) from Hi-C data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cellular-genomics/python-FIREcaller",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires='>=3.6',
    install_requires=[
        "cooler>=0.8.7",
        "pandas>=1.0.1",
        "statsmodels>=0.11.0",
        "dask>=2.10.1",
        "fsspec>=0.6.2",
        "tables>=3.6.1",        
    ],
    entry_points = {
        'console_scripts': ['FIREcaller=FIREcaller.FIREcaller:main'],
    }    
)