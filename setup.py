from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='Sheriff',
    version='1.1.3',
    author='Brad Balderson, Michael Lorenzini, Aaron Ho',
    author_email='bbalderson@salk.edu',
    packages=find_packages(),
    license='BSD',
    long_description_content_type="text/markdown",
    long_description=long_description,
    scripts=['bin/sheriff'],
    install_requires=[
        'pysam>=0.19.0',
        'gtfparse==2.5.0',
        'pyranges>=0.0.111',
        'biopython==1.81',
        'numpy==1.26.4',
        'pandas>=1.5.0',
        'polars>=0.18.0',
        'scipy>=1.10.0',
        'numba>=0.56.0',
        'faiss-cpu==1.10.0',
        'typer>=0.9.0',
        'typing_extensions>=4.5.0',
    ],
    entry_points={
        'console_scripts': [
        'sheriff=sheriff.__main__:main',
        'Sheriff=sheriff.__main__:main'
        ]
    },
    python_requires='>=3.10',
    description=("Sheriff calls CRISPR/cas9 edit sites, and quantifies the number of edited alleles at each site along "
                 "with gene expression in single cells using Superb-seq data."),
    keywords=['Sheriff', 'sheriff', 'Superb-seq', 'superb-seq', 'CRISPR', 'single-cell', 'bioinformatics'],
    url='https://github.com/BradBalderson/Sheriff',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
    ],
)
