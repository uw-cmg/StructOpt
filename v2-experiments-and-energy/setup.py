#!/usr/bin/env python

import sys

from setuptools import setup, find_packages


install_requires = [
    "ase",
    "numpy",
    "scipy",
    #"mpi4py"
]

setup(
    name="StructOpt",
    version="0.5",
    author="UW Madison, Materials Science",
    author_email="maldonis@wisc.edu",
    url="https://github.com/uw-cmg/StructOpt_modular",
    description="Atomic Stucture Optimization Framework",
    keywords=["Atomic", "Stucture", "Optimization", "Materials Science"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Environment :: Web Environment",
        "Operating System :: Linux",
        "Intended Audience :: Developers",
    ],
    packages=find_packages(),
    zip_safe=True,
    install_requires=install_requires
)
