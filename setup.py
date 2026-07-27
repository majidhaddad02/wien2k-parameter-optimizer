# -*- coding: utf-8 -*-
"""Install WIEN2k Parameter Optimizer."""

import sys
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_desc = fh.read()

setup(
    name="opt-wien2k",
    version="1.1.0",
    author="open2code",
    description="Automatic WIEN2k parameter optimizer — RMT, RKMAX, GMAX, LMAX, k-mesh, mixing, core/valence",
    long_description=long_desc,
    long_description_content_type="text/markdown",
    url="https://github.com/majidhaddad02/wien2k-parameter-optimizer",
    packages=find_packages(),
    py_modules=["optimize_wien2k"],
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "opt_wien2k = optimize_wien2k:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
    ],
)
