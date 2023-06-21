# Author: Kobie Kirven
# Davenport Lab - Penn State University
# Date: 9-2-2021

import setuptools
import sys
from pathlib import Path

CURRENT_DIR = Path(__file__).parent
sys.path.insert(0, str(CURRENT_DIR))

setuptools.setup(
    name="scims",
    version="0.01",
    packages=["scims", "scims"],
    include_package_data = True,
    entry_points={"console_scripts": ["scims=scims.__main__:main",],},
    description="SCiMS: Sex Calling for Metagenomic Sequences",
    install_requires=["setuptools", "biopython", "pandas", "matplotlib", "scipy", "numpy"],
    python_requires=">=3.6",
)