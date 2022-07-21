<img src="static/scims_logo.png" align=right width="25%">


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
 ``SCiMS``: **S**ex **C**alling from **M**etagenomic **S**equences        

#### `SCiMS` is a tool for identifying the sex of a host organism based on the alignment of metagenomic sequences. 

----
## Installation 

`SCiMS` can be easily installed with the pip package manager

```
pip3 install git+https://github.com/Kobie-Kirven/SCIMS-V1.0
```
 
To confirm that the instillation was successful, run:
```
scims -h
```
---
## Usage
`SCiMS` can be used on any alignment data, regardless of the platform used for sequencing or the aligner that generated the alignment file. 

```
usage: scims [-h] [-v] [--i INPUT] [--x HET] [--w WINDOW] [--dir OUT_DIR]
             [--pre PREFIX] [--scaffolds SCAFF]

Sex Calling from Metagenomic Sequences

optional arguments:
  -h, --help         show this help message and exit
  -v, --version      show program's version number and exit
  --i INPUT          Input SAM or BAM file
  --x HET            ID of heterogametic sex chromosome (ex. X)
  --w WINDOW         Window size for the coverage calculation
                     (default=10,000)
  --dir OUT_DIR      Ouput directory to hold the restuls
  --pre PREFIX       Prefix for the output files
  --scaffolds SCAFF  Scaffolds IDs to use in the analysis
```