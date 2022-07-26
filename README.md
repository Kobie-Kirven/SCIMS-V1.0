[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<img src="static/scims_logo.png" align=right width="25%">

<h1>
 <strong>SCiMS</strong>: <strong><u>S</u></strong>ex <strong><u>C</u></strong>alling from <strong><u>M</u></strong>etagenomic <strong><u>S</u></strong>equences   
 </h1>     

 
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


## Input Files
### Scaffolds.txt
Since most assemblies include scaffolds representing other DNA than simply genomic (ex. mitochondrial), it is necessary to define what scaffolds we are interested in using for our analysis. This can be specified with a ```scaffolds.txt`` file. This is a single-column text file where each row is a scaffold ID. Here is an example, 
```
NC_000001.11
NC_000002.12
NC_000003.12
NC_000004.12
NC_000005.10
NC_000006.12
NC_000007.14
NC_000008.11
...
``` 
Several pre-compiled scaffolds lists are already available in SCiMS, including ```GRCh38```. 

## Output Files

The main output of the SCiMS program is an report like the one shown below.

<img src="static/scims_report_output.png">