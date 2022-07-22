#######################################################
# SCiMS: Sex calling for Metagenomic Sequences
#
# Author: Kobie Kirven 
#
# Davenport Lab
# The Pennsylvania State University
#######################################################

###################
# -- Arguments -- #
###################
import argparse
from .determine_sex import *
from pathlib import Path
import sys 

def scims():
    """
    Main function for the SCiMS program
    """

    parser = argparse.ArgumentParser(
        description="Sex Calling from Metagenomic Sequences"
    )

    parser.add_argument("-v", "--version", action="version", version="scims 0.0.1")

    ######################################
    # determine-sex
    ######################################
    parser.add_argument(
        "--i",
        dest="input",
        help="Input SAM or BAM file",
    )

    parser.add_argument(
        "--x",
        dest="het",
        help="ID of heterogametic sex chromosome (ex. X)",
    )

    parser.add_argument(
        "--w",
        dest="window",
        help="Window size for the coverage calculation (default=10,000)",
    )

    parser.add_argument(
        "--dir",
        dest="out_dir",
        help="Ouput directory to hold the results",
    )

    parser.add_argument(
        "--pre",
        dest="prefix",
        help="Prefix for the output files",
    )
    parser.add_argument(
        "--scaffolds",
        dest="scaff",
        help="Scaffolds IDs to use in the analysis",
    )



    args = parser.parse_args()

    #################
    # -- Pipeline --#
    #################
    print("\n")
    if args.scaff == "GRCh38":
        scaffolds = read_scaffolds(f"{Path(__file__).parent}/static/GRCh38_scaffolds.txt")
    else:
        try:
            scaffolds = read_scaffolds(args.scaff)
        except:
            print("Error: The scaffolds file does not exist or is invalid!")
            sys.exit(1)

    # Read in the alignment file
    handle = get_alignment_handle(args.input)

    # Get the chromosome lengths
    chrom_lengths = chrom_len_from_sam(handle)

    #Initiate chromosome coverage dict
    chrom_cov_dict = build_chrom_coverage_dict(chrom_lengths,scaffolds,  int(args.window))

    # Get the coverages
    coverages = get_chrom_windows_coverage(handle, chrom_cov_dict, int(args.window))

    # Split coverages into homogametic and heterogametic
    hom, het = get_hom_het_lists(coverages, args.het.split(","))

    # Get the p-value
    p_value = get_ks_stat(hom, het)[1]

    # Create the outputs
    create_results_directory(f"{args.out_dir}/{args.prefix}" )
    make_html_output(f"{args.out_dir}/{args.prefix}/{args.prefix}", p_value)
    plot_coverage_hist(hom, het, f"{args.out_dir}/{args.prefix}/images/hist.jpg")


    # Print Messages
    print(f"Results are in {args.out_dir}/{args.prefix}")
    print("Thank you for using SCiMS!")










