#######################################################
# SCiMS: Sex calling for Metagenomic Sequences
#
# Author: Kobie Kirven 
#
# Davenport Lab
# The Pennsylvania State University
#######################################################

#################
# Imports
#################
import pysam 
from os.path import exists
from .errors import *
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os
from pathlib import Path
import sys

def get_alignment_handle(file_name):
    """
    Get the alignment handle from a SAM or BAM formatted file

    Parameters:
        - file_name(str): Path to the SAM or BAM file

    Returns:
        - handle(pysam object) 
    """
    if not exists(file_name):
        print("Error: The alignment file does not exist")
        sys.exit(1)
    else:
        try:
            handle = pysam.AlignmentFile(file_name, "rb", require_index=False)
        except:
            try:
                handle = pysam.AlignmentFile(file_name, "r", require_index=False)
            finally:
                print("Error: This file does not look like a BAM file!")
                sys.exit(1)
        return handle

def read_scaffolds(file_prefix):
    try:
        with open(file_prefix) as fn:
            return [line.strip("\n") for line in fn.readlines() if (len(line.strip("\n")) > 1)]
            
    except:
        print("Error: Something is wrong with the scaffolds file")
        sys.exit(1)

def chrom_len_from_sam(handle):
    """
    Get the lengths of the chromosomes from the SAM header

    Parameters:
        - handle(pysam handle): output from get_alignment_handle
    
    Returns:
        - chrom_lengths(dict): dictionary where keys are chrom names
                                and values are lengths
    """
    chrom_lengths = {}
    for i in range(len(handle.references)):
        chrom_lengths[handle.references[i]] = int(handle.lengths[i])
    return chrom_lengths

def build_chrom_coverage_array(length, window_size=10000):
    """
    Build a dictionary to hold future coverage stats

    Parameters:
        - length(int): Length of the chromosome
        - window_size(int): window size for the coverage
    """
    if type(window_size) != int:
        print("Error: Window size needs to be an integer!")
        sys.exit(1)
    
    elif window_size <= 0:
        print("Error: You can not have a window size of 0")
        sys.exit(1)
    
    if length <= window_size:
        print("Error: Window size is too large, consider going smaller")
        sys.exit(1)
    
    return np.zeros(length // window_size)

def build_chrom_coverage_dict(chrom_lengths_dict, scaffolds, window_size=10000):
    """
    Build a dictionary to hold the coverage information
    """
    chrom_dict = {}
    for key in chrom_lengths_dict:
        if key in scaffolds:
            chrom_dict[key] = build_chrom_coverage_array(chrom_lengths_dict[key], window_size)
    return chrom_dict


def add_to_coverage_dict(coverage_array, start, align_len, window_size):
    """
    Add the number of nucleotides that cover a specific region
    """
    coverage_array = list(coverage_array)
    index = start // window_size
    second = (start + align_len) // window_size

    if second == index:
        coverage_array[int(index)] += align_len

    elif (second == index + 1) and second < len(coverage_array):
        for i in range(0, align_len):
            new_index = ((start + i) // window_size)
            if new_index > index:
                coverage_array[int(second)] += 1
            else:
                coverage_array[int(index)] +=  1

    elif second > (index + 1):
        print("Exit: The window size is greater than the read length. Please increase the window size")
        sys.exit(1)

    return np.asarray(coverage_array)

def decompose_sam_flag(flag):
    """
    Decompose SAM flag into its component parts
    Parameters:
        flag(int): Sam flag
    Returns:
        (list): Elements of sam flag
    """
    out_list = []
    flag_list = ["PAIRED", "PROPER_PAIR", "UNMAP", "MUNMAP",
                "REVERSE", "MREVERSE", "READ1", "READ2",
                "SECONDARY", "QCFAIL", "DUP", "SUPPLEMENTARY"]

    binary = str(f"{flag:b}"[::-1])
    for i in range(len(binary)):
        if binary[i] == "1":
            out_list.append(flag_list[i])
    return out_list

def check_alignment(alignment):
    """
    Check to make sure that the alignment passes the quality filters
    """
    bad = ["SECONDARY", "UNMAP", "MUNMAP", "SUPPLEMENTARY"]
    if not any(item in bad for item in decompose_sam_flag(alignment.flag)):
        return alignment.reference_name, alignment.reference_start, alignment.query_alignment_length

def get_ks_stat(array1, array2):
    """
    Compute the two-sample Kolmogorov-Smirnov test between two arrays
    """
    return st.ks_2samp(array1, array2)


def get_chrom_windows_coverage(handle, chrom_dict, window_size):
    """
    Get the coverage of each window across each chromosome

    Parameters:
        - handle(pysam handle): Handle of SAM or BAM file
        - chrom_dict(dict): output from  build_chrom_coverage_dict
        - window_size(int): size of the window to use

    """
    total = 0
    for rec in handle:
        if rec.reference_name in chrom_dict:
            align = check_alignment(rec)
            total += 1
            if align:
                chrom_dict[align[0]] = add_to_coverage_dict(chrom_dict[align[0]],align[1], align[2], window_size)
    return chrom_dict

def get_hom_het_lists(coverage_dict, heterogamtic_ids):
    """
    Get the window coverages of the homogametic and heterogamtic elements 
    """
    homogametic_coverage = []
    heterogametic_coverage = []

    for id in heterogamtic_ids:
        if id not in coverage_dict:
            print(f"Error: The heterogametic ID {id} is not valid")
            sys.exit(1)

    for key in coverage_dict:
        if key in heterogamtic_ids:
            heterogametic_coverage += list(coverage_dict[key])
        else:
            homogametic_coverage += list(coverage_dict[key])
    
    return homogametic_coverage, heterogametic_coverage

def plot_coverage_hist(homogametic, heterogametic, output):
    """
    Plot the window coverages for the homogametic and heterogametic elements
    """
    # Make it look pretty 
    plt.style.use('seaborn-colorblind')
    fig, ax = plt.subplots()
    plt.setp(ax.spines.values(), linewidth=2)
    ax.xaxis.set_tick_params(width=2)
    ax.yaxis.set_tick_params(width=2)
    plt.rc('font', weight='bold')
    plt.rc('xtick.major', size=5, pad=7)
    plt.rc('xtick', labelsize=15)

    # Make the plot
    plt.hist(homogametic, label="Homogametic", bins=100, density=True)
    plt.hist(heterogametic, label="Heterogametic", bins=100, density=True)
    plt.ylabel("Density", fontweight='bold')
    plt.xlabel("Coverage", fontweight='bold')
    plt.legend()
    plt.rcdefaults()
    plt.savefig(output, bbox_inches="tight")

def make_html_output(output, p_value):
    """
    Create the HTML for the report
    """
    with open(output + ".html", "w") as fn:
        fn.write('<img src="images/scims_logo.png" width="100">')
        fn.write('<h2>Sex Determination:</h2>')
        if float(p_value) < 0.05:
            sig = "is"
        else:
            sig = "is not"
        fn.write(f'<strong>Conclusion:</strong> Homogametic coverage <strong>{sig} significantly different</strong> to heterogametic coverage')
        fn.write(f'<p><strong>P-Value: {p_value}</strong></p>')
        fn.write('<h2>Coverage </h2>')
        fn.write('<img src="images/hist.jpg">')

def create_results_directory(location):
    try:
        os.mkdir(location)
        os.mkdir(location + "/images")
    
    except:
        print("Error: Can't make results directory")
        sys.exit(1)

    os.system(f"cp {Path(__file__).parent}/static/scims_logo.png {location}/images")