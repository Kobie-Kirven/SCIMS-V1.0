#######################################################
# SCiMS: Sex calling for Metagenomic Sequences
#
#######################################################

#################
# Imports
#################
import pysam 
from os.path import exists
from errors import *
import numpy as np
import scipy.stats as st

def get_alignment_handle(file_name):
    """
    Get the alignment handle from a SAM or BAM formatted file

    Parameters:
        - file_name(str): Path to the SAM or BAM file

    Returns:
        - handle(pysam object) 
    """
    if not exists(file_name):
        raise DoesNotExist("The input file does not exist")
    else:
        try:
            handle = pysam.AlignmentFile(file_name, "rb", require_index=False)
        except:
            try:
                handle = pysam.AlignmentFile(file_name, "r", require_index=False)
            finally:
                raise WrongFile("This file does not look like a BAM file")
        return handle

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
        raise ValueError("Window size needs to be an integer")
    
    elif window_size <= 0:
        raise ValueError("You can not have a window size of 0")
    
    if length <= window_size:
        raise ValueError("Window size is too large, consider going smaller")
    
    return np.zeros(length // window_size)

def build_chrom_coverage_dict(chrom_names, lengths, window_size=10000):
    """
    Build a dictionary to hold the coverage information
    """
    if len(chrom_names) != len(lengths):
        raise ValueError("The number of chromosomes is not equal to the number of lengths")

    chrom_dict = {}
    for i in range(len(chrom_names)):
        chrom_dict[chrom_names[i]] = build_chrom_coverage_array(lengths[i], window_size)
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
        raise ValueError("The window size is greater than the read length. Please increase the window size")

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

