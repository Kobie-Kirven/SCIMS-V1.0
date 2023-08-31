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
import numpy as np
import pandas as pd
import statsmodels.api as sm

# Function to calculate Rt values
def calculate_Rt(idxstats, total_ref, total_map):
    return (idxstats[:, 1] / total_map) / (idxstats[:, 0] / total_ref)

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Sex Assignment Script')
    parser.add_argument('--i',dest= "idxtats", type=str, help='idxstats file')
    parser.add_argument('--s',dest="scaffold_ids_file", type=str, help='File containing scaffold IDs of interest')
    parser.add_argument('--X-id',dest="x_id", type=str, help='Scaffold ID for X chromosome')
    parser.add_argument('--Y-id',dest="y_id",type=str, help='Scaffold ID for Y chromosome')
    args = parser.parse_args()

    # Read scaffold IDs from file
    with open(args.scaffold_ids_file, 'r') as file:
        scaffold_ids = file.read().splitlines()

    # Read idxstats file
    idxstats = pd.read_table(args.idxtats, header=None, index_col=0)
    idxstats = idxstats.loc[scaffold_ids]

    c1 = np.array(idxstats.iloc[:, 0], dtype=float)
    c2 = np.array(idxstats.iloc[:, 1], dtype=float)
    total_ref = np.sum(c1)
    total_map = np.sum(c2)

    # Fit linear regression model
    X = sm.add_constant(c2)
    LM = sm.OLS(c1, X).fit()
    # print(LM.summary())

    print("=================================================")
    print("""
                                                    
    _|_|_|    _|_|_|  _|_|_|  _|      _|    _|_|_|  
    _|        _|          _|    _|_|  _|_|  _|        
    _|_|    _|          _|    _|  _|  _|    _|_|    
        _|  _|          _|    _|      _|        _|  
    _|_|_|      _|_|_|  _|_|_|  _|      _|  _|_|_|    
    """)
    print("=================================================")

    # Calculate Rt values
    Rt_values = calculate_Rt(idxstats.values, total_ref, total_map)

    # Calculate Rx
    copy = list(set(list(np.delete(Rt_values, scaffold_ids.index(args.x_id))) + list(np.delete(Rt_values, scaffold_ids.index(args.y_id)))))
    tot = Rt_values[scaffold_ids.index(args.x_id)] / copy
    Rx = np.mean(tot)
    print("Rx:", np.round(Rx, 3))

    # Calculate confidence interval for Rx
    conf_interval = (np.std(tot) / np.sqrt(22)) * 1.96
    CI1 = Rx - conf_interval
    CI2 = Rx + conf_interval
    print("95% CI:", np.round(CI1, 3), np.round(CI2, 3))

    # Calculate Ry
    tot = Rt_values[scaffold_ids.index(args.y_id)] / copy
    Ry = np.mean(tot)
    print("Ry:", np.round(Ry, 3))

    # Calculate confidence interval for Ry
    conf_interval = (np.std(tot) / np.sqrt(22)) * 1.96
    CI1 = Ry - conf_interval
    CI2 = Ry + conf_interval
    print("95% CI:", np.round(CI1, 3), np.round(CI2, 3))












