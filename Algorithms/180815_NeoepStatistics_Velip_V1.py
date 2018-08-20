#/usr/bin/python
import sys
import importlib

import pandas as pd

# sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/180815_NeoStream/Algorithms" )
from SVSv7 import NeoepStatistics
from mokhaPy import mokhaPy

#Constants - directories
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
# DIR_CURR = DIR_PROJ + "/180203_GenomAltToNeoepV3"
DIR_CURR = DIR_PROJ + "/180815_NeoStream"
DIR_DATA = DIR_CURR + "/Data"
DIR_RESULTS = DIR_CURR + "/Results"
DIR_RESULTS_VELIP = DIR_RESULTS + "/180203_Velip_V1"

DIR_DATA_VELIP = DIR_PROJ + "/170304_NeoantigenExtract/Data/Velip"
#get mapped reads
DIR_RNASEQ = DIR_PROJ + "/150802_TophatSamples"
#Constants - gene expression
DIR_EXPRESSION = DIR_PROJ + "/161001_SJExonDiff/Results"
DIR_GENOME = DIR_PROJ + '/ArchiveData/hg19.fa'      #directory for samtool-indexed genome

def create_hash_filter( frame_stat, neoep_stat, bool_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq ):
    """
    Creates a hash filter for the function quantify_all_alterations_v2()
    Args:
        -frame_stat = integer with the following status: (NOTE: Default = 0)
            -0 = will only consider both the alteration type (mutation, indel) (column = 'my_change_type')
            -1 = will only consider the reading frame (column = 'in_frame')
            -2 = will consider both the alteration type (mutation, indel) & frame-preservation (column = 'alt_readframe_type')
        -neoep_stat = integer that will perform a specific type of analysis with neoepitopes. Default = 0
            -0 = will consider ALL genomic alterations that produce any neoepitope
            -1 = will only consider high-affinity MHC neoepitopes that result from genomic alterations
            -2 = will only consider high-efficacious neoepitopes that result from genomic alterations
        -bool_NMD = boolean that considers Nonsense-mediate Decay (NMD) -> default is False
            -True = will only count events (genomic alterations) where NMD should not occur ()
            -False = will consider all events regardless if it is susceptible to NMD or not
        -thres_gene_exp = number (integer or float) that serves as gene expression threshold, but only if the expression level is greater than -1. Default = -1
        -thres_prot = number (integer or float) that serves as the proteasome percentile threshold, but only if the expression level is greater than -1. Default = -1
        -thres_tap = number (integer or float) that serves as the TAP percentile threshold, but only if the expression level is greater than -1. Default = -1
        -thres_endogenous_freq = the limit of the number times the peptide is found as an endogenous peptide -> I will usually use 0 for this to find peptides that are unique. Default = -1
    """
    hash_filters = {
        "frame_stat": frame_stat,
        "neoep_stat": neoep_stat,
        "bool_NMD": bool_NMD,
        "thres_percent_exp": thres_gene_exp,
        "thres_percent_prot": thres_prot,
        "thres_percent_tap": thres_tap,
        "thres_endogenous_freq": thres_endogenous_freq
    }

    return hash_filters

def count_alterations_sample( ):
    """
    quantify the number of alterations per sample
    """

    pass

print "------------ Algorithm: 180203_NeoepStatistics_Velip_V1.py ------------"
"""
Algorithm: Similar to algorithm 180125_NeoepStatistics_PD1_V2.py
Associated Notes:
    -NOTE: same as 171207_NeoepStatistics_PD1_V1.py but uses new function quantify_all_alterations_v2() instead of quantify_all_alterations() from module "171207_SuiteNeoepAlgs.py"
    -NOTE: this algorithm is exactly the same as 171010_NeoepStatistics_V1.py, just a difference in columns so this computes the output from 171031_NeoantigenMHC_Eval_V5.py whereas 171010_NeoepStatistics_V1.py computes output from 171010_NeoantigenMHC_Eval_V4.py.
        -column difference: "aff_cat_x_1" in 171010 version, "aff_cat_x_1" in this version
        -column difference: "aff_cat_x_2" in 171010 version, "aff_cat_x_12" in this version
Protocol:
    -open file of interest
    -
"""

# suite_neoeps = importlib.import_module( "171207_SuiteNeoepAlgs" )
# suite_neoeps = importlib.import_module( "180128_SuiteNeoepAlgs_V2" )
# suite_neoeps = importlib.import_module( "180130_SuiteNeoepClass" )


df = pd.read_csv( DIR_RESULTS_VELIP + "/180203_Thres0_Velip_PC_NeoepCompare_V1.txt", sep = '\t' )
# df = pd.read_csv( DIR_RESULTS_VELIP + "/180203_Thres0_Velip_DM_NeoepCompare_V1.txt", sep = '\t' )




print "------- Quantify Number of Alterations producing High-affinity for each MHC -------"

print "------- Quantify Neoep 2a: Only Gene Expression + NMD-irresistant -------"
# create_hash_filter( frame_stat, neoep_stat, bool_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
hash_filter = create_hash_filter( 0, 0, True, 50, -1, -1, -1 )
obj_filters = NeoepStatistics( df, hash_filter )

df_alt_count = obj_filters.retrieve_mutation_count_all_samples()
df_neoep_count = obj_filters.quant_neoeps_samples()


##TEST: print 
print "df_alt_count:\n", df_alt_count
print "df_neoep_count:\n", df_neoep_count


# suite_neoeps.quantify_all_alterations_v2( df_neoep, file_summary, hash_filter )

print "------------ Algorithm Completed: 180203_NeoepStatistics_Velip_V1.py ------------"