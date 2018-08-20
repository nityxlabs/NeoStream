#/usr/bin/python

import sys
import time

import importlib

import scipy
from scipy import stats         #can use percentileofscore, rankdata
from scipy.stats import ttest_ind       #this performs t-tests

import pandas as pd
import numpy as np

#Graphing library
import matplotlib
# matplotlib.use("Agg")     #comment out this line if I want the plot to show, else use this line if I want to save to a file
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns       #new matplotlib library for aesthetically-pleasing graphs

#PATH ON: Atlas, Pleione, Morpheus
sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv7 import NeoepStatistics
from mokhaPy import mokhaPy


#Constants - directories
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/180203_GenomAltToNeoepV3"
DIR_DATA = DIR_CURR + "/Data"
DIR_RESULTS = DIR_CURR + "/Results"
# DIR_RESULTS_FOLDER = DIR_RESULTS + "/171204_NeoepProcess_V15"
# DIR_RESULTS_FOLDER = DIR_RESULTS + "/180203_Velip_V1"
# DIR_RESULTS_FOLDER = DIR_RESULTS + "/180403_Velip_V2"
DIR_RESULTS_FOLDER = DIR_RESULTS + "/180531_Velip_V3"



def combine_patient_count( series_sample_neoep, version = 1 ):
    """
    Args:
        -series_sample_neoep = Pandas Series where key = ['case_id', 'irRECIST'], & value is count of specific value (e.g. # of neoepitopes, # of mutations, # of nonsilent SNVs, etc.) Examples of commands to retrieve this series:
            -series_sample_neoep = df_filtered.groupby(['case_id', 'irRECIST'])['peptide_2'].count()        #counts all neoepitopes, even if it is double-counted
            -series_sample_neoep = df_filtered.drop_duplicates('peptide_2').groupby(['case_id', 'irRECIST'])['peptide_2'].count()      #number of unique neoepitopes
        -version = integer that should be the type of information that should be recorded in each element
            -1 = corresponds neoepitope count for each patient ID, irRECIST (treatment response)
            -2 = corresponds neoepitope count for each patient ID, irRECIST, HLA allele group with strongest affinity
    Returns:
        -an array of arrays, where each inner array has the following information: [0] = Patient ID, [1] = irRECIST (basically response type), [2] = the value of interest (e.g. # of neoepitopes for sample, # of nonsilent SNVs) -> I actually don't use the Patient ID as of now.
    """
    list_condition_neoep_count = []     #array of arrays, where each array [0] = patient ID, [1] = the response type (e.g. Partial Response, Complete Response, Progressive Disease), & [2] = the # of neoepitopes found in each sample
    for i_sample, v in series_sample_neoep.iteritems():        #i = sample ID, v = count of neoeps in sample ID
        print "Series: i = ", i_sample, " & v = ", v, " type = ", type( v )
        
        if version == 1:        #usually used with file 180130_Thres0_SampAllExpress_NeoepCompare_irRECIST_V1.txt
            patient_id = i_sample
            list_condition_neoep_count.append( [ patient_id, v ] )
        elif version == 2:      #for file 180225_Thres0_SampAllExpress_NeoepCompare_irRECIST_V2.txt
            patient_id = i_sample[0]
            hla_subgroup = i_sample[1]
            list_condition_neoep_count.append( [ patient_id, hla_subgroup, v ] )

    return list_condition_neoep_count


def combine_patient_count_addlabel( series_sample_neoep, list_values, label, version = 1):
    """
    This is the same as def "combine_irRECIST_patient_count()", but also adds label for each element & returns the array that records all results (as to potentially add more elements later) - this is intended to combine all information needed into an array of arrays 
    Args:
        -series_sample_neoep = Pandas Series where key = ['case_id', 'irRECIST'], & value is count of specific value (e.g. # of neoepitopes, # of mutations, # of nonsilent SNVs, etc.) Examples of commands to retrieve this series:
            -series_sample_neoep = df_filtered.groupby(['case_id', 'irRECIST'])['peptide_2'].count()        #counts all neoepitopes, even if it is double-counted
            -series_sample_neoep = df_filtered.drop_duplicates('peptide_2').groupby(['case_id', 'irRECIST'])['peptide_2'].count()      #number of unique neoepitopes
        -list_values = an array that will record the neoepitope with specific conditions associated with the count. This will replace "list_condition_neoep_count" variable used in def "combine_irRECIST_patient_count()"
            -array of arrays, where each array [0] = patient ID, [1] = the response type (e.g. Partial Response, Complete Response, Progressive Disease), & [2] = the # of neoepitopes found in each sample
        -label = added label to keep track of specific conditions, does not have to be a string, can be string, integer, float, etc.
        -version = integer that should be the type of information that should be recorded in each element
            -1 = corresponds neoepitope count for each patient ID, irRECIST (treatment response)
            -2 = corresponds neoepitope count for each patient ID, irRECIST, HLA allele group with strongest affinity
    Returns:
        -an array of arrays, where each inner array has the following information: [0] = Patient ID, [1] = irRECIST (basically response type), [2] = the value of interest (e.g. # of neoepitopes for sample, # of nonsilent SNVs) -> I actually don't use the Patient ID as of now.
    """
    hash_conditions_TEST = {}        #k = irRECIST recovery status, v = array of counts for each sample
    list_condition_neoep_count = []     #array of arrays, where each array [0] = patient ID, [1] = the response type (e.g. Partial Response, Complete Response, Progressive Disease), & [2] = the # of neoepitopes found in each sample
    for i_sample, v in series_sample_neoep.iteritems():        #i = sample ID, v = count of neoeps in sample ID
        print "Series: i = ", i_sample, " & v = ", v, " type = ", type( v )

        
        if version == 1:        #usually used with file 180130_Thres0_SampAllExpress_NeoepCompare_irRECIST_V1.txt
            patient_id = i_sample
            list_values.append( [ patient_id, label, v ] )
        elif version == 2:      #for file 180225_Thres0_SampAllExpress_NeoepCompare_irRECIST_V2.txt
            patient_id = i_sample[0]
            hla_subgroup = i_sample[1]
            list_values.append( [ patient_id, irRECIST_type, hla_subgroup, label, v ] )

    #     ##TEST::
    #     if not irRECIST_type in hash_conditions_TEST:
    #         hash_conditions_TEST[irRECIST_type] = []
    #     hash_conditions_TEST[irRECIST_type].append( v )

    # ##TEST::
    # print "To see read counts for each irRECIST - hash_conditions_TEST:\n", hash_conditions_TEST

    return list_values

def combine_patient_count_addlabel_v2( series_sample_neoep, hash_sample_label, list_values, label, version = 1):
    """
    This is the same as def "combine_patient_count_addlabel()", but also adds label for each sample subgroup (using hash_sample_label) - this is intended to separate each trial for each treatment condition subgroup
    Args:
        -series_sample_neoep = Pandas Series where key = ['case_id', 'irRECIST'], & value is count of specific value (e.g. # of neoepitopes, # of mutations, # of nonsilent SNVs, etc.) Examples of commands to retrieve this series:
            -series_sample_neoep = df_filtered.groupby(['case_id', 'irRECIST'])['peptide_2'].count()        #counts all neoepitopes, even if it is double-counted
            -series_sample_neoep = df_filtered.drop_duplicates('peptide_2').groupby(['case_id', 'irRECIST'])['peptide_2'].count()      #number of unique neoepitopes
        -hash_sample_label = hash where the key = sample ID in "series_sample_neoep" and the value = the associated group label
        -list_values = an array that will record the neoepitope with specific conditions associated with the count. This will replace "list_condition_neoep_count" variable used in def "combine_irRECIST_patient_count()"
            -array of arrays, where each array [0] = patient ID, [1] = the response type (e.g. Partial Response, Complete Response, Progressive Disease), & [2] = the # of neoepitopes found in each sample
        -label = added label to keep track of specific conditions, does not have to be a string, can be string, integer, float, etc.
        -version = integer that should be the type of information that should be recorded in each element
            -1 = corresponds neoepitope count for each patient ID, irRECIST (treatment response)
            -2 = corresponds neoepitope count for each patient ID, irRECIST, HLA allele group with strongest affinity
    Returns:
        -an array of arrays, where each inner array has the following information: [0] = Patient ID, [1] = irRECIST (basically response type), [2] = the value of interest (e.g. # of neoepitopes for sample, # of nonsilent SNVs) -> I actually don't use the Patient ID as of now.
    """
    hash_conditions_TEST = {}        #k = irRECIST recovery status, v = array of counts for each sample
    list_condition_neoep_count = []     #array of arrays, where each array [0] = patient ID, [1] = the response type (e.g. Partial Response, Complete Response, Progressive Disease), & [2] = the # of neoepitopes found in each sample
    for i_sample, v in series_sample_neoep.iteritems():        #i = sample ID, v = count of neoeps in sample ID
        print "Series: i = ", i_sample, " & v = ", v, " type = ", type( v )

        
        if version == 1:        #usually used with file 180130_Thres0_SampAllExpress_NeoepCompare_irRECIST_V1.txt
            patient_id = i_sample
            list_values.append( [ patient_id, hash_sample_label[patient_id], label, v ] )
        elif version == 2:      #for file 180225_Thres0_SampAllExpress_NeoepCompare_irRECIST_V2.txt
            patient_id = i_sample[0]
            hla_subgroup = i_sample[1]
            list_values.append( [ patient_id, hash_sample_label[patient_id], irRECIST_type, hla_subgroup, label, v ] )

    #     ##TEST::
    #     if not irRECIST_type in hash_conditions_TEST:
    #         hash_conditions_TEST[irRECIST_type] = []
    #     hash_conditions_TEST[irRECIST_type].append( v )

    # ##TEST::
    # print "To see read counts for each irRECIST - hash_conditions_TEST:\n", hash_conditions_TEST

    return list_values

def combine_patient_count_addlabel_v3( series_sample_neoep, list_values, label, version = 1):
    """
    NOTE: This function is similar to def "combine_patient_count_addlabel()", but considers 
    This is the same as def "combine_irRECIST_patient_count()", but also adds label for each element & returns the array that records all results (as to potentially add more elements later) - this is intended to combine all information needed into an array of arrays 
    Args:
        -series_sample_neoep = Pandas Series where key = ['case_id', 'irRECIST'], & value is count of specific value (e.g. # of neoepitopes, # of mutations, # of nonsilent SNVs, etc.) Examples of commands to retrieve this series:
            -series_sample_neoep = df_filtered.groupby(['case_id', 'irRECIST'])['peptide_2'].count()        #counts all neoepitopes, even if it is double-counted
            -series_sample_neoep = df_filtered.drop_duplicates('peptide_2').groupby(['case_id', 'irRECIST'])['peptide_2'].count()      #number of unique neoepitopes
        -list_values = an array that will record the neoepitope with specific conditions associated with the count. This will replace "list_condition_neoep_count" variable used in def "combine_irRECIST_patient_count()"
            -array of arrays, where each array [0] = patient ID, [1] = the response type (e.g. Partial Response, Complete Response, Progressive Disease), & [2] = the # of neoepitopes found in each sample
        -label = added label to keep track of specific conditions, does not have to be a string, can be string, integer, float, etc.
        -version = integer that should be the type of information that should be recorded in each element
            -1 = corresponds neoepitope count for each patient ID, irRECIST (treatment response)
            -2 = corresponds neoepitope count for each patient ID, irRECIST, HLA allele group with strongest affinity
    Returns:
        -an array of arrays, where each inner array has the following information: [0] = Patient ID, [1] = irRECIST (basically response type), [2] = the value of interest (e.g. # of neoepitopes for sample, # of nonsilent SNVs) -> I actually don't use the Patient ID as of now.
    """
    hash_conditions_TEST = {}        #k = irRECIST recovery status, v = array of counts for each sample
    list_condition_neoep_count = []     #array of arrays, where each array [0] = patient ID, [1] = the response type (e.g. Partial Response, Complete Response, Progressive Disease), & [2] = the # of neoepitopes found in each sample
    for i_sample, v in series_sample_neoep.iteritems():        #i = sample ID, v = count of neoeps in sample ID
        print "Series: i = ", i_sample, " & v = ", v, " type = ", type( v )

        
        if version == 1:        #usually used with file 180130_Thres0_SampAllExpress_NeoepCompare_irRECIST_V1.txt
            patient_id = i_sample[0]
            patient_category = i_sample[1]
            list_values.append( [ patient_id, patient_category, label, v ] )
        elif version == 2:      #for file 180225_Thres0_SampAllExpress_NeoepCompare_irRECIST_V2.txt
            patient_id = i_sample[0]
            patient_category = i_sample[1]
            hla_subgroup = i_sample[2]
            list_values.append( [ patient_id, patient_category, hla_subgroup, label, v ] )

    #     ##TEST::
    #     if not irRECIST_type in hash_conditions_TEST:
    #         hash_conditions_TEST[irRECIST_type] = []
    #     hash_conditions_TEST[irRECIST_type].append( v )

    # ##TEST::
    # print "To see read counts for each irRECIST - hash_conditions_TEST:\n", hash_conditions_TEST

    return list_values

def combine_patient_count_addlabel_v4( df, col_oi, list_values, str_label, version ):
    """
    This is similar to function combine_patient_count_addlabel_v3(), but the input is a pandas dataframe (most likely filtered) and will be used to quantify
    Args:
        -df = pandas dataframe, it could filtered form of the dataframe
        -col_oi = string that is the column of interest (e.g. "ic_score_2", "mhc_score_2")
        -list_values = an array that will record the neoepitope with specific conditions associated with the count. This will replace "list_condition_neoep_count" variable used in def "combine_irRECIST_patient_count()"
        -array of arrays, where each array [0] = patient ID, [1] = the response type (e.g. Partial Response, Complete Response, Progressive Disease), & [2] = the # of neoepitopes found in each sample
        -label = added label to keep track of specific conditions, does not have to be a string, can be string, integer, float, etc.
        -version = integer that should be the type of information that should be recorded in each element
            -1 = corresponds neoepitope count for each patient ID, irRECIST (treatment response)
            -2 = corresponds neoepitope count for each patient ID, irRECIST, HLA allele group with strongest affinity
    Returns:
        -an array of arrays, where each inner array has the following information: [0] = Patient ID, [1] = irRECIST (basically response type), [2] = the value of interest (e.g. # of neoepitopes for sample, # of nonsilent SNVs) -> I actually don't use the Patient ID as of now.
    """
    hash_conditions_TEST = {}        #k = irRECIST recovery status, v = array of counts for each sample
    list_condition_neoep_count = []     #array of arrays, where each array [0] = patient ID, [1] = the response type (e.g. Partial Response, Complete Response, Progressive Disease), & [2] = the # of neoepitopes found in each sample

    for i_row, row in df.iterrows():
        if version == 1:
            list_values.append( [row['case_id'], row['case_category'], str_label, row[col_oi]] )
        elif version == 2:
            list_values.append( [row['case_id'], row['case_category'], row['preferred_HLA_subgroup'], str_label, row[col_oi] ] )

    return list_values


def calc_ttest_btwn_conds( df_condition_neoep_mut, cond_1_label, cond_2_label, col_oi ):
    """
    Performs a t-test between conditions "cond_1_label" & "cond_2_label"
    Args:
        -df_condition_neoep_mut = pandas DataFrame that contains information about neoepitopes across different samples (e.g. # of neoepitopes in a condition like pcVelip or BRDM,)
    Returns:
        outputs a string of the t-test values
    """
    all_processing_bins = df_condition_neoep_mut['bin_window'].unique()     #retrieve the binning categories, in this case each bin refers to a processing efficiency (e.g. 0-10 is 0 to 10th percentile of neoepitope processing efficiency)
    hash_bin_ttest = {}      #hash where key = processing bin value, value = t-statistic & p-value
    for each_bin in all_processing_bins:
        bin_cond_1 = df_condition_neoep_mut[ (df_condition_neoep_mut["case_category"].str.contains(cond_1_label, case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ][col_oi]
        bin_cond_2 = df_condition_neoep_mut[ (df_condition_neoep_mut["case_category"].str.contains(cond_2_label, case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ][col_oi]

        # ##TEST::
        # print each_bin, " - T-TEST: bin_cond_1:\n", bin_cond_1
        # print each_bin, " - T-TEST: bin_cond_2:\n", bin_cond_2

        t_stat, p_val = ttest_ind( bin_cond_1.values, bin_cond_2.values )
        hash_bin_ttest[each_bin] = (t_stat, p_val)

        # ##TEST:: make sure I'm extracting the correct groups
        # bin_responder = df_condition_neoep_mut[ (~df_condition_neoep_mut["case_category"].str.contains("Non-Respon", case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]
        # bin_nonresponder = df_condition_neoep_mut[ (df_condition_neoep_mut["case_category"].str.contains("Non-Respon", case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]
        # print each_bin, " - bin_responder =\n", bin_responder
        # print each_bin, " - bin_nonresponder =\n", bin_nonresponder
        # print "---------\n"


    # #show the t-test values from all bins
    # print "T-TEST STATISTICS!!! hash_bin_ttest - compare " + cond_1_label + " to " + cond_2_label + ":\n", hash_bin_ttest
    # # for k,v in hash_bin_ttest.iteritems():
    # #     print "bin = ", k, " & t-stat = ", v[0], " & p-value = ", v[1]

    print "T-TEST STATISTICS!!! hash_bin_ttest - compare " + cond_1_label + " to " + cond_2_label + ":"
    all_processing_bins = sorted( df_condition_neoep_mut['bin_window'].unique() )
    for each_bin in all_processing_bins:
        v = hash_bin_ttest[each_bin]
        print "bin = ", each_bin, " & t-stat = ", v[0], " & p-value = ", v[1]

print "------------ Algorithm: 180404_Velip_Graph_V1.py ------------"


# ##Experiment 1 - look at specific conditions of Anti-PD1 - ["Progressive Disease", "Partial Response", "Complete Response"] and the # of neoepitopes associated with them. Very similar to Experiment 15B, but I look at low- & high-mutational load & also look at following comparisons: all mutations vs. high HLA affinity mutations, all mutations vs. high Proteasome score mutations, all mutations vs. high TAP score mutations for all-, low-, & high-mutational load 
# #ACTUAL FILE
# #Exome I dataset
# df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File1_NeoepCompare_V1.txt", sep = '\t' )
# #Exome II dataset
# df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File2_NeoepCompare_V1.txt", sep = '\t' )
# #Exome III dataset - pooled
# df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File3_NeoepCompare_V1.txt", sep = '\t' )
# #WGS I dataset
# df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File4_NeoepCompare_V1.txt", sep = '\t' )
# #WGS II dataset
# df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File5_NeoepCompare_V1.txt", sep = '\t' )

# #need to append "_WGS" (or something) to each row in "case_id" of WGS file or else it will be confused with Exome I & II datasets
# df_4["case_id"] = df_4["case_id"].astype(str) + "_WGS"
# df_5["case_id"] = df_5["case_id"].astype(str) + "_WGS"

# df = pd.concat( [df_1, df_2, df_3, df_4, df_5] )

# print "------------ Experiment 1 ------------"

# print "------------ Create filters for neoepitopes ------------"
# # hash_filters = create_hash_filter( frame_stat, neoep_stat, bool_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# #DIFFERENCE: "mhc_affinity" is different between experiment 12C & 12D
# mhc_affinity = 0        #looks at neoepitope affinities to HLA alleles (e.g. LOW, MID, HIGH), 2 = select high & mid affinity neoepitopes
# frame_stat = 0
# neoep_stat = 0      #looks at all neoeps, high/mid-affinity neoeps, or high-efficacy neoeps
# check_NMD = True
# thres_gene_exp = 0      #as of now, I do not have gene expression for Velip Data
# thres_prot = 0
# thres_tap = 0
# thres_endogenous_freq = -1      #as of now, I do not have endogenous frequency for Velip Data
# hash_filters = NeoepStatistics.create_hash_filter( mhc_affinity, frame_stat, neoep_stat, check_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# obj_neoep_stat = NeoepStatistics( df, hash_filters )

# print "------------ Apply filters for neoepitopes ------------"
# #applied the filters to the dataset
# df_filtered = obj_neoep_stat.filter_neoeps()
# # filter_label = obj_neoep_stat.filter_correspond_labels( hash_filters )

# series_sample_neoep_mut_total = df_filtered.drop_duplicates('genome_pos').groupby( ['case_id'] )['genome_pos'].count()      #counts all neoepitopes, even if it is double-counted
# # series_sample_neoep_mut_hiaff = df_hiaff.drop_duplicates('genome_pos').groupby( ['case_id', 'irRECIST' ] )['genome_pos'].count()      #counts all neoepitopes, even if it is double-counted

# #All Mutations: when mhc_affinity = 0, check_NMD = False
# graph_title = "All Mutations Per Treatment Cohort"
# #High-Affinity Mutations: when mhc_affinity = 2, check_NMD = True
# # graph_title = "High-Affinity Mutations Per Treatment Cohort"
# #High-Efficacy Mutations: when mhc_affinity = 2, check_NMD = True, thres_prot = 90, thres_tap = 90
# # graph_title = "High-Efficacy Mutations Per Treatment Cohort"



# # #count the number of mutations per sample
# # version = 1
# # list_condition_neoep_count = []
# # list_condition_neoep_count = combine_patient_count( series_sample_neoep_mut_total, version )

# # #NOTE: pcMV (this is pcVelip in pooled) & BMV (this is BRVelip in pooled) are from the pooled samples
# # items_Velip = ['pcVelip1', 'pcVelip2', 'pcVelip3', 'pcVelip1_WGS', 'BRVelip1', 'BRVelip2', 'BRVelip3', 'BRVelip1_WGS', 'pcMV', 'BMV']
# # items_DMSO = ['pcDM1', 'pcDM2', 'pcDM3', 'pcDM1_WGS', 'BRDM1', 'BRDM2', 'BRDM3', 'BRDM1_WGS', 'pcMD', 'BMD']

# # #Attempt 1: graph Velip vs. DMSO
# # version = 1
# # list_condition_neoep_count = []
# # #Get Velip conditions
# # series_Velip = series_sample_neoep_mut_total.filter( items = items_Velip, axis = 0)
# # list_condition_neoep_count = combine_patient_count_addlabel( series_Velip, list_condition_neoep_count, "Velip", version )
# # #Get DMSO
# # series_DMSO = series_sample_neoep_mut_total.filter( items = items_DMSO, axis = 0)
# # list_condition_neoep_count = combine_patient_count_addlabel( series_DMSO, list_condition_neoep_count, "DMSO", version )

# # df_condition_neoep_count = pd.DataFrame( list_condition_neoep_count, columns = ['sample_id', 'exp_cond', 'neoep_mut_count'] )
# # #need to sort labels by therapy responses so I can see a trend - ["Progressive Disease", "Partial Response", "Complete Response"]
# # df_condition_neoep_count.sort_values( by = ['exp_cond'], ascending = True, inplace = True )


# #Attempt 2: graph pc vs. br & Velip vs. DMSO
# items_pc_Velip = ['pcVelip1', 'pcVelip2', 'pcVelip3', 'pcVelip1_WGS', 'pcMV']
# items_br_Velip = ['BRVelip1', 'BRVelip2', 'BRVelip3', 'BRVelip1_WGS', 'BMV']
# items_pc_DMSO = ['pcDM1', 'pcDM2', 'pcDM3', 'pcDM1_WGS', 'pcMD']
# items_br_DMSO = ['BRDM1', 'BRDM2', 'BRDM3', 'BRDM1_WGS', 'BMD']
# list_hue_order = items_pc_Velip + items_br_Velip + items_pc_DMSO + items_br_DMSO        #for sns.barplot

# version = 1
# list_condition_neoep_count = []
# #Get pc + Velip conditions
# series_pc_Velip = series_sample_neoep_mut_total.filter( items = items_pc_Velip, axis = 0)
# list_condition_neoep_count = combine_patient_count_addlabel( series_pc_Velip, list_condition_neoep_count, "pcVelip", version )
# #Get br + Velip conditions
# series_br_Velip = series_sample_neoep_mut_total.filter( items = items_br_Velip, axis = 0)
# list_condition_neoep_count = combine_patient_count_addlabel( series_br_Velip, list_condition_neoep_count, "brVelip", version )
# #Get pc + DMSO conditions
# series_pc_DMSO = series_sample_neoep_mut_total.filter( items = items_pc_DMSO, axis = 0)
# list_condition_neoep_count = combine_patient_count_addlabel( series_pc_DMSO, list_condition_neoep_count, "pcDMSO", version )
# #Get br + DMSO conditions
# series_br_DMSO = series_sample_neoep_mut_total.filter( items = items_br_DMSO, axis = 0)
# list_condition_neoep_count = combine_patient_count_addlabel( series_br_DMSO, list_condition_neoep_count, "brDMSO", version )

# df_condition_neoep_count = pd.DataFrame( list_condition_neoep_count, columns = ['sample_id', 'exp_cond', 'neoep_mut_count'] )
# #need to sort labels by therapy responses so I can see a trend between the treatment conditions
# df_condition_neoep_count.sort_values( by = ['exp_cond'], ascending = True, inplace = True )


# # list_condition_neoep_count = combine_patient_count_addlabel( series_sample_neoep_mut_hiaff, list_condition_neoep_count, label_filter, version )


# df_condition_neoep_count = pd.DataFrame( list_condition_neoep_count, columns = ['sample_id', 'grouped_trials', 'exp_cond', 'neoep_mut_count'] )
# #need to sort labels by therapy responses so I can see a trend - ["Progressive Disease", "Partial Response", "Complete Response"]
# df_condition_neoep_count.sort_values( by = ['exp_cond'], ascending = True, inplace = True )


# ##TEST::
# print "df_condition_neoep_count:\n", df_condition_neoep_count


# #graph the results
# # sns.set( font_scale = 1.2 )
# # sns.set()
# # sns.set_style( "darkgrid" )
# # sns.set_style( "whitegrid" )
# # sns.set_style( "dark" )
# sns.set_style( "white" )

# # sns.set_context( "paper" )        #smallest features ( thinnest lines, smallest font for ticks and labels )
# # sns.set_context( "notebook" )     #3rd largest features
# sns.set_context( "talk" )     #2nd largest features
# # sns.set_context( "poster" )     #largest features ( thickest lines, largest font for ticks and labels )

# #graph Seaborn pyplot try to create a box-whisker plot for each Proteasome
# # ax = sns.boxplot( x = 'response_type', y = 'neoep_mut_count', hue = 'mut_type', data = df_condition_neoep_count, hue_order = ['all', label_filter], showfliers = False, palette = "Blues" )
# # ax = sns.swarmplot( x = 'exp_cond', y = 'neoep_mut_count', hue = 'sample_id', data = df_condition_neoep_count, hue_order = ['all', label_filter], palette = "Blues", color = 0.25 )

# # ax = sns.factorplot( x = 'exp_cond', y = 'neoep_mut_count', hue = 'sample_id', data = df_condition_neoep_count, kind = "bar", palette = "Blues" )
# ax = sns.barplot( x = 'exp_cond', y = 'neoep_mut_count', hue = 'sample_id', data = df_condition_neoep_count, hue_order = list_hue_order, palette = "Blues" )
# # sns.boxplot( x = 'exp_cond', y = 'neoep_mut_count', hue = 'sample_id', data = df_condition_neoep_count, palette = "Blues" )
# # sns.swarmplot( x = 'exp_cond', y = 'neoep_mut_count', hue = 'sample_id', data = df_condition_neoep_count, palette = "Reds" )
# # sns.boxplot( x = 'response_type', y = 'neoep_mut_count', data = df_condition_neoep_count, palette = "Blues" )
# sns.despine()
# # pyplot.title( "Neoepitope Count Per irRECIST Condition: Unique Neoepitopes" )
# # pyplot.title( "Exp 12: Neoepitope Count for HLA Subgroup Per irRECIST Condition: " + annex_title, fontsize = 10 )
# # graph_title = "Expressed Mutations Per Treatment Cohort"    #when thres_gene_exp > 0 (usu. thres_gene_exp = 1)
# # graph_title = "High-Affinity Mutations Per Treatment Cohort"    #when thres_gene_exp > 0 (usu. thres_gene_exp = 1), mhc_affinity = 2 (looks for high/mid-affinity neoepitopes), & thres_endogenous_freq = 0 (neoepitope sequence is novel, no matching endogenous AA sequence )
# pyplot.title( graph_title, fontsize = 30 )
# # pyplot.title( graph_title, fontsize = 14 )
# # pyplot.xlabel( "irRECIST" )
# pyplot.xlabel( "" )
# pyplot.ylabel( "Number of Mutations", fontsize = 18 )
# # ax.legend_.remove()       #use this to remove the legend
# # ax.legend().set_title( "" )     #This is how to change the legend title
# # pyplot.ylim( (0, 1600) )

# pyplot.show()









# ##Experiment 2 - Look at the different clone lines between
# print "------------ Experiment 1 ------------"
# #ACTUAL FILE
# #Exome I dataset
# df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File1_NeoepCompare_V1.txt", sep = '\t' )
# #Exome II dataset
# df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File2_NeoepCompare_V1.txt", sep = '\t' )
# #Exome III dataset - pooled
# df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File3_NeoepCompare_V1.txt", sep = '\t' )
# #WGS I dataset
# df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File4_NeoepCompare_V1.txt", sep = '\t' )
# #WGS II dataset
# df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File5_NeoepCompare_V1.txt", sep = '\t' )

# #need to append "_WGS" (or something) to each row in "case_id" of WGS file or else it will be confused with Exome I & II datasets
# df_4["case_id"] = df_4["case_id"].astype(str) + "_WGS"
# df_5["case_id"] = df_5["case_id"].astype(str) + "_WGS"

# df = pd.concat( [df_1, df_2, df_3, df_4, df_5] )



# print "------------ Create filters for neoepitopes ------------"
# # hash_filters = create_hash_filter( frame_stat, neoep_stat, bool_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# #DIFFERENCE: "mhc_affinity" is different between experiment 12C & 12D
# mhc_affinity = 2        #looks at neoepitope affinities to HLA alleles (e.g. LOW, MID, HIGH), 2 = select high & mid affinity neoepitopes
# frame_stat = 0
# neoep_stat = 0      #looks at all neoeps, high/mid-affinity neoeps, or high-efficacy neoeps
# check_NMD = True
# thres_gene_exp = 0      #as of now, I do not have gene expression for Velip Data
# thres_prot = 90
# thres_tap = 90
# thres_endogenous_freq = -1      #as of now, I do not have endogenous frequency for Velip Data
# hash_filters = NeoepStatistics.create_hash_filter( mhc_affinity, frame_stat, neoep_stat, check_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# obj_neoep_stat = NeoepStatistics( df, hash_filters )

# print "------------ Apply filters for neoepitopes ------------"
# #applied the filters to the dataset
# df_filtered = obj_neoep_stat.filter_neoeps()
# # filter_label = obj_neoep_stat.filter_correspond_labels( hash_filters )

# series_sample_neoep_mut_total = df_filtered.drop_duplicates('genome_pos').groupby( ['case_id'] )['genome_pos'].count()      #counts all neoepitopes, even if it is double-counted
# # series_sample_neoep_mut_hiaff = df_hiaff.drop_duplicates('genome_pos').groupby( ['case_id', 'irRECIST' ] )['genome_pos'].count()      #counts all neoepitopes, even if it is double-counted

# # #All Mutations: when mhc_affinity = 0, check_NMD = False
# # graph_title = "All Mutations Per Treatment Cohort"
# # # High-Affinity Mutations: when mhc_affinity = 2, check_NMD = True
# # graph_title = "High-Affinity Mutations Per Treatment Cohort"
# #High-Efficacy Mutations: when mhc_affinity = 2, check_NMD = True, thres_prot = 90, thres_tap = 90
# graph_title = "High-Efficacy Mutations Per Treatment Cohort"

# items_pc_Velip = ['pcVelip1', 'pcVelip2', 'pcVelip3', 'pcVelip1_WGS', 'pcMV']
# items_br_Velip = ['BRVelip1', 'BRVelip2', 'BRVelip3', 'BRVelip1_WGS', 'BMV']
# items_pc_DMSO = ['pcDM1', 'pcDM2', 'pcDM3', 'pcDM1_WGS', 'pcMD']
# items_br_DMSO = ['BRDM1', 'BRDM2', 'BRDM3', 'BRDM1_WGS', 'BMD']
# list_hue_order = items_pc_Velip + items_br_Velip + items_pc_DMSO + items_br_DMSO        #for sns.barplot

# #need to create hash the corresponds clones, WGS, & pooled
# hash_clone_groups = {}      #hash where keys = clone groups labels, values = each sample
# list_clone_groups = ["clone1", "clone2", "clone3", "WGS", "pooled"]
# for c, clone_group in enumerate( list_clone_groups ):
#     hash_clone_groups[clone_group] = [ items_pc_Velip[c], items_br_Velip[c], items_pc_DMSO[c], items_br_DMSO[c] ]
# #if the hash is 
# reverse_hash_clone_groups = {}      #since this is the reverse of "hash_clone_groups", k = each sample name, v = the group it belongs to
# for k,v in hash_clone_groups.iteritems():       #k = labeled clone group, v = array of samples that fit in this group
#     for each_v in v:
#         reverse_hash_clone_groups[each_v] = k

# print "hash_clone_groups: ", hash_clone_groups
# print "reverse_hash_clone_groups: ", reverse_hash_clone_groups

# version = 1
# list_condition_neoep_count = []
# #Get pc + Velip conditions
# series_pc_Velip = series_sample_neoep_mut_total.filter( items = items_pc_Velip, axis = 0)
# list_condition_neoep_count = combine_patient_count_addlabel_v2( series_pc_Velip, reverse_hash_clone_groups, list_condition_neoep_count, "pcVelip", version )
# #Get br + Velip conditions
# series_br_Velip = series_sample_neoep_mut_total.filter( items = items_br_Velip, axis = 0)
# list_condition_neoep_count = combine_patient_count_addlabel_v2( series_br_Velip, reverse_hash_clone_groups, list_condition_neoep_count, "brVelip", version )
# #Get pc + DMSO conditions
# series_pc_DMSO = series_sample_neoep_mut_total.filter( items = items_pc_DMSO, axis = 0)
# list_condition_neoep_count = combine_patient_count_addlabel_v2( series_pc_DMSO, reverse_hash_clone_groups, list_condition_neoep_count, "pcDMSO", version )
# #Get br + DMSO conditions
# series_br_DMSO = series_sample_neoep_mut_total.filter( items = items_br_DMSO, axis = 0)
# list_condition_neoep_count = combine_patient_count_addlabel_v2( series_br_DMSO, reverse_hash_clone_groups, list_condition_neoep_count, "brDMSO", version )

# # list_condition_neoep_count = combine_patient_count_addlabel( series_sample_neoep_mut_hiaff, list_condition_neoep_count, label_filter, version )


# df_condition_neoep_count = pd.DataFrame( list_condition_neoep_count, columns = ['sample_id', 'clone_group', 'exp_cond', 'neoep_mut_count'] )
# #need to sort labels by therapy responses so I can see a trend - ["Progressive Disease", "Partial Response", "Complete Response"]
# df_condition_neoep_count.sort_values( by = ['exp_cond'], ascending = True, inplace = True )


# ##TEST::
# print "df_condition_neoep_count:\n", df_condition_neoep_count


# #graph the results
# # sns.set( font_scale = 1.2 )
# # sns.set()
# # sns.set_style( "darkgrid" )
# # sns.set_style( "whitegrid" )
# # sns.set_style( "dark" )
# sns.set_style( "white" )

# # sns.set_context( "paper" )        #smallest features ( thinnest lines, smallest font for ticks and labels )
# # sns.set_context( "notebook" )     #3rd largest features
# sns.set_context( "talk" )     #2nd largest features
# # sns.set_context( "poster" )     #largest features ( thickest lines, largest font for ticks and labels )

# #graph Seaborn pyplot try to create a box-whisker plot for each Proteasome
# # ax = sns.boxplot( x = 'response_type', y = 'neoep_mut_count', hue = 'mut_type', data = df_condition_neoep_count, hue_order = ['all', label_filter], showfliers = False, palette = "Blues" )
# # ax = sns.swarmplot( x = 'exp_cond', y = 'neoep_mut_count', hue = 'sample_id', data = df_condition_neoep_count, hue_order = ['all', label_filter], palette = "Blues", color = 0.25 )

# # ax = sns.factorplot( x = 'exp_cond', y = 'neoep_mut_count', hue = 'sample_id', data = df_condition_neoep_count, kind = "bar", palette = "Blues" )
# ax = sns.barplot( x = 'exp_cond', y = 'neoep_mut_count', hue = 'clone_group', data = df_condition_neoep_count, hue_order = list_clone_groups, palette = "Blues" )
# # sns.boxplot( x = 'exp_cond', y = 'neoep_mut_count', hue = 'sample_id', data = df_condition_neoep_count, palette = "Blues" )
# # sns.swarmplot( x = 'exp_cond', y = 'neoep_mut_count', hue = 'sample_id', data = df_condition_neoep_count, palette = "Reds" )
# # sns.boxplot( x = 'response_type', y = 'neoep_mut_count', data = df_condition_neoep_count, palette = "Blues" )
# sns.despine()
# # pyplot.title( "Neoepitope Count Per irRECIST Condition: Unique Neoepitopes" )
# # pyplot.title( "Exp 12: Neoepitope Count for HLA Subgroup Per irRECIST Condition: " + annex_title, fontsize = 10 )
# # graph_title = "Expressed Mutations Per Treatment Cohort"    #when thres_gene_exp > 0 (usu. thres_gene_exp = 1)
# # graph_title = "High-Affinity Mutations Per Treatment Cohort"    #when thres_gene_exp > 0 (usu. thres_gene_exp = 1), mhc_affinity = 2 (looks for high/mid-affinity neoepitopes), & thres_endogenous_freq = 0 (neoepitope sequence is novel, no matching endogenous AA sequence )
# pyplot.title( graph_title, fontsize = 30 )
# # pyplot.title( graph_title, fontsize = 14 )
# # pyplot.xlabel( "irRECIST" )
# pyplot.xlabel( "" )
# pyplot.ylabel( "Number of Mutations", fontsize = 18 )
# # ax.legend_.remove()       #use this to remove the legend
# # ax.legend().set_title( "" )     #This is how to change the legend title
# # pyplot.ylim( (0, 1600) )

# pyplot.show()








# ##Experiment 3 - Retrieve the highest MHC affinity for each mutation and then graph the average MHC affinity (in nM) for each condition
# def get_mean_std_experconditions( df, series_std, exp_cond, clone_order ):
#     """
#     Retrieves the corresponding average value & standard deviation for each experimental condition (e.g. pc, br + Velip, DMSO, so examples are "pcVelip", "pcDMSO", "brVelip", "brDMSO")
#     Args:
#         -df_oi = Pandas Dataframe that is filtered for specific experimental condition
#         -series_std = Pandas Series that has the standard deviations for each sample ID 
#         -exp_cond = string that is the experimental condition (i.e. "pcVelip", "pcDMSO", "brVelip", "brDMSO")
#         -clone_order = array that is the order of the clones -> ['clone1', 'clone2', 'clone3', 'pooled', 'WGS']
#     """
#     clone_order = ['clone1', 'clone2', 'clone3', 'pooled', 'WGS']
#     df_oi = df[ df['exp_cond'] == exp_cond ]

#     mean_cond = []
#     std_cond = []
#     # mean_cond = [df_oi[df_oi["clone_group"] == x]["avg_IC50"].values[0] for x in clone_order]
#     for each_c in clone_order:
#         mean_cond.append( df_oi[df_oi["clone_group"] == each_c]["avg_IC50"].values[0] )
#         std_cond.append( series_std[ df_oi[df_oi["clone_group"] == each_c]["sample_id"].values[0] ] )

#     #returns 2 arrays, 1 with the average values & the other of 
#     return [mean_cond, std_cond]

# def get_mean_std_clones( df, series_std, clone_group, exp_cond_order ):
#     """
#     Retrieves the corresponding average value & standard deviation for each experimental condition (e.g. pc, br + Velip, DMSO, so examples are "pcVelip", "pcDMSO", "brVelip", "brDMSO")
#     Args:
#         -df_oi = Pandas Dataframe that is filtered for specific experimental condition
#         -series_std = Pandas Series that has the standard deviations for each sample ID 
#         -clone_group = string that is the clone group (i.e. "clone1", "clone2", "clone3", "pooled", "WGS")
#         -exp_cond_order = array of the order of experimental conditions -> "['brDMSO', 'brVelip', 'pcDMSO', 'pcVelip']"
#     """
#     exp_cond_order = ['brDMSO', 'brVelip', 'pcDMSO', 'pcVelip']
#     df_oi = df[ df['clone_group'] == clone_group ]

#     mean_cond = []
#     std_cond = []
#     # mean_cond = [df_oi[df_oi["clone_group"] == x]["avg_IC50"].values[0] for x in clone_order]
#     for each_c in exp_cond_order:
#         mean_cond.append( df_oi[df_oi["exp_cond"] == each_c]["avg_IC50"].values[0] )
#         std_cond.append( series_std[ df_oi[df_oi["exp_cond"] == each_c]["sample_id"].values[0] ] )

#     #returns 2 arrays, 1 with the average values & the other of 
#     return [mean_cond, std_cond]

# print "------------ Experiment 3 ------------"
# # #ACTUAL FILE - Files that DO NOT have MHC scores ranked for each mutation
# # #Exome I dataset
# # df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File1_NeoepCompare_V1.txt", sep = '\t' )
# # #Exome II dataset
# # df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File2_NeoepCompare_V1.txt", sep = '\t' )
# # #Exome III dataset - pooled
# # df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File3_NeoepCompare_V1.txt", sep = '\t' )
# # #WGS I dataset
# # df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File4_NeoepCompare_V1.txt", sep = '\t' )
# # #WGS II dataset
# # df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File5_NeoepCompare_V1.txt", sep = '\t' )

# #ACTUAL FILE - Files that do have MHC scores ranked for each mutation
# #Exome I dataset
# df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/180411_Thres0_File1_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #Exome II dataset
# df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/180411_Thres0_File2_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #Exome III dataset - pooled
# df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/180411_Thres0_File3_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #WGS I dataset
# df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/180411_Thres0_File4_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #WGS II dataset
# df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/180411_Thres0_File5_NeoepCompare_V1_Extended.txt", sep = '\t' )

# #need to append "_WGS" (or something) to each row in "case_id" of WGS file or else it will be confused with Exome I & II datasets
# df_4["case_id"] = df_4["case_id"].astype(str) + "_WGS"
# df_5["case_id"] = df_5["case_id"].astype(str) + "_WGS"

# df = pd.concat( [df_1, df_2, df_3, df_4, df_5] )

# ##TEST::
# print "TEST:: length before rank 1 filter = ", len( df )
# #filter for the highest MHC score for each mutation
# df = df[ df['mut_neoep_rank_1'] == 1 ]
# print "TEST:: length AFTER rank 1 filter = ", len( df )
# print "TEST - df ALL sample names = ",  np.unique( df['case_id'] )

# #this is just to write the file of all mutations assigned to a specific MHC allele based on highest affinity
# df.to_csv( DIR_RESULTS_FOLDER + "/180411_VelipMutNeoep.txt", sep = '\t' )


# #ERROR WITH THIS - for some reason, after applying this filter I am missing 'case_id' that were present in the original 'df' -> not sure why
# # #retrieve the highest MHC affinity for each mutated peptide (should have the highest MHC score / lowest IC50 value)
# # # df_mhc_max = df.iloc[ df.groupby(['case_id', 'genome_pos'])['mhc_score_2'].agg( pd.Series.idxmax ) ]
# # df_mhc_max = df.iloc[ df.groupby(['case_id', 'genome_pos'])['ic50_score_2'].agg( pd.Series.idxmin ) ]

# # print "df_mhc_max:\n", df_mhc_max

# # df_mhc_max.to_csv( DIR_RESULTS_FOLDER + "/180411_MaxMHCTest.txt", sep = '\t' )

# # ##TEST::
# # print "TEST 2- df_mhc_max sample names = ",  np.unique( df_mhc_max['case_id'] )


# print "------------ Create filters for neoepitopes ------------"
# # hash_filters = create_hash_filter( frame_stat, neoep_stat, bool_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# #DIFFERENCE: "mhc_affinity" is different between experiment 12C & 12D
# mhc_affinity = 2        #looks at neoepitope affinities to HLA alleles (e.g. LOW, MID, HIGH), 2 = select high & mid affinity neoepitopes
# frame_stat = 0
# neoep_stat = 0      #looks at all neoeps, high/mid-affinity neoeps, or high-efficacy neoeps
# check_NMD = True
# thres_gene_exp = 0      #as of now, I do not have gene expression for Velip Data
# thres_prot = 0
# thres_tap = 0
# thres_endogenous_freq = -1      #as of now, I do not have endogenous frequency for Velip Data
# hash_filters = NeoepStatistics.create_hash_filter( mhc_affinity, frame_stat, neoep_stat, check_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# obj_neoep_stat = NeoepStatistics( df, hash_filters )

# print "------------ Apply filters for neoepitopes ------------"
# #applied the filters to the dataset
# df_filtered = obj_neoep_stat.filter_neoeps()
# # filter_label = obj_neoep_stat.filter_correspond_labels( hash_filters )

# if mhc_affinity == 0:
#     graph_title = "Average Binding Affinity of All Neoepitopes"
# elif mhc_affinity == 2:
#     graph_title = "Average Binding Affinity of High-Affinity Neoepitopes"

# # #Count the number of mutations for each case_id
# # series_sample_neoep_mut_total = df_filtered.drop_duplicates('genome_pos').groupby( ['case_id'] )['genome_pos'].count()      #counts all neoepitopes, even if it is double-counted
# # # series_sample_neoep_mut_hiaff = df_hiaff.drop_duplicates('genome_pos').groupby( ['case_id', 'irRECIST' ] )['genome_pos'].count()      #counts all neoepitopes, even if it is double-counted


# series_sample_neoep_mhcAff_mean = df_filtered.drop_duplicates('genome_pos').groupby( ['case_id'] )['ic50_score_2'].mean()      #counts all neoepitopes, even if it is double-counted
# series_sample_neoep_mhcAff_std = df_filtered.drop_duplicates('genome_pos').groupby( ['case_id'] )['ic50_score_2'].std()      #counts all neoepitopes, even if it is double-counted



# ##TEST::
# print "series_sample_neoep_mhcAff_mean - MEAN:\n", series_sample_neoep_mhcAff_mean
# print "series_sample_neoep_mhcAff_std - STD:\n", series_sample_neoep_mhcAff_std
# print "series_sample_neoep_mhcAff_std - STD - BRDM1_WGS:\n", series_sample_neoep_mhcAff_std['BRDM1_WGS']


# items_pc_Velip = ['pcVelip1', 'pcVelip2', 'pcVelip3', 'pcVelip1_WGS', 'pcMV']
# items_br_Velip = ['BRVelip1', 'BRVelip2', 'BRVelip3', 'BRVelip1_WGS', 'BMV']
# items_pc_DMSO = ['pcDM1', 'pcDM2', 'pcDM3', 'pcDM1_WGS', 'pcMD']
# items_br_DMSO = ['BRDM1', 'BRDM2', 'BRDM3', 'BRDM1_WGS', 'BMD']
# list_hue_order = items_pc_Velip + items_br_Velip + items_pc_DMSO + items_br_DMSO        #for sns.barplot

# #need to create hash the corresponds clones, WGS, & pooled
# hash_clone_groups = {}      #hash where keys = clone groups labels, values = each sample
# list_clone_groups = ["clone1", "clone2", "clone3", "WGS", "pooled"]
# for c, clone_group in enumerate( list_clone_groups ):
#     hash_clone_groups[clone_group] = [ items_pc_Velip[c], items_br_Velip[c], items_pc_DMSO[c], items_br_DMSO[c] ]
# #if the hash is 
# reverse_hash_clone_groups = {}      #since this is the reverse of "hash_clone_groups", k = each sample name, v = the group it belongs to
# for k,v in hash_clone_groups.iteritems():       #k = labeled clone group, v = array of samples that fit in this group
#     for each_v in v:
#         reverse_hash_clone_groups[each_v] = k

# print "hash_clone_groups: ", hash_clone_groups
# print "reverse_hash_clone_groups: ", reverse_hash_clone_groups

# version = 1
# list_condition_neoep_count = []
# #Get pc + Velip conditions
# series_pc_Velip = series_sample_neoep_mhcAff_mean.filter( items = items_pc_Velip, axis = 0)
# list_condition_neoep_count = combine_patient_count_addlabel_v2( series_pc_Velip, reverse_hash_clone_groups, list_condition_neoep_count, "pcVelip", version )
# #Get br + Velip conditions
# series_br_Velip = series_sample_neoep_mhcAff_mean.filter( items = items_br_Velip, axis = 0)
# list_condition_neoep_count = combine_patient_count_addlabel_v2( series_br_Velip, reverse_hash_clone_groups, list_condition_neoep_count, "brVelip", version )
# #Get pc + DMSO conditions
# series_pc_DMSO = series_sample_neoep_mhcAff_mean.filter( items = items_pc_DMSO, axis = 0)
# list_condition_neoep_count = combine_patient_count_addlabel_v2( series_pc_DMSO, reverse_hash_clone_groups, list_condition_neoep_count, "pcDMSO", version )
# #Get br + DMSO conditions
# series_br_DMSO = series_sample_neoep_mhcAff_mean.filter( items = items_br_DMSO, axis = 0)
# list_condition_neoep_count = combine_patient_count_addlabel_v2( series_br_DMSO, reverse_hash_clone_groups, list_condition_neoep_count, "brDMSO", version )

# # list_condition_neoep_count = combine_patient_count_addlabel( series_sample_neoep_mut_hiaff, list_condition_neoep_count, label_filter, version )


# df_condition_neoep_count = pd.DataFrame( list_condition_neoep_count, columns = ['sample_id', 'clone_group', 'exp_cond', 'avg_IC50'] )
# #need to sort labels by therapy responses so I can see a trend - ["Progressive Disease", "Partial Response", "Complete Response"]
# df_condition_neoep_count.sort_values( by = ['exp_cond'], ascending = True, inplace = True )


# ##TEST::
# print "df_condition_neoep_count:\n", df_condition_neoep_count


# #Generate the graphs for all experimental conditions & clones
# # bool_errorbars = True
# bool_errorbars = False

# sns.set_style( "white" )

# # sns.set_context( "paper" )        #smallest features ( thinnest lines, smallest font for ticks and labels )
# # sns.set_context( "notebook" )     #3rd largest features
# sns.set_context( "talk" )     #2nd largest features
# # sns.set_context( "poster" )     #largest features ( thickest lines, largest font for ticks and labels )

# #custom color palette
# # current_palette = sns.color_palette( "husl", 4 )        #color-wheel palette, for now skip colors between indices 1 to 3
# current_palette = sns.color_palette( "Blues", 5 )
# # # sns.set_palette( [ current_palette[0] ] )
# sns.set_palette( current_palette )

# width = 0.15
# exp_cond_order = ['brDMSO', 'brVelip', 'pcDMSO', 'pcVelip']

# # clone_order = ['clone1', 'clone2', 'clone3', 'pooled', 'WGS']
# # [mean_brDMSO, std_brDMSO] = get_mean_std_experconditions( df_condition_neoep_count, series_sample_neoep_mhcAff_std, 'brDMSO', clone_order )
# [mean_clone1, std_clone1] = get_mean_std_clones( df_condition_neoep_count, series_sample_neoep_mhcAff_std, 'clone1', exp_cond_order )
# idx = np.arange( len(mean_clone1) )
# # mean_brDMSO = [df_brDMSO[df_brDMSO["clone_group"] == x] for x in clone_order] 
# print "mean_clone1 = ", mean_clone1
# print "std_clone1 = ", std_clone1
# if bool_errorbars:
#     pyplot.bar( idx, mean_clone1, width, yerr = std_clone1, label = 'clone1')
# else:
#     pyplot.bar( idx, mean_clone1, width, label = 'clone1')      #no error bars

# [mean_clone2, std_clone2] = get_mean_std_clones( df_condition_neoep_count, series_sample_neoep_mhcAff_std, 'clone2', exp_cond_order )
# idx = np.arange( len(mean_clone2) )
# # mean_brDMSO = [df_brDMSO[df_brDMSO["clone_group"] == x] for x in clone_order] 
# print "mean_clone2 = ", mean_clone2
# print "std_clone2 = ", std_clone2
# if bool_errorbars:
#     pyplot.bar( idx + width, mean_clone2, width, yerr = std_clone2, label = 'clone2')
# else:
#     pyplot.bar( idx + width, mean_clone2, width, label = 'clone2')      #no error bars

# [mean_clone3, std_clone3] = get_mean_std_clones( df_condition_neoep_count, series_sample_neoep_mhcAff_std, 'clone3', exp_cond_order )
# idx = np.arange( len(mean_clone3) )
# # mean_brDMSO = [df_brDMSO[df_brDMSO["clone_group"] == x] for x in clone_order] 
# print "mean_clone3 = ", mean_clone3
# print "std_clone3 = ", std_clone3
# if bool_errorbars:
#     pyplot.bar( idx + 2*width, mean_clone3, width, yerr = std_clone3, label = 'clone3')
# else:
#     pyplot.bar( idx + 2*width, mean_clone3, width, label = 'clone3')              #no error bars

# [mean_WGS, std_WGS] = get_mean_std_clones( df_condition_neoep_count, series_sample_neoep_mhcAff_std, 'WGS', exp_cond_order )
# idx = np.arange( len(mean_WGS) )
# # mean_brDMSO = [df_brDMSO[df_brDMSO["clone_group"] == x] for x in clone_order] 
# print "mean_WGS = ", mean_WGS
# print "std_WGS = ", std_WGS
# if bool_errorbars:
#     pyplot.bar( idx + 3*width, mean_WGS, width, yerr = std_WGS, label = 'WGS')
# else:
#     pyplot.bar( idx + 3*width, mean_WGS, width, label = 'WGS')      #no error bars

# [mean_pooled, std_pooled] = get_mean_std_clones( df_condition_neoep_count, series_sample_neoep_mhcAff_std, 'pooled', exp_cond_order )
# idx = np.arange( len(mean_pooled) )
# # mean_brDMSO = [df_brDMSO[df_brDMSO["clone_group"] == x] for x in clone_order] 
# print "mean_pooled = ", mean_pooled
# print "std_pooled = ", std_pooled
# if bool_errorbars:
#     pyplot.bar( idx + 4*width, mean_pooled, width, yerr = std_pooled, label = 'pooled')
# else:
#     pyplot.bar( idx + 4*width, mean_pooled, width, label = 'pooled')        #no error bars

# sns.despine()

# pyplot.title( graph_title, fontsize = 30 )
# # pyplot.title( graph_title, fontsize = 14 )
# # pyplot.xlabel( "irRECIST" )
# pyplot.xlabel( "" )
# pyplot.ylabel( "MHC affinity (nM)", fontsize = 18 )

# pyplot.xticks(idx + 2 * width, exp_cond_order, rotation = 0 )
# pyplot.legend()

# pyplot.show()













# ##Experiment 4 - Calculate the fold change in MHC affinity before & after alteration
# print "------------ Experiment 4 ------------"
# # input_file_date = "180411"      #for directory "/180403_Velip_V2"
# # input_file_date = "180605"      #for directory "/180531_Velip_V3"
# input_file_date = "180614"      #for directory "/180531_Velip_V3"


# # #ACTUAL FILE - Files that DO NOT have MHC scores ranked for each mutation
# # #Exome I dataset
# # df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File1_NeoepCompare_V1.txt", sep = '\t' )
# # #Exome II dataset
# # df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File2_NeoepCompare_V1.txt", sep = '\t' )
# # #Exome III dataset - pooled
# # df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File3_NeoepCompare_V1.txt", sep = '\t' )
# # #WGS I dataset
# # df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File4_NeoepCompare_V1.txt", sep = '\t' )
# # #WGS II dataset
# # df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File5_NeoepCompare_V1.txt", sep = '\t' )

# #ACTUAL FILE - Files that do have MHC scores ranked for each mutation
# #Exome I dataset
# df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File1_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #Exome II dataset
# df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File2_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #Exome III dataset - pooled
# df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File3_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #WGS I dataset
# df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File4_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #WGS II dataset
# df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File5_NeoepCompare_V1_Extended.txt", sep = '\t' )

# # #need to append "_WGS" (or something) to each row in "case_id" of WGS file or else it will be confused with Exome I & II datasets
# # df_4["case_id"] = df_4["case_id"].astype(str) + "_WGS"
# # df_5["case_id"] = df_5["case_id"].astype(str) + "_WGS"


# # df = pd.concat( [df_1, df_2, df_3, df_4, df_5] )
# #NOTE: this excludes the pooled samples contained in "df_3"
# df = pd.concat( [df_1, df_2, df_4, df_5] )      #this excludes the "pooled" dataset

# # #ERROR WITH ASSIGNING MAX "mhc_score_2" to each mutation - for some reason, after applying this filter I am missing 'case_id' that were present in the original 'df' -> not sure why
# # #assign neoepitope with highest affinity to each mutation
# # df_mhc_max = df.iloc[ df.groupby(['case_id', 'genome_pos'])['mhc_score_2'].agg( pd.Series.idxmax ) ]
# # # df_mhc_max = df.iloc[ df.groupby(['case_id', 'genome_pos'])['ic50_score_2'].agg( pd.Series.idxmin ) ]


# ##TEST:: BEFORE
# print "TEST:: length before rank 1 filter = ", len( df )
# #filter for the highest MHC score for each mutation
# df_mhc_max = df[ df['mut_neoep_rank_1'] == 1 ]
# df_mhc_max.drop_duplicates('genome_pos', inplace = True)

# ##TEST:: AFTER
# print "TEST:: length AFTER rank 1 filter = ", len( df_mhc_max )
# print "TEST - df ALL sample names = ",  np.unique( df_mhc_max['case_id'] )

# ##TEST::
# print "df_mhc_max.dtypes:\n", df_mhc_max.dtypes

# # ##TEST:: Want to see if the altered neoepitope with the maximum MHC affinity is assigned to each genomic mutation - column "mut_neoep_rank_1" should be 1 for all mutations
# # print "df_mhc_max:\n", df_mhc_max

# print "------------ Create filters for neoepitopes ------------"
# # hash_filters = create_hash_filter( frame_stat, neoep_stat, bool_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# #DIFFERENCE: "mhc_affinity" is different between experiment 12C & 12D
# mhc_affinity = 0        #looks at neoepitope affinities to HLA alleles (e.g. LOW, MID, HIGH), 2 = select high & mid affinity neoepitopes
# frame_stat = 0
# neoep_stat = 0      #looks at all neoeps, high/mid-affinity neoeps, or high-efficacy neoeps
# check_NMD = True
# thres_gene_exp = 0      #as of now, I do not have gene expression for Velip Data
# thres_prot = 0
# thres_tap = 0
# thres_endogenous_freq = -1      #as of now, I do not have endogenous frequency for Velip Data
# hash_filters = NeoepStatistics.create_hash_filter( mhc_affinity, frame_stat, neoep_stat, check_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# # obj_neoep_stat = NeoepStatistics( df, hash_filters )
# obj_neoep_stat = NeoepStatistics( df_mhc_max, hash_filters )

# print "------------ Apply filters for neoepitopes ------------"
# #applied the filters to the dataset
# df_filtered = obj_neoep_stat.filter_neoeps()
# # filter_label = obj_neoep_stat.filter_correspond_labels( hash_filters )

# #calculate the fold change in the MHC affinity
# # def calc_affinity_fold_change( row ):
# #     row['mhc_fold_change'] = float( row['mhc_score_1'] ) / row['mhc_score_2']
# #     return row
# # df_filtered = df_filtered.apply( calc_affinity_fold_change, axis = 1 )

# # df_filtered['mhc_fold_change'] = df_filtered.apply( lambda r: (float( r['mhc_score_1'] ) / r['mhc_score_2']), axis = 1 )
# df_filtered['ic50_fold_change'] = df_filtered.apply( lambda r: (float( r['ic50_score_1'] ) / r['ic50_score_2']), axis = 1 )


# df_filtered.reset_index( inplace = True )
# ##retrieve thresholds for high-affinity, mid-affinity, & low affinity
# # df_mhc_hi = df_filtered.iloc[ df_filtered[ df_filtered['ic50_score'] >= 50 ]['ic50_score'].idxmin() ]
# # df_mhc_mid = df_filtered.iloc[ df_filtered[ df_filtered['ic50_score'] >= 500 ]['ic50_score'].idxmin() ]
# # df_mhc_low = df_filtered.iloc[ df_filtered[ df_filtered['ic50_score'] >= 5000 ]['ic50_score'].idxmin() ]
# df_mhc_hi = df_filtered.iloc[ df_filtered[ df_filtered['aff_cat_2'] == "HIGH" ]['ic50_score_2'].idxmax() ]
# df_mhc_mid = df_filtered.iloc[ df_filtered[ df_filtered['aff_cat_2'] == "MID" ]['ic50_score_2'].idxmax() ]
# df_mhc_low = df_filtered.iloc[ df_filtered[ df_filtered['aff_cat_2'] == "LOW" ]['ic50_score_2'].idxmax() ]


# ##TEST::
# print "df_mhc_hi:\n", df_mhc_hi['aff_cat_2'], " & ", df_mhc_hi['ic50_score_2'], " & ", df_mhc_hi['mhc_score_2']
# print "df_mhc_mid:\n", df_mhc_mid['aff_cat_2'], " & ", df_mhc_mid['ic50_score_2'], " & ", df_mhc_mid['mhc_score_2']
# print "df_mhc_low:\n", df_mhc_low['aff_cat_2'], " & ", df_mhc_low['ic50_score_2'], " & ", df_mhc_low['mhc_score_2']

# #Assign the general experimental condition based on the specific case ID
# items_pc_Velip = ['pcVelip1', 'pcVelip2', 'pcVelip3', 'pcVelip1_WGS', 'pcMV']
# items_br_Velip = ['BRVelip1', 'BRVelip2', 'BRVelip3', 'BRVelip1_WGS', 'BMV']
# items_pc_DMSO = ['pcDM1', 'pcDM2', 'pcDM3', 'pcDM1_WGS', 'pcMD']
# items_br_DMSO = ['BRDM1', 'BRDM2', 'BRDM3', 'BRDM1_WGS', 'BMD']
# hash_general_cond = {}
# hash_general_cond.update( {k:'pcVelip' for k in items_pc_Velip} )
# hash_general_cond.update( {k:'brVelip' for k in items_br_Velip} )
# hash_general_cond.update( {k:'pcDMSO' for k in items_pc_DMSO} )
# hash_general_cond.update( {k:'brDMSO' for k in items_br_DMSO} )
# df_filtered['gen_case_id'] = df_filtered.apply( lambda r: hash_general_cond[r['case_id']], axis = 1 )

# ##NOT SURE IF I NEED THIS YET
# list_hue_order = items_pc_Velip + items_br_Velip + items_pc_DMSO + items_br_DMSO        #for sns.barplot


# #PARAMETER 1: consider neoepitope binders - NOTE: make sure to assign 1 neoepitope to each mutation before performing this parameter
# select_neoep_type = -1
# if select_neoep_type == 1:
#     df_filtered = df_filtered[ df_filtered['aff_cat_2'].isin( ["HIGH", "MID"] ) ]
#     graph_title_append = "MHC Binders Only"
# elif select_neoep_type == 0:       #consider neoepitope nonbinders
#     df_filtered = df_filtered[ df_filtered['aff_cat_2'].isin( ["LOW", "none"] ) ]
#     graph_title_append = "MHC Non-Binders Only"
# elif select_neoep_type == -1:
#     # df_filtered = df_filtered
#     graph_title_append = "All Mutations"



# ##calculate fold change differences
# print "Fold Change < 1:\n", df_filtered[ df_filtered['ic50_fold_change'] < 1 ].groupby( ['gen_case_id'] )['gen_case_id'].count()
# print "Fold Change == 1:\n", df_filtered[ df_filtered['ic50_fold_change'] == 1 ].groupby( ['gen_case_id'] )['gen_case_id'].count()
# print "Fold Change > 1:\n", df_filtered[ df_filtered['ic50_fold_change'] > 1 ].groupby( ['gen_case_id'] )['gen_case_id'].count()
# print "Fold Change > 100:\n", df_filtered[ df_filtered['ic50_fold_change'] > 100 ].groupby( ['gen_case_id'] )['gen_case_id'].count()
# print "Fold Change > 500:\n", df_filtered[ df_filtered['ic50_fold_change'] > 500 ].groupby( ['gen_case_id'] )['gen_case_id'].count()

# print "Total # of mutations:\n", df_filtered.groupby( ['gen_case_id'] )['genome_pos'].count()

# #Begin graphing 
# #graph the results
# # sns.set( font_scale = 1.2 )
# # sns.set()
# # sns.set_style( "darkgrid" )
# # sns.set_style( "whitegrid" )
# # sns.set_style( "dark" )
# sns.set_style( "white" )

# # sns.set_context( "paper" )        #smallest features ( thinnest lines, smallest font for ticks and labels )
# # sns.set_context( "notebook" )     #3rd largest features
# sns.set_context( "talk" )     #2nd largest features
# # sns.set_context( "poster" )     #largest features ( thickest lines, largest font for ticks and labels )

# get_palette = sns.color_palette( "Blues", 16 )
# # sns.set_palette( get_palette )
# # sns.palplot( get_palette )
# #retrieve a subset of Blues to make it easier to see the graphs
# palette_subset = [get_palette[8], get_palette[10], get_palette[12], get_palette[14]]

# x_axis_order = ['brDMSO', 'brVelip', 'pcDMSO', 'pcVelip']
# dot_size = 5

# # ##GRAPH 1: Dot plot of MHC score affinities
# # # #Graph 1A - graph the the MHC scores for the neoepitopes from IEDB
# # # #Graph 1A - plot the MHC affinity score for the neoepitopes (i.e. altered neoepitopes)
# # # sns.swarmplot( x = 'gen_case_id', y = 'mhc_score_2', data = df_filtered, palette = palette_subset, order = x_axis_order, size = dot_size )
# # # graph_title = "Neoepitope Affinity Across Different Conditions"
# # # graph_ylabel = "MHC Affinity Score"
# # # #draw the threshold lines for affinity
# # # pyplot.axhline( y = df_mhc_hi['mhc_score_2'], color = get_palette[15] )
# # # pyplot.axhline( y = df_mhc_mid['mhc_score_2'], color = get_palette[10] )
# # # pyplot.axhline( y = df_mhc_low['mhc_score_2'], color = get_palette[6] )

# # #graph 1b - graph IC50 scores
# # sns.swarmplot( x = 'gen_case_id', y = 'ic50_score_2', data = df_filtered, palette = palette_subset, order = x_axis_order, size = dot_size )
# # graph_title = "Neoepitope Affinity Across Different Conditions"
# # graph_ylabel = "IC50 Affinity Score (IC50)"
# # #draw the threshold lines for affinity
# # pyplot.axhline( y = df_mhc_hi['ic50_score_2'], color = get_palette[15] )
# # pyplot.axhline( y = df_mhc_mid['ic50_score_2'], color = get_palette[10] )
# # pyplot.axhline( y = df_mhc_low['ic50_score_2'], color = get_palette[6] )



# #GRAPH 2: Dot plot of MHC affinity fold change
# #plot the MHC fold change for each sample
# sns.swarmplot( x = 'gen_case_id', y = 'ic50_fold_change', data = df_filtered, palette = palette_subset, order = x_axis_order, size = dot_size )
# graph_title = "MHC Affinity Fold Change Across Conditions (" + graph_title_append + ")"
# graph_ylabel = "MHC Affinity Fold Change (original IC50 / altered IC50)"

# sns.despine()

# pyplot.title( graph_title, fontsize = 26)
# # pyplot.title( graph_title, fontsize = 14 )
# # pyplot.xlabel( "irRECIST" )
# pyplot.xlabel( "" )
# pyplot.ylabel( graph_ylabel, fontsize = 14 )

# # pyplot.xticks(idx + 2 * width, exp_cond_order, rotation = 0 )
# pyplot.legend()

# pyplot.show()








# ##Experiment 5 - Compare the number of MHC binders between drug treatment conditions
# print "------------ Experiment 5 ------------"
# # input_file_date = "180411"      #for directory "/180403_Velip_V2"
# # input_file_date = "180605"      #for directory "/180531_Velip_V3"
# input_file_date = "180614"      #for directory "/180531_Velip_V3"


# # #ACTUAL FILE - Files that DO NOT have MHC scores ranked for each mutation
# # #Exome I dataset
# # df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File1_NeoepCompare_V1.txt", sep = '\t' )
# # #Exome II dataset
# # df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File2_NeoepCompare_V1.txt", sep = '\t' )
# # #Exome III dataset - pooled
# # df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File3_NeoepCompare_V1.txt", sep = '\t' )
# # #WGS I dataset
# # df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File4_NeoepCompare_V1.txt", sep = '\t' )
# # #WGS II dataset
# # df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File5_NeoepCompare_V1.txt", sep = '\t' )

# #ACTUAL FILE - Files that do have MHC scores ranked for each mutation
# #Exome I dataset
# df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File1_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #Exome II dataset
# df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File2_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #Exome III dataset - pooled
# df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File3_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #WGS I dataset
# df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File4_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #WGS II dataset
# df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File5_NeoepCompare_V1_Extended.txt", sep = '\t' )

# #need to append "_WGS" (or something) to each row in "case_id" of WGS file or else it will be confused with Exome I & II datasets
# df_4["case_id"] = df_4["case_id"].astype(str) + "_WGS"
# df_5["case_id"] = df_5["case_id"].astype(str) + "_WGS"

# # df = pd.concat( [df_1, df_2, df_3, df_4, df_5] )
# df = pd.concat( [df_1, df_2, df_4, df_5] )      #this excludes the "pooled" dataset

# ##TEST::
# print "df unique mutations = ", df['genome_pos'].nunique()

# ##TEST:: BEFORE
# print "TEST:: length before rank 1 filter = ", len( df )
# #filter for the highest MHC score for each mutation
# df_mhc_max = df[ df['mut_neoep_rank_1'] == 1 ]
# df_mhc_max.drop_duplicates('genome_pos', inplace = True)

# print "len( df_mhc_max ) [should be the same as df unique mutations] = ", len( df_mhc_max )

# ##TEST:: AFTER
# print "TEST:: length AFTER rank 1 filter = ", len( df_mhc_max )
# print "TEST - df ALL sample names = ",  np.unique( df_mhc_max['case_id'] )

# ##TEST::
# print "df_mhc_max.dtypes:\n", df_mhc_max.dtypes
# print "df_mhc_max:\n", df_mhc_max




# print "------------ Create filters for neoepitopes ------------"
# # hash_filters = create_hash_filter( frame_stat, neoep_stat, bool_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# #DIFFERENCE: "mhc_affinity" is different between experiment 12C & 12D
# mhc_affinity = 0        #looks at neoepitope affinities to HLA alleles (e.g. LOW, MID, HIGH), 2 = select high & mid affinity neoepitopes
# frame_stat = 0
# neoep_stat = 0      #looks at all neoeps, high/mid-affinity neoeps, or high-efficacy neoeps
# check_NMD = True
# thres_gene_exp = 0      #as of now, I do not have gene expression for Velip Data
# thres_prot = 0
# thres_tap = 0
# thres_endogenous_freq = -1      #as of now, I do not have endogenous frequency for Velip Data
# hash_filters = NeoepStatistics.create_hash_filter( mhc_affinity, frame_stat, neoep_stat, check_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# # obj_neoep_stat = NeoepStatistics( df, hash_filters )
# obj_neoep_stat = NeoepStatistics( df_mhc_max, hash_filters )

# print "------------ Apply filters for neoepitopes ------------"
# #applied the filters to the dataset
# df_filtered = obj_neoep_stat.filter_neoeps()
# # filter_label = obj_neoep_stat.filter_correspond_labels( hash_filters )


# #Attempt 2: graph pc vs. br & Velip vs. DMSO
# items_pc_Velip = ['pcVelip1', 'pcVelip2', 'pcVelip3', 'pcVelip1_WGS', 'pcMV']
# items_br_Velip = ['BRVelip1', 'BRVelip2', 'BRVelip3', 'BRVelip1_WGS', 'BMV']
# items_pc_DMSO = ['pcDM1', 'pcDM2', 'pcDM3', 'pcDM1_WGS', 'pcMD']
# items_br_DMSO = ['BRDM1', 'BRDM2', 'BRDM3', 'BRDM1_WGS', 'BMD']

# hash_case_labels = {}
# hash_case_labels.update( {x:"pcVelip" for x in items_pc_Velip} )
# hash_case_labels.update( {x:"BRVelip" for x in items_br_Velip} )
# hash_case_labels.update( {x:"pcDM" for x in items_pc_DMSO} )
# hash_case_labels.update( {x:"BRDM" for x in items_br_DMSO} )

# print "hash_case_labels = ", hash_case_labels

# #label samples either Responders or Non-Responders based on irRECIST
# # df_filtered['case_group'] = df_filtered.apply( lambda row: "Responder" if not 'disease' in row['irRECIST'].lower() else "Non-Responder", axis = 1 )
# df_filtered['case_category'] = df_filtered.apply( lambda row: hash_case_labels[row['case_id']], axis = 1 )


# #PARAMETER 1: consider neoepitope binders - NOTE: make sure to assign 1 neoepitope to each mutation before performing this parameter
# select_neoep_type = 1
# if select_neoep_type == 1:
#     df_filtered = df_filtered[ df_filtered['aff_cat_2'].isin( ["HIGH", "MID"] ) ]
#     graph_title = "Veliparib: Mutations Producing MHC Binders Per Threshold Window" 
# elif select_neoep_type == 0:       #consider neoepitope nonbinders
#     df_filtered = df_filtered[ df_filtered['aff_cat_2'].isin( ["LOW", "none"] ) ]
#     graph_title = "Veliparib: Mutations Producing MHC Non-Binders Per Threshold Window"
# elif select_neoep_type == -1:
#     # df_filtered = df_filtered
#     graph_title = "Veliparib: All Mutations Per Threshold Window"


# #PARAMETER 2: calculated # of mutations in each processing efficiency bin
# #IN ER: for neoepitopes considered "in the ER space"
# list_values = []
# window_hue_order = []       #this will be used for "hue_order" property to order boxes in seaborn.boxplot
# # step_val = 10
# # step_val = 20
# # step_val = 25
# step_val = 50
# for i_thres in range( 0, 100, step_val ):
#     thres_prot_ER_in = i_thres
#     thres_TAP_ER_in = i_thres
#     df_filtered_temp = df_filtered[ 
#     (df_filtered['proteasome_percentile_2'] >= thres_prot_ER_in) & 
#     (df_filtered['proteasome_percentile_2'] < (thres_prot_ER_in + step_val)) &
#     (df_filtered['tap_percentile_2'] >= thres_TAP_ER_in) & 
#     (df_filtered['tap_percentile_2'] < (thres_TAP_ER_in + step_val)) ]

#     series_sample_neoep_mut = df_filtered_temp.drop_duplicates('genome_pos').groupby( ['case_id', 'case_category'] )['genome_pos'].count()      #counts all neoepitopes, even if it is double-counted

#     # str_label = ">=" + str( i_thres )
#     str_label = str(i_thres) + "-" + str(i_thres + step_val)
#     window_hue_order.append( str_label )
#     version = 1
#     list_values = combine_patient_count_addlabel_v3( series_sample_neoep_mut, list_values, str_label, version )


#     ##TEST::
#     print "str_label = ", str_label, " & series_sample_neoep_mut:\n", series_sample_neoep_mut

# ##TEST::
# print "list_values of mutations from bracketed:\n"
# print list_values

# #create a dataframe for the number of neoepitopes assigned to a given sample, response type (e.g. Complete Response, Progressive Disease), the HLA subgroup (the HLA allele assigned to the neoepitope with the highest affinity)
# df_condition_neoep_mut = pd.DataFrame( list_values, columns = ['sample_id', 'case_category', 'bin_window', 'neoep_mut_count'] )


# #PARAMETER 3: Select Experiment to Run
# #JUST FOR REFERENCE: the different case_categories
# # list_case_cat = ["pcVelip", "BRVelip", "pcDM", "BRDM"]
# select_experiment = 5
# if select_experiment == 1:      #Experiment 1: pcVelip vs. Control
#     cond_1_label = "pcVelip"
#     cond_2_label = "BRDM"
#     sort_ascending = False
#     experiment_label = "(Veliparib + BRCA vs. Control)"
# elif select_experiment == 2:      #Experiment 2: brVelip vs. pcDM
#     cond_1_label = "BRVelip"
#     cond_2_label = "pcDM"
#     sort_ascending = True
#     experiment_label = "(Veliparib Only vs. BRCA)"
# elif select_experiment == 3:      #Experiment 3: pcVelip vs. pcDM
#     cond_1_label = "pcVelip"
#     cond_2_label = "pcDM"
#     sort_ascending = False
#     experiment_label = "(Veliparib + BRCA vs. BRCA)"
# elif select_experiment == 4:      #Experiment 4: brVelip vs. BRDM
#     cond_1_label = "BRVelip"
#     cond_2_label = "BRDM"
#     sort_ascending = False
#     experiment_label = "(BRCA + Veliparib vs. Control)"
# elif select_experiment == 5:      #Experiment 4: brVelip vs. BRDM
#     cond_1_label = "pcDM"
#     cond_2_label = "BRDM"
#     sort_ascending = False
#     experiment_label = "(BRCA Only vs. Control)"

# #PARAMETER 4: select the conditions of interest for the graph
# bool_show_exp_cond = False       #if this is False, then this will show all conditions in 1 graph, if True then will only show conditions of interest (in cond_1_label & cond_2_label)
# # bool_show_exp_cond = True
# if bool_show_exp_cond:
#     list_conds = [cond_1_label, cond_2_label]
#     df_condition_neoep_mut = df_condition_neoep_mut[ df_condition_neoep_mut['case_category'].isin(list_conds) ]
#     graph_title += " " + experiment_label

# #need to sort labels by therapy responses so I can see a trend - ["Progressive Disease", "Partial Response", "Complete Response"]
# df_condition_neoep_mut.sort_values( by = ['case_category'], ascending = sort_ascending, inplace = True )


# ##TEST::
# print "df_condition_neoep_mut:\n", df_condition_neoep_mut


# all_processing_bins = df_condition_neoep_mut['bin_window'].unique()
# hash_bin_ttest = {}      #hash where key = processing bin value, value = t-statistic & p-value
# for each_bin in all_processing_bins:
#     bin_cond_1 = df_condition_neoep_mut[ (df_condition_neoep_mut["case_category"].str.contains(cond_1_label, case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]['neoep_mut_count']
#     bin_cond_2 = df_condition_neoep_mut[ (df_condition_neoep_mut["case_category"].str.contains(cond_2_label, case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]['neoep_mut_count']


#     # ##TEST::
#     # print each_bin, " - T-TEST: bin_cond_1:\n", bin_cond_1
#     # print each_bin, " - T-TEST: bin_cond_2:\n", bin_cond_2

#     t_stat, p_val = ttest_ind( bin_cond_1.values, bin_cond_2.values )
#     hash_bin_ttest[each_bin] = (t_stat, p_val)

#     # ##TEST:: make sure I'm extracting the correct groups
#     # bin_responder = df_condition_neoep_mut[ (~df_condition_neoep_mut["case_category"].str.contains("Non-Respon", case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]
#     # bin_nonresponder = df_condition_neoep_mut[ (df_condition_neoep_mut["case_category"].str.contains("Non-Respon", case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]
#     # print each_bin, " - bin_responder =\n", bin_responder
#     # print each_bin, " - bin_nonresponder =\n", bin_nonresponder
#     # print "---------\n"


# #show the t-test values from all bins
# print "T-TEST STATISTICS!!! hash_bin_ttest - compare " + cond_1_label + " to " + cond_2_label + ":\n", hash_bin_ttest
# # for k,v in hash_bin_ttest.iteritems():
# #     print "bin = ", k, " & t-stat = ", v[0], " & p-value = ", v[1]

# all_processing_bins = sorted( df_condition_neoep_mut['bin_window'].unique() )
# for each_bin in all_processing_bins:
#     v = hash_bin_ttest[each_bin]
#     print "bin = ", each_bin, " & t-stat = ", v[0], " & p-value = ", v[1]


# # sns.set()
# # sns.set_style( "darkgrid" )
# # sns.set_context( "paper" )        #smallest features ( thinnest lines, smallest font for ticks and labels )
# # sns.set_context( "notebook" )     #3rd largest features
# sns.set_context( "talk" )     #2nd largest features
# # sns.set_context( "poster" )     #largest features ( thickest lines, largest font for ticks and labels )
# #custom color palette
# # current_palette = sns.color_palette( "husl", 12 )        #color-wheel palette, for now skip colors between indices 1 to 3
# # # sns.set_palette( [ current_palette[4] ] )
# current_palette = sns.color_palette( "Blues" )
# sns.set_palette( current_palette )

# #graph the results
# #graph Seaborn pyplot try to create a box-whisker plot for each Proteasome
# ax = sns.boxplot( x = 'case_category', y = 'neoep_mut_count', hue = 'bin_window', data = df_condition_neoep_mut, showfliers = False, hue_order = window_hue_order, palette = "Blues" )     #this produces a rainbow of pastel colors by default

# # #graph both boxplot + swarmplot
# # sns.boxplot( x = col_response, y = 'neoep_mut_count', hue = 'bin_window', data = df_condition_neoep_mut, showfliers = False, hue_order = window_hue_order, palette = "Blues" )     #this produces a rainbow of pastel colors by default
# # sns.swarmplot( x = col_response, y = 'neoep_mut_count', hue = 'bin_window', data = df_condition_neoep_mut, hue_order = window_hue_order, color = ".25" )      #Use this if I want to show the dots
# # ax = sns.boxplot( x = 'response_type', y = 'neoep_mut_count', hue = 'bin_window', data = df_condition_neoep_mut, showfliers = False, hue_order = window_hue_order, palette = "Blues" )
# sns.despine()
# # pyplot.title( "Neoepitope Count Per irRECIST Condition: Unique Neoepitopes" )
# # pyplot.title( "Exp 12: Neoepitope Count for HLA Subgroup Per irRECIST Condition: " + annex_title, fontsize = 10 )
# # pyplot.title( graph_title )
# graph_title_final = "Processing of MHC Binders: Veliparib & BRCA"
# pyplot.title( graph_title_final, fontsize = 30 )
# # pyplot.title( graph_title, fontsize = 10 )
# # pyplot.xlabel( "irRECIST" )
# pyplot.xlabel( "" )
# pyplot.ylabel( "Number of Mutations", fontsize = 24 )
# pyplot.yticks(np.arange(0, 21, 5))

# # ax.legend().set_title( "Processing\nEfficiency" )     #This is how to change the legend title
# # pyplot.legend( fontsize = 14 )
# # # pyplot.legend( loc = "upper left", fontsize = 14 )

# ax.legend().set_title( "Processing\nEfficiency" )     #This is how to change the legend title
# # ax.legend().set_title( "Window" )     #This is how to change the legend title
# # pyplot.legend( fontsize = 14 )
# pyplot.setp(ax.get_legend().get_texts(), fontsize = '14') # for legend text

# pyplot.show()











# ##Experiment 6 - Compare the number of MHC binders between drug treatment conditions
# print "------------ Experiment 6 ------------"
# # input_file_date = "180411"      #for directory "/180403_Velip_V2"
# # input_file_date = "180605"      #for directory "/180531_Velip_V3"
# input_file_date = "180614"      #for directory "/180531_Velip_V3"


# # #ACTUAL FILE - Files that DO NOT have MHC scores ranked for each mutation
# # #Exome I dataset
# # df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File1_NeoepCompare_V1.txt", sep = '\t' )
# # #Exome II dataset
# # df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File2_NeoepCompare_V1.txt", sep = '\t' )
# # #Exome III dataset - pooled
# # df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File3_NeoepCompare_V1.txt", sep = '\t' )
# # #WGS I dataset
# # df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File4_NeoepCompare_V1.txt", sep = '\t' )
# # #WGS II dataset
# # df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File5_NeoepCompare_V1.txt", sep = '\t' )

# #ACTUAL FILE - Files that do have MHC scores ranked for each mutation
# #Exome I dataset
# df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File1_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #Exome II dataset
# df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File2_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #Exome III dataset - pooled
# df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File3_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #WGS I dataset
# df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File4_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #WGS II dataset
# df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File5_NeoepCompare_V1_Extended.txt", sep = '\t' )


# ##NOTE: I do not need this since I already labeled files with "_WGS" labels
# # #need to append "_WGS" (or something) to each row in "case_id" of WGS file or else it will be confused with Exome I & II datasets
# # df_4["case_id"] = df_4["case_id"].astype(str) + "_WGS"
# # df_5["case_id"] = df_5["case_id"].astype(str) + "_WGS"

# # df = pd.concat( [df_1, df_2, df_3, df_4, df_5] )
# #NOTE: this excludes the pooled samples contained in "df_3"
# df = pd.concat( [df_1, df_2, df_4, df_5] )      #this excludes the "pooled" dataset

# ##TEST::
# print "df unique mutations = ", df['genome_pos'].nunique()

# ##TEST:: BEFORE
# print "TEST:: length before rank 1 filter = ", len( df )
# #filter for the highest MHC score for each mutation
# df_mhc_max = df[ df['mut_neoep_rank_1'] == 1 ]
# df_mhc_max.drop_duplicates('genome_pos', inplace = True)

# print "len( df_mhc_max ) [should be the same as df unique mutations] = ", len( df_mhc_max )

# ##TEST:: AFTER
# print "TEST:: length AFTER rank 1 filter = ", len( df_mhc_max )
# print "TEST - df ALL sample names = ",  np.unique( df_mhc_max['case_id'] )

# ##TEST::
# print "df_mhc_max.dtypes:\n", df_mhc_max.dtypes
# print "df_mhc_max:\n", df_mhc_max




# print "------------ Create filters for neoepitopes ------------"
# # hash_filters = create_hash_filter( frame_stat, neoep_stat, bool_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# #DIFFERENCE: "mhc_affinity" is different between experiment 12C & 12D
# mhc_affinity = 0        #looks at neoepitope affinities to HLA alleles (e.g. LOW, MID, HIGH), 2 = select high & mid affinity neoepitopes
# frame_stat = 0
# neoep_stat = 0      #looks at all neoeps, high/mid-affinity neoeps, or high-efficacy neoeps
# check_NMD = True
# thres_gene_exp = 0      #as of now, I do not have gene expression for Velip Data
# thres_prot = 0
# thres_tap = 0
# thres_endogenous_freq = -1      #as of now, I do not have endogenous frequency for Velip Data
# hash_filters = NeoepStatistics.create_hash_filter( mhc_affinity, frame_stat, neoep_stat, check_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# # obj_neoep_stat = NeoepStatistics( df, hash_filters )
# obj_neoep_stat = NeoepStatistics( df_mhc_max, hash_filters )

# print "------------ Apply filters for neoepitopes ------------"
# #applied the filters to the dataset
# df_filtered = obj_neoep_stat.filter_neoeps()
# # filter_label = obj_neoep_stat.filter_correspond_labels( hash_filters )
# df_filtered['ic50_fold_change'] = df_filtered.apply( lambda r: (float( r['ic50_score_1'] ) / r['ic50_score_2']), axis = 1 )

# df_filtered.reset_index( inplace = True )

# #Attempt 2: graph pc vs. br & Velip vs. DMSO
# items_pc_Velip = ['pcVelip1', 'pcVelip2', 'pcVelip3', 'pcVelip1_WGS', 'pcMV']
# items_br_Velip = ['BRVelip1', 'BRVelip2', 'BRVelip3', 'BRVelip1_WGS', 'BMV']
# items_pc_DMSO = ['pcDM1', 'pcDM2', 'pcDM3', 'pcDM1_WGS', 'pcMD']
# items_br_DMSO = ['BRDM1', 'BRDM2', 'BRDM3', 'BRDM1_WGS', 'BMD']

# hash_case_labels = {}
# hash_case_labels.update( {x:"pcVelip" for x in items_pc_Velip} )
# hash_case_labels.update( {x:"BRVelip" for x in items_br_Velip} )
# hash_case_labels.update( {x:"pcDM" for x in items_pc_DMSO} )
# hash_case_labels.update( {x:"BRDM" for x in items_br_DMSO} )

# print "hash_case_labels = ", hash_case_labels

# #label samples either Responders or Non-Responders based on irRECIST
# # df_filtered['case_group'] = df_filtered.apply( lambda row: "Responder" if not 'disease' in row['irRECIST'].lower() else "Non-Responder", axis = 1 )
# df_filtered['case_category'] = df_filtered.apply( lambda row: hash_case_labels[row['case_id']], axis = 1 )


# #PARAMETER 1: consider neoepitope binders - NOTE: make sure to assign 1 neoepitope to each mutation before performing this parameter
# select_neoep_type = 1
# if select_neoep_type == 1:
#     df_filtered = df_filtered[ df_filtered['aff_cat_2'].isin( ["HIGH", "MID"] ) ]
#     graph_title_append = "MHC Binders Only"
# elif select_neoep_type == 0:       #consider neoepitope nonbinders
#     df_filtered = df_filtered[ df_filtered['aff_cat_2'].isin( ["LOW", "none"] ) ]
#     graph_title_append = "MHC Non-Binders Only"
# elif select_neoep_type == -1:
#     # df_filtered = df_filtered
#     graph_title_append = "All Mutations"


# #PARAMETER 2: calculated # of mutations in each processing efficiency bin
# #IN ER: for neoepitopes considered "in the ER space"
# #COLUMN 1: IC50 score
# col_oi_df = "ic50_score_2";
# y_axis_label = "MHC-Neoepitope Affinity (IC50)"
# # #COLUMN 2: calculated MHC affinity
# # col_oi_df = "mhc_score_2";
# # y_axis_label = "MHC Affinity Score"
# # #COLUMN 3: MHC fold change
# # col_oi_df = "ic50_fold_change";
# # y_axis_label = "Affinity Fold Change (IC50)"
# list_values = []
# window_hue_order = []       #this will be used for "hue_order" property to order boxes in seaborn.boxplot
# # step_val = 10
# # step_val = 20
# # step_val = 25
# step_val = 50
# for i_thres in range( 0, 100, step_val ):
#     thres_prot_ER_in = i_thres
#     thres_TAP_ER_in = i_thres
#     df_filtered_temp = df_filtered[ 
#     (df_filtered['proteasome_percentile_2'] >= thres_prot_ER_in) & 
#     (df_filtered['proteasome_percentile_2'] < (thres_prot_ER_in + step_val)) &
#     (df_filtered['tap_percentile_2'] >= thres_TAP_ER_in) & 
#     (df_filtered['tap_percentile_2'] < (thres_TAP_ER_in + step_val)) ]
#     df_filtered_temp = df_filtered_temp.drop_duplicates('genome_pos');

#     #record the array of values for a "col_oi" (column of interest) based on the filtered DataFrame
#     # str_label = ">=" + str( i_thres )
#     str_label = str(i_thres) + "-" + str(i_thres + step_val)
#     window_hue_order.append( str_label )
#     version = 1
#     # list_values = combine_patient_count_addlabel_v3( series_sample_neoep_mut, list_values, str_label, version )
#     list_values = combine_patient_count_addlabel_v4( df_filtered_temp, col_oi_df, list_values, str_label, version )


# ##TEST::
# print "list_values of mutations from bracketed:\n"
# print list_values

# #create a dataframe for the number of neoepitopes assigned to a given sample, response type (e.g. Complete Response, Progressive Disease), the HLA subgroup (the HLA allele assigned to the neoepitope with the highest affinity)
# df_condition_neoep_mut = pd.DataFrame( list_values, columns = ['sample_id', 'case_category', 'bin_window', 'neoep_affinity'] )


# #PARAMETER 3: Select Experiment to Run
# #JUST FOR REFERENCE: the different case_categories
# # list_case_cat = ["pcVelip", "BRVelip", "pcDM", "BRDM"]
# select_experiment = 5
# if select_experiment == 1:      #Experiment 1: pcVelip vs. Control
#     cond_1_label = "pcVelip"
#     cond_2_label = "BRDM"
#     sort_ascending = False
#     experiment_label = "(Veliparib + BRCA vs. Control)"
# elif select_experiment == 2:      #Experiment 2: brVelip vs. pcDM
#     cond_1_label = "BRVelip"
#     cond_2_label = "pcDM"
#     sort_ascending = True
#     experiment_label = "(Veliparib Only vs. BRCA)"
# elif select_experiment == 3:      #Experiment 3: pcVelip vs. pcDM
#     cond_1_label = "pcVelip"
#     cond_2_label = "pcDM"
#     sort_ascending = False
#     experiment_label = "(Veliparib + BRCA vs. BRCA)"
# elif select_experiment == 4:      #Experiment 4: brVelip vs. BRDM
#     cond_1_label = "BRVelip"
#     cond_2_label = "BRDM"
#     sort_ascending = False
#     experiment_label = "(BRCA + Veliparib vs. Control)"
# elif select_experiment == 5:      #Experiment 4: brVelip vs. BRDM
#     cond_1_label = "pcDM"
#     cond_2_label = "BRDM"
#     sort_ascending = False
#     experiment_label = "(BRCA Only vs. Control)"

# #PARAMETER 4: select the conditions of interest for the graph
# bool_show_exp_cond = False       #if this is False, then this will show all conditions in 1 graph, if True then will only show conditions of interest (in cond_1_label & cond_2_label)
# # bool_show_exp_cond = True
# if bool_show_exp_cond:
#     list_conds = [cond_1_label, cond_2_label]
#     df_condition_neoep_mut = df_condition_neoep_mut[ df_condition_neoep_mut['case_category'].isin(list_conds) ]
#     graph_title += " " + experiment_label

# #need to sort labels by therapy responses so I can see a trend - ["Progressive Disease", "Partial Response", "Complete Response"]
# df_condition_neoep_mut.sort_values( by = ['case_category'], ascending = sort_ascending, inplace = True )


# ##TEST::
# print "df_condition_neoep_mut:\n", df_condition_neoep_mut


# all_processing_bins = df_condition_neoep_mut['bin_window'].unique()
# hash_bin_ttest = {}      #hash where key = processing bin value, value = t-statistic & p-value
# for each_bin in all_processing_bins:
#     bin_cond_1 = df_condition_neoep_mut[ (df_condition_neoep_mut["case_category"].str.contains(cond_1_label, case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]['neoep_affinity']
#     bin_cond_2 = df_condition_neoep_mut[ (df_condition_neoep_mut["case_category"].str.contains(cond_2_label, case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]['neoep_affinity']


#     # ##TEST::
#     # print each_bin, " - T-TEST: bin_cond_1:\n", bin_cond_1
#     # print each_bin, " - T-TEST: bin_cond_2:\n", bin_cond_2

#     t_stat, p_val = ttest_ind( bin_cond_1.values, bin_cond_2.values )
#     hash_bin_ttest[each_bin] = (t_stat, p_val)

#     # ##TEST:: make sure I'm extracting the correct groups
#     # bin_responder = df_condition_neoep_mut[ (~df_condition_neoep_mut["case_category"].str.contains("Non-Respon", case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]
#     # bin_nonresponder = df_condition_neoep_mut[ (df_condition_neoep_mut["case_category"].str.contains("Non-Respon", case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]
#     # print each_bin, " - bin_responder =\n", bin_responder
#     # print each_bin, " - bin_nonresponder =\n", bin_nonresponder
#     # print "---------\n"


# #show the t-test values from all bins
# print "T-TEST STATISTICS!!! hash_bin_ttest - compare " + cond_1_label + " to " + cond_2_label + ":\n", hash_bin_ttest
# # for k,v in hash_bin_ttest.iteritems():
# #     print "bin = ", k, " & t-stat = ", v[0], " & p-value = ", v[1]

# all_processing_bins = sorted( df_condition_neoep_mut['bin_window'].unique() )
# for each_bin in all_processing_bins:
#     v = hash_bin_ttest[each_bin]
#     print "bin = ", each_bin, " & t-stat = ", v[0], " & p-value = ", v[1]


# # sns.set()
# # sns.set_style( "darkgrid" )
# # sns.set_context( "paper" )        #smallest features ( thinnest lines, smallest font for ticks and labels )
# # sns.set_context( "notebook" )     #3rd largest features
# sns.set_context( "talk" )     #2nd largest features
# # sns.set_context( "poster" )     #largest features ( thickest lines, largest font for ticks and labels )
# #custom color palette
# # current_palette = sns.color_palette( "husl", 12 )        #color-wheel palette, for now skip colors between indices 1 to 3
# # # sns.set_palette( [ current_palette[4] ] )
# # current_palette = sns.color_palette( "Blues" )
# current_palette = ["#008cff", "#ff0048"]
# sns.set_palette( current_palette )


# x_axis_order = ['BRDM', 'BRVelip', 'pcDM', 'pcVelip']
# # x_axis_order = ['brDMSO', 'brVelip', 'pcDMSO', 'pcVelip']

# dot_size = 7

# #graph the results
# #GRAPH: Boxplot
# #graph Seaborn pyplot try to create a box-whisker plot for each Proteasome
# # # ax = sns.boxplot( x = 'case_category', y = 'neoep_affinity', hue = 'bin_window', data = df_condition_neoep_mut, showfliers = False, hue_order = window_hue_order, palette = "Blues" )
# # ax = sns.boxplot( x = 'case_category', y = 'neoep_affinity', hue = 'bin_window', data = df_condition_neoep_mut, showfliers = False, order = x_axis_order, hue_order = window_hue_order, palette = "Blues" )     #This is the same as the graph above, but this has a defined x-axis order

# #GRAPH: Swarm Plot
# ax = sns.swarmplot( x = 'case_category', y = 'neoep_affinity', hue = 'bin_window', data = df_condition_neoep_mut, order = x_axis_order, hue_order = window_hue_order, size = dot_size )      #Use this if I want to show the dots

# # #GRAPH: Boxplot + Swarm Plot
# # sns.boxplot( x = 'case_category', y = 'neoep_affinity', hue = 'bin_window', data = df_condition_neoep_mut, showfliers = False, hue_order = window_hue_order, palette = "Blues" )     #this produces a rainbow of pastel colors by default
# # sns.swarmplot( x = 'case_category', y = 'neoep_affinity', hue = 'bin_window', data = df_condition_neoep_mut, hue_order = window_hue_order, color = ".25" )      #Use this if I want to show the dots
# # # ax = sns.boxplot( x = 'response_type', y = 'neoep_affinity', hue = 'bin_window', data = df_condition_neoep_mut, showfliers = False, hue_order = window_hue_order, palette = "Blues" )
# sns.despine()
# # pyplot.title( "Neoepitope Count Per irRECIST Condition: Unique Neoepitopes" )
# # pyplot.title( "Exp 12: Neoepitope Count for HLA Subgroup Per irRECIST Condition: " + annex_title, fontsize = 10 )
# # pyplot.title( graph_title )
# graph_title_final = "Processing and MHC Affinity: Veliparib & BRCA (" + graph_title_append + ")"
# pyplot.title( graph_title_final, fontsize = 24 )
# # pyplot.title( graph_title, fontsize = 10 )
# # pyplot.xlabel( "irRECIST" )
# pyplot.xlabel( "" )
# pyplot.ylabel( y_axis_label, fontsize = 20 )
# # pyplot.yticks(np.arange(0, 800, 100))           #Scale for: IC50
# # pyplot.yticks(np.arange(-3, 1, 0.5))          #Scale for: MHC affinity score

# # ax.legend().set_title( "Processing\nEfficiency" )     #This is how to change the legend title
# # pyplot.legend( fontsize = 14 )
# # # pyplot.legend( loc = "upper left", fontsize = 14 )

# ax.legend().set_title( "Processing\nEfficiency" )     #This is how to change the legend title
# # ax.legend().set_title( "Window" )     #This is how to change the legend title
# # pyplot.legend( fontsize = 14 )
# pyplot.setp(ax.get_legend().get_texts(), fontsize = '14') # for legend text

# pyplot.show()









# ##Experiment 7 - Compare the number of MHC binders between drug treatment conditions. This is the same as Experiment 6 but this will also incorporates gene expression
# print "------------ Experiment 7 ------------"
# # input_file_date = "180411"      #for directory "/180403_Velip_V2"
# # input_file_date = "180605"      #for directory "/180531_Velip_V3"
# input_file_date = "180614"      #for directory "/180531_Velip_V3"


# # #ACTUAL FILE - Files that DO NOT have MHC scores ranked for each mutation
# # #Exome I dataset
# # df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File1_NeoepCompare_V1.txt", sep = '\t' )
# # #Exome II dataset
# # df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File2_NeoepCompare_V1.txt", sep = '\t' )
# # #Exome III dataset - pooled
# # df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File3_NeoepCompare_V1.txt", sep = '\t' )
# # #WGS I dataset
# # df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File4_NeoepCompare_V1.txt", sep = '\t' )
# # #WGS II dataset
# # df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File5_NeoepCompare_V1.txt", sep = '\t' )

# #ACTUAL FILE - Files that do have MHC scores ranked for each mutation
# #Exome I dataset
# df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File1_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #Exome II dataset
# df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File2_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #Exome III dataset - pooled
# df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File3_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #WGS I dataset
# df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File4_NeoepCompare_V1_Extended.txt", sep = '\t' )
# #WGS II dataset
# df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File5_NeoepCompare_V1_Extended.txt", sep = '\t' )


# ##NOTE: I do not need this since I already labeled files with "_WGS" labels
# # #need to append "_WGS" (or something) to each row in "case_id" of WGS file or else it will be confused with Exome I & II datasets
# # df_4["case_id"] = df_4["case_id"].astype(str) + "_WGS"
# # df_5["case_id"] = df_5["case_id"].astype(str) + "_WGS"

# # df = pd.concat( [df_1, df_2, df_3, df_4, df_5] )
# #NOTE: this excludes the pooled samples contained in "df_3"
# df = pd.concat( [df_1, df_2, df_4, df_5] )      #this excludes the "pooled" dataset

# ##TEST::
# print "df unique mutations = ", df['genome_pos'].nunique()

# ##TEST:: BEFORE
# print "TEST:: length before rank 1 filter = ", len( df )
# #filter for the highest MHC score for each mutation
# df_mhc_max = df[ df['mut_neoep_rank_1'] == 1 ]
# df_mhc_max.drop_duplicates('genome_pos', inplace = True)

# print "len( df_mhc_max ) [should be the same as df unique mutations] = ", len( df_mhc_max )

# ##TEST:: AFTER
# print "TEST:: length AFTER rank 1 filter = ", len( df_mhc_max )
# print "TEST - df ALL sample names = ",  np.unique( df_mhc_max['case_id'] )

# ##TEST::
# print "df_mhc_max.dtypes:\n", df_mhc_max.dtypes
# print "df_mhc_max:\n", df_mhc_max




# print "------------ Create filters for neoepitopes ------------"
# # hash_filters = create_hash_filter( frame_stat, neoep_stat, bool_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# #DIFFERENCE: "mhc_affinity" is different between experiment 12C & 12D
# mhc_affinity = 0        #looks at neoepitope affinities to HLA alleles (e.g. LOW, MID, HIGH), 2 = select high & mid affinity neoepitopes
# frame_stat = 0
# neoep_stat = 0      #looks at all neoeps, high/mid-affinity neoeps, or high-efficacy neoeps
# check_NMD = True
# thres_gene_exp = 0      #as of now, I do not have gene expression for Velip Data
# thres_prot = 0
# thres_tap = 0
# thres_endogenous_freq = -1      #as of now, I do not have endogenous frequency for Velip Data
# hash_filters = NeoepStatistics.create_hash_filter( mhc_affinity, frame_stat, neoep_stat, check_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# # obj_neoep_stat = NeoepStatistics( df, hash_filters )
# obj_neoep_stat = NeoepStatistics( df_mhc_max, hash_filters )

# print "------------ Apply filters for neoepitopes ------------"
# #applied the filters to the dataset
# df_filtered = obj_neoep_stat.filter_neoeps()
# # filter_label = obj_neoep_stat.filter_correspond_labels( hash_filters )
# df_filtered['ic50_fold_change'] = df_filtered.apply( lambda r: (float( r['ic50_score_1'] ) / r['ic50_score_2']), axis = 1 )

# df_filtered.reset_index( inplace = True )

# #Attempt 2: graph pc vs. br & Velip vs. DMSO
# items_pc_Velip = ['pcVelip1', 'pcVelip2', 'pcVelip3', 'pcVelip1_WGS', 'pcMV']
# items_br_Velip = ['BRVelip1', 'BRVelip2', 'BRVelip3', 'BRVelip1_WGS', 'BMV']
# items_pc_DMSO = ['pcDM1', 'pcDM2', 'pcDM3', 'pcDM1_WGS', 'pcMD']
# items_br_DMSO = ['BRDM1', 'BRDM2', 'BRDM3', 'BRDM1_WGS', 'BMD']

# hash_case_labels = {}
# hash_case_labels.update( {x:"pcVelip" for x in items_pc_Velip} )
# hash_case_labels.update( {x:"BRVelip" for x in items_br_Velip} )
# hash_case_labels.update( {x:"pcDM" for x in items_pc_DMSO} )
# hash_case_labels.update( {x:"BRDM" for x in items_br_DMSO} )

# print "hash_case_labels = ", hash_case_labels

# #label samples either Responders or Non-Responders based on irRECIST
# # df_filtered['case_group'] = df_filtered.apply( lambda row: "Responder" if not 'disease' in row['irRECIST'].lower() else "Non-Responder", axis = 1 )
# df_filtered['case_category'] = df_filtered.apply( lambda row: hash_case_labels[row['case_id']], axis = 1 )


# #PARAMETER 1: consider neoepitope binders - NOTE: make sure to assign 1 neoepitope to each mutation before performing this parameter
# select_neoep_type = 1
# if select_neoep_type == 1:
#     df_filtered = df_filtered[ df_filtered['aff_cat_2'].isin( ["HIGH", "MID"] ) ]
#     graph_title_append = "MHC Binders Only"
# elif select_neoep_type == 0:       #consider neoepitope nonbinders
#     df_filtered = df_filtered[ df_filtered['aff_cat_2'].isin( ["LOW", "none"] ) ]
#     graph_title_append = "MHC Non-Binders Only"
# elif select_neoep_type == -1:
#     # df_filtered = df_filtered
#     graph_title_append = "All Mutations"

# #PARAMETER 1B: gene expression
# filter_gene_express = 1
# if filter_gene_express == 1:
#     df_filtered = df_filtered[ df_filtered['gene_percentile'] > 0]
#     graph_title_append += ", Expression > 0"
# else:
#     graph_title_append += ", Expression Not Considered"


# #PARAMETER 2: calculated # of mutations in each processing efficiency bin
# #IN ER: for neoepitopes considered "in the ER space"
# # #COLUMN 1: IC50 score
# # col_oi_df = "ic50_score_2";
# # y_axis_label = "MHC-Neoepitope Affinity (IC50)"
# # #COLUMN 2: calculated MHC affinity
# # col_oi_df = "mhc_score_2";
# # y_axis_label = "MHC Affinity Score"
# #COLUMN 3: MHC fold change
# col_oi_df = "ic50_fold_change";
# y_axis_label = "Affinity Fold Change (IC50)"
# list_values = []
# window_hue_order = []       #this will be used for "hue_order" property to order boxes in seaborn.boxplot
# # step_val = 10
# # step_val = 20
# # step_val = 25
# step_val = 50
# for i_thres in range( 0, 100, step_val ):
#     thres_prot_ER_in = i_thres
#     thres_TAP_ER_in = i_thres
#     df_filtered_temp = df_filtered[ 
#     (df_filtered['proteasome_percentile_2'] >= thres_prot_ER_in) & 
#     (df_filtered['proteasome_percentile_2'] < (thres_prot_ER_in + step_val)) &
#     (df_filtered['tap_percentile_2'] >= thres_TAP_ER_in) & 
#     (df_filtered['tap_percentile_2'] < (thres_TAP_ER_in + step_val)) ]
#     df_filtered_temp = df_filtered_temp.drop_duplicates('genome_pos');

#     #record the array of values for a "col_oi" (column of interest) based on the filtered DataFrame
#     # str_label = ">=" + str( i_thres )
#     str_label = str(i_thres) + "-" + str(i_thres + step_val)
#     window_hue_order.append( str_label )
#     version = 1
#     # list_values = combine_patient_count_addlabel_v3( series_sample_neoep_mut, list_values, str_label, version )
#     list_values = combine_patient_count_addlabel_v4( df_filtered_temp, col_oi_df, list_values, str_label, version )


# ##TEST:: "list_values" records the counts associated with a specific combination of conditions (e.g. [#, "brVelip", "50-100"])
# print "list_values of mutations from bracketed:\n"
# print list_values

# #create a dataframe for the number of neoepitopes assigned to a given sample, response type (e.g. Complete Response, Progressive Disease), the HLA subgroup (the HLA allele assigned to the neoepitope with the highest affinity)
# df_condition_neoep_mut = pd.DataFrame( list_values, columns = ['sample_id', 'case_category', 'bin_window', 'neoep_affinity'] )


# #PARAMETER 3: Select Experiment to Run
# #JUST FOR REFERENCE: the different case_categories
# # list_case_cat = ["pcVelip", "BRVelip", "pcDM", "BRDM"]
# select_experiment = 5
# if select_experiment == 1:      #Experiment 1: pcVelip vs. Control
#     cond_1_label = "pcVelip"
#     cond_2_label = "BRDM"
#     sort_ascending = False
#     experiment_label = "(Veliparib + BRCA vs. Control)"
# elif select_experiment == 2:      #Experiment 2: brVelip vs. pcDM
#     cond_1_label = "BRVelip"
#     cond_2_label = "pcDM"
#     sort_ascending = True
#     experiment_label = "(Veliparib Only vs. BRCA)"
# elif select_experiment == 3:      #Experiment 3: pcVelip vs. pcDM
#     cond_1_label = "pcVelip"
#     cond_2_label = "pcDM"
#     sort_ascending = False
#     experiment_label = "(Veliparib + BRCA vs. BRCA)"
# elif select_experiment == 4:      #Experiment 4: brVelip vs. BRDM
#     cond_1_label = "BRVelip"
#     cond_2_label = "BRDM"
#     sort_ascending = False
#     experiment_label = "(BRCA + Veliparib vs. Control)"
# elif select_experiment == 5:      #Experiment 4: brVelip vs. BRDM
#     cond_1_label = "pcDM"
#     cond_2_label = "BRDM"
#     sort_ascending = False
#     experiment_label = "(BRCA Only vs. Control)"

# #PARAMETER 4: select the conditions of interest for the graph
# bool_show_exp_cond = False       #if this is False, then this will show all conditions in 1 graph, if True then will only show conditions of interest (in cond_1_label & cond_2_label)
# # bool_show_exp_cond = True
# if bool_show_exp_cond:
#     list_conds = [cond_1_label, cond_2_label]
#     df_condition_neoep_mut = df_condition_neoep_mut[ df_condition_neoep_mut['case_category'].isin(list_conds) ]
#     graph_title += " " + experiment_label

# #need to sort labels by therapy responses so I can see a trend - ["Progressive Disease", "Partial Response", "Complete Response"]
# df_condition_neoep_mut.sort_values( by = ['case_category'], ascending = sort_ascending, inplace = True )


# ##TEST:: print "df_condition_neoep_mut:\n", df_condition_neoep_mut


# #Perform t-test
# print "NOTE: this is comparing the neoepitope affinity score (CONJ: perhaps the average?) between conditions, not the number of neoepitopes in each condition"

# #Perform T-test VERSION 2
# # calc_ttest_btwn_conds( df_condition_neoep_mut, cond_1_label, cond_2_label, "neoep_affinity" )
# calc_ttest_btwn_conds( df_condition_neoep_mut, "BRVelip", "BRDM", "neoep_affinity" )
# calc_ttest_btwn_conds( df_condition_neoep_mut, "pcVelip", "pcDM", "neoep_affinity" )

# # #Perform T-test VERSION 1
# # all_processing_bins = df_condition_neoep_mut['bin_window'].unique()
# # hash_bin_ttest = {}      #hash where key = processing bin value, value = t-statistic & p-value
# # for each_bin in all_processing_bins:
# #     bin_cond_1 = df_condition_neoep_mut[ (df_condition_neoep_mut["case_category"].str.contains(cond_1_label, case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]['neoep_affinity']
# #     bin_cond_2 = df_condition_neoep_mut[ (df_condition_neoep_mut["case_category"].str.contains(cond_2_label, case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]['neoep_affinity']


# #     # ##TEST::
# #     # print each_bin, " - T-TEST: bin_cond_1:\n", bin_cond_1
# #     # print each_bin, " - T-TEST: bin_cond_2:\n", bin_cond_2

# #     t_stat, p_val = ttest_ind( bin_cond_1.values, bin_cond_2.values )
# #     hash_bin_ttest[each_bin] = (t_stat, p_val)

# #     # ##TEST:: make sure I'm extracting the correct groups
# #     # bin_responder = df_condition_neoep_mut[ (~df_condition_neoep_mut["case_category"].str.contains("Non-Respon", case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]
# #     # bin_nonresponder = df_condition_neoep_mut[ (df_condition_neoep_mut["case_category"].str.contains("Non-Respon", case = False) ) & (df_condition_neoep_mut["bin_window"] == each_bin) ]
# #     # print each_bin, " - bin_responder =\n", bin_responder
# #     # print each_bin, " - bin_nonresponder =\n", bin_nonresponder
# #     # print "---------\n"


# # #show the t-test values from all bins
# # print "T-TEST STATISTICS!!! hash_bin_ttest - compare " + cond_1_label + " to " + cond_2_label + ":\n", hash_bin_ttest
# # # for k,v in hash_bin_ttest.iteritems():
# # #     print "bin = ", k, " & t-stat = ", v[0], " & p-value = ", v[1]

# # all_processing_bins = sorted( df_condition_neoep_mut['bin_window'].unique() )
# # for each_bin in all_processing_bins:
# #     v = hash_bin_ttest[each_bin]
# #     print "bin = ", each_bin, " & t-stat = ", v[0], " & p-value = ", v[1]


# # sns.set()
# # sns.set_style( "darkgrid" )
# # sns.set_context( "paper" )        #smallest features ( thinnest lines, smallest font for ticks and labels )
# # sns.set_context( "notebook" )     #3rd largest features
# sns.set_context( "talk" )     #2nd largest features
# # sns.set_context( "poster" )     #largest features ( thickest lines, largest font for ticks and labels )
# #custom color palette
# # current_palette = sns.color_palette( "husl", 12 )        #color-wheel palette, for now skip colors between indices 1 to 3
# # # sns.set_palette( [ current_palette[4] ] )
# # current_palette = sns.color_palette( "Blues" )
# current_palette = ["#008cff", "#ff0048"]
# sns.set_palette( current_palette )


# x_axis_order = ['BRDM', 'BRVelip', 'pcDM', 'pcVelip']
# # x_axis_order = ['brDMSO', 'brVelip', 'pcDMSO', 'pcVelip']

# dot_size = 7

# #graph the results
# #GRAPH: Boxplot
# #graph Seaborn pyplot try to create a box-whisker plot for each Proteasome
# # # ax = sns.boxplot( x = 'case_category', y = 'neoep_affinity', hue = 'bin_window', data = df_condition_neoep_mut, showfliers = False, hue_order = window_hue_order, palette = "Blues" )
# # ax = sns.boxplot( x = 'case_category', y = 'neoep_affinity', hue = 'bin_window', data = df_condition_neoep_mut, showfliers = False, order = x_axis_order, hue_order = window_hue_order, palette = "Blues" )     #This is the same as the graph above, but this has a defined x-axis order

# #GRAPH: Swarm Plot
# ax = sns.swarmplot( x = 'case_category', y = 'neoep_affinity', hue = 'bin_window', data = df_condition_neoep_mut, order = x_axis_order, hue_order = window_hue_order, size = dot_size )      #Use this if I want to show the dots

# # #GRAPH: Boxplot + Swarm Plot
# # sns.boxplot( x = 'case_category', y = 'neoep_affinity', hue = 'bin_window', data = df_condition_neoep_mut, showfliers = False, hue_order = window_hue_order, palette = "Blues" )     #this produces a rainbow of pastel colors by default
# # sns.swarmplot( x = 'case_category', y = 'neoep_affinity', hue = 'bin_window', data = df_condition_neoep_mut, hue_order = window_hue_order, color = ".25" )      #Use this if I want to show the dots
# # # ax = sns.boxplot( x = 'response_type', y = 'neoep_affinity', hue = 'bin_window', data = df_condition_neoep_mut, showfliers = False, hue_order = window_hue_order, palette = "Blues" )
# sns.despine()
# # pyplot.title( "Neoepitope Count Per irRECIST Condition: Unique Neoepitopes" )
# # pyplot.title( "Exp 12: Neoepitope Count for HLA Subgroup Per irRECIST Condition: " + annex_title, fontsize = 10 )
# # pyplot.title( graph_title )
# graph_title_final = "Processing and MHC Affinity: Veliparib & BRCA\n(" + graph_title_append + ")"
# pyplot.title( graph_title_final, fontsize = 30 )
# # pyplot.title( graph_title, fontsize = 10 )
# # pyplot.xlabel( "irRECIST" )
# pyplot.xlabel( "" )
# pyplot.ylabel( y_axis_label, fontsize = 20 )
# # pyplot.yticks(np.arange(0, 800, 100))           #Scale for: IC50
# # pyplot.yticks(np.arange(-3, 1, 0.5))          #Scale for: MHC affinity score

# # ax.legend().set_title( "Processing\nEfficiency" )     #This is how to change the legend title
# # pyplot.legend( fontsize = 14 )
# # # pyplot.legend( loc = "upper left", fontsize = 14 )

# ax.legend().set_title( "Processing\nEfficiency" )     #This is how to change the legend title
# # ax.legend().set_title( "Window" )     #This is how to change the legend title
# # pyplot.legend( fontsize = 14 )
# pyplot.setp(ax.get_legend().get_texts(), fontsize = '14') # for legend text

# pyplot.show()















##Experiment 8 - Count the number of neoepitopes per condition
print "------------ Experiment 8 ------------"
# input_file_date = "180411"      #for directory "/180403_Velip_V2"
# input_file_date = "180605"      #for directory "/180531_Velip_V3"
input_file_date = "180614"      #for directory "/180531_Velip_V3"


# #ACTUAL FILE - Files that DO NOT have MHC scores ranked for each mutation
# #Exome I dataset
# df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File1_NeoepCompare_V1.txt", sep = '\t' )
# #Exome II dataset
# df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File2_NeoepCompare_V1.txt", sep = '\t' )
# #Exome III dataset - pooled
# df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File3_NeoepCompare_V1.txt", sep = '\t' )
# #WGS I dataset
# df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File4_NeoepCompare_V1.txt", sep = '\t' )
# #WGS II dataset
# df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/180403_Thres0_File5_NeoepCompare_V1.txt", sep = '\t' )

#ACTUAL FILE - Files that do have MHC scores ranked for each mutation & remove genome position duplicates
#Exome I dataset
df_1 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File1_NeoepCompare_V1_Extended.txt", sep = '\t' )
#Exome II dataset
df_2 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File2_NeoepCompare_V1_Extended.txt", sep = '\t' )
#Exome III dataset - pooled
df_3 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File3_NeoepCompare_V1_Extended.txt", sep = '\t' )
#WGS I dataset
df_4 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File4_NeoepCompare_V1_Extended.txt", sep = '\t' )
#WGS II dataset
df_5 = pd.read_csv( DIR_RESULTS_FOLDER + "/" + input_file_date + "_Thres0_File5_NeoepCompare_V1_Extended.txt", sep = '\t' )


# #remove all genome mutation duplicates
df_1.drop_duplicates('genome_pos', inplace = True)
df_2.drop_duplicates('genome_pos', inplace = True) 
df_3.drop_duplicates('genome_pos', inplace = True) 
df_4.drop_duplicates('genome_pos', inplace = True) 
df_5.drop_duplicates('genome_pos', inplace = True)


##NOTE: I do not need this since I already labeled files with "_WGS" labels
# #need to append "_WGS" (or something) to each row in "case_id" of WGS file or else it will be confused with Exome I & II datasets
# df_4["case_id"] = df_4["case_id"].astype(str) + "_WGS"
# df_5["case_id"] = df_5["case_id"].astype(str) + "_WGS"

# df = pd.concat( [df_1, df_2, df_3, df_4, df_5] )
#NOTE: this excludes the pooled samples contained in "df_3"
df = pd.concat( [df_1, df_2, df_4, df_5] )      #this excludes the "pooled" dataset

##TEST::
print "df unique mutations = ", df['genome_pos'].nunique()

##TEST:: BEFORE
print "TEST:: length before rank 1 filter = ", len( df )
#filter for the highest MHC score for each mutation
df_mhc_max = df[ df['mut_neoep_rank_1'] == 1 ]
# df_mhc_max.drop_duplicates('genome_pos', inplace = True)      #I DO NOT NEED TO DROP MUTATION DUPLICATES SINCE I DID THIS EARLIER

print "len( df_mhc_max ) [should be the same as df unique mutations] = ", len( df_mhc_max )

##TEST:: AFTER
print "TEST:: length AFTER rank 1 filter = ", len( df_mhc_max )
print "TEST - df ALL sample names = ",  np.unique( df_mhc_max['case_id'] )

##TEST::
print "df_mhc_max.dtypes:\n", df_mhc_max.dtypes
print "df_mhc_max:\n", df_mhc_max


print "------------ Create filters for neoepitopes ------------"
# hash_filters = create_hash_filter( frame_stat, neoep_stat, bool_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
#DIFFERENCE: "mhc_affinity" is different between experiment 12C & 12D
mhc_affinity = 0        #looks at neoepitope affinities to HLA alleles (e.g. LOW, MID, HIGH), 2 = select high & mid affinity neoepitopes
frame_stat = 0
neoep_stat = 0      #looks at all neoeps, high/mid-affinity neoeps, or high-efficacy neoeps
check_NMD = True
thres_gene_exp = 0      #as of now, I do not have gene expression for Velip Data
thres_prot = 0
thres_tap = 0
thres_endogenous_freq = -1      #as of now, I do not have endogenous frequency for Velip Data
hash_filters = NeoepStatistics.create_hash_filter( mhc_affinity, frame_stat, neoep_stat, check_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq )
# obj_neoep_stat = NeoepStatistics( df, hash_filters )
obj_neoep_stat = NeoepStatistics( df_mhc_max, hash_filters )

print "------------ Apply filters for neoepitopes ------------"
#applied the filters to the dataset
df_filtered = obj_neoep_stat.filter_neoeps()
# filter_label = obj_neoep_stat.filter_correspond_labels( hash_filters )
df_filtered['ic50_fold_change'] = df_filtered.apply( lambda r: (float( r['ic50_score_1'] ) / r['ic50_score_2']), axis = 1 )

df_filtered.reset_index( inplace = True )

#append a case category for each experiment
items_pc_Velip = ['pcVelip1', 'pcVelip2', 'pcVelip3', 'pcVelip1_WGS', 'pcMV']
items_br_Velip = ['BRVelip1', 'BRVelip2', 'BRVelip3', 'BRVelip1_WGS', 'BMV']
items_pc_DMSO = ['pcDM1', 'pcDM2', 'pcDM3', 'pcDM1_WGS', 'pcMD']
items_br_DMSO = ['BRDM1', 'BRDM2', 'BRDM3', 'BRDM1_WGS', 'BMD']

hash_case_labels = {}
hash_case_labels.update( {x:"pcVelip" for x in items_pc_Velip} )
hash_case_labels.update( {x:"BRVelip" for x in items_br_Velip} )
hash_case_labels.update( {x:"pcDM" for x in items_pc_DMSO} )
hash_case_labels.update( {x:"BRDM" for x in items_br_DMSO} )

print "hash_case_labels = ", hash_case_labels

#label samples either Responders or Non-Responders based on irRECIST
# df_filtered['case_group'] = df_filtered.apply( lambda row: "Responder" if not 'disease' in row['irRECIST'].lower() else "Non-Responder", axis = 1 )
df_filtered['case_category'] = df_filtered.apply( lambda row: hash_case_labels[row['case_id']], axis = 1 )      #CONJ: I'm pretty sure "axis=1" means add new column, whereas "axis=0" means add new row


df_filtered.to_csv( DIR_RESULTS_FOLDER + "/180723_TEST_NEOEP_COUNT.txt", sep = '\t' )

#calculate the number of neoepitopes for each experiment
select_cat_clone = 2        #need to select either counting based on condition (1) or on clone (2)

if select_cat_clone == 1:
    #-------COUNT SET 1 - count the number of neoepitopes per condition (e.g. pcVelip, BRVelip, pcDM, etc.)
    #count total number of all neoepitopes per condition
    # series_1 = df_filtered.groupby( ['case_id', 'case_category'] )['genome_pos'].count()
    series_1 = df_filtered.groupby( 'case_category' )['genome_pos'].count()

    #count the total number of high-affinity neoepitopes
    df_filtered = df_filtered[ df_filtered['aff_cat_2'].isin( ["HIGH", "MID"] ) ]
    # series_2 = df_filtered.groupby( ['case_id', 'case_category'] )['genome_pos'].count()
    series_2 = df_filtered.groupby( 'case_category' )['genome_pos'].count()

    #count the number of high-efficacy neoepitopes
    thres_prot_ER_in = 50
    thres_TAP_ER_in = 50
    df_filtered = df_filtered[ 
    (df_filtered['proteasome_percentile_2'] >= thres_prot_ER_in) & 
    (df_filtered['tap_percentile_2'] >= thres_TAP_ER_in) ]
    # series_3 = df_filtered.groupby( ['case_id', 'case_category'] )['genome_pos'].count()
    series_3 = df_filtered.groupby( 'case_category' )['genome_pos'].count()

    #count the number of high-efficiacy neoepitopes 
    df_filtered = df_filtered[ df_filtered['gene_percentile'] > 0]
    # series_4 = df_filtered.groupby( ['case_id', 'case_category'] )['genome_pos'].count()
    series_4 = df_filtered.groupby( 'case_category' )['genome_pos'].count()

    print "COUNT BY CONDITION: count the number of neoepitopes per condition"
    print "CONDITION: # of neoepitopes - series_1:\n", series_1
    print "CONDITION: # of high-affinity neoepitopes - series_2:\n", series_2
    print "CONDITION: # of high-efficacy neoepitopes - series_3:\n", series_3
    print "CONDITION: # of high-efficacy neoepitopes that are expressed - series_4:\n", series_4

elif select_cat_clone == 2:
    #------COUNT SET 2 - count the number of neoepitopes per clone ("case_id" e.g. pcVelip, BRVelip, pcDM, etc.)
    #count total number of all neoepitopes per condition
    # series_1 = df_filtered.groupby( ['case_id', 'case_category'] )['genome_pos'].count()
    series_B1 = df_filtered.groupby( 'case_id' )['genome_pos'].count()

    #count the total number of high-affinity neoepitopes
    df_filtered = df_filtered[ df_filtered['aff_cat_2'].isin( ["HIGH", "MID"] ) ]
    # series_2 = df_filtered.groupby( ['case_id', 'case_category'] )['genome_pos'].count()
    series_B2 = df_filtered.groupby( 'case_id' )['genome_pos'].count()

    #count the number of high-efficacy neoepitopes
    thres_prot_ER_in = 50
    thres_TAP_ER_in = 50
    df_filtered = df_filtered[ 
    (df_filtered['proteasome_percentile_2'] >= thres_prot_ER_in) & 
    (df_filtered['tap_percentile_2'] >= thres_TAP_ER_in) ]
    # series_3 = df_filtered.groupby( ['case_id', 'case_category'] )['genome_pos'].count()
    series_B3 = df_filtered.groupby( 'case_id' )['genome_pos'].count()

    #count the number of high-efficiacy neoepitopes 
    df_filtered = df_filtered[ df_filtered['gene_percentile'] > 0]
    # series_4 = df_filtered.groupby( ['case_id', 'case_category'] )['genome_pos'].count()
    series_B4 = df_filtered.groupby( 'case_id' )['genome_pos'].count()

    print "\nCOUNT BY CLONES: count the number of neoepitopes per clone"
    print "CLONE: # of neoepitopes - series_1:\n", series_B1
    print "CLONE: # of high-affinity neoepitopes - series_2:\n", series_B2
    print "CLONE: # of high-efficacy neoepitopes - series_3:\n", series_B3
    print "CLONE: # of high-efficacy neoepitopes that are expressed - series_4:\n", series_B4


print "------------ Algorithm Completed: 180404_Velip_Graph_V1.py ------------"