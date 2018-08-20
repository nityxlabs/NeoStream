# #/usr/bin/python
# import sys
# import time

# import scipy
# from scipy import stats
# import pandas as pd
# import numpy as np

# from cruzdb import Genome

# sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
# from SVSv7 import Isoform, EnsemblVEP, SimpleNeoepitopeAllV2, SimpleNeoepitopeIsoformV2, TranscribeTranscript, TranslateTranscript
# from mokhaPy import mokhaPy


#/usr/bin/python
import sys
import time

import pandas as pd


#Constants - directories
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/180203_GenomAltToNeoepV3"
DIR_DATA = DIR_CURR + "/Data"
DIR_RESULTS = DIR_CURR + "/Results"
DIR_RESULTS_FOLDER = DIR_RESULTS + "/180531_Velip_V3"

DIR_DATA_VELIP = DIR_PROJ + "/170304_NeoantigenExtract/Data/Velip"
#get mapped reads
DIR_RNASEQ = DIR_PROJ + "/150802_TophatSamples"
#Constants - gene expression
DIR_EXPRESSION = DIR_PROJ + "/161001_SJExonDiff/Results"
DIR_GENOME = DIR_PROJ + '/ArchiveData/hg19.fa'      #directory for samtool-indexed genome

#Constants - values for algorithm
SEC_BREAK = 2       #number of seconds to take break
AA_LEN = 9          #number of amino acids to retrieve around an altered AA


print "------------ Algorithm: 180611_FilterNeoepSweasy_V1.py ------------"
"""
This filters neoepitopes and rearranges them for Joann Sweasy's specifications

PROTOCOL:
-retrieve file using pandas dataframe
-filter using the following
    -bool_NMD
    -make sure the original & altered sequence are not the same
        -it would be helpful to have the endogenous frequency recorded for each peptide sequence
-filter neoepitopes rank for each mutation
-
"""
input_files = []
date_input = "180614"
date_output = "180625"

# list_filenum = range(1,6)       #include all files
list_filenum = [1,2,4,5]        #exclude the pooled files

for num_file in list_filenum:
    input_files.append( date_input + "_Thres0_File" + str(num_file) + "_NeoepCompare_V1_Extended.txt" )

#hash_correspond_col = hash table where the "key" is the column name of the outputted file & "value" is the column for the inputted file
# hash_correspond_col = {"case_id": "case_id",
#     "gene_symbol": "gene_symbol",
#     "allele": "allele",
#     "genome_pos": "genome_pos",
#     "myVEP_aa_orig_alt": "myVEP_aa_orig_alt",
#     "my_change_type": "my_change_type",
#     "in_frame": "in_frame",
#     "alt_readframe_type": "alt_readframe_type",
#     "original_peptide": "peptide_1",
#     "ic50_score_orig": "ic50_score_1",
#     "aff_cat_orig": "aff_cat_1",
#     "altered_peptide": "peptide_2",
#     "ic50_score_alt": "ic50_score_2",
#     "aff_cat_alt": "aff_cat_2",
#     "bool_NMD": "bool_NMD",
#     "Fold Change (orig/alt)": "fold_change" }

list_hash_col_order = [("case_id", "case_id"),
    ("gene_symbol", "gene_symbol"),
    ("allele", "allele"),
    ("genome_pos", "genome_pos"),
    ("myVEP_aa_orig_alt", "myVEP_aa_orig_alt"),
    ("my_change_type", "my_change_type"),
    ("in_frame", "in_frame"),
    ("alt_readframe_type", "alt_readframe_type"),
    ("original_peptide", "peptide_1"),
    ("ic50_score_orig", "ic50_score_1"),
    ("aff_cat_orig", "aff_cat_1"),
    ("proteasome_score_orig", "proteasome_score_1"),
    ("tap_score_orig", "tap_score_1"),
    ("proteasome_percentile_orig", "proteasome_percentile_1"),
    ("tap_percentile_orig", "tap_percentile_1"),
    ("altered_peptide", "peptide_2"),
    ("ic50_score_alt", "ic50_score_2"),
    ("aff_cat_alt", "aff_cat_2"),
    ("proteasome_score_alt", "proteasome_score_2"),
    ("tap_score_alt", "tap_score_2"),
    ("proteasome_percentile_alt", "proteasome_percentile_2"),
    ("tap_percentile_alt", "tap_percentile_2"),
    ("bool_NMD", "bool_NMD"),
    ("Fold Change (orig/alt)", "fold_change"),
    ("mut_neoep_rank_1", "mut_neoep_rank_1")]
hash_correspond_col = {x[0]:x[1] for x in list_hash_col_order}


list_columns = ["sample_row"] + [x[0] for x in list_hash_col_order]
# list_columns = hash_correspond_col.values()
# list_columns.append( "row_index" )
df_neoep_filtered = pd.DataFrame( data=[], columns=list_columns )


##TEST::
print "df_neoep_filtered = ", df_neoep_filtered

for i_file, each_file in enumerate(input_files):

    df_neoep = pd.read_csv( DIR_RESULTS_FOLDER + "/" + each_file, sep = "\t", header = 0 )


    ##TEST::
    print "columns:\n", df_neoep.columns.values
    print "neoep dtypes:\n", df_neoep.dtypes

    #filter each neoepitope based on specific parameters
    #Filter based on NMD
    # df_neoep = df_neoep[ ~df_neoep['bool_NMD'].str.contains("tru", case=False) ]        #remove all transcripts susceptible to NMD degradation
    df_neoep = df_neoep[ df_neoep['bool_NMD'].str.contains("fal", case=False) ]        #remove all transcripts susceptible to NMD degradation
    #only record the neoepitope with the highest MHC affinity for each mutation (rank = 1)
    df_neoep = df_neoep[ df_neoep['mut_neoep_rank_1'] == 1 ]        #select the neoepitopes with the highest MHC affinity for a given mutation (NOTE: this rank is for altered peptides)
    df_neoep = df_neoep[ df_neoep['peptide_1'] != df_neoep['peptide_2'] ]   #filter neoepitopes that match
    #remove any instances where the AA change is "None"
    df_neoep = df_neoep[ ~df_neoep['myVEP_aa_orig_alt'].str.contains("non", case=False) ]

    #calculate the MHC affinity fold change before & after 
    df_neoep['fold_change'] = df_neoep.apply( lambda row: float(row['ic50_score_1']) / row['ic50_score_2'], axis=1 )


    #go through each row in the file & record the rows of interest
    for i, (i_row, row) in enumerate( df_neoep.iterrows() ):

        ##TEST::
        print "\tTEST row ", i, " - ", row

        hash_extract_info = {}      #this records the necessary info from list of neoepitopes
        for k,v in hash_correspond_col.iteritems():
            hash_extract_info[k] = row[v]

        sample_row = str(i_file) + "_" + str(i)
        hash_extract_info['sample_row'] = sample_row

        ##TEST::
        print "\tTEST row ", i, " - ", hash_extract_info, " & len(df_neoep_filtered) = ", len(df_neoep_filtered)


        df_neoep_filtered = df_neoep_filtered.append( [hash_extract_info], ignore_index=True )



#write output of filtered neoepitopes
df_neoep_filtered.to_csv( DIR_RESULTS_FOLDER + "/" + date_output + "_NeoepCompare_Rank1_V1.txt", sep="\t", index_label="row", index=True)

#write output of neoepitopes
df_neoep_filtered_hiaff = df_neoep_filtered[ ~df_neoep_filtered['aff_cat_alt'].str.contains("none", case=False) ]
df_neoep_filtered_hiaff.to_csv( DIR_RESULTS_FOLDER + "/" + date_output + "_NeoepCompare_Rank1_HiAff_V1.txt", sep="\t", index_label="row", index=True )


##TEST::
print "list_columns = ", list_columns

print "------------ Algorithm Completed: 180611_FilterNeoepSweasy_V1.py ------------"