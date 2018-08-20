#/usr/bin/python
import os
import sys
import time

import scipy
from scipy import stats
import pandas as pd

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
# from SVSv5 import Isoform, MultiIsoform, IsoformSJ, TranslateTranscript
from mokhaPy import mokhaPy

#Constants - directories
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/170304_NeoantigenExtract"
DIR_DATA = DIR_CURR + "/Data"
DIR_DATA_VELIP = DIR_DATA + "/Velip"
DIR_RESULTS = DIR_CURR + "/Results"
# DIR_RESULTS_VELIP = DIR_RESULTS + "/170614_VelipNeoepitopes"
# DIR_RESULTS_VELIP = DIR_RESULTS + "/170712_VelipNeoepitopes"
DIR_RESULTS_VELIP = DIR_RESULTS + "/170724_VelipNeoepitopes"

#expression type - the gene expression profile for all genes - these are names of different conditions in date_Velip_GeneCount.txt (e.g. 170602_Velip_GeneCounts.txt & 170712_Velip_GeneCount.txt)
#for "Velip" condition
# EXP_COL_PC = 'pcVelip2'
# EXP_COL_BR = 'BRVelip2'
#for "DM" condition
EXP_COL_PC = 'pcDM2'
EXP_COL_BR = 'BRDM2'

def calculate_percentileofscore( list_gene_exp, curr_exp ):
    """
    Args:
        list_gene_exp = an array (list or numpy array) of gene expression values. As of now, this should not contain any duplicate values
        curr_exp = float value that is the gene expression value
    Function: determine the 
    """
    #Method 1
    return scipy.stats.percentileofscore( list_gene_exp, curr_exp, 'strict' )

def percentile_all_genes( df_express ):
    """
    Calculates the percentile expression for the sample
    """
    #calculate the percentile expression ]
    # list_express_pc = df_express[EXP_COL_PC].tolist()
    # list_express_br = df_express[EXP_COL_BR].tolist()
    list_express_pc = df_express[df_express[EXP_COL_PC] > 0][EXP_COL_PC].tolist()
    list_express_br = df_express[df_express[EXP_COL_BR] > 0][EXP_COL_BR].tolist()
    hash_pe_pc = {}     #hash_pe_pc = hash Percentile Expression for 'pc condition', where k = gene_symbol, & v = percentile
    hash_pe_br = {}     #hash_pe_br = hash Percentile Expression for 'pc condition'


    total_rows = len( df_express )
    update = 0
    for i, (i_row, row) in enumerate( df_express.iterrows() ):
        ##TEST::
        # if i > 100:
        #     break

        calc_percentile_pc = calculate_percentileofscore( list_express_pc, row[EXP_COL_PC] )
        hash_pe_pc[ row['GeneName'] ] = {'GeneName': row['GeneName'], 'express': row[EXP_COL_PC], 'percentile': calc_percentile_pc}

        calc_percentile_br = calculate_percentileofscore( list_express_br, row[EXP_COL_BR] )
        hash_pe_br[ row['GeneName'] ] = {'GeneName': row['GeneName'], 'express': row[EXP_COL_BR], 'percentile': calc_percentile_br}

        update = mokhaPy.dataMeasure_rowsPeriodicUpdate( total_rows, i, update, 0.002 )


    #create dataframe for each
    df_percentile_pc = pd.DataFrame( hash_pe_pc.values(), columns = ['GeneName', 'express', 'percentile']  )
    df_percentile_br = pd.DataFrame( hash_pe_br.values(), columns = ['GeneName', 'express', 'percentile']  )

    return [df_percentile_pc, df_percentile_br]

def mutation_genes_list( df_mut, df_percentile_pc, df_percentile_br, thres_percentile, stat_mutate = 1 ):
    """
    counts the number of genes present with mutations

    Args:
        df_mut = pandas dataframe that contains mutations of interest
        df_percentile_pc = dataframe that is the gene expression in condition 'pc' and the associated percentile
        df_percentile_br = dataframe that is the gene expression in condition 'pc' and the associated percentile
        thres_percentile = float that is the threshold value
        stat_mutate = integer that is the type of defect I am interested in:
            -1 = find missense mutations
            -2 = find frameshift mutations
    
    Returns:
        hash that is the # of genes that contains mutations 

    """
    if stat_mutate == 1:
        df_mut_missense = df_mut[ df_mut['Consequence'].str.contains('missense', case = False) ]
        df_oi = df_mut_missense
    else:
        df_mut_frameshift = df_mut[ df_mut['Consequence'].str.contains('frameshift', case = False) ]
        df_oi = df_mut_frameshift

    #separate the conditions 'pc' & 'br'
    df_pc = df_oi[ df_oi['Tumor_Sample_Barcode'].str.contains('pc', case = False) ]
    df_br = df_oi[ df_oi['Tumor_Sample_Barcode'].str.contains('BR', case = False) ]

    genes_pc = df_percentile_pc[ df_percentile_pc['percentile'] >= thres_percentile ]['GeneName'].tolist()
    genes_br = df_percentile_br[ df_percentile_br['percentile'] >= thres_percentile ]['GeneName'].tolist()

    ##TEST::
    label = 'missense' if stat_mutate == 1 else 'frameshift'
    print "MGL: pc ", label," - df count = ", df_oi[ df_oi['Hugo_Symbol'].isin(genes_pc) ]['Hugo_Symbol'].count()
    print df_oi[ df_oi['Hugo_Symbol'].isin(genes_pc) ]['Hugo_Symbol']
    print "MGL: br ", label," -  df count = ", df_oi[ df_oi['Hugo_Symbol'].isin(genes_br) ]['Hugo_Symbol'].count()
    print df_oi[ df_oi['Hugo_Symbol'].isin(genes_br) ]['Hugo_Symbol']
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"

    return { 'pc_count': df_pc[ df_pc['Hugo_Symbol'].isin(genes_pc) ]['Hugo_Symbol'].count(), 'br_count': df_br[ df_br['Hugo_Symbol'].isin(genes_br) ]['Hugo_Symbol'].count() }

def find_mutation_express( df_mut, df_express, stat_mutate = 1 ):
    """
    determines the frequency of mutations in each condition

    Args:
        df_mut = dataframe that is the file of the spreadsheet
        df_express = dataframe that is the expression for genes in specific conditions.
        stat_mutate = integer that is the type of defect I am interested in:
            -1 = find missense mutations
            -2 = find frameshift mutations

    Output:
    """
    #extract all missense mutations and frameshift events (frameshifts can be caused by insertions/deletions)
    if stat_mutate == 1:
        df_mut_missense = df_mut[ df_mut['Consequence'].str.contains('missense', case = False) ]
        df_oi = df_mut_missense
    else:
        df_mut_frameshift = df_mut[ df_mut['Consequence'].str.contains('frameshift', case = False) ]
        df_oi = df_mut_frameshift

    #count the total number of defects per event
    #NOTE: so in the experiment, a drug has been used to induce mutation for the sake of generating neoepitopes (I'll called it mutation-inducing agent. 'pc' is the "control" (hasn't been treated by mutation-induced agent), & 'BR' is the "control"
    df_pc = df_oi[ df_oi['Tumor_Sample_Barcode'].str.contains('pc', case = False) ]
    df_br = df_oi[ df_oi['Tumor_Sample_Barcode'].str.contains('BR', case = False) ]

    ##TEST::
    # print "len( df_pc ) = ", len( df_pc )
    # print "df_pc:"
    # print df_pc
    # print "len( df_br ) = ", len( df_br )
    # print "df_br:"
    # print df_br

    #determine which mutations are in expressed genes
    expressed_genes_pc = df_express[ df_express[EXP_COL_PC] > 0 ]['GeneName'].tolist()
    expressed_genes_br = df_express[ df_express[EXP_COL_BR] > 0 ]['GeneName'].tolist()
    df_pc_express = df_pc[ df_pc['Hugo_Symbol'].isin(expressed_genes_pc) ]
    df_br_express = df_br[ df_br['Hugo_Symbol'].isin(expressed_genes_br) ]

    #TEST:: calculate number of mutations associated
    # print "df_pc.count = ", df_pc['Variant_Type'].count()
    # print "df_pc.groupby('Variant_Classification').count = ", df_pc.groupby('Variant_Classification')['Variant_Classification'].count()
    # print "df_pc.groupby('Variant_Type').count = ", df_pc.groupby('Variant_Type')['Variant_Type'].count()

    ##TEST:: print genes that contain mutation
    label = 'missense' if stat_mutate == 1 else 'frameshift'
    print "ALL MUTATIONS - ", label
    print "df_pc all mutations: ", label
    print df_pc['Hugo_Symbol']
    print "df_br all mutations: ", label
    print df_br['Hugo_Symbol']
    print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n"


    return {'mut_total_pc': df_pc['Variant_Type'].count(), 'mut_total_br': df_br['Variant_Type'].count(), 'mut_expressed_pc': df_pc_express['Variant_Type'].count(), 'mut_expressed_br': df_br_express['Variant_Type'].count() }


# print "------------ Algorithm: 170712_Velip_SNV_Type_Express.py ------------"
print "------------ Algorithm: 180201_Velip_SNV_Type_Express.py --> I HAVE NOT MADE THE ADJUSTMENTS YET ------------"
"""
Algorithm: identify SNVs that can generate neoantigens for a sample
PROTOCOL: retrieve file that contains SNV & gene expression level -> determine the genes & isoforms the SNV belongs to -> retrieve the transcript, incorporate the mutation --> translate the sequence, retrieve N sequences before & after the mutation as well as the mutation itself
Notes on algorithm:
    -for each mutation, retrieve the mutation from the plus strand regardless if gene is on + or - strand. If on - strand, then will need to find the base complements
"""

start_time = time.time()

# path_file_express = DIR_DATA_VELIP + "/170712_Velip_GeneCount.txt"
path_file_express = DIR_DATA_VELIP + "/170724_Velip_GeneCount.txt"
df_express = pd.read_csv( path_file_express, sep = '\t')

# df_muts_exome = pd.read_csv( DIR_DATA_VELIP + "/170602_Velip_Exome.txt", skiprows = 1, sep = '\t' )      #whole exome sequencing
# df_muts_genome = pd.read_csv( DIR_DATA_VELIP + "/170712_Velip_WGS.txt", skiprows = 1,  sep = '\t' )      #whole genome sequencing
# path_file_oi = DIR_DATA_VELIP + "/170712_Velip_WGS.txt"
path_file_oi = DIR_DATA_VELIP + "/170724_Velip_WES_DM.txt"
df_muts_genome = pd.read_csv( path_file_oi, skiprows = 1,  sep = '\t' )      #whole genome sequencing



#extract only somatic mutations (exclude germline mutations)
# df_muts_exome = df_muts_exome[ df_muts_exome['Mutation_Status'].str.contains('Somatic', case = False) ]
df_muts_genome = df_muts_genome[ df_muts_genome['Mutation_Status'].str.contains('Somatic', case = False) ]


##TEST::
# print "df_express.index: "
# print df_express.index.tolist()
print "df_express['GeneName']: "
print df_express['GeneName']


# ##TEST::
# # print "columns:"
# # print df_mut.columns.values
# # print df_mut.dtypes

# #extract all missense mutations and frameshift events (frameshifts can be caused by insertions/deletions)
# df_mut_missense = df_mut[ df_mut['Consequence'].str.contains('missense', case = False) ]
# df_mut_frameshift = df_mut[ df_mut['Consequence'].str.contains('frameshift', case = False) ]

# #count the total number of defects per event
# #NOTE: so in the experiment, a drug has been used to induce mutation for the sake of generating neoepitopes (I'll called it mutation-inducing agent. 'pc' is the "control" (hasn't been treated by mutation-induced agent), & 'BR' is the "control"
# df_pc = df_mut_missense[ df_mut_missense['Tumor_Sample_Barcode'].str.contains('pc', case = False) ]
# df_br = df_mut_missense[ df_mut_missense['Tumor_Sample_Barcode'].str.contains('BR', case = False) ]

# for i, (i_row, row) in enumerate( df_pc.iterrows() ):
#     if i > 10:
#         break

#     print "PC ", i, ": i_row = ", i_row, " & row = ", row

# print "df_pc.count = ", df_pc['Variant_Type'].count()
# print "df_pc.groupby('Variant_Classification').count = ", df_pc.groupby('Variant_Classification')['Variant_Classification'].count()
# print "df_pc.groupby('Variant_Type').count = ", df_pc.groupby('Variant_Type')['Variant_Type'].count()

#calculate gene expression percentile
[df_percentile_pc, df_percentile_br] = percentile_all_genes( df_express )
##TEST::
# print "df_percentile_pc = ", df_percentile_pc
# print "df_percentile_br = ", df_percentile_br

#find missense mutations
# hash_express_exome_missense = find_mutation_express( df_muts_exome, df_express, 1 )
hash_express_genome_missense = find_mutation_express( df_muts_genome, df_express, 1 )
# print "hash_express_exome_missense = ", hash_express_exome_missense
print "hash_express_genome_missense = ", hash_express_genome_missense

#find frameshift mutations
# hash_express_exome_frameshift = find_mutation_express( df_muts_exome, df_express, 2 )
hash_express_genome_frameshift = find_mutation_express( df_muts_genome, df_express, 2 )
# print "hash_express_exome_frameshift = ", hash_express_exome_frameshift
print "hash_express_genome_frameshift = ", hash_express_genome_frameshift

#missense mutations
#get top 25% of expressed genes
# mut_missense_exome_25 = mutation_genes_list( df_muts_exome, df_percentile_pc, df_percentile_br, 25, 1 )
mut_missense_genome_25 = mutation_genes_list( df_muts_genome, df_percentile_pc, df_percentile_br, 25, 1 )
#get top 50% of expressed genes
# mut_missense_exome_50 = mutation_genes_list( df_muts_exome, df_percentile_pc, df_percentile_br, 50, 1 )
mut_missense_genome_50 = mutation_genes_list( df_muts_genome, df_percentile_pc, df_percentile_br, 50, 1 )
#get top 75% of expressed genes
# mut_missense_exome_75 = mutation_genes_list( df_muts_exome, df_percentile_pc, df_percentile_br, 75, 1 )
mut_missense_genome_75 = mutation_genes_list( df_muts_genome, df_percentile_pc, df_percentile_br, 75, 1 )

# print "mut_missense_exome_25 = ", mut_missense_exome_25
# print "mut_missense_exome_50 = ", mut_missense_exome_50
# print "mut_missense_exome_75 = ", mut_missense_exome_75

print "mut_missense_genome_25 = ", mut_missense_genome_25
print "mut_missense_genome_50 = ", mut_missense_genome_50
print "mut_missense_genome_75 = ", mut_missense_genome_75

#frameshift mutations
#get top 25% of expressed genes
# mut_frameshift_exome_25 = mutation_genes_list( df_muts_exome, df_percentile_pc, df_percentile_br, 25, 2 )
mut_frameshift_genome_25 = mutation_genes_list( df_muts_genome, df_percentile_pc, df_percentile_br, 25, 2 )
#get top 50% of expressed genes
# mut_frameshift_exome_50 = mutation_genes_list( df_muts_exome, df_percentile_pc, df_percentile_br, 50, 2 )
mut_frameshift_genome_50 = mutation_genes_list( df_muts_genome, df_percentile_pc, df_percentile_br, 50, 2 )
#get top 75% of expressed genes
# mut_frameshift_exome_75 = mutation_genes_list( df_muts_exome, df_percentile_pc, df_percentile_br, 75, 2 )
mut_frameshift_genome_75 = mutation_genes_list( df_muts_genome, df_percentile_pc, df_percentile_br, 75, 2 )

# print "mut_frameshift_exome_25 = ", mut_frameshift_exome_25
# print "mut_frameshift_exome_50 = ", mut_frameshift_exome_50
# print "mut_frameshift_exome_75 = ", mut_frameshift_exome_75

print "mut_frameshift_genome_25 = ", mut_frameshift_genome_25
print "mut_frameshift_genome_50 = ", mut_frameshift_genome_50
print "mut_frameshift_genome_75 = ", mut_frameshift_genome_75


rows_df = ['genome_pc', 'genome_br', 'exome_pc', 'exome_br']
columns_df = ['missense_total', 'missense_expressed', 'missense>=25_percentile', 'missense>=50_percentile', 'missense>=75_percentile']
columns_df+= ['frameshift_total', 'frameshift_expressed', 'frameshift>=25_percentile', 'frameshift>=50_percentile',  'frameshift>=75_percentile']


print "------------ Write Summary Table ------------"
df_write_file = pd.DataFrame( [], index = rows_df, columns = columns_df )


#write columns to files - 'genome_pc'
#missense information
df_write_file['missense_total']['genome_pc'] = hash_express_genome_missense['mut_total_pc']
df_write_file['missense_total']['genome_br'] = hash_express_genome_missense['mut_total_br']
df_write_file['missense_expressed']['genome_pc'] = hash_express_genome_missense['mut_expressed_pc']
df_write_file['missense_expressed']['genome_br'] = hash_express_genome_missense['mut_expressed_br']

df_write_file['frameshift_total']['genome_pc'] = hash_express_genome_frameshift['mut_total_pc']
df_write_file['frameshift_total']['genome_br'] = hash_express_genome_frameshift['mut_total_br']
df_write_file['frameshift_expressed']['genome_pc'] = hash_express_genome_frameshift['mut_expressed_pc']
df_write_file['frameshift_expressed']['genome_br'] = hash_express_genome_frameshift['mut_expressed_br']

#missense mutation - pc & br
df_write_file['missense>=25_percentile']['genome_pc'] = mut_missense_genome_25['pc_count']
df_write_file['missense>=50_percentile']['genome_pc'] = mut_missense_genome_50['pc_count']
df_write_file['missense>=75_percentile']['genome_pc'] = mut_missense_genome_75['pc_count']

df_write_file['missense>=25_percentile']['genome_br'] = mut_missense_genome_25['br_count']
df_write_file['missense>=50_percentile']['genome_br'] = mut_missense_genome_50['br_count']
df_write_file['missense>=75_percentile']['genome_br'] = mut_missense_genome_75['br_count']

#frameshift - pc & br
df_write_file['frameshift>=25_percentile']['genome_pc'] = mut_frameshift_genome_25['pc_count']
df_write_file['frameshift>=50_percentile']['genome_pc'] = mut_frameshift_genome_50['pc_count']
df_write_file['frameshift>=75_percentile']['genome_pc'] = mut_frameshift_genome_75['pc_count']

df_write_file['frameshift>=25_percentile']['genome_br'] = mut_frameshift_genome_25['br_count']
df_write_file['frameshift>=50_percentile']['genome_br'] = mut_frameshift_genome_50['br_count']
df_write_file['frameshift>=75_percentile']['genome_br'] = mut_frameshift_genome_75['br_count']

print "df_write_file:"
print df_write_file


df_write_file.to_csv( DIR_RESULTS_VELIP + "/170724_Velip_Mut_Summary.txt", sep = '\t' )

elapse_time = time.time() - start_time
mokhaPy.timeElapse_convertToHMS( elapse_time )
print "------------ Algorithm Completed: 180201_Velip_SNV_Type_Express.py ------------"