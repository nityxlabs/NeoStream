#/usr/bin/python
import sys
import time

import scipy
from scipy import stats
import pandas as pd
import numpy as np

from cruzdb import Genome

# sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/180815_NeoStream/Algorithms" )
from SVSv7 import Isoform, EnsemblVEP, SimpleNeoepitopeAllV2, SimpleNeoepitopeIsoformV2, TranscribeTranscript, TranslateTranscript
from mokhaPy import mokhaPy

#Constants - directories
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
# DIR_CURR = DIR_PROJ + "/180203_GenomAltToNeoepV3"
DIR_CURR = DIR_PROJ + "/180815_NeoStream"
DIR_DATA = DIR_CURR + "/Data"
DIR_RESULTS = DIR_CURR + "/Results"

DIR_DATA_VELIP = DIR_PROJ + "/170304_NeoantigenExtract/Data/Velip"
#get mapped reads
DIR_RNASEQ = DIR_PROJ + "/150802_TophatSamples"
#Constants - gene expression
DIR_EXPRESSION = DIR_PROJ + "/161001_SJExonDiff/Results"
DIR_GENOME = DIR_PROJ + '/ArchiveData/hg19.fa'      #directory for samtool-indexed genome

#Constants - values for algorithm
SEC_BREAK = 2       #number of seconds to take break
AA_LEN = 9          #number of amino acids to retrieve around an altered AA


#Columns for Veliparib dataset
COL_GENE_SYM ='Hugo_Symbol'

#Columns for Anti-PD1 dataset
# COL_GENE_SYM ='gene_symbol'


def format_df_expression( df ):     #I haven't used this function in this script yet...
    """
    Args:
        df = Dataframe where each row is a gene
    Function: formats dataframe to remove certain elements
    """
    #remove '-' to 0 so the column can be converted to float
    df.replace('-', 0, inplace = True)

    #preserve gene expression genes
    # df = df[ df.index.str.contains('NM_') ]
    #OR remove all non-coding
    # df = df[ ~df.index.str.contains('NR_') ]

    #convert column to float
    for each_col in df.columns.values:      #each_col = each column
        df[[each_col]] = df[[each_col]].astype(float)

    ##TEST:: print "AFTER PROCESSING: len( df ) = ", len( df )

    return df

def format_df_v2( df ):
    """
    Args:
        df = Dataframe where each row is a gene
    Function: formats dataframe to remove certain elements - NOTE: this won't work because columns 'bool_NMD' & 'bool_NSD' do not have just "True", "False", but also "-"
    """
    df = df[ ~df['isoform_id'].str.contains('NR_') ]       #remove all non-coding instances

    #make columns into boolean
    df[['bool_NMD']] = df[['bool_NMD']].astype(bool)
    df[['bool_NSD']] = df[['bool_NSD']].astype(bool)

    return df

def format_df_v3( df ):
    """
    Args:
        df = Dataframe where each row is a gene
    Function: formats dataframe to remove certain elements - NOTE: this won't work because columns 'bool_NMD' & 'bool_NSD' do not have just "True", "False", but also "-"
    """
    #make columns into boolean
    df[['Chromosome']] = df[['Chromosome']].astype(str)

    return df


"""
Functions: evaluate gene expression
"""
def calculate_percentileofscore( list_gene_exp, curr_exp ):
    """
    Args:
        list_gene_exp = an array (list or numpy array) of gene expression values. As of now, this should not contain any duplicate values
        curr_exp = float value that is the gene expression value
    Function: determine the 
    """
    #Method 1
    return scipy.stats.percentileofscore( list_gene_exp, curr_exp, 'strict' )

def percentile_all_genes( df_express, pick_condition ):
    """
    Calculates the percentile expression for the sample

    Args:
        df_express = dataframe that reads in the expression 
        pick_condition = string that is a label for which condition to look into

    Returns:
        returns dataframe that contains information for genes, gene expression, and gene expression percentile with respect to all EXPRESSED genes (does not consider non-expressed genes)
    """
    # condition_label = EXP_COL_PC if condition_oi == 'pc' else EXP_COL_BR
    # EXP_COL_PC = 'pcVelip2'
    # EXP_COL_BR = 'BRVelip2'
    # #for "DM" condition - this is the control samples ("DM" means DMSO)
    # # EXP_COL_PC = 'pcDM2'
    # # EXP_COL_BR = 'BRDM2'

    ##CAUTION: This does not apply to new datasets (e.g. does not apply to 180403_Velip)
    if "pcvelip" in pick_condition.lower():
        condition_label = "pcVelip2"
    elif "brvelip" in pick_condition.lower():
        condition_label = "BRVelip2"
    elif "pcdm" in pick_condition.lower():
        condition_label = "pcDM2"
    elif "brdm" in pick_condition.lower():
        condition_label = "BRDM2"
    ##CAUTION: This is incorrect, but this is just extra just in case to handle extra datasets (esp. 180403_VelipData/ExomeIII)
    elif "pcmv" in pick_condition.lower():
        condition_label = "pcVelip2"
    elif "bmv" in pick_condition.lower():
        condition_label = "BRVelip2"
    elif "pcmd" in pick_condition.lower():
        condition_label = "pcDM2"
    elif "bmd" in pick_condition.lower():
        condition_label = "BRDM2"
    ##CAUTION: AGAIN, This is incorrect, but this is just extra just in case to handle extra datasets
    elif "p" in pick_condition.lower():
        condition_label = "pcVelip2"
    elif "b" in pick_condition.lower():
        condition_label = "BRVelip2"

    #extract list of expression levels that are above 0 (i.e. that are expressed)
    list_express = df_express[df_express[condition_label] > 0][condition_label].tolist()
    hash_pe = {}     #hash_pe = hash Percentile Expression

    total_rows = len( df_express )
    update = 0
    for i, (i_row, row) in enumerate( df_express.iterrows() ):
        ##TEST::
        # if i > 1000:
        #     break

        calc_percentile = calculate_percentileofscore( list_express, row[condition_label] )
        hash_pe[ row['GeneName'] ] = {'GeneName': row['GeneName'], 'express': row[condition_label], 'percentile': calc_percentile}

        update = mokhaPy.dataMeasure_rowsPeriodicUpdate( total_rows, i, update, 0.05 )

    #create dataframe for each
    df_percentile = pd.DataFrame( hash_pe.values(), index = hash_pe.keys(), columns = ['GeneName', 'express', 'percentile']  )

    return df_percentile


def mutation_genes_list( df_mut, df_percentile, thres_percentile, stat_mutate = 1 ):
    """
    counts the number of genes present with mutations

    Args:
        df_mut = pandas dataframe that contains mutations of interest
        df_percentile_pc = dataframe that is the gene expression for a specific condition (in this case, condition could be 'pc' or 'br')
        thres_percentile = float that is the threshold value
        stat_mutate = integer that is the type of defect I am interested in:
            -1 = find missense mutations
            -2 = find frameshift mutations
            -any other number = include all genomic alterations
    
    Returns:
        Dataframe that contains all mutations in genes that are highly-expressed

    """
    if stat_mutate == 1:        #select only missense mutations
        df_mut_missense = df_mut[ df_mut['Consequence'].str.contains('missense', case = False) ]
        df_oi = df_mut_missense
    elif stat_mutate == 2:      #select only alterations that cause frameshift
        df_mut_frameshift = df_mut[ df_mut['Consequence'].str.contains('frameshift', case = False) ]
        df_oi = df_mut_frameshift
    else:               #select all genomic alterations (mutation, indels / in-frame, frameshift)
        df_oi = df_mut

    #extract genes that are above expression
    genes_thres = df_percentile[ df_percentile['percentile'] >= thres_percentile ]['GeneName'].tolist()

    ##TEST::
    # label = 'missense' if stat_mutate == 1 else 'frameshift'
    # print "MGL: pc ", label," - df count = ", df_oi[ df_oi['Hugo_Symbol'].isin(genes_thres) ]['Hugo_Symbol'].count()
    # print df_oi[ df_oi['Hugo_Symbol'].isin(genes_thres) ]['Hugo_Symbol']
    # print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"

    return df_oi[ df_oi['Hugo_Symbol'].isin(genes_thres) ]

"""
FUNCTIONS: Dealing with general script stuff (REST API breaks, error messages)
"""

def showError( r_stat, row, i, header ):
    """
    creates a line for errors when the request doesn't go well
    Args:
        -r_stat = integer that is the status of the request, from class EnsemblVEP.get_vep_request()
        -row = information in file
        -i = integer of row
        -header = header for file, just counts the # of columns in file
    """
    subtract_col = 5
    line = str(i) + '\t'
    line+= str(row['Chromosome']) + '\t'
    line+= str(row['Start_Position']) + '\t'
    line+= str(row['End_Position']) + '\t'
    if r_stat == -1:
        line+= "TIMEOUT_ERROR\t"
    else:
        line+= "Error_no_info\t"

    # num_cols_to_fill = len( header.split('\t') ) - subtract_col
    num_cols_to_fill = 4
    for x in range(0, num_cols_to_fill):
        line+= "??\t"
    line+= "\n"

    return line

def take_request_break( i_counter, i, sec_break ):
    if i_counter % 5 == 0:
        print "Body: at row ", i, " - taking a ", sec_break, " second break (Apparently Ensembl REST API only allows 15 requests per second)"
        time.sleep( sec_break )

print "------------ Algorithm: 180825_ProcessGenomeAlt_V2.py ------------"
"""
Algorithm: Determine the consequence of a genomic alteration & records it in file. This also retrieves the nucleotides before & after the position of the genomic alteration

PROTOCOL: open file that contains genomic alterations -> calculate the consequence of each mutation by performing a VEP request -> record the results of the request (position, reading frame, relative CDS position, )
#LATER: retrieve the AA sequence before & after the genomic alteration
"""


start_time = time.time()

g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )


"""
mode_calc_exp -> will be used to determine if the gene expression should be calculated or not, where 0 = do not calculate gene expression, 1 = calculate gene expression.
"""
mode_calc_exp = 0     #use this if I do not want to consider thresholding by gene expression percentile. Also, this may include all occurrences of mutations, including non-coding? 
# mode_calc_exp = 1     #use this to consider thresholding by gene expression percentile


#retrieve user-inputted parameters
arg_date_output = sys.argv[1]
thres_express = int( sys.argv[2] )      #threshold for the gene expression percentile to accept
path_velip_file = sys.argv[3]
output_dir = sys.argv[4]
is_seq_WGS = sys.argv[5]


##MAY DELETE, I DO NOT NEED THIS BECAUSE OF "path_velip_file"
# # df_muts_orig = pd.read_csv( DIR_DATA_VELIP + "/" + input_velip_file, skiprows = 1, sep = '\t' )
# ##TEST DATA
# df_muts_orig = pd.read_csv( DIR_DATA_VELIP + "/180203_Test1.txt", skiprows = 1, sep = '\t' )
# ##ACTUAL DATA
# if "_PC" in arg_date_output:
#     df_muts_orig = pd.read_csv( DIR_DATA_VELIP + "/180203_Velip_PC.txt", skiprows = 1, sep = '\t' )
# elif "_DM" in arg_date_output:
#     df_muts_orig = pd.read_csv( DIR_DATA_VELIP + "/180203_Velip_DM.txt", skiprows = 1, sep = '\t' )
# else:
#     raw_input( "Error: press ctrl + C as there is no _PC or _DM in the file name" )

df_muts_orig = pd.read_csv( path_velip_file, skiprows = 1, sep = '\t' )


# ##NS MODIFY - NEED TO REMOVE THIS....
# #Filter 1: extract specific mutations, in this case "somatic" mutations - I want to exclude the "Germline" mutations as I want to see "Somatic mutations" arise from the treatment conditions (e.g. somatic mutations from Veliparib drug, pc vector that doesn't contain rescue BRCA1, br vector that does contain rescue BRCA1)
# df_muts_orig = df_muts_orig[ df_muts_orig['Mutation_Status'].str.contains('Somatic', case = False) ]



##TEST::
print "FILTER: mutations after somatic filter = ", len( df_muts_orig )


#Filter 2: conveys how often the mutated allele is observed betweent he tumor & the normal tissue. The lower the "n_alt_count", the less frequently the mutation is observed in normal and the more frequently it is observed in tumor
df_muts_orig = df_muts_orig[ df_muts_orig['n_alt_count'] < 3 ]



##TEST::
print "FILTER: mutations after 'n_alt_count' = ", len( df_muts_orig )



# # this only retrieves missense mutations
# df_muts_orig = df_muts_orig[ df_muts_orig['Consequence'].str.contains('missense', case = False) ]
#remove all mutations that do not cause an amino acid change
df_muts_orig['Amino_acids'].replace( '', np.nan, inplace = True )
df_muts_orig.dropna( subset = ['Amino_acids'], inplace = True )
df_muts_orig = df_muts_orig[ df_muts_orig['Amino_acids'] != '' ]        #I do not think I need this, this is a bit extra


#sort by sample ID (there are several conditions to consider)
df_muts_orig.sort_values( by = ['Tumor_Sample_Barcode'], inplace = True )



##METHOD 1 for analyzing data: retrieve genes that are above a specific threshold expression
if mode_calc_exp == 1:
    df_express = pd.read_csv( DIR_DATA_VELIP + "/170724_Velip_GeneCount.txt", sep = '\t', header = 0 )
    #calculate percentile expression for genes
    # df_percentile = percentile_all_genes( df_express, pick_condition )
    
    ##QUES: DO I NEED THIS "stat_mutate = 9"??
    #calculate genes that contain mutation that is above threshold
    # thres_express = 50
    # stat_mutate = 1         #stat_mutate -> selects the type of alteration (missense, frameshift, etc.)
    # stat_mutate = 9         #stat_mutate -> selects the type of alteration (missense, frameshift, etc.)
    # df_muts = mutation_genes_list( df_muts_orig, df_percentile, thres_express, stat_mutate )
    # count_mut_thres = df_mut_thres.count()        #for counting the total number of genes that contain a mutation that is above the threshold
##METHOD 2 for analyzing data: retrieve all entries
elif mode_calc_exp == 0:
    #extract the condition of interest
    # df_muts = df_muts_orig[ df_muts_orig['Tumor_Sample_Barcode'].str.contains(condition_oi, case = False) ]
    # df_muts = df_muts_orig
    pass

#format specific columns of file
df_muts = format_df_v3( df_muts_orig )


##TEST::
print "Show column types:"
print df_muts.dtypes





###----------OLD WAY TO CALCULATE--------

###-------DELETE ZONE 1 START-------

# arg_date_output = sys.argv[1]
# file_oi = int(sys.argv[2])          #selects between DM (0) & Exp (1 - this is using the Veliparib drug) 
# input_condition_type = int(sys.argv[3])
# thres_express = int( sys.argv[4] )      #threshold for the gene expression percentile to accept
# output_dir = sys.argv[5]
# file_seq_platform = sys.argv[6]


# #Possible combinations for "pick_condition" -> 'pc_velip', 'br_velip', 'pc_dm', 'br_dm'
# if input_condition_type == 0:
#     pick_condition = "pc_"
# elif input_condition_type == 1:
#     pick_condition = "br_"

# date_input = "170724"
# if file_oi == 0:
#     input_velip_file = date_input + "_Velip_" + file_seq_platform + "_DM.txt"
#     pick_condition += "dm"
# elif file_oi == 1:
#     input_velip_file = date_input + "_Velip_"  + file_seq_platform + "_Exp.txt"
#     pick_condition += "velip"

# condition_oi = pick_condition.split('_')[0]

# df_muts_orig = pd.read_csv( DIR_DATA_VELIP + "/" + input_velip_file, skiprows = 1, sep = '\t' )

# ##TEST::
# print "current file = ", pick_condition
# print "TOTAL: df_muts_orig length = ", len( df_muts_orig )

# #extract the condition of interest
# df_muts_orig = df_muts_orig[ df_muts_orig['Tumor_Sample_Barcode'].str.contains(condition_oi, case = False) ]
# #extract specific mutations, in this case "somatic" mutations - I want to exclude the "Germline" mutations as I want to see "Somatic mutations" arise from the treatment conditions (e.g. somatic mutations from Veliparib drug, pc vector that doesn't contain rescue BRCA1, br vector that does contain rescue BRCA1)
# df_muts_orig = df_muts_orig[ df_muts_orig['Mutation_Status'].str.contains('Somatic', case = False) ]
# # # this only retrieves missense mutations
# # df_muts_orig = df_muts_orig[ df_muts_orig['Consequence'].str.contains('missense', case = False) ]

# ##TEST::
# print "Filter 1: df_muts_orig length = ", len( df_muts_orig )

# ##METHOD 1 for analyzing data: retrieve genes that are above a specific threshold expression
# if mode_calc_exp == 1:
#     df_express = pd.read_csv( DIR_DATA_VELIP + "/" + date_input + "_Velip_GeneCount.txt", sep = '\t')
#     #calculate percentile expression for genes
#     df_percentile = percentile_all_genes( df_express, pick_condition )
    
#     #calculate genes that contain mutation that is above threshold
#     # thres_express = 50
#     # stat_mutate = 1         #stat_mutate -> selects the type of alteration (missense, frameshift, etc.)
#     stat_mutate = 9         #stat_mutate -> selects the type of alteration (missense, frameshift, etc.)
#     df_muts = mutation_genes_list( df_muts_orig, df_percentile, thres_express, stat_mutate )
#     # count_mut_thres = df_mut_thres.count()        #for counting the total number of genes that contain a mutation that is above the threshold
# ##METHOD 2 for analyzing data: retrieve all entries
# elif mode_calc_exp == 0:
#     #extract the condition of interest
#     # df_muts = df_muts_orig[ df_muts_orig['Tumor_Sample_Barcode'].str.contains(condition_oi, case = False) ]
#     df_muts = df_muts_orig



# #CHANGE for NonsilentSNV - only select events that change amino acid (may not all be single SNV but whatever)
# # df_muts = df_muts[ df_muts["HGVSp_Short"].str.contains("p") == True ]

# ##TEST::
# print "Filter 2: df_muts length = ", len( df_muts )

# ##TEST::
# print "columns in df_muts: ", df_muts.columns.values
# print "columns in df_mut_thres: ", df_muts.columns.values
# print "arg_date_output = ", arg_date_output

# #get the input file


###-------DELETE ZONE 1 ENDS-------

#write the file that will contain the output
file_write = open( output_dir + "/" + arg_date_output + "_ProcessGenomeAlt_V1.txt", 'w' )
# list_columns = []       #records all the headers that will be record with the 
header = ""

# #find first column that is a missense mutation
# [list_columns_gen, list_columns_transcript] = write_first_missense( df_muts )
# list_columns_start = ['row', 'chrom', 'start_pos', 'end_pos']
# header = '\t'.join( list_columns_start ) + '\t'
# header+= '\t'.join( list_columns_gen ) + '\t'
# header+= '\t'.join( list_columns_transcript ) + '\t'

# list_columns_start = ['row', 'sample_barcode', 'chrom', 'start_pos', 'end_pos']
list_columns_start = ['row', 'case_id', 'chrom', 'start_pos', 'end_pos']        #'case_id' is the same as 'sample_barcode'
list_columns_start+= ['Hugo_Symbol', 'Variant_Classification', 'Variant_Type']
list_columns_start+= ['file_hgvsc', 'file_hgvsp', 'file_short_hgvsp']
list_columns_start+= ['file_codon', 'file_amino_acid', 'file_protein_mrna_pos']
header = '\t'.join( list_columns_start ) + '\t'


##OLD HEADER VERSION
# header+= "Express_GeneName\tExpress_exp\tExpress_percentile\t"

# header+= "gene_symbol\tisoform_id\tgene_strand\t"
# header+= "allele_string\tnuc_orig\tnuc_alt\tcodon_change\t"
# header+= "in_frame\tvariant_class\tmy_change_type\t"
# header+= "alt_genome_range\tcds_start\tcds_end\t"
# header+= "alt_consequence\tbool_alteration\t"
# header+= "my_HGVSc\tmy_HGVSp\t"
# header+= "aa_change\taa_start\taa_end\t"
# header+= "codon_orig_alt\taa_orig_alt\t"
# header+= "myVEP_codon_orig_alt\tmyVEP_aa_orig_alt\t"
# header+= "mRNA_orig\tmRNA_alt\t"
# header+= "AA_orig_full\tAA_alt_full\tAA_orig_STOP\tAA_alt_STOP\t"
# header+= "possible_neoep_orig\tpossible_neoep_alt\tmatch_possible_neoep\t"

# #columns for double-checking columns
# header+= "compare_aa\tmatch_aa\t"
# header+= "compare_codon\tmatch_codon\t"
# header+= "compare_geneID\tmatch_geneID\t"
# header+= "compare_transcriptID\tmatch_transcriptID\t"

# #information about alteration
# header+= "exon_with_alt_start\texon_with_alt_end\t"
# header+= "penultimate_exon\talt_btwn_tp_penulti\t"
# #informatino about NMD & NSD
# header+= "bool_NMD\tbool_NSD\t"


##NEW HEADER VERSION
header+= "gene_symbol\tisoform_id\tgene_strand\t"
header+= "express_count\texpress_percentile\t"
header+= "allele_string\tnuc_orig\tnuc_alt\tcodon_change\t"
header+= "in_frame\tvariant_class\tmy_change_type\t"
header+= "alt_genome_range\tcds_start\tcds_end\t"
header+= "alt_consequence\tbool_alteration\t"
header+= "my_HGVSc\tmy_HGVSp\t"
header+= "aa_change\taa_start\taa_end\t"
header+= "codon_orig_alt\taa_orig_alt\t"
header+= "myVEP_codon_orig_alt\tmyVEP_aa_orig_alt\t"
header+= "mRNA_orig\tmRNA_alt\t"
header+= "AA_orig_full\tAA_alt_full\tAA_orig_STOP\tAA_alt_STOP\t"
header+= "possible_neoep_orig\tpossible_neoep_alt\tmatch_possible_neoep\t"

#columns for double-checking columns
header+= "compare_aa\tmatch_aa\t"
header+= "compare_codon\tmatch_codon\t"
header+= "compare_geneID\tmatch_geneID\t"
header+= "compare_transcriptID\tmatch_transcriptID\t"

#information about alteration
header+= "exon_with_alt_start\texon_with_alt_end\t"
header+= "penultimate_exon\talt_btwn_tp_penulti\t"
#informatino about NMD & NSD
header+= "bool_NMD\tbool_NSD\t"

#MISSING: can't calculate this just yet
# header+= "NMD_SENS_pos\tNMD_SENS_nuc_seq\tNMD_SENS_aa\t"
# header+= "NMD_IRREL_pos\tNMD_IRREL_nuc_seq\tNMD_IRREL_aa\t"
# header+= "bool_NMD\tbool_NSD\t"
#MISSING: can't calculate this just yet
# header+= "feature_start\tfeature_end\t"

header+= "\n"
file_write.write( header )


print "\n\n------------------------------\n\n"

##TEST::
print "header = ", header

print "\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n"


#keep track of transitioning to different sample types
current_case_id = None
counter_case_id = 1     #counts the number of case IDs the file has analyzed

total_rows = len( df_muts )
update = 0
i_counter = -1      #this counter is for determining when to take a break from making calls to VEP, also this is needed for the "retry_counter"
for i, (i_row, row) in enumerate( df_muts.iterrows() ):
    i_counter += 1

    # ##TEST::
    # if i > 30:
    #     break
    ##TEST:: to make the script faster
    # if i < 11:
    #     continue
    # if i > 12:
    #     break
    # if i < 2492:      #test row 2492 - isoform ENST00000366958, deletion position 1:214170534-214170535, nuc = AG/G
    #     continue
    # if i > 2492:
    #     break
    # if i < 88:
    #     continue

    #Filters
    #filter 1: remove any chromosome value that is not a number or
    if not row['Chromosome'].isdigit():
        if not row['Chromosome'] in ["X", "Y"]:
            continue

    #check if analyzing a new sample (i.e. new case ID) - this is a perfect place to determine gene expression for new sample
    if current_case_id != row['Tumor_Sample_Barcode']:
        current_case_id = row['Tumor_Sample_Barcode']
        counter_case_id += 1

        print "*****Patient ", counter_case_id, ": ", current_case_id


        ##QUES: I think I can delete all the non-commented stuff here...
        # #retrieve expression profile for this sample
        # hash_patient_exceptions = {'Pt16': 'Pt16.OnTx', 'Pt27': 'Pt27A.baseline'}
        # col_patient = current_case_id + '.baseline' if not current_case_id in hash_patient_exceptions else hash_patient_exceptions[current_case_id]
        
        # #calculate percentile expression for genes
        # if col_patient in df_express:
        #     df_percentile = percentile_all_genes( df_express, col_patient )
        # else:
        #     patient_no_express_data.append( row['Tumor_Sample_Barcode'] )

        #consider gene expression
        if mode_calc_exp == 1:
            df_percentile = percentile_all_genes( df_express, current_case_id )

    ##TEST::
    print "#######row ", i, " & i_row = ", i_row, " - start = ", row["Start_Position"], "|| end = ", row["End_Position"], " & gene = ", row["Hugo_Symbol"]

    #see if algorithm needs to take a break from making requests - only X requests allowed per second
    take_request_break( i_counter, i, SEC_BREAK )

    genome_pos = str(row['Chromosome']) + ':' + str(row['Start_Position']) + '-' + str(row['End_Position'])
    nuc_orig = row['Reference_Allele'] if row['Reference_Allele'] == "-" else None
    alt = row['Tumor_Seq_Allele2']
    strand = 1
    # opt_param = "variant_class=1&hgvs=1&refseq=1
    opt_param = "variant_class=1&hgvs=1"
    build_hg38 = False
    #determine if I want to only analyze the specific isoform of interest or not
    # curr_isoform_id = None
    curr_isoform_id = str(row["Transcript_ID"])

    ##TEST::
    print "&&&&&&&&row ", i, " & i_row = ", i_row, " || opt_param = ", opt_param


    obj_sna = SimpleNeoepitopeAllV2( genome_pos, nuc_orig, alt, DIR_GENOME, strand, opt_param, build_hg38, curr_isoform_id )

    retry_counter = 0
    if obj_sna.r_stat < 1:
        while retry_counter < 3 and obj_sna.r_stat < 1:
            obj_sna = SimpleNeoepitopeAllV2( genome_pos, nuc_orig, alt, DIR_GENOME, strand, opt_param )
            retry_counter += 1
            i_counter += 1
            print "ERROR at row ", i, ": Retry Ensembl VEP call - retry_counter = ", retry_counter, " & i_counter = ", i_counter
            #see if algorithm needs to take a break from making requests - only X requests allowed per second
            take_request_break( i_counter, i, SEC_BREAK ) 
        if obj_sna.r_stat < 1:
            error_line = showError( obj_sna.r_stat, row, i, header )
            print "STOP RETRY - ERROR at row ", i, ": FAILED!!! VEP call - retry_counter = ", retry_counter, " & i_counter = ", i_counter  
            file_write.write( error_line )
            #see if algorithm needs to take a break from making requests - only X requests allowed per second
            retry_counter = 0
            i_counter += 1
            take_request_break( i_counter, i, SEC_BREAK )
            continue

    ##TEST:: print "row ", i, ": list_info = ", list_info
    print "row ", i, " - genome_pos = ", genome_pos, " & alt = ", alt

    #this is the genomic position information for the genomic alteration
    sample_id = row['Tumor_Sample_Barcode']
    #append any information I may need at the end of sample_id (e.g. if Whole-Genome Seq, then "WGS")
    if "Y" in is_seq_WGS.upper():     #either "Y" for yes or "N" for no
        sample_id += "_WGS"

    line_prepend = str(i) + '\t'
    line_prepend+= sample_id + '\t'
    line_prepend+= str(row['Chromosome']) + '\t'
    line_prepend+= str(row['Start_Position']) + '\t'
    line_prepend+= str(row['End_Position']) + '\t'

    line_prepend+= str(row['Hugo_Symbol']) + '\t'
    line_prepend+= str(row['Variant_Classification']) + '\t'
    line_prepend+= str(row['Variant_Type']) + '\t'

    line_prepend+= str(row['HGVSc']) + '\t'
    line_prepend+= str(row['HGVSp']) + '\t'
    line_prepend+= str(row['HGVSp_Short']) + '\t'

    line_prepend+= str(row['Codons']) + '\t'
    line_prepend+= str(row['Amino_acids']) + '\t'
    line_prepend+= str(row['Protein_position']) + '\t'

    ##
    for i_iso, each_iso in enumerate(obj_sna.list_isoforms):        #each_iso = SimpleNeoepitopeIsoformV2 instance
        print "mRNA: ----------", i_iso, "----------"
        print "string = ", each_iso, " & len = ", len( each_iso.mRNA )

        mrna_orig = each_iso.mRNA
        mrna_alt = each_iso.mRNA_alt

        # [subseq_orig, subseq_alt] = each_iso.retrieve_mRNA_neoep_subseq( mrna_orig, mrna_alt, 1 )
        # # print "LEN subseq_org = ", len(subseq_orig)
        # # print "LEN subseq_alt = ", len(subseq_alt)
        # print "str subseq_org = ", ''.join(subseq_orig)
        # print "str subseq_alt = ", ''.join(subseq_alt)


        if mode_calc_exp == 1:
            find_express_count = df_percentile[df_percentile['GeneName'] == row[COL_GENE_SYM]]
            find_express_percentile = df_percentile[df_percentile['GeneName'] == row[COL_GENE_SYM]]
            gene_express_count = 0 if find_express_count.empty else find_express_count.iloc[0]['express']
            gene_express_percentile = 0 if find_express_percentile.empty else find_express_percentile.iloc[0]['percentile']
        else:
            gene_express_count = 0
            gene_express_percentile = 0

        #determine if alteration is in-frame or frameshifting
        nuc_change = each_iso.allele_string
        [count_change, in_frame, change_type] = TranscribeTranscript.determine_alt_rf( nuc_change )

        #determine change in amino acid sequence
        mrna_orig = each_iso.mRNA
        mrna_alt = each_iso.mRNA_alt

        #retrieve the codons & amino acids, original & altered, for each event
        [codon_orig, codon_alt] = each_iso.compare_codon_orig_alt( True )
        [aa_orig, aa_alt] = each_iso.compare_aa_orig_alt()
        str_combine_codon = str(codon_orig) + "/" + str(codon_alt)
        str_combine_aa = str(aa_orig) + "/" + str(aa_alt)

        #retrieve codon & amino acid change AGAIN (Yeah I know I know), but this time present it Ensembl VEP style
        [codon_orig_vep, codon_alt_vep] = each_iso.extract_changed_codon_self()
        [aa_orig_vep, aa_alt_vep] = each_iso.extract_changed_aa_self()
        str_combine_codon_vep = str(codon_orig_vep) + "/" + str(codon_alt_vep)
        str_combine_aa_vep = str(aa_orig_vep) + "/" + str(aa_alt_vep)

        #retrieve the nucleotide sequence that differs between both the original & mutated nucleotide sequence
        # #version 1: retrieve differing nucleotide sequence between original & altered
        # [subseq_mrna_orig, subseq_mrna_alt] = each_iso.retrieve_mRNA_neoep_subseq_v2( mrna_orig, mrna_alt, AA_LEN )
        #version 2: retrieve differing nucleotide sequence between original & altered
        [subseq_mrna_orig, subseq_mrna_alt] = each_iso.retrieve_orig_alt_neoeps_v3( mrna_orig, mrna_alt, AA_LEN )
        """
        For the test output:
        -aa_orig = from mrna_orig, the full amino acid sequence that is the original, non-mutated sequence
        -aa_alt = from mrna_alt, the full alterated amino acid sequence
        -aa_orig_2 = same as aa_orig, but only the subsequence before the first stop signal identified in "aa_alt_2"
        -aa_alt_2 = same as aa_alt, but only the subsequence before the first stop signal identified in "aa_alt_2"
        """

        ##TEST:: print "i:i_iso = ", i, ":", i_iso, " - mrna_orig = ", mrna_orig, " & mrna_alt = ", mrna_alt

        #version 1: to retrieve original & altered nucleotide seq / AA sequence
        # [aa_orig, aa_alt, aa_orig_2, aa_alt_2]  = each_iso.retrieve_comparative_neoeps( mrna_orig, mrna_alt, AA_LEN )
        #version 2: to retrieve original & altered nucleotide seq / AA sequence
        [aa_orig, aa_alt, aa_orig_2, aa_alt_2]  = each_iso.retrieve_comparative_neoeps_v2( mrna_orig, mrna_alt, AA_LEN )

        line_transcript = ""
        line_transcript+= each_iso.gene_symbol + "\t"
        line_transcript+= each_iso.isoform_id + "\t"
        line_transcript+= str( each_iso.gene_strand ) + "\t"

        line_transcript+= '%.4f' % gene_express_count + "\t"
        line_transcript+= '%.4f' % gene_express_percentile + "\t"

        line_transcript+= each_iso.allele_string + "\t"
        line_transcript+= each_iso.nuc_orig + "\t"
        line_transcript+= each_iso.nuc_alt + "\t"
        line_transcript+= each_iso.codon_change + "\t"

        hash_readframe_status = { True: "inframe", False: "frameshift" }
        hash_alt_type = { 0: "mutation", 1: "insertion", 2: "deletion" }
        line_transcript+= hash_readframe_status.get( in_frame, 'Unknown' ) + "\t"       #writes "Unknown" if in_frame is neither True nor False
        line_transcript+= str(each_iso.variant_class) + "\t"
        line_transcript+= hash_alt_type[change_type] + "\t"

        alt_genome_range = each_iso.chrom + ":" + str( each_iso.genome_start ) + "-" + str( each_iso.genome_end )
        line_transcript+= alt_genome_range + "\t"
        line_transcript+= str( each_iso.cds_start ) + "\t"
        line_transcript+= str( each_iso.cds_end ) + "\t"

        line_transcript+= str( each_iso.consequence ) + "\t"
        line_transcript+= str( each_iso.bool_alteration ) + "\t"

        line_transcript+= str( each_iso.hgvsc ) + "\t"
        line_transcript+= str( each_iso.hgvsp ) + "\t"

        line_transcript+= each_iso.aa_change + "\t"
        line_transcript+= str( each_iso.aa_start ) + "\t"
        line_transcript+= str( each_iso.aa_end ) + "\t"

        line_transcript+= str_combine_codon + "\t"
        line_transcript+= str_combine_aa + "\t"

        line_transcript+= str_combine_codon_vep + "\t"
        line_transcript+= str_combine_aa_vep + "\t"

        line_transcript+= ''.join( mrna_orig ) + "\t"
        line_transcript+= ''.join( mrna_alt ) + "\t"

        if len( aa_orig_2 ) == 0:
            num_potential_neoep_orig = 0
            num_potential_neoep_alt = 0
            potential_neoep_len = "-"
        else:
            num_potential_neoep_orig = len( aa_orig_2 ) - AA_LEN + 1
            num_potential_neoep_alt = len( aa_alt_2 ) - AA_LEN + 1
            potential_neoep_len = (num_potential_neoep_orig == num_potential_neoep_alt)
        line_transcript+= aa_orig + "\t"
        line_transcript+= aa_alt + "\t"
        line_transcript+= aa_orig_2 + "\t"
        line_transcript+= aa_alt_2 + "\t"
        line_transcript+= str( num_potential_neoep_orig ) + "\t"      
        line_transcript+= str( num_potential_neoep_alt ) + "\t"      
        line_transcript+= str( potential_neoep_len ) + "\t"      

        ##TEST::
        print "row ", i, " - str_combine_codon = ", str_combine_codon, " --> row[Codons] = ", row["Codons"], " & HGVSc = ", each_iso.hgvsc, " && HGVSp ", each_iso.hgvsp

        #VERSION 2: comparing codons & amino acids to file's version 
        line_transcript+= str_combine_codon_vep + " | " + str(row["Codons"]) + "\t"
        if pd.isnull( row["Codons"] ) or "None" in str_combine_codon_vep:
            line_transcript+= "-codon_empty-\t"
        else:
            line_transcript+= str( str_combine_codon_vep == row["Codons"] ) + "\t"    
        
        line_transcript+= str_combine_aa_vep + " | " + str(row["Amino_acids"]) + "\t"
        if pd.isnull( row["Amino_acids"] ) or "None" in str_combine_aa_vep:
            line_transcript+= "-aa_empty-\t"
        else:
            line_transcript+= str( str_combine_aa_vep == row["Amino_acids"] ) + "\t"
        
        #compare gene symbols & transcript IDs
        line_transcript+= each_iso.gene_symbol + " | " + str(row["Hugo_Symbol"]) + "\t"
        line_transcript+= str( each_iso.gene_symbol == row["Hugo_Symbol"] ) + "\t"
        line_transcript+= each_iso.isoform_id + " | " + str(row["Transcript_ID"]) + "\t"
        line_transcript+= str( each_iso.isoform_id == row["Transcript_ID"] ) + "\t"

        #find the 2nd to last exon
        line_transcript+= str( each_iso.alt_feat_start ) + "\t"
        line_transcript+= str( each_iso.alt_feat_end ) + "\t"

        if not each_iso.genome_penulti_three:
            penulti_feat = None
        else:
            pos_penult_pos = str(each_iso.chrom) + ":" + str(each_iso.genome_penulti_three) + "-" + str(each_iso.genome_penulti_three)
            penulti_feat = each_iso.find_pos_features( pos_penult_pos, 1)
        line_transcript+= str( penulti_feat ) + "\t"
        line_transcript+= str( each_iso.is_alt_between_tp_penulti_end() ) + "\t"

        line_transcript+= str( each_iso.determine_nmd_susceptible() ) + "\t"
        line_transcript+= str( each_iso.determine_nsd_susceptible() ) + "\t"

        #VERSION 1: comparing codons & amino acids to file's version 
        # line_transcript+= str_combine_codon + " | " + str(row["Codons"]) + "\t"
        # line_transcript+= str( str(str_combine_codon).upper() == str(row["Codons"]).upper() ) + "\t"
        # line_transcript+= str_combine_aa + " | " + str(row["Amino_acids"]) + "\t"
        # line_transcript+= str( str(str_combine_aa).upper() == str(row["Amino_acids"]).upper() ) + "\t"
        # line_transcript+= each_iso.gene_symbol + " | " + str(row["Hugo_Symbol"]) + "\t"
        # line_transcript+= str( each_iso.gene_symbol == row["Hugo_Symbol"] ) + "\t"
        # line_transcript+= each_iso.isoform_id + " | " + str(row["Transcript_ID"]) + "\t"
        # line_transcript+= str( each_iso.isoform_id == row["Transcript_ID"] ) + "\t"

        line = line_prepend + line_transcript + '\n'
        file_write.write( line )

    update = mokhaPy.dataMeasure_rowsPeriodicUpdate( total_rows, i, update, 0.0001 )

file_write.close()


elapse_time = time.time() - start_time
mokhaPy.timeElapse_convertToHMS( elapse_time )
print "------------ Algorithm Completed: 180825_ProcessGenomeAlt_V2.py ------------"