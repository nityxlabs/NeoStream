#/usr/bin/python
import os
import sys
import time
import pandas as pd

# sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
# sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/180815_NeoStream/Algorithms" )
sys.path.insert( 0, "../Algorithms" )
from SVSv7 import NeoepitopeMHC, MHC_IEDB_V2
from mokhaPy import mokhaPy


#Constants - directories
# DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = "../Algorithms"
DIR_DATA =  "../Data"
DIR_RESULTS = "../Results"

DIR_TAP_PROT_FILE = DIR_DATA + "/180106_Proteasome_TAP_Chart.txt"

SEC_BREAK = 5       #number of seconds to take break
#these are the constants for checking if a file has a header recorded for it
HEADFILE_RANK = True
HEADFILE_COMPARE = True


def neoep_rank_sort( df_neoep ):
    """
    In DataFrame of neoepitopes 'df_neoep', applies rank to different scores for antigen-processing, sorts the ranks, and then returns the dataframe
    Args:
        -df_neoep = DataFrame that contains all neoepitopes
    """

    #apply rankings to neoepitopes - neoepitope affinity, proteasomal processing, TAP transport, 
    df_neoep['rank_IC50'] = df_neoep['ic50_score'].rank( ascending = 1 )      #sort from least to greatest, where least gets rank 1
    #NOTE: mhc_score maybe a bit redundant because I am ranking IC50, but lets use it to double check
    df_neoep['rank_mhc'] = df_neoep['mhc_score'].rank( ascending = 0 )        #sort from greatest to least, where greatest gets rank 1
    df_neoep['rank_proteasome'] = df_neoep['proteasome_score'].rank( ascending = 0 )        #sort from greatest to least, where greatest gets rank 1
    df_neoep['rank_tap'] = df_neoep['tap_score'].rank( ascending = 0 )      #sort from greatest to least, where greatest gets rank 1
    #processing_score = for altered neoepitope, this sums the scores "proteasome_score + tap_score"
    df_neoep['rank_processing'] = df_neoep['processing_score'].rank( ascending = 0 )      #sort from greatest to least, where greatest gets rank 1
    #total_score = for altered neoepitope, this sums the scores "proteasome_score + tap_score + mhc_score"
    df_neoep['rank_total'] = df_neoep['total_score'].rank( ascending = 0 )      #sort from greatest to least, where greatest gets rank 1

    #create a "total sum of ranks" column
    col_list = ['rank_mhc', 'rank_proteasome', 'rank_tap']        #only chose these columns as 'rank_IC50' = 'rank_mhc', 'rank_processing' is a combination of 'rank_proteasome' & 'rank_tap', and 'rank_total' is a combination of 'rank_mhc', 'rank_proteasome', & 'rank_tap'
    df_neoep['sum_all_ranks'] = df_neoep[col_list].sum( axis = 1 )       #this sums for each column value for a given row (axis = 1)

    #sort by rankings
    df_neoep = df_neoep.sort_values( by = ['sum_all_ranks', 'rank_mhc', 'rank_proteasome', 'rank_tap', 'rank_processing', 'rank_total'], ascending = True )      #sort in descending order

    return df_neoep

def take_request_break( i_counter, i_row, sec_break ):
    if (i_counter % 1) == 0:
        print "Body: at row ", i, " & counter = ", i_counter, " - taking a ", sec_break, " second break (I don't know IEDB REST API's request limit.... - they suggested 1 request per second)"
        print "Body: at row ", i, " & counter = ", i_counter, " - taking a ", sec_break, " second break (Guideline for a limit to request information is 1 request per second)"
        time.sleep( sec_break )

def showError( i, row, neoep_seq_orig, neoep_seq_mut, r_stat ):
    """
    prints an error line for the file of interest
    """
    error_line = str( i ) + '\t'
    error_line+= row['gene_symbol'] + '\t'
    error_line+= row['isoform_id'] + '\t'
    error_line+= row['compare_codon'] + '\t'
    error_line+= row['compare_aa'] + '\t'
    error_line+= neoep_seq_orig + '\t'
    error_line+= neoep_seq_mut + '\t'
    error_line+= str( r_stat ) + '\n'

    return error_line


def calc_pep_processing( obj_mhc, neoep_seq, TEST_COUNTER, i ):
    """
    retrieve the analysis scores for peptide processing
    Args:
        -neoep_seq = string that is the peptide sequence
        -TEST_COUNTER = integer that is used for counting through the loop, MAINLY used for determining when to give IEDB web services a break
        -i = the integer that is the row
    """
    # #Method 1 to calculating peptide affinity to MHC - use both "MHC Only" analysis & combined proteasome + TAP + MHC analysis
    # [r_stat, df_pep] = obj_mhc.retrieve_epitope_analysis_all( neoep_seq )
    #Method 2 to calculating peptide affinity to MHC - only use combined proteasome + TAP + MHC analysis
    analysis_type = 2       #analysis_type = 2 -> perform analysis of peptide processing with proteasome, TAP, & MHC
    pep_match_freq = 0      #DOES NOT determine how frequently the peptide sequence occurs endogenously
    # pep_match_freq = 1      #determine how frequently the peptide sequence occurs endogenously
    calc_pt_perc = 3        #calculate percentile for TAP & proteasome scores
    [r_stat, df_pep] = obj_mhc.retrieve_epitope_analysis( neoep_seq, analysis_type, pep_match_freq, calc_pt_perc )


    retry_counter = 0
    if r_stat < 1:
        while retry_counter < 2 and r_stat < 1:
            #METHOD 1 to calculating peptide affinity to MHC:
            # [r_stat, df_pep] = obj_mhc.retrieve_epitope_analysis_all( neoep_seq )
            #METHOD 2 to calculating peptide affinity to MHC: 
            [r_stat, df_pep] = obj_mhc.retrieve_epitope_analysis( neoep_seq, analysis_type, pep_match_freq, calc_pt_perc )
            retry_counter += 1
            TEST_COUNTER += 1
            print "ERROR at row ", i, ": Retry Ensembl IEDB call - retry_counter = ", retry_counter, " & TEST_COUNTER = ", TEST_COUNTER
            #see if algorithm needs to take a break from making requests - only X requests allowed per second
            take_request_break( TEST_COUNTER, i, SEC_BREAK ) 
        if r_stat < 1:
            print "STOP RETRY - ERROR at row ", i, ": FAILED!!! IEDB call - retry_counter = ", retry_counter, " & TEST_COUNTER = ", TEST_COUNTER  
            error_line = showError( i, row, neoep_seq_orig, neoep_seq_mut, r_stat )
            file_error.write( error_line )
            #see if algorithm needs to take a break from making requests - only X requests allowed per second
            retry_counter = 0
            TEST_COUNTER += 1
            take_request_break( TEST_COUNTER, i, SEC_BREAK )
            df_pep = pd.DataFrame({'A' : []})       #create an empty dataframe

    return [df_pep, TEST_COUNTER, r_stat]


def write_to_file( df_pep, row, file_path, header_stat = False ):
    """
    writes the peptide information to the file
    """
    #append mutation information to MHC affinity
    genome_pos = "chr" + str(row['chrom']) + ":" + str(row['start_pos']) + "-" + str(row['end_pos'])

    df_pep["genome_pos"] = genome_pos
    # list_col = ["amino_acids", "gene_id", "transcript_id", "consequence_terms", "compare_aa", "match_aa", "compare_codon", "match_codon", "bool_NMD", "bool_NSD"]
    list_col = ["case_id", "cancer_type", "gene_symbol", "isoform_id", "express_count", "express_percentile", "my_change_type", "in_frame", "my_HGVSc", "my_HGVSp", "myVEP_codon_orig_alt", "myVEP_aa_orig_alt", "bool_NMD", "bool_NSD"]
    for each_col in list_col:
        if each_col in row:    
            df_pep[each_col] = row[each_col]
        else:
            df_pep[each_col] = "~_~"

    #record the genomic alteration and if it perturbs the reading frame or not
    def combine_alt( x ):
        return str( x['my_change_type'] ) + "_" + str( x['in_frame'] )
    df_pep['alt_readframe_type'] = df_pep.apply( combine_alt, axis = 1 )

    df_pep.to_csv( file_path, mode = 'a', header = header_stat, sep = '\t' )


def merge_compare_neoeps( df_orig, df_mut ):
    """
    A simple comparison betweeen peptide sequence 1 & peptide sequence 2 for antigen processing (proteasome, TAP, & MHC_1 binding)
    Args:
        -pep1 & pep2 = string that is amino acid sequence to be evaluated for affinity towards MHC alleles of interest. Usually pep1 is the original epitope & pep2 is the mutated epitope
        -analysis_type = integer with the following meaning
            -1 = perform analysis of peptide affinity to MHC
            -2 = perform analysis of peptide processing with proteasome, TAP, & MHC
    Returns: returns an array where [0] = integer that is status of request & [1] = a pandas dataframe that contains the peptide-binding affinity to a given proteasome & TAP. It has the following columns: allele, seq_num, start, end, length, peptide, ic50, percentile, rank
        -[0] values:
            - -1 = request error
            -0 = request not ok (r.ok is False)
            -1 = request is ok 
    """
    ic50 = "ic50_score"
    df_all = pd.merge( df_orig, df_mut, on = ['allele', 'start', 'end', 'length'], suffixes = ("_1", "_2") )
    #create a new column to compare the affinity between peptide 1 & peptide 2. delta_aff = change in affinity from pep1 to pep2, where a lower value in 'ic50' will be a stronger affinity, and a lower 'percentile' (percentile_rank) means a higher affinity for MHC allele

    ##TEST:: show columns
    print "df_orig.columns = ", df_orig.columns.values
    print "df_mut.columns = ", df_mut.columns.values
    print "df_orig.empty?? - ", df_orig.empty
    print "df_mut.empty?? - ", df_mut.empty
    print "TEST - merge_compare_neoeps: df_all.column values:"
    print df_all.columns.values
    print "df_all column datatypes:"
    print df_all.dtypes


    def calc_delta_aff( row, col_name, bool_inverse ):
        """
        determine if the new epitope '_2' (usually the mutated epitope) has a higher affinity to MHC allele than '_1' (usually the original epitope)
        Args:
            -row = each row from the dataframe
            -col_name = the prefix of the column name that will be compared
            -bool_inverse = boolean used for comparison between col_2 & col_1 (a comparison between numbers that are usually conveying some sort of affinity score). Note that '_2' is usually the after-effect (e.g. mutation) & '_1' is the original state (e.g. original nucleotide sequenced)
                -True = means if col_2 > col_1, then the after-state affinity is lower than the original state & if col_2 < col_1, then the after-state affinity is higher than the original.
                    -This is true for IC50 or a percentile rank where the lower the percentile the higher the affinity
                -False = means if col_2 > col_1, then the after-state affinity is higher than the original state & if col_2 < col_1, then the after-state affinity is higher than the original.
                    -This is true for -1*IC50 (e.g. -log(IC50))
        """
        col_1 = col_name + '_1'
        col_2 = col_name + '_2'
        if row[col_2] > row[col_1]:
            return "lower" if bool_inverse else "higher"
        elif row[col_2] < row[col_1]:
            return "higher" if bool_inverse else "lower"
        else:
            return "equal"

    #compare the scores between the original & altered peptide
    #bool_inverse = False for all of these because if col_2 > col_1 value (col_2 = altered peptide, col_1 = non-altered peptide), this means higher concentration (proteasome) or higher affinity (TAP & MHC).
    #NOTE: according to IEDB API (bottom of link: http://tools.iedb.org/processing/help/), the scores (proteasome, TAP, MHC) have been calculated so the higher the score, the higher the concentration (proteasome) or affinity (TAP & MHC, both use -log(IC50)). 
    #As they use MHC = -log( IC50 ): higher the peptide affinity for MHC allele -> higher the -log( IC50 ) -> the lower the IC50
    #As they use MHC = -log( IC50 ): lower the peptide affinity for MHC allele -> lower the -log( IC50 ) -> the higher the IC50
    df_all['delta_mhc_score'] = df_all.apply( calc_delta_aff, axis = 1, args = ('mhc_score',False,) )
    df_all['delta_proteasome_score'] = df_all.apply( calc_delta_aff, axis = 1, args = ('proteasome_score',False,) )
    df_all['delta_TAP_score'] = df_all.apply( calc_delta_aff, axis = 1, args = ('tap_score',False,) )

    return df_all

def find_last_index( df_getlastrow, df_neoep ):
    """
    If this program crashes during analysis, this function finds the last row of the list of genomic alterations in 'df_neoep' so it can continue analysis where it left off
    Args:
        -df_getlastrow = pandas dataframe, the file that is currently being written to. This file will be used to find the last genomic position record so I know where to restart (in this case this file is the list of neoepitopes)
        -df_neoep = DataFrame that is the list of genomic alterations
    """
    #retrieve the last row from the neoepitope file
    df_tail = df_getlastrow.tail(1)
    tail_genome_pos = df_tail.iloc[0]['genome_pos'].replace('chr', '')      #the format of the genomic position string is different between the neoepitope file and the mutation file (the neoepitope file has a 'chr' whereas the mutation file does not)

    ##TEST::
    # print "df_tail All = ", df_tail
    # print "df_tail genome_pos = ", df_tail['genome_pos']
    # print "df_tail genome_pos = ", tail_genome_pos
    # print "df_neoep.tail = ", df_neoep.tail(5).iloc[0]['alt_genome_range']

    #retrieve the first row that contains this genomic position from the mutation file
    df_gp = df_neoep[ df_neoep['alt_genome_range'] == tail_genome_pos ]       #df_gp = dataframe genomic position
    
    if not df_gp.empty:
        row_start = int( df_gp.index[0] )

        ##TEST::
        print "Method 1 - find by index!!!"
    else:       #else if genomic position not found, then just look for unique genomic position
        row_start = len( df_getlastrow['genome_pos'].unique() )

        ##TEST::
        print "Method 2 - unique genome_pos!!!"

    ##TEST::
    test_row = len( df_getlastrow['genome_pos'].unique() ) 
    print "TEST:: row_start = ", row_start, " & unique genome pos = ", test_row
    
    return row_start

def create_new_MHC_IEDB_V2( df_patient_mhc, patient, path_prot_tap ):
    """
    creates an MHC_IEDB_V2 instance based on patient alleles
    """
    info = df_patient_mhc[ df_patient_mhc['Patient ID'] == patient ].iloc[0]
    hla_a = ','.join( [('HLA-' + x) for x in info['HLA-A.1'].split(',')] )
    hla_b = ','.join( [('HLA-' + x) for x in info['HLA-B.1'].split(',')] )
    hla_c = ','.join( [('HLA-' + x) for x in info['HLA-C.1'].split(',')] )
    all_hla = hla_a + ',' + hla_b + ',' + hla_c
    pep_length = ','.join( [ '9' for i in range( 0, len(all_hla.split(',')) ) ] )

    ##TEST::
    print "all_hla = ", all_hla
    print "pep_length = ", pep_length

    #create an MHC IEDB instance
    # pred_method = "ann"
    pred_method = "netmhcpan"
    obj_mhc = MHC_IEDB_V2( all_hla, pep_length, pred_method, path_prot_tap )

    return obj_mhc

# print "------------ Algorithm: 180112_NeoantigenMHC_Eval_PD1_V2.py ------------"
print "------------ Algorithm: 180203_NeoantigenMHC_Eval_Velip_V1.py ------------"
"""
NOTE: this is the same as "180112_NeoantigenMHC_Eval_PD1_V2.py", just using this for Veliparib dataset
NOTE for 180112_NeoantigenMHC_Eval_PD1_V2.py: this algorithm is the same as "171205_NeoantigenMHC_Eval_PD1_V1.py", but uses MHC_IEDB_V2 instead of MHC_IEDB --> MHC_IEDB_V2 allows to calculate frequency of peptide sequence in endogenous proteins & percentiles for both proteasome & TAP
Algorithm: This algorithm will extract information neoepitope information
PROTOCOL:
    -open file that contains neoantigen sequences
    -extract each column with neoantigen sequences
    -will use 8 amino acids that are canonical, NMD-sensitive, and NMD-irrelevant
        -check to see if there is an early stop signal in NMD-sensitive or NMD-irrelevant --> if early stop signal in canonical then do not use
    -header information: sample name, sj ID, sj read count, SJ ratio value, gene symbol, isoform ID, gene expression, gene expression percentile, neoepitope window
        -neoepitope window = format 0:8, 1:9, 2:10, etc. -> size of window depends on constant NEOEPITOPE_LEN
"""

start_time = time.time()

date_input = sys.argv[1]
output_dir = sys.argv[2]
version_num_genome_alt = "V1"
version_num_neoep_calc = "V1"


#read-in file
df_neoep = pd.read_csv( output_dir + "/" + date_input + "_ProcessGenomeAlt_" + version_num_genome_alt + ".txt", sep = '\t', header = 0 )

# df_neoep = pd.read_csv( DIR_RESULTS + "/170830_ProcessGenomeAlt_V2_2.txt", sep = '\t' )

##TEST::
print "BEFORE LEN: df_neoep = ", len( df_neoep )

##TEST::
# print "df.index:"
# print df_neoep.index.tolist()
# print "column values:"
# print df_neoep.columns.values
print "column datatypes:"
print df_neoep.dtypes

#create files that will record information
##NOTE: other columns to consider - sample prevalence, control prevalence, read count, sj ratio, 
# file_fasta_write = open( DIR_RESULTS + "/170906_FASTANeoantigenSeq_V1.faa", "w" )
# path_neoep_compare = DIR_RESULTS + "/170910_NeoepMHC_V1.txt"

# #METHOD 1 for creating output: Just make fresh files each time
# path_neoep_compare = output_dir + "/" + date_input + "_NeoepCompare_" + version_num_neoep_calc + ".txt"
# file_neoep_compare = open( path_neoep_compare, 'w' ) 
# path_neoep_rank = output_dir + "/" + date_input + "_AllNeoepRanked_" + version_num_neoep_calc + ".txt"
# file_neoep_rank = open( path_neoep_rank, 'w' ) 

#METHOD 2 for creating output: As the program sometimes crashes or halts, need to pick up where I left off with the neoepitope analysis
path_neoep_compare = output_dir + "/" + date_input + "_NeoepCompare_" + version_num_neoep_calc + ".txt"
path_neoep_rank = output_dir + "/" + date_input + "_AllNeoepRanked_" + version_num_neoep_calc + ".txt"
path_error = output_dir + "/" + date_input + "_ERRORS_NeoepCompare_" + version_num_neoep_calc + ".txt"
if not os.path.exists( path_neoep_compare ) or os.path.getsize( path_neoep_compare ) == 0:        #create new files since the files do not already exist
    file_neoep_compare = open( path_neoep_compare, 'w' )
    file_neoep_rank = open( path_neoep_rank, 'w' )
    file_error = open( path_error, 'w' )
    row_start = 0

    ##TEST::
    print "DOES NOT FIND FILE: ", path_neoep_compare
else:       #else pick up where I left off with the peptide analysis
    file_neoep_compare = open( path_neoep_compare, 'a' )
    file_neoep_rank = open( path_neoep_rank, 'a' )
    file_error = open( path_error, 'a' )

    ##TEST::
    print "File does exist: ", path_neoep_compare
    
    #find the last row to continue from
    df_getlastrow = pd.read_csv( path_neoep_compare, sep = '\t', header = 0 )
    # df_getlastrow[['row']] = df_getlastrow[['row']].astype( int )
    # row_start = len( df_getlastrow['genome_pos'].unique() ) + 1      #retrieve the last row and continue from there

    row_start = find_last_index( df_getlastrow, df_neoep )
    row_start += 1      #add +1 because needs to start at next index of list of genomic alterations

    ##TEST:: just to see where the file left off
    print "row_start = ", row_start, " & BEFORE LEN: df_neoep = ", len( df_neoep ), " & percent = ", float( row_start ) / len( df_neoep )


#this will keep track of the # of errors
header_error = "row\tgene_id\ttranscript_id\t"
header_error+= "compare_codon\tcompare_aa\t"
header_error+= "neoep_orig\tneoep_alt\t"
header_error+= "req_stat\n"
file_error.write( header_error )

#may need to remove instances where NMD/NSD-susceptible and MAYBE ones that do not match the isoform IDs
# df_neoep = df_neoep[ (df_neoep['AA_orig_STOP'] != "") & (df_neoep['AA_alt_STOP'] != "") ]
df_neoep = df_neoep.dropna( subset = ['AA_orig_STOP', 'AA_alt_STOP'] )
# df_neoep = df_neoep[ (df_neoep['bool_NMD'] == True) & (df_neoep['bool_NSD'] == True) & (df_neoep['match_transcriptID'] == True) & (df_neoep['before_alt_aa'] != '-') & (df_neoep['after_alt_aa'] != '-') ]

##TEST::
print "AFTER LEN: df_neoep = ", len( df_neoep )

# #MHC-peptide prediction - NOT SURE WHICH MHC allele TO USE
# # mhc_alleles = "HLA-A*23:01,HLA-B*07:02,HLA-B*40:01,HLA-A*24:02"
# # length = "9,9,9,9"
# # mhc_alleles = "HLA-A*23:01"
# # length = "9"
# #these are the most common alleles in caucasians (only use for TCGA!), look at paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3459602/ ("HLA-A, B and DRB1 allele and haplotype frequencies in volunteer bone marrow donors from the north of Parana State"/Bardi, Jarduli)
# mhc_alleles = "HLA-A*01:01,HLA-B*35:01"
# length = "9,9"
# # pred_method = "ann"
# pred_method = "netmhcpan"
# obj_mhc = MHC_IEDB_V2( mhc_alleles, length, pred_method )

# df_patient_mhc = pd.read_csv( DIR_DATA + "/171204_PatientMHCAllele.txt", sep = '\t', header = 0 )

#MHC-peptide prediction
mhc_alleles = "HLA-A*23:01,HLA-B*07:02,HLA-B*40:01,HLA-A*24:02"
length = "9,9,9,9"
# mhc_alleles = "HLA-A*23:01"
# length = "9"
# pred_method = "ann"
pred_method = "netmhcpan"       #I think I'll only use "netmhcpan" for Veliparib dataset, not for TCGA... (I don't know why yet...)
obj_mhc = MHC_IEDB_V2( mhc_alleles, length, pred_method, DIR_TAP_PROT_FILE )


len_aa_radius = 9
obj_neoep = NeoepitopeMHC( len_aa_radius )
# window_size = len_aa_radius       #I DON'T THINK I NEED THIS VARIABLE?...
header_stat = True
TEST_COUNTER = 0

total_rows = len( df_neoep )
update = 0
current_case_id = None
for i, (i_row, row) in enumerate( df_neoep.iterrows() ):
    print "Loop i = ", i, " & i_row = ", i_row

    #use this especially if restarting analysis of file that has been halted or terminated early
    if i < row_start:
        continue

    #if new patient, then need to change the instance of the MHC alleles - use this only for Anti-PD1, not Veliparib (only use this line when the MHC alleles change for each patient)
    # if current_case_id != row['case_id']:
    #     obj_mhc = create_new_MHC_IEDB_V2( df_patient_mhc, row['case_id'], DIR_TAP_PROT_FILE )


    ##TEST::
    # if i > 3:
    #     break
    # if i < 58:
    #     continue

    ##TEST::
    # print "TEST i - ", i, " i_row = ", i_row, " & row\n"
    # print row
    # print "row.index = ", row.index.tolist()
    # print "dir( row ) = ", dir( row )

    #if row['bool_alteration'] is False, this means this alteration does not change the amino acid sequence or there is an error, else if True this means the genomic alteration has changed the amino acid change 
    if not row['bool_alteration']:
        continue

    #retrieve the unaltered peptide (orig) & altered peptide (mut) for analysis against IEDB algorithms
    neoep_seq_orig = row['AA_orig_STOP']
    neoep_seq_mut = row['AA_alt_STOP']

    
    ##TEST::
    print "TEST: neoep_seq_orig = ", neoep_seq_orig, " & obj_mhc = ", obj_mhc, " & TEST_COUNTER = ", TEST_COUNTER, " & i = ", i



    #analyze neoepitopes individually
    [df_orig, TEST_COUNTER, r_stat_orig] = calc_pep_processing( obj_mhc, neoep_seq_orig, TEST_COUNTER, i )
    if not df_orig.empty:
        df_orig['peptide_status'] = 'orig'
        write_to_file( df_orig.copy( deep = True ), row, path_neoep_rank, HEADFILE_RANK )
        if HEADFILE_RANK:
            HEADFILE_RANK = False
    [df_mut, TEST_COUNTER, r_stat_mut] = calc_pep_processing( obj_mhc, neoep_seq_mut, TEST_COUNTER, i )
    if not df_mut.empty:
        df_mut['peptide_status'] = 'alt'
        write_to_file( df_mut.copy( deep = True ), row, path_neoep_rank, HEADFILE_RANK )
        if HEADFILE_RANK:
            HEADFILE_RANK = False

    #need to compare the original and altered peptide
    if df_orig.empty or df_mut.empty:
        error_line = showError( i, row, neoep_seq_orig, neoep_seq_mut, r_stat_orig )
        file_error.write( error_line )
        continue

    #determine the comparisons between the original & mutated neoepitope
    # df_neoep_compare = obj_mhc.compare_epitope_affinity( neoep_seq_orig, neoep_seq_mut )
    # df_neoep_compare = obj_mhc.compare_epitope_processing( neoep_seq_orig, neoep_seq_mut, 1 )
    # [r_stat, df_neoep_compare] = obj_mhc.compare_epitope_processing_all( neoep_seq_orig, neoep_seq_mut )
    # df_neoep_compare = pd.merge( df_orig, df_mut, on = ['allele', 'start', 'end', 'length'], suffixes = ("_1", "_2") )
    df_neoep_compare = merge_compare_neoeps( df_orig, df_mut )

    write_to_file( df_neoep_compare, row, path_neoep_compare, HEADFILE_COMPARE )
    if HEADFILE_COMPARE:
        HEADFILE_COMPARE = False
    
    # #append mutation information to MHC affinity
    # genome_pos = "chr" + str(row['chrom']) + ":" + str(row['start_pos']) + "-" + str(row['end_pos'])
    # df_neoep_compare["genome_pos"] = genome_pos
    # # list_col = ["amino_acids", "gene_id", "transcript_id", "consequence_terms", "compare_aa", "match_aa", "compare_codon", "match_codon", "bool_NMD", "bool_NSD"]
    # list_col = ["gene_symbol", "isoform_id", "Express_exp", "Express_percentile", "my_change_type", "in_frame", "my_HGVSc", "my_HGVSp", "myVEP_codon_orig_alt", "myVEP_aa_orig_alt", "bool_NMD", "bool_NSD"]
    # for each_col in list_col:
    #     if each_col in row:    
    #         df_neoep_compare[each_col] = row[each_col]
    #     else:
    #         df_neoep_compare[each_col] = "~_~"
    # df_neoep_compare.to_csv( path_neoep_compare, mode = 'a', header = HEADFILE_RANK, sep = '\t' )
    # if HEADFILE_RANK:
    #     HEADFILE_RANK = False
    
    TEST_COUNTER += 1
    take_request_break( TEST_COUNTER, i, SEC_BREAK )

    ##TEST:: just keep track of time
    print "Timestamp: ", time.strftime("%c")
    

    update = mokhaPy.dataMeasure_rowsPeriodicUpdate( total_rows, i, update, 0.0001 )


file_neoep_compare.close()
file_neoep_rank.close()


#need to add rank to file "file_neoep_rank"
df_rank = pd.read_csv( path_neoep_rank, header = 0, sep = '\t' )
df_rank = neoep_rank_sort( df_rank )
df_rank.to_csv( path_neoep_rank , mode = 'w', header = True, sep = '\t' )


elapse_time = time.time() - start_time
mokhaPy.timeElapse_convertToHMS( elapse_time )

# print "------------ Algorithm Completed: 180112_NeoantigenMHC_Eval_PD1_V2.py ------------"
print "------------ Algorithm Completed: 180203_NeoantigenMHC_Eval_Velip_V1.py ------------"