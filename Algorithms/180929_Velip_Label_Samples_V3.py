#/usr/bin/python
import os
import sys
import time

import pandas as pd

import scipy
from scipy import stats

# sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
sys.path.insert( 0, "../Algorithms" )
from mokhaPy import mokhaPy

#Constants - directories
# DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
# DIR_CURR = DIR_PROJ + "/180203_GenomAltToNeoepV3"
DIR_CURR = "../Algorithms"
DIR_DATA = "../Data"
DIR_RESULTS = "../Results"


#Constants 
HEADFILE_STAT = True

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

    if df_express == None:
        return [None, None]

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

def retrieve_condition_id_express( str_cond ):
    """
    METHOD: this retrieves the condition ID corresponding the gene expression file
    NOTE: I think 'bmd', 'bmv', 'pcmd', & 'pcmv' are from the pooled dataset
    """
    if any( x in str_cond.lower() for x in ['brdm', 'bmd', 'bmv'] ):
        return "BRDM_avg"
    elif 'brvelip' in str_cond.lower():
        return "BRVelip_avg"
    elif any( x in str_cond.lower() for x in ['pcdm', 'pcmd', 'pcmv'] ):
        return "pcDM_avg"
    elif 'pcvelip' in str_cond.lower():
        return "pcVelip_avg"
    else:
        return None



# print "------------ Algorithm: 180224_AntiPD1_Label_Samples_V2.py ------------"
# print "------------ Algorithm: 180325_AntiPD1_Label_Samples_V3.py ------------"
# print "------------ Algorithm: 180411_Velip_Label_Samples_V1.py ------------"
# print "------------ Algorithm: 180815_Velip_Label_Samples_V1.py ------------"
# print "------------ Algorithm: 180825_Velip_Label_Samples_V2.py ------------"
print "------------ Algorithm: 180929_Velip_Label_Samples_V3.py ------------"
"""
This algorithm will label patients with the correct response type - ["Progressive Disease", "Partial Response", "Complete Response"]
-NOTE: This algorithm is the same as 180825_Velip_Label_Samples_V2.py but is designed to be more general as to be applied to a dataset more generally (e.g. no gene expression file)
-NOTE: This is the same as 180224_AntiPD1_Label_Samples_V2.py, but also assigns a specific neoepitope to mutation so I only count 1 neoepitope per mutation. I assign neoepitope with highest predicted HLA affinity to given mutatin
-NOTE: this algorithm is similar to "171209_NeoantigenMHC_AppendInfo.py" in terms of appending new information to a file, however I like the approach applied to this file better.
"""

start_time = time.time()

##BE OPEN TO ACCEPTING GENE EXPRESSION DATA
# path_file_express = DIR_DATA_VELIP + "/170712_Velip_GeneCount.txt"
# path_file_express = DIR_DATA_VELIP + "/170724_Velip_GeneCount.txt"
# path_file_express = DIR_DATA_VELIP + "/180711_Velip_GeneCount.txt"
# df_express = pd.read_csv( path_file_express, sep = '\t')
# df_express = df_express.set_index('GeneName')     #perhaps hold off on make the gene symbol the index

df_express = None               #this is the default for gene expression as of now as I am not expecting gene expression

input_filename = sys.argv[1]
output_dir = sys.argv[2]

for ( dirpath, dirnames, filenames ) in os.walk( DIR_RESULTS ):
    break
#filter files of interest
all_files_oi = [x for x in filenames if input_filename in x and "_NeoepCompare_" in x]


for file_num in all_files_oi:
    HEADFILE_STAT = True        #NOTE: Only place "HEADFILE_STAT" here for loops

    # path_write = DIR_RESULT_ANTIPD1 + "/180130_Thres0_SampAllExpress_NeoepCompare_irRECIST_V1.txt"
    # path_write = DIR_RESULT_ANTIPD1 + "/180224_Thres0_SampAllExpress_NeoepCompare_irRECIST_V2.txt"
    # path_write = DIR_RESULT_ANTIPD1 + "/180225_Thres0_SampAllExpress_NeoepCompare_irRECIST_V2.txt"
    # path_write = DIR_RESULT_ANTIPD1 + "/180325_Thres0_SampAllExpress_NeoepCompare_irRECIST_V3.txt"
    # path_write = DIR_RESULTS_FOLDER + "/180411_Thres0_File" + str(file_num) + "_NeoepCompare_V1_Extended.txt"       #NOTE: for directory "/180403_Velip_V2"
    # path_write = DIR_RESULTS_FOLDER + "/" + input_filename + "_Thres0_File" + str(file_num) + "_NeoepCompare_V1_Extended.txt"       #NOTE: for directory "/180531_Velip_V3"
    path_write = output_dir + "/" + input_filename + "_NeoepCompare_V1_Extended.txt"       #NOTE: for directory "/180531_Velip_V3"
    file_write = open( path_write, 'w' )


    print "----- Read File -----"
    ##ACTUAL FILE
    # path_read_actual = DIR_RESULT_ANTIPD1 + "/180130_Thres0_SampAllExpress_NeoepCompare_V1.txt"
    # path_read_actual = DIR_RESULTS_FOLDER + "/180403_Thres0_File" + str(file_num) + "_NeoepCompare_V1.txt"  #NOTE: for directory "/180403_Velip_V2"
    # path_read_actual = DIR_RESULTS_FOLDER + "/" + input_filename + "_Thres0_File" + str(file_num) + "_NeoepCompare_V1.txt"  #NOTE: for directory ""/180531_Velip_V3"
    path_read_actual = output_dir + "/" + input_filename + "_NeoepCompare_V1.txt"  #NOTE: for directory ""/180531_Velip_V3"
    df = pd.read_csv( path_read_actual, sep = '\t' )
    ##TEST FILE
    # df = pd.read_csv( path_read_actual, sep = '\t', nrows = 200 )

    #calculate the MHC affinity fold change before & after 
    print "----- Calculate the MHC affinity fold change before & after mutation -----"
    df['fold_change'] = df.apply( lambda row: float(row['ic50_score_1']) / row['ic50_score_2'], axis=1 )

    #make sure the neoepitope sequence considered before & after are not the same peptide sequence     
    print "----- Determine if the neoepitope before & after alteration are the same sequence -----"
    df['doublecheck_pep_seq'] = df.apply( lambda row: row['peptide_1'] == row['peptide_2'], axis=1 )
    

    print "----- Find Preferred HLA for each altered peptide -----"
    #retrieve the max HLA alleles for each altered peptide
    series_pep_mhc_max = df.groupby('peptide_2')['mhc_score_2'].max()
    hash_pep_hla = {}       #key = altered peptide sequence ('peptide_2'), v = another hash with the following keys "preferred_HLA" & "preferred_HLA_subgroup" string that are the HLA alleles (HLA are the human version of MHC)
    total_rows = len( series_pep_mhc_max )
    update = 0
    for i, (k,v) in enumerate( series_pep_mhc_max.iteritems() ):      #k = peptide_2, v = mhc_score_2
        ##TEST::
        # print "print k = ", k, " & v = ", v

        hla_alleles = df[ (df['peptide_2'] == k) & (df['mhc_score_2'] == v) ]['allele'].tolist()
        hash_pep_hla[k] = {}
        hash_pep_hla[k]['preferred_HLA'] = ','.join( hla_alleles )
        hash_pep_hla[k]['preferred_HLA_subgroup'] = ','.join( list( set([x.split('*')[0] for x in hla_alleles]) ) )
        # hash_pep_hla[k] = ','.join( hla_alleles )

        update = mokhaPy.dataMeasure_rowsPeriodicUpdate( total_rows, i, update, 0.001 )


    #DIFFERENCE: This is where 180325_AntiPD1_Label_Samples_V3.p & 180224_AntiPD1_Label_Samples_V2.py are different, as 180224_AntiPD1_Label_Samples_V2.py does not have this
    #assign a specific neoepitope to a given mutation - do this by ranking neoepitopes for a given mutation based on MHC affinity value
    print "----- Rank Neoepitopes for each mutation based on MHC affinity -----"
    """
    As each mutation generates multiple neoepitopes & each neoepitope is evaluated against multiple MHC alleles, this ranks each neoepitope for a given mutation by the MHC affinity (highest affinity), makes it easier to assign a specific neoepitope to a mutation (just pick neoepitope with highest affinity)
    -NOTE: as 'mhc_score_2' is a -log(), the least negative number corresponds to the highest affinity. I could have used 'ic50_score_2' with "ascending = True" instead
    """
    df['mut_neoep_rank_1'] = df.groupby( ['case_id', 'genome_pos'] )['mhc_score_2'].rank( ascending = False )
    # df.sort_values( ['case_id', 'genome_pos', 'mut_neoep_rank_1'], inplace = True )


    print "----- Begin Writing File -----"
    total_rows = len( df )
    update = 0
    for i, (i_row, row) in enumerate( df.iterrows() ):
        #create new columns with desired information
        row['preferred_HLA'] = hash_pep_hla[ row['peptide_2'] ]['preferred_HLA']
        row['preferred_HLA_subgroup'] = hash_pep_hla[ row['peptide_2'] ]['preferred_HLA_subgroup']
        
        #retrieve the gene expression
        gene_sym = row['gene_symbol']
        express_cond_id = retrieve_condition_id_express( row['case_id'] )
        if express_cond_id and df_express != None:
            find_express_count = df_express[df_express['GeneName'] == gene_sym]


            ##TEST:: make sure correct column of gene expression is extracted
            print "TEST:: express_cond_id = ", express_cond_id, " & find_express_count = ", find_express_count, " & row['case_id'] = ", row['case_id']

            gene_express_count = 0 if find_express_count.empty else find_express_count.iloc[0][express_cond_id]
            
            #calculate the gene percentile
            find_express_percentile = calculate_percentileofscore( df_express[express_cond_id], gene_express_count )
        else:
            gene_express_count = 0
            find_express_percentile = 0;

        row['gene_express'] = gene_express_count
        row['gene_percentile'] = find_express_percentile


        #need to transpose (to.frame().T) the row else it will not be written horizontally to match up with the columns
        row.to_frame().T.to_csv( path_write, sep = '\t', mode = 'a', header = HEADFILE_STAT, index = False )
        if HEADFILE_STAT:
            HEADFILE_STAT = False

        update = mokhaPy.dataMeasure_rowsPeriodicUpdate( total_rows, i, update, 0.001 )

    print "----- Writing File Complete -----"


    ##------- I MAY BE ABLE TO DELETE THIS ------
    # #modify actual file
    # #read in file
    # #ACTUAL FILE
    # # path_read_actual = DIR_RESULT_ANTIPD1 + "/180130_Thres0_SampAllExpress_NeoepCompare_V1.txt"
    # # df = pd.read_csv( path_read_actual, sep = '\t', chunksize = 5000 )
    # # #TEST FILE
    # # df = pd.read_csv( path_read_actual, sep = '\t', nrows = 200000, , chunksize = 5000 )
    # # #TEST FILE
    # path_read_test = DIR_RESULT_ANTIPD1 + "/180210_NeoepTEST.txt"
    # df = pd.read_csv( path_read_test, sep = '\t', chunksize = 100 )

    # for i_chunk, chunk in enumerate( df ):

    #     # print "current chunk: ", i_chunk, " & len( chunk ) = ", len( chunk ), " & type( chunk ) = ", type( chunk )
    #     print "current chunk: ", i_chunk, " & len( chunk ) = ", len( chunk )


    #     #group chunk by altered peptide sequence, find the HLA subgroup (HLA-A, HLA-B, HLA-C)
    #     chunk_group = chunk.groupby('peptide_2')
    #     for name, group in chunk_group:
    #         pass

    #     #retrieve the max HLA alleles for each altered peptide
    #     series_pep_mhc_max = chunk.groupby('peptide_2')['mhc_score_2'].max()
    #     hash_pep_hla = {}       #key = altered peptide sequence ('peptide_2'), v = string that are the HLA alleles (HLA are the human version of MHC)
    #     for k,v in series_pep_mhc_max.iteritems():      #k = peptide_2, v = mhc_score_2
    #         ##TEST::
    #         print "print k = ", k, " & v = ", v

    #         hla_alleles = df[ (df['peptide_2'] == k) & (df['mhc_score_2'] == v) ]['allele'].tolist()
    #         hash_pep_hla[k] = ','.join( hla_alleles )


    #     total_rows = len( chunk )
    #     update = 0
    #     for i, (i_row, row) in enumerate( chunk.iterrows() ):
    #         #create new columns with desired information
    #         row['irRECIST'] = hash_sample_irRECIST[ row['case_id'] ]
    #         row['preferred_HLA'] = hash_pep_hla[ row['peptide_2'] ]

    #         #need to transpose (to.frame().T) the row else it will not be written horizontally to match up with the columns
    #         row.to_frame().T.to_csv( path_write, sep = '\t', mode = 'a', header = HEADFILE_STAT, index = False )
    #         if HEADFILE_STAT:
    #             HEADFILE_STAT = False

    #         update = mokhaPy.dataMeasure_rowsPeriodicUpdate( total_rows, i, update, 0.05 )

    file_write.close()


elapse_time = time.time() - start_time
mokhaPy.timeElapse_convertToHMS( elapse_time )
# print "------------ Algorithm Completed: 180224_AntiPD1_Label_Samples_V2.py ------------"
# print "------------ Algorithm Completed: 180325_AntiPD1_Label_Samples_V3.py ------------"
# print "------------ Algorithm Completed: 180815_Velip_Label_Samples_V1.py ------------"
# print "------------ Algorithm Completed: 180825_Velip_Label_Samples_V2.py ------------"
print "------------ Algorithm Completed: 180929_Velip_Label_Samples_V3.py ------------"

