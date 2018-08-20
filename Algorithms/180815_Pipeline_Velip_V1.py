#/usr/bin/python
import os
import sys
import subprocess
import time

# sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/180815_NeoStream/Algorithms" )
from mokhaPy import mokhaPy

#Constants - directories
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/180815_NeoStream"
DIR_DATA = DIR_CURR + "/Data"
DIR_RESULTS = DIR_CURR + "/Results"
# DIR_RESULTS_FOLDER = DIR_RESULTS + "/171204_NeoepProcess_V15"
# DIR_RESULTS_FOLDER = DIR_RESULTS + "/180203_Velip_V1"
# DIR_RESULTS_FOLDER = DIR_RESULTS + "/180403_Velip_V2"
# DIR_RESULTS_FOLDER = DIR_RESULTS + "/180531_Velip_V3"
DIR_RESULTS_FOLDER = DIR_RESULTS + "/180815_Velip_V1"

DIR_DATA_VELIP = DIR_PROJ + "/170304_NeoantigenExtract/Data/Velip"


print "------------ Algorithm: 180815_Pipeline_Velip_V1.py ------------"
"""
Example of python call:
    "python 171005_ProcessGenomeAlt_V5.py " + date_analysis + " " + str(velip_file_oi) + " " + pick_condition
    -python 171005_Pipeline_Neoep_MHC_V2.py 171005 0 0 - DM & pc
    -python 171005_Pipeline_Neoep_MHC_V2.py 171005 1 0 - Velip & pc
    -python 171005_Pipeline_Neoep_MHC_V2.py 171005 0 1 - DM & br
    -python 171005_Pipeline_Neoep_MHC_V2.py 171005 1 1 - Velip & br
"""

start_time = time.time()

#check if output directory exists
check_dir = os.path.exists( DIR_RESULTS_FOLDER )
if not check_dir:
    wait = input( "ERROR: directory " + DIR_RESULTS_FOLDER + " does not exist! -> Ctrl+C to leave script" )

date_analysis = sys.argv[1]
# select_file = int( sys.argv[2] )
#enter range of samples: if I only want to analyze sample 3, then file_range_1 & file_range_2 = 3
# file_range_1 = int( sys.argv[2] )
# file_range_2 = int( sys.argv[3] )
file_range_1 = 1
file_range_2 = 5        #ignore file 6
# for select_file in range(1, 5):
for select_file in range( file_range_1, file_range_2 + 1 ):
    # #DATASET 1 - 180403_VelipData
    # if select_file == 1:
    #     path_velip_file = DIR_DATA_VELIP + "/180403_VelipData/Exome/ExomeI/SNV_WES.maf"
    # elif select_file == 2:
    #     path_velip_file = DIR_DATA_VELIP + "/180403_VelipData/Exome/ExomeII/SNV_WES_DM.maf"
    # elif select_file == 3:
    #     path_velip_file = DIR_DATA_VELIP + "/180403_VelipData/Exome/ExomeIII/SNV_WES_pooled.maf"
    # elif select_file == 4:      #this is the whole-genome sequencing sample
    #     path_velip_file = DIR_DATA_VELIP + "/180403_VelipData/WGS/WGSI/SNV_WGS_somatic.maf.txt"
    # elif select_file == 5:      #this is the whole-genome sequencing sample
    #     path_velip_file = DIR_DATA_VELIP + "/180403_VelipData/WGS/WGSII/SNV_WGS_somatic_DM.maf.txt"

    #DATASET 2 - 180531_VelipData
    if select_file == 1:
        path_velip_file = DIR_DATA_VELIP + "/180531_VelipData/SNV_WES_Velip.maf"
        is_seq_WGS = "N"
    elif select_file == 2:
        path_velip_file = DIR_DATA_VELIP + "/180531_VelipData/SNV_WES_DM.maf"
        is_seq_WGS = "N"
    elif select_file == 3:
        path_velip_file = DIR_DATA_VELIP + "/180531_VelipData/SNV_WES_pooled.maf"
        is_seq_WGS = "N"
    elif select_file == 4:      #this is the whole-genome sequencing sample
        path_velip_file = DIR_DATA_VELIP + "/180531_VelipData/SNV_WGS_Velip_somatic.maf"
        is_seq_WGS = "Y"
    elif select_file == 5:      #this is the whole-genome sequencing sample
        path_velip_file = DIR_DATA_VELIP + "/180531_VelipData/SNV_WGS_DM_somatic.maf"
        is_seq_WGS = "Y"
    elif select_file == 6:      #IGNORE THIS FILE FOR NOW - it doesn't follow the same output structure
        path_velip_file = DIR_DATA_VELIP + "/180531_VelipData/SNV_WGS_OrigHCC_somatic.maf"
        is_seq_WGS = "Y"

    thres_percentile_gene_exp = 0    #this is gene expression percentile threshold - select 0 so I can record all genes & their respective expression level & percentile

    # input_pc_dm = int( sys.argv[2] )
    # # thres_percentile_gene_exp = int( sys.argv[2] )     #this is gene expression percentile threshold
    # ##DO NOT NEED THIS
    # # extra_label = "Velip_TEST"
    # if input_pc_dm == 1:
    #     extra_label = "Velip_PC"
    # else:
    #     extra_label = "Velip_DM"

    # thres_percentile_gene_exp = 50     #this is gene expression percentile threshold
    # velip_file_oi = int( sys.argv[2] )      #0 for "_DM" (control) or 1 for Velip (experimental)
    # pick_condition = int( sys.argv[3] )     #look at script "171005_ProcessGenomeAlt_V5.py", where 0 = "pc_" as this is the empty DNA vector, 1 = "br_" as this is the vector contain normal BRCA1

    list_inputs = [1]     #the inputs for parameters "velip_file_oi" & "pick_condition"

    for i, dumb_test in enumerate( list_inputs ):

        curr_file_name = date_analysis    #this will record the outputted file name
        curr_file_name += "_Thres" + str( thres_percentile_gene_exp )
        curr_file_name += "_File" + str( select_file )

        print "------ Analysis: ", i, " curr_file_name = ", curr_file_name


        # ##VERSION 1
        # #first evaluate each genomic alteration
        # os.system( "python 170829_ProcessGenomeAlt_V2.py " + curr_file_name )
        # #second - determine the neoepitope affinity for specific MHC alleles
        # os.system( "python 170909_NeoantigenMHC_Eval.py " + curr_file_name )


        # ##VERSION 2
        #first evaluate each genomic alteration
        #170925_ProcessGenomeAlt_V4.py - used for analyzing all genomic alterations
        # os.system( "python 170925_ProcessGenomeAlt_V4.py " + curr_file_name + " " + file_oi )
        #170925_ProcessGenomeAlt_V4_NonsilentSNV.py - used for analyzing only single point mutations, though it may look at insertions & deletions, but for now just only use for single point mutations
        # os.system( "python 170925_ProcessGenomeAlt_V4_NonsilentSNV.py " + curr_file_name + " " + file_oi )

        # ##VERSION 3
        # #Step 1 - evaluate each genomic alteration
        # os.system( "python 171010_ProcessGenomeAlt_V6.py " + curr_file_name + " " + str( velip_file_oi ) + " " + str(pick_condition ) + " " + str( thres_percentile_gene_exp ) + " " + DIR_RESULTS_FOLDER )
        # #Step 2 - determine the neoepitope affinity for specific MHC alleles
        # os.system( "python 171010_NeoantigenMHC_Eval_V4.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # #Step 3 - summarize the statistics for the neoepitopes
        # os.system( "python 171010_NeoepStatistics_V1.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )

        # ##VERSION 4 --> CURRENTLY THIS IS THE SAFEST VERSION
        # #Step 1 - evaluate each genomic alteration
        # p1 = os.system( "python 171010_ProcessGenomeAlt_V6.py " + curr_file_name + " " + str( velip_file_oi ) + " " + str( pick_condition ) + " " + str( thres_percentile_gene_exp ) + " " + DIR_RESULTS_FOLDER + " " + seq_platform )
        # if p1 != 0:
        #     status_loop = False
        #     print "FULL PROGRAM CEASE: Error in Algorithm 1: 171010_ProcessGenomeAlt_V6.py"
        #     break
        
        # #Step 2 - determine the neoepitope affinity for specific MHC alleles
        # p2 = os.system( "python 171031_NeoantigenMHC_Eval_V5.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # if p2 != 0:
        #     status_loop = False
        #     print "FULL PROGRAM CEASE: Error in Algorithm 2: 171031_NeoantigenMHC_Eval_V5.py"
        #     break
        
        # #Step 3 - summarize the statistics for the neoepitopes
        # # p3 = os.system( "python 171031_NeoepStatistics_V1.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # p3 = os.system( "python 171124_NeoepStatistics_V2.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # if p3 != 0:
        #     status_loop = False
        #     print "FULL PROGRAM CEASE: Error in Algorithm 3: 171031_NeoepStatistics_V1.py"
        #     break

        
        # ##VERSION 5 --> EXPERIMENTAL VERSION (can perform NMD/NSD analysis)
        # ##NOTE: COMMENTING OUT THIS FOR NOW, JUST NEED PERFORM IEDB ANALYSIS --> NEED IT UNCOMMENT LATER FOR FULL ANALYSIS
        # # #Step 1 - evaluate each genomic alteration
        # # p1 = os.system( "python 171205_ProcessGenomeAlt_AntiPD1_V1.py " + curr_file_name + " " + str(thres_percentile_gene_exp) + " " + DIR_RESULTS_FOLDER )
        # # if p1 != 0:
        # #     status_loop = False
        # #     print "FULL PROGRAM CEASE: Error in Algorithm 1: 171205_ProcessGenomeAlt_AntiPD1_V1.py"
        # #     break
        

        # # #Step 2 - determine the neoepitope affinity for specific MHC alleles
        # # # p2 = os.system( "python 171205_NeoantigenMHC_Eval_PD1_V1.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # # p2 = os.system( "python 180112_NeoantigenMHC_Eval_PD1_V2.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # # if p2 != 0:
        # #     status_loop = False
        # #     print "FULL PROGRAM CEASE: Error in Algorithm 2: 171205_NeoantigenMHC_Eval_PD1_V1.py"
        # #     break
        
        # #NOTE: COMMENTING OUT THIS FOR NOW, JUST NEED PERFORM IEDB ANALYSIS --> NEED IT UNCOMMENT LATER FOR FULL ANALYSIS
        # #Step 3 - summarize the statistics for the neoepitopes
        # # p3 = os.system( "python 171124_NeoepStatistics_V2.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # # p3 = os.system( "python 171129_NeoepStatistics_V2.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # # p3 = os.system( "python 171204_NeoepStatistics_V3.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # p3 = os.system( "python 180125_NeoepStatistics_PD1_V2.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # if p3 != 0:
        #     status_loop = False
        #     print "FULL PROGRAM CEASE: Error in Algorithm 3: 180125_NeoepStatistics_PD1_V2.py"
        #     break

        

        ##VERSION 6
        #Step 1 - evaluate each genomic alteration
        # p1 = os.system( "python 171205_ProcessGenomeAlt_AntiPD1_V1.py " + curr_file_name + " " + str(thres_percentile_gene_exp) + " " + DIR_RESULTS_FOLDER )
        # run_alg_1 = "180203_ProcessGenomeAlt_V1.py"
        run_alg_1 = "180815_ProcessGenomeAlt_V1.py"
        p1 = os.system( "python " + run_alg_1 + " " + curr_file_name + " " + str(thres_percentile_gene_exp) + " " + path_velip_file + " " + DIR_RESULTS_FOLDER + " " + is_seq_WGS )
        if p1 != 0:
            status_loop = False
            print "FULL PROGRAM CEASE: Error in Algorithm 1: " + run_alg_1
            break

        
        #Step 2 - determine the neoepitope affinity for specific MHC alleles
        # run_alg_2 = "171205_NeoantigenMHC_Eval_PD1_V1.py.py"
        # run_alg_2 = "180112_NeoantigenMHC_Eval_PD1_V2.py"
        # run_alg_2 = "180203_NeoantigenMHC_Eval_Velip_V1.py"
        run_alg_2 = "180815_NeoantigenMHC_Eval_Velip_V1.py"
        p2 = os.system( "python " + run_alg_2 + " " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        if p2 != 0:
            status_loop = False
            print "FULL PROGRAM CEASE: Error in Algorithm 2: " + run_alg_2
            break

        #Step 2B - need to label the neoepitope data with information like ranking neoepitopes per mutation based on MHC affinity, patient data, etc.
        # run_alg_2b = "180411_Velip_Label_Samples_V1.py"
        run_alg_2b = "180815_Velip_Label_Samples_V1.py"
        p2b = os.system( "python " + run_alg_2b + " " + str(select_file) )
        if p2b != 0:
            status_loop = False
            print "FULL PROGRAM CEASE: Error in Algorithm 2B: " + run_alg_2b
            break

                                        
        
        #NOTE: COMMENTING OUT THIS FOR NOW, JUST NEED PERFORM IEDB ANALYSIS --> NEED IT UNCOMMENT LATER FOR FULL ANALYSIS
        #Step 3 - summarize the statistics for the neoepitopes
        # p3 = os.system( "python 171124_NeoepStatistics_V2.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # p3 = os.system( "python 171129_NeoepStatistics_V2.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # p3 = os.system( "python 171204_NeoepStatistics_V3.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # p3 = os.system( "python 180125_NeoepStatistics_PD1_V2.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        # p3 = os.system( "python 180203_NeoepStatistics_Velip_V1.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        p3 = os.system( "python 180815_NeoepStatistics_PD1_V2.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
        if p3 != 0:
            status_loop = False
            print "FULL PROGRAM CEASE: Error in Algorithm 3: 180125_NeoepStatistics_PD1_V2.py"
            break
 

# #Step 2C - Joann Sweasy output program where I select neoepitopes of interest. NOTE: I THINK I SHOULD RUN THIS SEPARATELY FROM THE PIPELINE - as of now the algorithm "180611_FilterNeoepSweasy_V1.py" is not suitable for the loop section of this pipeline
# run_alg_4 = "180611_FilterNeoepSweasy_V1.py"
# p4 = os.system( "python " + run_alg_4 + " " + select_file )
# if p4 != 0:
#     status_loop = False
#     print "FULL PROGRAM CEASE: Error in Algorithm 4: " + run_alg_4
#     break


elapse_time = time.time() - start_time
mokhaPy.timeElapse_convertToHMS( elapse_time )
print "------------ Algorithm Completed: 180815_Pipeline_Velip_V1.py ------------"