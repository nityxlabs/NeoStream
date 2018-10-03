#/usr/bin/python
import os
import sys
import subprocess
import time


# sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
sys.path.insert( 0, "../Algorithms" )
from mokhaPy import mokhaPy

#Constants - directories
# DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = "../Algorithms"
DIR_DATA =  "../Data"
DIR_RESULTS = "../Results"
# DIR_RESULTS_FOLDER = DIR_RESULTS + "/171204_NeoepProcess_V15"
# DIR_RESULTS_FOLDER = DIR_RESULTS + "/180203_Velip_V1"
# DIR_RESULTS_FOLDER = DIR_RESULTS + "/180403_Velip_V2"
# DIR_RESULTS_FOLDER = DIR_RESULTS + "/180531_Velip_V3"
# DIR_RESULTS_FOLDER = DIR_RESULTS + "/180815_Velip_V1"
DIR_RESULTS_FOLDER = DIR_RESULTS

# DIR_DATA_VELIP = DIR_PROJ + "/170304_NeoantigenExtract/Data/Velip"


print "------------ Algorithm: 180820_Pipeline_NeoStream_V2.py ------------"
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

#enter the inputs
select_file = sys.argv[1]       #enter the file 
date_analysis = sys.argv[2]     #enter a desirable file prefix

#create path to file being analyzed
path_velip_file = DIR_DATA + "/" + select_file
check_dir = os.path.exists( path_velip_file )
if not check_dir:
    wait = input( "ERROR: directory " + path_velip_file + " does not exist! -> Ctrl+C to leave script" )


thres_percentile_gene_exp = 0    #this is gene expression percentile threshold - select 0 so I can record all genes & their respective expression level & percentile

#create file name
curr_file_name = date_analysis    #this will record the outputted file name
curr_file_name += "_Thres" + str( thres_percentile_gene_exp )
curr_file_name += "_File" + str( select_file )

print "------ Analysis: curr_file_name = ", curr_file_name


##VERSION 6
#Step 1 - evaluate each genomic alteration
# p1 = os.system( "python 171205_ProcessGenomeAlt_AntiPD1_V1.py " + curr_file_name + " " + str(thres_percentile_gene_exp) + " " + DIR_RESULTS_FOLDER )
# run_alg_1 = "180203_ProcessGenomeAlt_V1.py"
# run_alg_1 = "180815_ProcessGenomeAlt_V1.py"
run_alg_1 = "180825_ProcessGenomeAlt_V2.py"
is_seq_WGS = "N"        #this is just used to append suffix to condition
p1 = os.system( "python " + run_alg_1 + " " + curr_file_name + " " + str(thres_percentile_gene_exp) + " " + path_velip_file + " " + DIR_RESULTS_FOLDER + " " + is_seq_WGS )
if p1 != 0:
    status_loop = False
    print "FULL PROGRAM CEASE: Error in Algorithm 1: " + run_alg_1
    # break
    #Use the "input" instead to stop script (can only use "break" in loop)
    wait = input( "FULL PROGRAM CEASE: Error in Algorithm 1: " + run_alg_1 + " -> Ctrl+C to leave script" )


#Step 2 - determine the neoepitope affinity for specific MHC alleles
# run_alg_2 = "171205_NeoantigenMHC_Eval_PD1_V1.py.py"
# run_alg_2 = "180112_NeoantigenMHC_Eval_PD1_V2.py"
# run_alg_2 = "180203_NeoantigenMHC_Eval_Velip_V1.py"
run_alg_2 = "180815_NeoantigenMHC_Eval_Velip_V1.py"
p2 = os.system( "python " + run_alg_2 + " " + curr_file_name + " " + DIR_RESULTS_FOLDER )
if p2 != 0:
    status_loop = False
    print "FULL PROGRAM CEASE: Error in Algorithm 2: " + run_alg_2
    # break
    #Use the "input" instead to stop script (can only use "break" in loop) 
    wait = input( "FULL PROGRAM CEASE: Error in Algorithm 2: " + run_alg_2 + " -> Ctrl+C to leave script" )

#Step 2B - need to label the neoepitope data with information like ranking neoepitopes per mutation based on MHC affinity, patient data, etc.
# run_alg_2b = "180411_Velip_Label_Samples_V1.py"
# run_alg_2b = "180815_Velip_Label_Samples_V1.py"
# run_alg_2b = "180825_Velip_Label_Samples_V2.py"
run_alg_2b = "180929_Velip_Label_Samples_V3.py"
p2b = os.system( "python " + run_alg_2b + " " + curr_file_name + " " + DIR_RESULTS_FOLDER )
if p2b != 0:
    status_loop = False
    print "FULL PROGRAM CEASE: Error in Algorithm 2B: " + run_alg_2b
    # break
    #Use the "input" instead to stop script (can only use "break" in loop) 
    wait = input( "FULL PROGRAM CEASE: Error in Algorithm 2B: " + run_alg_2b + " -> Ctrl+C to leave script" )

                                

# #NOTE: COMMENTING OUT THIS FOR NOW, JUST NEED PERFORM IEDB ANALYSIS --> NEED IT UNCOMMENT LATER FOR FULL ANALYSIS
# #Step 3 - summarize the statistics for the neoepitopes
# # p3 = os.system( "python 171124_NeoepStatistics_V2.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
# # p3 = os.system( "python 171129_NeoepStatistics_V2.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
# # p3 = os.system( "python 171204_NeoepStatistics_V3.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
# # p3 = os.system( "python 180125_NeoepStatistics_PD1_V2.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
# # p3 = os.system( "python 180203_NeoepStatistics_Velip_V1.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
# p3 = os.system( "python 180815_NeoepStatistics_PD1_V2.py " + curr_file_name + " " + DIR_RESULTS_FOLDER )
# if p3 != 0:
#     status_loop = False
#     print "FULL PROGRAM CEASE: Error in Algorithm 3: 180125_NeoepStatistics_PD1_V2.py"
#     # break
#     #Use the "input" instead to stop script (can only use "break" in loop) 
#     wait = input( "FULL PROGRAM CEASE: Error in Algorithm 3: " + run_alg_3 + " -> Ctrl+C to leave script" )
 

# #Step 2C - Joann Sweasy output program where I select neoepitopes of interest. NOTE: I THINK I SHOULD RUN THIS SEPARATELY FROM THE PIPELINE - as of now the algorithm "180611_FilterNeoepSweasy_V1.py" is not suitable for the loop section of this pipeline
# run_alg_4 = "180611_FilterNeoepSweasy_V1.py"
# p4 = os.system( "python " + run_alg_4 + " " + select_file )
# if p4 != 0:
#     status_loop = False
#     print "FULL PROGRAM CEASE: Error in Algorithm 4: " + run_alg_4
#     break
#     #Use the "input" instead to stop script (can only use "break" in loop) 
#     wait = input( "FULL PROGRAM CEASE: Error in Algorithm 4: " + run_alg_4 + " -> Ctrl+C to leave script" )






elapse_time = time.time() - start_time
mokhaPy.timeElapse_convertToHMS( elapse_time )
print "------------ Algorithm Completed: 180820_Pipeline_NeoStream_V2.py ------------"