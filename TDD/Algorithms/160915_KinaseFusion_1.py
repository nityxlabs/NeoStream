#/usr/bin/python
import os
import sys
import time
import re

import pandas as pd

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv5 import KinaseFusion
from mokhaPy import mokhaPy

DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv5"
DIR_DATA = DIR_CURR + "/TestData"
DIR_RESULTS = DIR_CURR + "/TestResults"
# DIR_RESULTS = DIR_CURR + "/Results/160729_Analyze_KF"
# DIR_RESULTS = DIR_CURR + "/Results/160731_Analyze_KF"
# DIR_RESULTS = DIR_CURR + "/Results/160909_Analyze_KF"

DIR_FUSION = DIR_PROJ + "/160510_GeneFusions"

#columns from kinase fusion (this is my random format, not an establish format)
COL_SAMPLE = 0
COL_LINE_ID = 2
COL_CHROM = 3
COL_POS_START = 4
COL_POS_END = 5
COL_FUSION_ID = 9
COL_ORIENTATION = 10
COL_COMPATIBLE = 11
COL_KINASE_POS = 12     #0 = no kinase in fusion, 1 = first gene is kinase, 2 = 2nd gene is kinase, 3 = both are kinases
COL_FUSION_NAME = 13
COL_GENE1 = 14
COL_GENE2 = 20
COL_ISO1 = 15
COL_ISO2 = 21
COL_STRAND1 = 16
COL_STRAND2 = 22
COL_KINASE_DOMAIN_1 = 19
COL_KINASE_DOMAIN_2 = 25

def get_control_fusions( path_fusion ):
    """
    Function: record all fusions that occur within the control samples, yunbmel & yuamp38f
    """
    control_sample = ['yunbmel', 'yuamp38f']

    control_fusions = []
    with open( path_fusion ) as f_f:
        for i, row in enumerate( f_f ):
            arr_rc = row.split( '\t' )

            if arr_rc[COL_SAMPLE] in control_sample and not arr_rc[COL_FUSION_NAME] in control_fusions:
                control_fusions.append( arr_rc[COL_FUSION_NAME] )

    return control_fusions

"""
Note on kinase domain inclusion:
-when denoting relative position, I use -1 (nucleotide position before kinase domain), 0 (nucleotide position within kinase domain), & 1 (nucleotide position after kinase domain)
-regardless of gene strand sign, the -1 = lower nucleotide position before kinase domain & 1 = higher nucleotide position after kinase domain
-kinase_first(): if the gene is on minus strand, then want to m
"""

def kinase_first( arr_rc ):
    """
    Function: looks at the strand sign & relative kinase domain orientation when the 1st gene in the fusion is a kinase
    NOTE: -1 means before kinase domain, 0 means within kinase domain, and 1 means after kinase domain
    """
    #if not a number (e.g. '?' or 'na')
    if not arr_rc[COL_KINASE_DOMAIN_1].isdigit():
        return False

    if list( arr_rc[COL_ORIENTATION] )[0].lower() == 'f':
         if int( arr_rc[COL_KINASE_DOMAIN_1] ) == 1 or int( arr_rc[COL_KINASE_DOMAIN_1] ) == 0:
            return True
    elif list( arr_rc[COL_ORIENTATION] )[0].lower() == 'r':
        if int( arr_rc[COL_KINASE_DOMAIN_1] ) == -1 or int( arr_rc[COL_KINASE_DOMAIN_1] ) == 0:
            return True

    return False


def kinase_second( arr_rc ):
    """
    Function: looks at the strand sign & relative kinase domain orientation when the 2nd gene in the fusion is a kinase
    NOTE: -1 means before kinase domain, 0 means within kinase domain, and 1 means after kinase domain
    """
    if not arr_rc[COL_KINASE_DOMAIN_2].isdigit():
        return False

    if list( arr_rc[COL_ORIENTATION] )[1].lower() == 'f':
         if int( arr_rc[COL_KINASE_DOMAIN_2] ) == -1 or int( arr_rc[COL_KINASE_DOMAIN_2] ) == 0:
            return True
    elif list( arr_rc[COL_ORIENTATION] )[1].lower() == 'r':
        if int( arr_rc[COL_KINASE_DOMAIN_2] ) == 1 or int( arr_rc[COL_KINASE_DOMAIN_2] ) == 0:
            return True

    return False
   

def kinase_pos( arr_rc ):
    """
    Function: checks if the kinase domain is kept in the gene fusion. Returns True if kinase is in fusion, else returns false
    """

    #if the column containing the kinase position is 0, means there is no kinase in this fusion
    if int( arr_rc[COL_KINASE_POS] ) < 1:
        return False

    #check the position the kinase is in in the fusion
    if int( arr_rc[COL_KINASE_POS] ) == 1:

        ##TEST::
        print "arr_rc = ", arr_rc

        stat = kinase_first( arr_rc )
    elif int( arr_rc[COL_KINASE_POS] ) == 2:

        ##TEST::
        print "arr_rc = ", arr_rc

        stat = kinase_second( arr_rc )
    elif int( arr_rc[COL_KINASE_POS] ) == 3:

        ##TEST::
        print "arr_rc = ", arr_rc

        if kinase_first( arr_rc ) or kinase_second( arr_rc ):
            stat = True
        else:
            stat = False

    return stat


print "------------ Algorithm: 160914_FilterKinaseFusions_V3.py ------------"
"""
Same as 160804_FilterKinaseFusions.py, however this also takes into account the orientation of the fusion (ff, rr, fr, rf)
"""


"""
-compatible
-not the same gene on both sides of the fusion
-can't be bound to non-coding RNA
-should I still keep things that are intronic??

Keeping Kinase domain in Fusion
kinase is plus strand OR forward??
-if kinase is 1st gene, then kinase domain = 1 (meaning the 2nd gene fuses after the kinase domain)
-if kinase is 2nd gene, then kinase domain = -1 (meaning the 1st gene fuses before the kinase domain)
kinase is minus strand OR reverse??
-if kinase is 1st gene, then kinase domain = -1 (meaning the 2nd gene fuses before the kinase domain)
-if kinase is 2nd gene, then kinase domain = 1 (meaning the 1st gene fuses after the kinase domain)
"""

start_time = time.time()

#set kinase gene annotation file
KinaseFusion.set_kinasefile( DIR_FUSION + "/Data/160910_KinaseAnnots_hg38_Final.txt" )

path_fusion = DIR_DATA + "/160914_KF_BRAF.txt"
# file_write = open( DIR_RESULTS + "/160914_KinaseFusion_Revised_V3.txt", 'w' )

#retrieve fusions that occur within the control sample
control_fusions = get_control_fusions( path_fusion )

#record all fusions as to not record duplicates - 2 measures to reduce redunancies
line_id_exists = []         #records the line ID to ensure no duplicates are recorded (measure 1)
fusion_id_exists = []       #records the fusion ID to ensure no duplicates are recorded (measure 2)

#read in kinase file
val_update = 0
total_rows = mokhaPy.dataMeasure_countAllLines( path_fusion )
with open( path_fusion ) as f_k:
    for i, row in enumerate( f_k ):

        ##TEST::
        if i != 11:
            continue

        #remove any new lines for each row with 
        row = re.sub( r'\r\n?', '\n', row )
        row = re.sub( r'\n?', '', row )

        arr_rc = row.split( '\t' )

        #save header for file
        if i == 0:
            #Write line in file
            # new_row = re.sub( r'\r\n?', '\t', row )
            # file_write.write( row + '\tkinase_domain_present\n' )
            #skip this row
            continue

        #should have integer in that column for each row, if not then skip
        if not arr_rc[COL_KINASE_POS].isdigit():
            continue

        #check if compatible
        # if not 'Yes' in arr_rc[COL_COMPATIBLE]:
        #     # print "compatible - ", arr_rc[COL_COMPATIBLE], " --> NO"
        #     continue
        #check to make sure both genes are not the same
        if arr_rc[COL_GENE1] == arr_rc[COL_GENE2]:
            # print "same gene - ", arr_rc[COL_GENE1], " & ", arr_rc[COL_GENE2] , " --> NO"
            continue
        #make sure no non-coding RNA
        # if 'NR_' in arr_rc[COL_ISO1] or 'NR_' in arr_rc[COL_ISO2]:
        #     # print "non-coding - ", arr_rc[COL_ISO1], " & ", arr_rc[COL_ISO2] , " --> NO"
        #     continue
        #should I check if intronic? - probably not now

        #check if fusion is within control
        # if arr_rc[COL_FUSION_NAME] in control_fusions:
        #     continue

        #check if fusion already exists, if so then skip, else record the fusion ID
        #NOTE: COL_FUSION_NAME is more stringent than COL_FUSION_ID, as it looks for gene fusion duplicates regardless of gene order in fusion
        # if arr_rc[COL_LINE_ID] in line_id_exists:
        #     continue
        # else:
        #     line_id_exists.append( arr_rc[COL_LINE_ID] )

        # if arr_rc[COL_FUSION_ID] in fusion_id_exists:
        #     continue
        # else:
        #     fusion_id_exists.append( arr_rc[COL_FUSION_ID] )

        

        print "YES I THRU!"

        hash_fusion = { "orientation": arr_rc[COL_ORIENTATION],
        "isoform_1": arr_rc[COL_ISO1],
        "isoform_2": arr_rc[COL_ISO2],
        "kinase_num": arr_rc[COL_KINASE_POS],
        "chrom_1": arr_rc[COL_CHROM].split('-')[0],
        "chrom_2": arr_rc[COL_CHROM].split('-')[1],
        "pos_1": arr_rc[COL_POS_START],
        "pos_2": arr_rc[COL_POS_END] }

        print "hash_fusion = ", hash_fusion
        
        obj_kf = KinaseFusion( hash_fusion )

        ##TEST:: see pandas Dataframe
        print "show df_1:\n", obj_kf.df_1
        print "show df_2:\n", obj_kf.df_2
        print "\n-------------------\n"

        ##TEST:: see exon range based on orientation of fusion
        print "exons for gene 1 = ", obj_kf.isoform_exon_range( 1 ), "\n"
        print "exons for gene 2 = ", obj_kf.isoform_exon_range( 2 ), "\n"

        print "number of exons coding for kinase = ", obj_kf.count_kinase_exons()

        
        # ##RETRIEVE INFORMATION ABOUT KINASE-CODING EXONS!!
        # #get number of exons encoding for kinase domain & total count of kinase-coding exons present
        # df_exon_present_1 = obj_kf.isoform_exon_range( 1 )
        # df_kinase_frac_1 = df_exon_present_1[ df_exon_present_1['exon_statKinaseDomain'].astype(int) == 0 ]
        # df_exon_present_2 = obj_kf.isoform_exon_range( 2 )
        # df_kinase_frac_2 = df_exon_present_2[ df_exon_present_2['exon_statKinaseDomain'].astype(int) == 0 ]


        # print "df_kinase_frac_1 = ", df_kinase_frac_1['exonNum']
        # print "df_kinase_frac_2 = ", df_kinase_frac_2['exonNum']

        # #quantify the number of exons that code for the kinsae domain
        # ( df_ake_1, df_ake_2 ) = obj_kf.retrieve_kinase_exons()     #df_ake = DataFrame All Kinase Exons

        # #get exons for kinase domain, count of number present, total number of kinase-coding exons, & percentage
        # #ke = kinase-coding exons
        # #kde = kinase domain exons
        # #ke_present = kinase-coding exons present
        # #ke_len_present = number of kinase-coding exons present
        # #ke_total = total number of kinase-coding exons in gene

        # ke_percent_present_1 = float( len( df_kinase_frac_1 ) ) / len( df_ake_1 ) if len( df_ake_1 ) > 0 else None
        # ke_percent_present_2 = float( len( df_kinase_frac_2 ) ) / len( df_ake_2 ) if len( df_ake_2 ) > 0 else None
        # hash_kde = {
        # "ke_present_1": list( df_kinase_frac_1['exonNum'] ),
        # "ke_len_present_1": len( df_kinase_frac_1 ),
        # "ke_total_1": list( df_ake_1['exonNum'] ),
        # "ke_total_len_1": len( df_ake_1 ),
        # "ke_percent_1": ke_percent_present_1,
        # "ke_bool_present_1": True if ke_percent_present_1 > 0 else False,
        # "ke_present_2": list( df_kinase_frac_2['exonNum'] ),
        # "ke_len_present_2": len( df_kinase_frac_2 ),
        # "ke_total_2": list( df_ake_2['exonNum'] ),
        # "ke_total_len_2": len( df_ake_2 ),
        # "ke_percent_2": ke_percent_present_2,
        # "ke_bool_present_2": True if ke_percent_present_2 > 0 else False,
        # }

        # #print information about kinase-coding exons
        # print "hash_kde: ", hash_kde



        # #check if kinase domain is present in gene
        # stat = kinase_pos( arr_rc )

        #write to file
        # new_row = re.sub( r'\r\n?', '', row )
        # file_write.write( row + '\t' + str( stat ) + '\n' )

        #update progress
        val_update = mokhaPy.dataMeasure_rowsPeriodicUpdate( total_rows, i, val_update, 0.01 )

# file_write.close()

elapse_time = time.time() - start_time
mokhaPy.timeElapse_convertToHMS( elapse_time )
print "------------ Algorithm Completed: 160914_FilterKinaseFusions_V3.py ------------"