#/usr/bin/python

#python libraries
import os
import sys
import time

import tabix
import pandas as pd
from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv5 import Exon, Isoform, MultiIsoform, SpliceJunction, IsoformSJ
from mokhaPy import mokhaPy

#Constants: columns for 160802_Library_ThresEstimates.txt
COL_SAMPLE = 0
COL_READCOUNT = 10

#Constants: file directory
# DIR_ALLSJ = DIR_PROJ + "/150802_TophatSamples"
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv5/TDD"
DIR_DATA = DIR_CURR + "/Data"
DIR_RESULTS = DIR_CURR + "/Results"
# DIR_RESULTS = DIR_CURR + "/Results/160729_Analyze_KF"
# DIR_RESULTS = DIR_CURR + "/Results/160731_Analyze_KF"
# DIR_RESULTS = DIR_CURR + "/Results/160909_Analyze_KF"

DIR_SVSv5_PLAY = DIR_PROJ + '/160704_SVSv5_Play'
DIR_TABIX = '/home/mokha/Documents/Krauthammer_Lab/160427_SJFalsePositives/Data/Tabix'
DIR_THRES = DIR_PROJ + "/160803_ReadCountToThreshold/Results"

def compile_sj_noncanon( hash_sj_noncanon ):
    """
    Args:
        hash_sj_noncanon = hash where k = splice junction ID & v = SpliceJunction instances where sj.canon is not True
    Function: compiles all non-canonical splice junctions into a string
    """
    str_sj = None
    for k,v in hash_sj_noncanon.iteritems():    #k = splice junction ID & v = SpliceJunction instances where sj.canon is not True
        if not str_sj:
            str_sj = v.str_aberrations()
        else:
            str_sj += '||' + v.str_aberrations()

    return str( str_sj )


print "------------ Algorithm: 160918_IsoformSJ_1.py ------------"
""" Reconstruct transcripts based on Splice Junctions """

start_time = time.time()

#assign all splice junctions to specific gene: go through cruzdb & find end points for each gene --> assign 
# g = Genome( 'sqlite:////tmp/hg19.db' )
g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )
#set the file that contains the prevalence count for all aberrant SJ
IsoformSJ.set_tabix_prevalence( DIR_SVSv5_PLAY + "/Results/sorted_160721_sj_prevalence.bed.gz" )

#retrieve file of thresholded SJ
df_thres = pd.read_csv( DIR_THRES + '/160802_Library_ThresEstimates.txt', sep = '\t', index_col = 0 )       #index_col = 0 assigns the index to the sample names (first column)

print "df_thres columns = ", df_thres.columns.values
print "df_thres rows = ", df_thres.index.tolist()


#total number of penalties allowed to make each transcript
max_penalty = int( sys.argv[1] )
canon_transcript_only = False if max_penalty > 0 else True

sj_thres_index = '-1' 

sample_start = 0
sample_end = 1

#retrieve all isoforms
# all_genes = g.refGene.filter_by().all()
##TEST:: genes to test
# all_genes = g.refGene.filter_by( name2 = 'TTN' ).first()
# all_genes = [all_genes]
all_genes = g.refGene.filter_by( name2 = 'AGRN' ).all()
# all_genes = g.refGene.filter_by( name2 = 'RPS8' ).all()
# all_genes = g.refGene.filter_by( name2 = 'AGRN' ).first()
# all_genes = g.refGene.filter_by( name2 = 'DIXDC1' ).all()
# all_genes = g.refGene.filter_by( name2 = 'BRAF' ).all()
# all_genes = g.refGene.filter_by( name2 = 'WLS' ).all()


#get all samples
for (dirpath, dirnames, filenames) in os.walk( DIR_TABIX ):
    break
all_tabix = [x for x in filenames if not '.tbi' in x]

print "Enumerated Sample List:"
for i, sample_tabix in enumerate( all_tabix ):
    print i, ": ", sample_tabix

#ACTUAL FILE
# file_write = open( DIR_RESULTS + '/160902_FFP_' + str( sample_start ) + '_' + str( sample_end ) + '_Penalty' + str(max_penalty) + '.txt', 'w' )
##TEST FILE
file_write = open( DIR_RESULTS + '/TRANSCRIPT_RECONSTRUCT_' + str( sample_start ) + '_' + str( sample_end ) + '_Penalty' + str(max_penalty) + '_YesThres_' + sj_thres_index + '.txt', 'w' )

header = "sample_name\tgene_name\tisoform_id\tstrand_sign\ttranscript\tsj_id\t"
header+= "aberrant_sj\tnum_aberrant_sj\taberrant_read_count\tsj_all_aberrant\t"
header+= "exonic\tintronic\tframeshift\texon_skips\t"
header+= "num_artificial\tprevalence_sample\tprevalence_control\tSpliceSightV1\n"
file_write.write( header )

list_isoforms = []         #records all genes that have been visited
errors_genes = []       #record genes that could not be retrieved

##TEST:: to shorten analysis of genes - all_genes = all_genes[0:100]

#retrieve all SJs for a specific sample, reconstruct transcripts, and write to file
total_rows = len( all_genes )
update = 0
for i, each_isoform in enumerate( all_genes ):       #each_isoform = cruzdb information about gene isoform
    isoform_name = each_isoform.name
    percent_complete = ( float(i) / total_rows ) * 100
    print "gene symbol = ", each_isoform.name2, " & isoform = ", each_isoform.name, " | percent = ", percent_complete

    #make sure gene is protein coding, else skip it
    if not each_isoform.is_coding:
        continue

    #make sure not gene name duplicate, and record gene name as to not repeat
    if isoform_name in list_isoforms:
        continue
    else:
        list_isoforms.append( isoform_name )

    isoform_range = str( each_isoform.chrom ) + ':' + str( each_isoform.txStart ) + '-' + str( each_isoform.txEnd )


    #go through each sample
    # total_rows_2 = len( all_tabix )
    # update_2 = 0
    for i2, each_sample in enumerate( all_tabix[sample_start:sample_end] ):
        #retrieve sample name, path to file, and SJ threshold
        sample_name = each_sample.split('_')[1]
        sample_path = DIR_TABIX + '/' + each_sample
        #retrieve the threshold
        # sample_sj_thres = df_thres[ df_thres['sample'] == sample_name.lower() ]['-1'].astype( float )
        sample_sj_thres = df_thres[ sj_thres_index ][ sample_name.lower() ].astype( float )
        # sample_sj_thres = 0

        ##TEST:: print "Sample ", sample_name, ": gene = ", each_isoform.name2, " & isoform = ", each_isoform.name, " & sample_sj_thres = ", sample_sj_thres

        obj_isoform_sj = IsoformSJ( each_isoform.name, sample_path, sample_sj_thres )


        #build transcript
        construct_transcripts = obj_isoform_sj.reconstruct_transcript( canon_transcript_only, max_penalty )

        for i3, each_transcript in enumerate( construct_transcripts ):

            ##TEST:: 
            print "Transcript ", i3, ": sample = ", sample_name, " & gene = ", each_isoform.name2, " & isoform = ", each_isoform.name

            transcript_pos = '>>'.join( [str( sj ) for sj in each_transcript] )
            transcript_sjid = ','.join( [sj.sj_id for sj in each_transcript] )
            # transcript_readcount = sum( [sj.read_count for sj in each_transcript] ) 
            transcript_sj_aberrant = [sj.sj_id for sj in each_transcript if not sj.canon]
            aberrant_readcount = sum( [sj.read_count for sj in each_transcript if not sj.canon] )
            transcript_effects = { sj.sj_id:sj.isoform_aberrants[ each_isoform.name ] for sj in each_transcript if each_isoform.name in sj.isoform_aberrants.keys() }


            # sj_intronic = [k for k,v in transcript_effects.iteritems() if 'intronic' in v]
            # sj_exonic = [k for k,v in transcript_effects.iteritems() if 'exonic' in v]

            #Version 1 - retrieve arrays for frameshift, exon skip
            # sj_frameshift = [k for k,v in transcript_effects.iteritems() if not v['frame_preserved']]
            # sj_exon_skip = [k for k,v in transcript_effects.iteritems() if v['exon_skip']]
            #Version 2 - retrieve hashes for frameshift, exon skip
            sj_noncanon = {sj.sj_id:sj for sj in each_transcript if not sj.canon and each_isoform.name in sj.isoform_aberrants.keys() }
            sj_all_noncanon = compile_sj_noncanon( sj_noncanon )
            sj_frameshift = {k:v for k,v in transcript_effects.iteritems() if not v['frame_preserved']}
            sj_exon_skip = {k:v for k,v in transcript_effects.iteritems() if v['exon_skip']}

            sj_simulant = [sj.sj_id for sj in each_transcript if sj.read_count < 0 ]

            #get sample prevalence & control prevalence counts
            #version 1 - look at splice junctions whose prevalence > 0
            sj_sample_prevalence = { sj.sj_id:sj.sample_prevalence for sj in each_transcript if sj.sample_prevalence > 0 }
            sj_control_prevalence = { sj.sj_id:sj.control_prevalence for sj in each_transcript if sj.control_prevalence > 0 }
            #version 2 - look at splice junctions that are not canonical
            # sj_sample_prevalence = { sj.sj_id:sj.sample_prevalence for sj in each_transcript if not sj.canon }
            # sj_control_prevalence = { sj.sj_id:sj.control_prevalence for sj in each_transcript if not sj.canon }

            #retrieve string that will help make SpliceSight string
            hash_json = {}
            hash_json['gene'] = obj_isoform_sj.json_isoform_info()
            hash_json['sj'] = obj_isoform_sj.json_allsj_info( each_transcript )

            header = "sample_name\tgene_name\tisoform_id\tstrand_sign\ttranscript\tsj_id\t"
            header+= "aberrant_sj\tnum_aberrant_sj\taberrant_read_count\tsj_all_aberrant\t"
            header+= "exonic\tintronic\tframeshift\texon_skips\t"
            header+= "num_artificial\tprevalence_sample\tprevalence_control\tSpliceSightV1\n"
            #write to file
            line = sample_name + '\t'
            line+= each_isoform.name2 + '\t'
            line+= each_isoform.name + '\t'
            line+= str( obj_isoform_sj.strand ) + '\t'
            line+= transcript_pos + '\t'
            line+= transcript_sjid + '\t'

            line+= ','.join( transcript_sj_aberrant ) + '\t'
            line+= str( len(transcript_sj_aberrant) ) + '\t'
            line+= str( aberrant_readcount ) + '\t'
            line+= sj_all_noncanon + '\t'
            
            # line+= ','.join( sj_exonic ) + '\t'
            # line+= ','.join( sj_intronic ) + '\t'
            line+= 'NA' + '\t'                      #currently not exonic status (all SJs should be exonic for now)
            line+= 'NA' + '\t'                      #currently not recording SJs that go intronic
            
            # line+= ','.join( sj_frameshift ) + '\t'
            # line+= ','.join( sj_exon_skip ) + '\t'
            line+= str( sj_frameshift ) + '\t'
            line+= str( sj_exon_skip ) + '\t'

            line+= str( len(sj_simulant) ) + '\t'
            line+= str( sj_sample_prevalence ) + '\t'
            line+= str( sj_control_prevalence ) + '\n'
            file_write.write( line )

            ##TEST::
            print each_isoform.name2, ":", each_isoform.name," - transcript ", i3, ": aberrant_readcount = ", aberrant_readcount, " & # sj in transcript = ", len(transcript_sj_aberrant) 

        # update_2 = mokhaPy.dataMeasure_rowsPeriodicUpdate( total_rows_2, i2, update_2, 0.1 )

    update = mokhaPy.dataMeasure_rowsPeriodicUpdate( total_rows, i, update, 0.002 )

file_write.close()


elapse_time = time.time() - start_time
mokhaPy.timeElapse_convertToHMS( elapse_time )
print "------------ Algorithm Completed: 160918_IsoformSJ_1.py ------------"