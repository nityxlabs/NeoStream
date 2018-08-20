#/usr/bin/python
import sys
import time

from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv5 import SpliceJunction, Isoform, IsoformSJ, TranscribeTranscript, TranslateTranscript
from mokhaPy import mokhaPy

#Constants: file paths
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv5"
DIR_DATA = DIR_CURR + "/TestData"
DIR_RESULTS = DIR_CURR + "/TestResults"
DIR_GENOME = DIR_PROJ + "/ArchiveData/hg19.fa"

def create_obj_sj( hash_sj_info ):
    """
    Args:
        hash_sj_info = a hash_sj_info from pandas Dataframe, where each hash_sj_info is indexed by the column labels
    Function: creates a SpliceJunction instance for the splice junction recorded in the file. Information about the splice junction is recorded in a hash_sj_info in the file contained in the variable "arr_rc"
    """
    sj_id = hash_sj_info['sj_id']
    hash_sj_pos = Isoform.split_genome_pos( hash_sj_info['sj_range'] )
    chrom = hash_sj_pos['chrom']
    start = hash_sj_pos['start']
    end = hash_sj_pos['end']
    strand = '-' if int( hash_sj_info['strand'] ) == -1 else '+'       #needs to be in string format
    read_count = hash_sj_info['read_count']
    gene_sym = hash_sj_info['gene_name']
    # isoform_id = hash_sj_info['isoform_id']
    isoform_id = None
    sample_prevalence = hash_sj_info['prevalence_all']
    control_prevalence = hash_sj_info['prevalence_control']
    bool_intronic = True
    # obj_sj = SpliceJunction( sj_id, chrom, start, end, strand, read_count, gene_sym, isoform_id = None, sample_prevalence = 0, control_prevalence = 0, bool_intronic = False )
    obj_sj = SpliceJunction( sj_id, chrom, start, end, strand, read_count, gene_sym, isoform_id, sample_prevalence, control_prevalence, bool_intronic )

    return obj_sj

print "------------ TDD: 170322_BuildTranscriptFromSingleSJ.py ------------"
start_time = time.time()

g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )


#Test SpliceJunction - 1
# sample_name = 'yukim'       #do not need this for hash_sj_info
# sj_id = 'JUNC00120493'
# sj_range = 'chr12:56355542-56359719'
# strand = -1
# read_count = 11
# gene_sym = 'PMEL'
# isoform_id = 'NM_001200053'
# sample_prevalence = 4
# control_prevalence = 1

#Test SJ 2 - this SJ splices outside of the isoform's gene structure
# sj_id = 'JUNC00090748'
# sj_range = 'chr12:56121123-56122061'
# strand = -1
# read_count = 13
# gene_sym = 'CD63'
# isoform_id = 'NM_001257400' 
# sample_prevalence = 0
# control_prevalence = 0

#Test SJ 3 - Having issues with this SJ & TranscribeTranslate_v4
# sj_id = 'JUNC00075691'
# sj_range = 'chr10:73587861-73588580'
# # sj_range = 'chr10:73587861-73594263'
# strand = -1
# read_count = 125
# gene_sym = 'PSAP'
# isoform_id = 'NM_001042466' 
# sample_prevalence = 0
# control_prevalence = 0

#TEST SJ 4
# line that contains information about aberrant SJ
# JUNC00025482, chr2:224866639-224903797, read = 46, strand sign = -1, canonical = False, assigned isoform   1036 s = NM_001136530 | NM_006216 | NR_073116 | NM_001136528
# sj_id = 'JUNC00025482'
# sj_range = 'chr2:224866639-224903797'
# strand = -1
# read_count = 46
# gene_sym = 'SERPINE2'
# isoform_id = 'NM_001136530' 
# sample_prevalence = 0
# control_prevalence = 0

#TEST SJ 5
# #yuame  - PASSED Filters: row  1547 : gene =  HLA-DPA1  & isoform_id =  NM_033554  & sj_id =  JUNC00050057  & sj canon =  False  & is it false?  False  & gene strand sign =
# #TTV4 - CEFSJ 0:  4  -  JUNC00050057, chr6:33010660-33053527, read = 124, strand sign = -1, canonical = False, assigned isoforms = NM_001242524 | NM_001242525 | NM_033554  && isoform IDs =  ['', '', '', '', u'NM_001242524 | NM_001242525 | NM_033554']
# sj_id = 'JUNC00050057'
# sj_range = 'chr6:33010660-33053527'
# strand = -1
# read_count = 124
# gene_sym = 'HLA-DPA1'
# isoform_id = 'NM_033554' 
# sample_prevalence = 0
# control_prevalence = 0

#TEST SJ 6
#yuame  - PASSED Filters: row  1547 : gene =  HLA-DPA1  & isoform_id =  NM_033554  & sj_id =  JUNC00050057  & sj canon =  False  & is it false?  False  & gene strand sign =
#TTV4 - CEFSJ 0:  4  -  JUNC00050057, chr6:33010660-33053527, read = 124, strand sign = -1, canonical = False, assigned isoforms = NM_001242524 | NM_001242525 | NM_033554  && isoform IDs =  ['', '', '', '', u'NM_001242524 | NM_001242525 | NM_033554']
# sj_id = 'JUNC00050057'
# sj_range = 'chr6:33010660-33053527'
# strand = -1
# read_count = 124
# gene_sym = 'HLA-DPA1'
# isoform_id = 'NM_033554' 
# sample_prevalence = 0
# control_prevalence = 0

##TEST SJ 7
#yuame  - PASSED Filters: row  35749 : gene =  TTN  & isoform_id =  NM_133437  & sj_id =  JUNC00023247  & sj canon =  False  & is it false?  False  & gene strand sign =  -1
# JUNC00023247, chr2:179407088-179436026, read = 26, strand sign = 1, canonical = False, assigned isoforms = NM_133432 | NM_133437 | NM_001256850 | NM_003319 | NM_133378 | NM_001267550

sj_id = 'JUNC00023247'
sj_range = 'chr2:179407088-179436026'
strand = -1
read_count = 26
gene_sym = 'TTN'
isoform_id = 'NM_133437' 
sample_prevalence = 0
control_prevalence = 0

hash_sj_info = {'sj_id': sj_id,
'sj_range': sj_range,
'strand': strand,
'read_count': read_count,
'gene_name': gene_sym,
# 'isoform_id': isoform_id,     #use this if want to assign SJ to a specific isoform
'isoform_id': None,       #use this if I just want to see which isoforms are assigned to SJ
'prevalence_all': sample_prevalence,
'prevalence_control': control_prevalence
}

obj_sj = create_obj_sj( hash_sj_info )
hash_sepr = obj_sj.spliced_elems_position_range( isoform_id )       #sepr = spliced_elems_position_range

##TEST::
print "obj_sj = ", obj_sj
print "obj_sj ligated exons canonical: "
for k,v in obj_sj.spliced_elems_position_range( isoform_id, False ).iteritems():
    print "k = ", k, " & v = ", v
print "##########################################\n"

##TEST:: this retrieves the 2 elements that are ligated by the SJ (elements could be exons, introns)
tuple_se = obj_sj.spliced_elems( isoform_id )    #se =spliced_elems
print "SJ.spliced_elems()"
print "start_x = ", tuple_se[0]
print "end_x = ", tuple_se[1]
print "~~~~~~~\n"

##TEST:: this retrieves the 2 exons that could be ligated 
print "SJ.spliced_elems_position_range()"
for k,v in hash_sepr.iteritems():
    print "k = ", k, " & v = ", v

print "\n---------- Reconstruct transcript from single SJ -----------\n"

##Experiment 1 - use IsoformSJ.reconstruct_transcript_single_sj() to reconstruct transcript with single SJ
list_sj = [obj_sj]
group_sj = 0            #do not perform any grouping of SJs
hash_pos = None

#reconstruct transcript around single SJ
print "create IsoformSJ instance..."
simulant_sj = False
iso_sj = IsoformSJ( isoform_id, list_sj, -10, hash_pos, simulant_sj, group_sj )

##EXPERIMENT 1 - create transcript just with single transcript
# transcript_1 = iso_sj.build_transcript_around_sj( obj_sj )
# list_transcripts = iso_sj.reconstruct_transcript_single_sj( obj_sj )


# #if want to reconstruct transcript the "old-fashioned" way
# simulant_sj = True
# iso_sj = IsoformSJ( isoform_id, list_sj, -10, hash_pos, simulant_sj, group_sj )
# canon_only = False
# max_penalty = 5
# list_transcripts = iso_sj.reconstruct_transcript( canon_only, max_penalty )


# ##TEST:: see each SJ that makes up transcript
# # print "Real transcript:"
# # for i, each_transcript in enumerate( list_transcripts ):
# #     print "transcript ", i, " - ", len( each_transcript )
# #     for i2, each_sj in enumerate( each_transcript ):
# #         print "real ", i2, " - ", each_sj




##EXPERIMENT 2 - use a different version of single sj transcript using IsoformSJ.reconstruct_transcript_single_sj_v2()
print "reconstruct transcript from single SJ..."
transcript_ssj = iso_sj.reconstruct_transcript_single_sj_v2( obj_sj )
# obj_tt = TranslateTranscript( transcript_ssj, iso_sj, DIR_GENOME, {} )
print "Test V2 - Transcript made with single SJ:"
for i, each_sj in enumerate( transcript_ssj ):
    # #display SJ & ligated exons
    # ligated_exons = each_sj.spliced_elems( isoform_id )    
    # print "real ", i, " - ", each_sj, " & exon 1 = ", ligated_exons[0], " & exon 2 = ", ligated_exons[1]
    print "real ", i, " - ", each_sj

canon_transcript = iso_sj.create_canon_transcript( False )
print "canonical SJ: "
for i, canon_sj in enumerate( canon_transcript ):
    print "canon ", i, " - ", canon_sj


##EXPERIMENT 3 - determine if making an IsoformSJ instance is possible
# gene_db = 1     #for gene_db, 1 = uses Isoform.obj_cruzdb.refGene, 2 = uses Isoform.obj_cruzdb.ensGene
# obj_poss = IsoformSJ.is_obj_possible( isoform_id, sj_range, gene_db )

# print "obj_poss = ", obj_poss

elapse_time = time.time() - start_time
mokhaPy.timeElapse_convertToHMS( elapse_time )
print "------------ TDD Completed: 170322_BuildTranscriptFromSingleSJ.py ------------"