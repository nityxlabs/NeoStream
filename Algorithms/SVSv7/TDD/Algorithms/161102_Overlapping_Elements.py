#/usr/bin/python

import sys

from cruzdb import Genome
import HTSeq

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv5 import IsoformSJ, Isoform, MultiIsoform, SpliceJunction

#Constants: file paths
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv5"
DIR_DATA = DIR_CURR + "/TestData"
DIR_RESULTS = DIR_CURR + "/TestResults"
#get mapped reads
DIR_RNASEQ = DIR_PROJ + "/150802_TophatSamples"

def count_list_sj( list_sj, obj_sj, bam_reader ):
    """
    Args:
        list_sj = array where each element is an Exon instance (in this case they are introns because I am looking for read counts mapping to SJ)
        obj_sj = instance of SpliceJunction, where the position is the aberrant SJ
        bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
    Function: returns a list of unique counts for each SJ in array "list_sj"
    """
    list_sj_count = []      #records all read counts for all SJs in "list_sj"
    for each_sj in list_sj:
        sj_range = each_sj.chrom + ':' + str(each_sj.start) + '-' + str(each_sj.end)        #NOTE: this SJ range is actually 
        overlap_count = obj_sj.sj_read_support( bam_reader, sj_range )
        list_sj_count.append( overlap_count['unique_count'] )

    return list_sj_count

print "------------ TDD: 161102_Overlapping_Elements.py ------------"
"""
Algorithm: this will test the functions in MultiIsoform -> find_overlapping_elements(), 
"""
g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )

sample_name = 'yuhimo'
path_bam = DIR_RNASEQ + '/tophat_sample_' + sample_name + '/accepted_hits.bam'
bam_reader = HTSeq.BAM_Reader( path_bam )

#find overlapping splice junctions
# gene_sym = 'CDK11B'
# sj_pos = 'chr1:1573952-1575753'
gene_sym = 'MLIP'
sj_pos = 'chr6:54025230-54034325'
hash_pos = MultiIsoform.split_genome_pos( sj_pos )
chrom = hash_pos['chrom']
start = hash_pos['start']
end = hash_pos['end']

#other parameters that aren't that important for this testing
strand = '+'
sj_id = "TEST"
read_count = 0
sample_prevalence = 0
control_prevalence = 0
bool_intronic = True

obj_mi = MultiIsoform( chrom, start, end, gene_sym )
obj_sj = SpliceJunction( sj_id, chrom, start, end, strand, read_count, sample_prevalence, control_prevalence, bool_intronic )

#get gene info
print "all isoforms = ", obj_mi.hash_isoforms.keys(), " & num iso = ", len( obj_mi.hash_isoforms.keys() )

#get overlapping elements
hash_overlapping_elem = obj_mi.find_overlapping_elements( sj_pos, gene_sym, True, True )

print "Elements overlapping ", sj_pos, ": "
for k,v in hash_overlapping_elem.iteritems():       #key = genomic range of element & value = Exon object
    print "k = ", k, " & v = ", v, " & Exon isoforms = ", v.arrAllIsoforms, " & num iso = ", len( v.arrAllIsoforms )

print "--- Unique read count for Overlapping SJ using HTSeq ---"
hash_overlap_sj = obj_sj.canonical_overlapping_sj( True )
overlap_sj_count = count_list_sj( hash_overlap_sj.values(), obj_sj, bam_reader )
print "overlap SJ count = ", overlap_sj_count

#Test exon skipping, full & partial
print "SJ exon skipping, full & partial - ", obj_sj.missing_exons()

#Test intronic retention]
print "any intronic retention? - ", obj_sj.intron_retained()


print "------------ TDD Completed: 161102_Overlapping_Elements.py ------------"