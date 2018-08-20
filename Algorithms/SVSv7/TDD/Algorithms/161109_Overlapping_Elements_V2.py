#/usr/bin/python

import os
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
DIR_TABIX = '/home/mokha/Documents/Krauthammer_Lab/160427_SJFalsePositives/Data/Tabix'
DIR_THRES = DIR_PROJ + "/160803_ReadCountToThreshold/Results"

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
        sj_range = each_sj.chrom + ':' + str(each_sj.exonPos.location.start) + '-' + str(each_sj.exonPos.location.end)        #NOTE: this SJ range is actually 
        overlap_count = obj_sj.sj_read_support( bam_reader, sj_range )
        list_sj_count.append( overlap_count['unique_count'] )

    return list_sj_count

def count_list_sj_v2( list_sj, obj_sj, bam_reader ):
    """
    Args:
        list_sj = array where each element is an Exon instance (in this case they are introns because I am looking for read counts mapping to SJ)
        obj_sj = instance of SpliceJunction, where the position is the aberrant SJ
        bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
    Function: returns a list of unique counts for each SJ in array "list_sj"
    """
    list_sj_count = {"all_count": [], "unique_count": [], "unique_count_50": []}      #records all read counts for all SJs in "list_sj"
    for each_sj in list_sj:
        sj_range = each_sj.chrom + ':' + str(each_sj.exonPos.location.start) + '-' + str(each_sj.exonPos.location.end)        #NOTE: this SJ range is actually 
        
        #assign all different types of reads 
        overlap_count = obj_sj.sj_read_support( bam_reader, sj_range )
        for k,v in overlap_count.iteritems():       #k = read count categories: "all_count", "unique_count", "unique_count_50", v = the read count for each category
            list_sj_count[k].append( v )

    return list_sj_count


def comparison_1( obj_sj, bam_reader ):
    """
    Args:
        obj_sj = instance of SpliceJunction, where the position is the aberrant SJ
        bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
    Function: compares the aberrant SJ read count to the overlapping canonical SJs, and outputs the ratio (aberrant SJ read count) / (max read count of overlapping, canonical SJs)
    
    Types of Comparisons to prove an aberrant SJ is real
    -Comparison 1: aberrant SJ to canonical SJ within the sample - see if the aberrant SJ is higher than the canonical SJs
    """
    #quantify number of unique reads that maps to aberrant SJ
    sj_count_aberrant = obj_sj.sj_read_support( bam_reader )

    #find all overlapping SJs to the aberrant SJ
    hash_overlap_sj = obj_sj.canonical_overlapping_sj( True )
    list_overlap_sj_count = count_list_sj_v2( hash_overlap_sj.values(), obj_sj, bam_reader )
    overlap_count_uniq = list_overlap_sj_count['unique_count'][:]
    #check if no read counts recorded for overlapping canonical
    if not overlap_count_uniq:
        overlap_count_uniq = [1]

    ##TEST::
    print "compare_1: aberrant_sj = ", obj_sj, " & isoform_aberrants = ", obj_sj.isoform_aberrants
    print "compare_1: sj_count_aberrant = ", sj_count_aberrant
    print "compare_1: hash_overlap_sj = ", hash_overlap_sj
    print "compare_1: overlap_count_uniq = ", overlap_count_uniq

    # aberrant_canon_ratio = float( sj_count_aberrant['unique_count'] ) / max( overlap_count_uniq )
    # aberrant_canon_ratio = float( sj_count_aberrant['unique_count'] ) / np.median( overlap_count_uniq )
    aberrant_canon_ratio = float( sj_count_aberrant['unique_count'] ) / np.mean( overlap_count_uniq )
    hash_c1 = {'aberrant_canon_ratio': aberrant_canon_ratio, 'sample_aberrant_sj': sj_count_aberrant['unique_count'], 'overlap_count_uniq': overlap_count_uniq, 'overlap_count_all': list_overlap_sj_count['all_count'] }

    return hash_c1


def comparison_2( obj_sj, bam_reader, gene_sym ):
    """
    Args:
        obj_sj = instance of SpliceJunction, where the position is the aberrant SJ
        bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
        gene_sym = string that is the gene symbol -> used to determine which elements (introns) are constitutive for that gene.
    Function: compares the overlapping canonical SJ & the aberrant SJ to the constitutive SJs in the gene

    Types of Comparisons to prove an aberrant SJ is real
    -Comparison 2: overlapping canonical SJ/constitutive SJs in sample -> compare to other samples (make sure other samples do not have aberrant SJ in that gene, else do not consider)
    -Comparison 2b: same as Comparison 2, but do this with aberrant SJ/constitutive SJs
    """

    #quantify number of unique reads that maps to aberrant SJ
    sj_count_aberrant = obj_sj.sj_read_support( bam_reader )

    #find all overlapping SJs to the aberrant SJ 
    hash_overlap_sj = obj_sj.canonical_overlapping_sj( True )
    list_overlap_sj_count = count_list_sj( hash_overlap_sj.values(), obj_sj, bam_reader )
    #check if no read counts recorded for overlapping canonical
    if not list_overlap_sj_count:
        list_overlap_sj_count = [0]     #0 because it will be in the numerator

    #find all constitutive SJs in the gene (these are basically the constitutive introns)
    constitutive_sj = obj_sj.find_constitutive_element( gene_sym, True )
    list_constitutive_sj_count = count_list_sj( constitutive_sj, obj_sj, bam_reader )
    #check if no read counts recorded for constitutive SJs (via constitutive introns)
    if not list_constitutive_sj_count:
        list_constitutive_sj_count = [1]        #1 because it will be in the numerator

    ratio_aberrant = float( sj_count_aberrant['unique_count'] ) / np.mean( list_constitutive_sj_count )
    # ratio_overlap_canonical = float( max( list_overlap_sj_count ) ) / np.median( list_constitutive_sj_count )
    # ratio_overlap_canonical = float( np.median( list_overlap_sj_count ) ) / np.median( list_constitutive_sj_count )
    ratio_overlap_canonical = float( np.mean( list_overlap_sj_count ) ) / np.mean( list_constitutive_sj_count )
    return { 'ratio_aberrant': ratio_aberrant, 'ratio_overlap_canonical': ratio_overlap_canonical, 'list_constitutive_sj_count': list_constitutive_sj_count }


def comparison_3( obj_sj, bam_reader, library_size, gene_sym ):
    """
    Args:
        obj_sj = instance of SpliceJunction, where the position is the aberrant SJ
        bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
        gene_sym = string that is the gene symbol -> used to determine which elements (introns) are constitutive for that gene.
        library_size = integer that is the number of mapped reads in the sample
    Function: compares exons that are affected by aberrant SJ (exon skip - full or partial) to constitutive exons in sample - Ratio = skipped exons/constitutive exons for sample affected vs. samples that are not affected

    Types of Comparisons to prove an aberrant SJ is real
    -Comparison 3: compare exons that are affected by aberrant SJ (exon skip - full or partial) to constitutive exons in sample - Ratio = skipped exons/constitutive exons for sample affected vs. samples that are not affected
    """
    #get skipped exons for each isoform -> calculate RPKM for those isoforms & take average -> get constitutive exons & divide

    #SAVE FOR LATER: this is just another way to find longest list of exons skipped
    # isoforms_exon_skips = [ v['exon_skip'] for k,v in obj_sj.isoform_aberrants.iteritems() ]
    # num_es_max = max( isoforms_exon_skips, key = lambda k: len(k) )     #get the maximum number of exons
    
    #METHOD 1 (get skipped exons): get the isoform with the maximum number of exon skips
    # isoform_max_es = max( obj_sj.isoform_aberrants, key = lambda k: len(obj_sj.isoform_aberrants[k]['exon_skip']) )
    # exons_skipped_obj = []      #records all positions of skipped exons
    # for x in obj_sj.isoform_aberrants[isoform_max_es]['exon_skip']:
    #     exons_skipped_obj.append( obj_sj.hash_isoforms[isoform_max_es].get_exon_num(x) )

    #get the expression of the constitutive exons for the gene 'gene_sym'
    constitutive_exons_express = get_constitutive_exon_expression( obj_sj, bam_reader, library_size, gene_sym )
    
    #METHOD 2 (get skipped exons): get all skipped exons, full skips and partial skips
    exons_skipped_obj = obj_sj.missing_exons()
    if not exons_skipped_obj:
        return {'ratio_skip_constitutive': [], 'exons_skipped_express': [], 'constitutive_exons_express': constitutive_exons_express}

    #calculate expression of skipped exons
    exons_skipped_express = all_exons_rpkm( exons_skipped_obj, bam_reader, library_size )

    # exon_skip_constitutive = float( max(exons_skipped_express) ) / np.median( constitutive_exons_express )
    # ratio_skip_constitutive = float( np.median(exons_skipped_express) ) / np.median( constitutive_exons_express )
    ratio_skip_constitutive = float( np.mean(exons_skipped_express) ) / np.mean( constitutive_exons_express )
    hash_c3 = {'ratio_skip_constitutive': ratio_skip_constitutive, 'exons_skipped_express': exons_skipped_express, 'constitutive_exons_express': constitutive_exons_express}

    return hash_c3


def comparison_4( obj_sj, bam_reader, library_size, gene_sym ):
    """
    Args:
        obj_sj = instance of SpliceJunction, where the position is the aberrant SJ
        bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
        gene_sym = string that is the gene symbol -> used to determine which elements (introns) are constitutive for that gene.
        library_size = integer that is the number of mapped reads in the sample
    Function: determines if the intron is retained in the aberrant SJ. Get introns that are retained -> calculate expression of that intronic region -> compare to constitutive exons
    """
    #get the expression of the constitutive exons for the gene 'gene_sym'
    constitutive_exons_express = get_constitutive_exon_expression( obj_sj, bam_reader, library_size, gene_sym )

    #get the potentially retained introns
    hash_introns_retained = obj_sj.intron_retained()
    if not hash_introns_retained:
        return {'ric_1': None, 'ric_2': None, 'intron_retained_1': None, 'intron_retained_2': None, 'constitutive_exons_express': constitutive_exons_express }

    #see if the start position of the SJ is intronic
    if hash_introns_retained['retained_intron_1']:
        express_intron_1 = all_exons_rpkm( hash_introns_retained['retained_intron_1'], bam_reader, library_size )
        # ric_1 = float( express_intron_1[0] ) / np.median( constitutive_exons_express )      #ric = ratio_intron_constitutive
        ric_1 = float( express_intron_1[0] ) / np.mean( constitutive_exons_express )      #ric = ratio_intron_constitutive
    else:
        express_intron_1 = None
        ric_1 = None        #ric = ratio_intron_constitutive

    #see if the end position of the SJ is intronic
    if hash_introns_retained['retained_intron_2']:
        express_intron_2 = all_exons_rpkm( hash_introns_retained['retained_intron_2'], bam_reader, library_size )
        # ric_2 = float( express_intron_2[0] ) / np.median( constitutive_exons_express )      #ric = ratio_intron_constitutive
        ric_2 = float( express_intron_2[0] ) / np.mean( constitutive_exons_express )      #ric = ratio_intron_constitutive
    else:
        express_intron_2 = None
        ric_2 = None    #ric = ratio_intron_constitutive

    #return intronic retention results
    intron_retain_info = {'ric_1': ric_1, 'ric_2': ric_2, 'intron_retained_1': express_intron_1[0], 'intron_retained_2': express_intron_2[0], 'constitutive_exons_express': constitutive_exons_express }
    return intron_retain_info


print "------------ TDD: 161109_Overlapping_Elements_V2.py ------------"
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
isoform_info = g.refGene.filter_by( name2 = gene_sym ).all()
isoform_id = isoform_info[0].name

sj_pos = 'chr6:54025230-54034325'
hash_pos = Isoform.split_genome_pos( sj_pos )
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
obj_sj = SpliceJunction( sj_id, chrom, start, end, strand, read_count, gene_sym, sample_prevalence, control_prevalence, bool_intronic )

#get all samples
for (dirpath, dirnames, filenames) in os.walk( DIR_TABIX ):
    break
all_tabix = [x for x in filenames if not '.tbi' in x and sample_name in x]
#create IsoformSJ instance - find all aberrant SJ for gene
path_tabix = dirpath + '/' + all_tabix[0]
print "sample tabix = ", all_tabix[0], " || path_tabix = ", path_tabix
iso_sj = IsoformSJ( isoform_id, path_tabix, 0 )

#show all aberrant SJ
for i, s in enumerate( iso_sj.list_sj ):
    print "sj ", i, " - ", s, " & aberrant = ", s.isoform_aberrants

#Comparison 1
c1 = comparison_1( obj_sj, bam_reader )

print "c1 = ", c1

# #get gene info
# print "all isoforms = ", obj_mi.hash_isoforms.keys(), " & num iso = ", len( obj_mi.hash_isoforms.keys() )

# #get overlapping elements
# hash_overlapping_elem = obj_mi.find_overlapping_elements( sj_pos, gene_sym, True, True )

# print "Elements overlapping ", sj_pos, ": "
# for k,v in hash_overlapping_elem.iteritems():       #key = genomic range of element & value = Exon object
#     print "k = ", k, " & v = ", v, " & Exon isoforms = ", v.arrAllIsoforms, " & num iso = ", len( v.arrAllIsoforms )

# print "--- Unique read count for Overlapping SJ using HTSeq ---"
# hash_overlap_sj = obj_sj.canonical_overlapping_sj( True )
# overlap_sj_count = count_list_sj( hash_overlap_sj.values(), obj_sj, bam_reader )
# print "overlap SJ count = ", overlap_sj_count

# #Test exon skipping, full & partial
# print "SJ exon skipping, full & partial - ", obj_sj.missing_exons()

# #Test intronic retention]
# print "any intronic retention? - ", obj_sj.intron_retained()

# #Comparison 1
# c1 = comparison_1( obj_sj, bam_reader )
# #Comparison 2
# c2 = comparison_2( obj_sj, bam_reader, gene_sym )
# #Comparison 3
# c3 = comparison_3( obj_sj, bam_reader, library_size, gene_sym )
# #Comparison 4
# c4 = comparison_4( obj_sj, bam_reader, library_size, gene_sym )


print "------------ TDD Completed: 161109_Overlapping_Elements_V2.py ------------"