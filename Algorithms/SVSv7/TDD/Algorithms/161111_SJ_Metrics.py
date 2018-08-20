#/usr/bin/python

import os
import sys

import numpy as np

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
DIR_TABIX = DIR_PROJ + '/160427_SJFalsePositives/Data/Tabix'
DIR_THRES = DIR_PROJ + "/160803_ReadCountToThreshold/Results"


def count_list_sj_v2( list_sj, obj_sj, bam_reader ):
    """
    Args:
        list_sj = array where each element is an Exon instance (in this case they are introns because I am looking for read counts mapping to SJ)
        obj_sj = instance of SpliceJunction, where the position is the aberrant SJ
        bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
    Function: returns a list of unique counts for each SJ in array "list_sj"
    """
    hash_sj_count = {"all_count": [], "unique_count": [], "unique_count_50": []}      #records all read counts for all SJs in "list_sj"
    for each_sj in list_sj:
        sj_range = each_sj.chrom + ':' + str(each_sj.exonPos.location.start) + '-' + str(each_sj.exonPos.location.end)        #NOTE: this SJ range is actually 
        
        #assign all different types of reads 
        overlap_count = obj_sj.sj_read_support( bam_reader, sj_range )
        for k,v in overlap_count.iteritems():       #k = read count categories: "all_count", "unique_count", "unique_count_50", v = the read count for each category
            hash_sj_count[k].append( v )

    return hash_sj_count

def all_exons_rpkm( list_exons, bam_reader, library_size ):
    """
    Args:
        list_exons = array of Exon objects
        bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
        gene_sym = string that is the gene symbol -> used to determine which elements (introns) are constitutive for that gene.
        library_size = integer that is the number of mapped reads in the sample
    Function: returns the RPKM for all exons in array 'list_exons'
    """
    #calculate expression of skipped exons
    list_exons_rpkm = []      #records RPKM for all exons
    for each_exon in list_exons:
        exon_count = Isoform.quant_genes_rpkm( bam_reader, each_exon.str_genomic_pos(), True )
        exon_rpkm = Isoform.calc_read_density( exon_count, each_exon.str_genomic_pos(), library_size )
        list_exons_rpkm.append( exon_rpkm )

    return list_exons_rpkm

def each_exon_rpkm( list_exons, bam_reader, library_size ):
    """
    Args:
        list_exons = array of Exon objects
        bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
        gene_sym = string that is the gene symbol -> used to determine which elements (introns) are constitutive for that gene.
        library_size = integer that is the number of mapped reads in the sample
    Function: returns the RPKM for all exons in array 'list_exons'
    """
    #calculate expression of skipped exons
    hash_exon_rpkm = {}     #k = string of position (format = chrom:start-end), v = hash that are 2 counts, RPKM from unique reads & all reads
    for each_exon in list_exons:
        exon_count_uniq = Isoform.quant_genes_rpkm( bam_reader, each_exon.str_genomic_pos(), True )
        exon_count_all = Isoform.quant_genes_rpkm( bam_reader, each_exon.str_genomic_pos(), True )

        str_exon_pos = each_exon.str_genomic_pos()
        hash_exon_rpkm[ str_exon_pos ] = {}
        hash_exon_rpkm[ str_exon_pos ]['rpkm_uniq'] = Isoform.calc_read_density( exon_count_uniq, str_exon_pos, library_size )
        hash_exon_rpkm[ str_exon_pos ]['rpkm_all'] = Isoform.calc_read_density( exon_count_all, str_exon_pos, library_size )

    return hash_exon_rpkm


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
    sj_overlapped = count_list_sj_v2( hash_overlap_sj.values(), obj_sj, bam_reader )

    #check if no read counts recorded for overlapping canonical
    # if not sj_overlapped['unique_count']:
    #     sj_overlapped['unique_count'] = [1]

    # if not sj_overlapped['all_count']:
    #     sj_overlapped['all_count'] = [1]

    #make sure the overlapped SJ count is not just 0
    if not sum( sj_overlapped['unique_count'] ) > 0:
        sj_overlapped['unique_count'] = [1]
    #make sure the overlapped SJ count is not just 0
    if not sum( sj_overlapped['all_count'] ) > 0:
        sj_overlapped['all_count'] = [1]


    # aberrant_canon_ratio = float( sj_count_aberrant['unique_count'] ) / max( sj_overlapped )
    # aberrant_canon_ratio = float( sj_count_aberrant['unique_count'] ) / np.median( sj_overlapped )
    aberrant_canon_ratio_uniq = float( sj_count_aberrant['unique_count'] ) / np.mean( sj_overlapped['unique_count'] )
    aberrant_canon_ratio_all = float( sj_count_aberrant['all_count'] ) / np.mean( sj_overlapped['all_count'] )
    hash_c1 = {
    'aberrant_canon_ratio_uniq': aberrant_canon_ratio_uniq,
    'sample_aberrant_sj_uniq': sj_count_aberrant['unique_count'],
    'sj_overlapped_uniq': sj_overlapped['unique_count'],
    'aberrant_canon_ratio_all': aberrant_canon_ratio_all,
    'sample_aberrant_sj_all': sj_count_aberrant['all_count'],
    'sj_overlapped_all': sj_overlapped['all_count']
    }

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
    #NOTE: don't need to check if overlapping SJs are zeroes because it will not be in the denominator in these calculations
    hash_overlap_sj = obj_sj.canonical_overlapping_sj( True )
    sj_overlapped = count_list_sj_v2( hash_overlap_sj.values(), obj_sj, bam_reader )

   
    #find all constitutive SJs in the gene (these are basically the constitutive introns)
    constitutive_sj = obj_sj.find_constitutive_element( gene_sym, True )        #hash where k = string that is genomic position of element, v = Exon instance
    sj_constitutive = count_list_sj_v2( constitutive_sj.values(), obj_sj, bam_reader )
    #check if no read counts recorded for constitutive SJs (via constitutive introns)

    ##TEST:: quantify expression of each constitutive exon
    c_sj_express = {}               #k = string element pos (chrom:start-end), v = hash that is different types of read counts
    for k,v in constitutive_sj.iteritems():     #k = string element pos (chrom:start-end), v = Exon object
        c_sj_express[k] = count_list_sj_v2( [v], obj_sj, bam_reader )

    alternative_sj = obj_sj.find_alternative_element( gene_sym, True )
    alt_sj_express = {}           #k = string element pos (chrom:start-end), v = hash that is different types of read counts
    for k,v in alternative_sj.iteritems():     #k = string element pos (chrom:start-end), v = Exon object
        alt_sj_express[k] = count_list_sj_v2( [v], obj_sj, bam_reader )

    # make sure the overlapped SJ count is not just 0
    if not sum( sj_constitutive['unique_count'] ) > 0:
        sj_constitutive['unique_count'] = [1]
    #make sure the overlapped SJ count is not just 0
    if not sum( sj_constitutive['all_count'] ) > 0:
        sj_constitutive['all_count'] = [1]

    #calculate the ratio between aberrant SJ/overlapping SJs & constitutive SJs
    ratio_aberrant_uniq = float( sj_count_aberrant['unique_count'] ) / np.mean( sj_constitutive['unique_count'] )
    ratio_aberrant_all = float( sj_count_aberrant['all_count'] ) / np.mean( sj_constitutive['all_count'] )
    # ratio_overlap_canonical = float( max( list_overlap_sj_count ) ) / np.median( list_constitutive_sj_count )
    # ratio_overlap_canonical = float( np.median( list_overlap_sj_count ) ) / np.median( list_constitutive_sj_count )
    ratio_overlap_canonical_uniq = float( np.mean( sj_overlapped['unique_count'] ) ) / np.mean( sj_constitutive['unique_count'] )
    ratio_overlap_canonical_all = float( np.mean( sj_overlapped['all_count'] ) ) / np.mean( sj_constitutive['all_count'] )
    
    return { 
    'ratio_aberrant_uniq': ratio_aberrant_uniq, 
    'ratio_aberrant_all': ratio_aberrant_all, 
    'ratio_overlap_canonical_uniq': ratio_overlap_canonical_uniq,
    'ratio_overlap_canonical_all': ratio_overlap_canonical_all,
    'sj_count_aberrant': sj_count_aberrant,
    'sj_constitutive': sj_constitutive,
    'c_sj_express': c_sj_express,
    'alt_sj_express': alt_sj_express
     }

def get_constitutive_exon_expression( obj_sj, bam_reader, library_size, gene_sym ):
    """
    Args:
        obj_sj = instance of SpliceJunction, where the position is the aberrant SJ
        bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
        gene_sym = string that is the gene symbol -> used to determine which elements (introns) are constitutive for that gene.
        library_size = integer that is the number of mapped reads in the sample
    Function: calculate the expression of all constitutive exons in the gene 'gene_sym'
    """

    #get the constitutive exons for the gene
    constitutive_exons = obj_sj.find_constitutive_element( gene_sym, False )     #hash where k = string that is genomic position of element, v = Exon instance
    constitutive_exons_express = all_exons_rpkm( constitutive_exons.values(), bam_reader, library_size )
    #check if no read counts recorded for constitutive exons
    if not constitutive_exons_express:
        constitutive_exons_express = [1]

    return constitutive_exons_express


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

    #calculate expression of skipped exons
    # hash_exon_skip = {}     #k = string of position (chrom:start-end), v = RPKM expression level
    # for each_exon in exons_skipped_obj:
    #     str_exon_pos = each_exon.str_genomic_pos()
    #     hash_exon_skip[str_exon_pos] = all_exons_rpkm( [each_exon], bam_reader, library_size )

    #Get expression of skipped exons
    skip_exons_rpkm = each_exon_rpkm( exons_skipped_obj, bam_reader, library_size )
    for v in exons_skipped_obj:
        k = v.str_genomic_pos()
        skip_exons_rpkm[k]['exon'] = 'exon ' + str( v.exonNum )

    #Get expression of constitutive exons
    constitutive_exons = obj_sj.find_constitutive_element( gene_sym, False )     #hash where k = string that is genomic position of element, v = Exon instance
    c_exons_rpkm = each_exon_rpkm( constitutive_exons.values(), bam_reader, library_size )
    for k,v in constitutive_exons.iteritems():
        c_exons_rpkm[k]['exon'] = 'exon ' + str( v.exonNum )

    #Get expression of constitutive exons
    alternative_exons = obj_sj.find_alternative_element( gene_sym, False )
    alt_exons_rpkm = each_exon_rpkm( alternative_exons.values(), bam_reader, library_size )
    for k,v in alternative_exons.iteritems():
        alt_exons_rpkm[k]['exon'] = 'exon ' + str( v.exonNum )



    hash_c3 = {
    'ratio_skip_constitutive': ratio_skip_constitutive,
    'exons_skipped_express': exons_skipped_express,
    'constitutive_exons_express': constitutive_exons_express,
    'skip_exons_rpkm': skip_exons_rpkm,
    'c_exons_rpkm': c_exons_rpkm,
    'alt_exons_rpkm': alt_exons_rpkm
    }

    return hash_c3



print "------------ TDD: 161111_SJ_Metrics.py ------------"

g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )

# sample_name = 'yuhimo'
# sample_name = 'yunige'
# sample_name = 'gopik'


#parameters for setting up classes
#example 1
# sample_name = 'yunige'
# gene_sym = 'MLIP'
# sj_pos = 'chr6:54025230-54034325'

#example 2
sample_name = 'yufulo'
gene_sym = 'MITF'
sj_pos = 'chr3:70014008-70014109'

isoform_info = g.refGene.filter_by( name2 = gene_sym ).all()
isoform_id = isoform_info[0].name
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


path_bam = DIR_RNASEQ + '/tophat_sample_' + sample_name + '/accepted_hits.bam'
bam_reader = HTSeq.BAM_Reader( path_bam )

##---------------------------------------------------------------------------------------------
#find all aberrant SJ for a gene
print "\n------ Find all aberrant SJ for gene ------"


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


##---------------------------------------------------------------------------------------------
#find canonical, overlapped elements
print "\n------ Find canonical, overlapped elements ------"

obj_mi = MultiIsoform( chrom, start, end, gene_sym )
obj_isoform = obj_mi.hash_isoforms[isoform_id]
obj_sj = SpliceJunction( sj_id, chrom, start, end, strand, read_count, gene_sym, sample_prevalence, control_prevalence, bool_intronic )
hash_overlap_sj = obj_sj.canonical_overlapping_sj( True )
print "sj of interest = ", sj_pos
# for k,v in obj_isoform.hashIntronList.iteritems():       #key = string that is exon range (chrom:start-end), value = Exon object
#     print "exon ", v.exonNum, ": ", k

for i in range( 1, len(obj_isoform.hashIntronList.values()) + 1 ):
    get_intron = obj_isoform.get_intron_num( i )
    print "intron ", get_intron.exonNum, ": ", get_intron


print "hash_overlap_sj = ", hash_overlap_sj
for k,v in hash_overlap_sj.iteritems():         #key = genomic range of element & value = Exon object
    print "overlapped intron ", v.exonNum, ": ", v


##---------------------------------------------------------------------------------------------
#compare aberrant SJ read count to constitutive SJ
print "\n------ Compare aberrant SJ to overlapped SJ ------"
hash_c1 = comparison_1( obj_sj, bam_reader )
print "hash_c1 = ", hash_c1

##---------------------------------------------------------------------------------------------
#compare aberrant SJ read count to constitutive SJ
print "\n------ Compare aberrant SJ to constitutive SJ ------"
hash_c2 = comparison_2( obj_sj, bam_reader, gene_sym )
for k,v in hash_c2.iteritems():
    if k == 'c_sj_express':
        print "expression of constitutive SJs: "
        for k2,v2 in v.iteritems():
            print "\t", k2, " - ", v2
    elif k == 'alt_sj_express':
        print "expression of Alternative SJs: "
        for k2,v2 in v.iteritems():
            print "\t", k2, " - ", v2
    else:
        print k, " - ", v

##---------------------------------------------------------------------------------------------
#compare exon expression of exons affected by aberrant SJ versus exons not affected by aberrant SJ
print "\n------ Compare exon expression ------"
library_size = Isoform.total_mapped_reads( path_bam )
hash_c3 = comparison_3( obj_sj, bam_reader, library_size, gene_sym )

for k,v in hash_c3.iteritems():
    if k == 'skip_exons_rpkm':
        print "expression of exons skipped: "
        for k2,v2 in v.iteritems():
            print "\t", k2, " - ", v2
    elif k == 'c_exons_rpkm':
        print "expression of Constitutive Exons: "
        for k2,v2 in v.iteritems():
            print "\t", k2, " - ", v2
    elif k == 'alt_exons_rpkm':
        print "expression of Alternative Exons: "
        for k2,v2 in v.iteritems():
            print "\t", k2, " - ", v2
    else:
        print k, " - ", v