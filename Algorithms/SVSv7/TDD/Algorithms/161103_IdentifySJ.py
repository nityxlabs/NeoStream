#/usr/bin/python

import sys

from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv5 import MultiIsoform, Isoform, SpliceJunction

#Constants: file paths
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv5"
DIR_DATA = DIR_CURR + "/TestData"
DIR_RESULTS = DIR_CURR + "/TestResults"
#get mapped reads
DIR_RNASEQ = DIR_PROJ + "/150802_TophatSamples"

print "------------ TDD: 161103_IdentifySJ.py ------------"
"""
Algorithm: this is meant to test SpliceJunction class to see if can determine a canonical SJ from an aberrant SJ
"""

g = Genome('sqlite:////tmp/hg19_v2.db')
Isoform.set_cruzdb( g )

# test_obj = Isoform.obj_cruzdb.knownToRefSeq.filter_by( value = 'NM_005112' ).first()

# print "test_obj = ", test_obj.name


# sj_pos = 'chr6:32410470-32410961'
# gene_sym = 'HLA-DRA'
# sj_pos = 'chr4:10084802-10086066'
# gene_sym = 'WDR1'
sj_pos = 'chr7:140500428-140501361'
gene_sym = 'BRAF'
isoform_id = 'NM_004333'

hash_sj = MultiIsoform.split_genome_pos( sj_pos )
chrom = hash_sj['chrom']
start = hash_sj['start']
end = hash_sj['end']

#other parameters that aren't that important for this testing
strand = '+'
sj_id = "TEST"
read_count = 0
sample_prevalence = 0
control_prevalence = 0
bool_intronic = True

obj_sj = SpliceJunction( sj_id, chrom, start, end, strand, read_count, gene_sym, sample_prevalence, control_prevalence, bool_intronic )

print "obj_sj = ", obj_sj
print "canonical = ", obj_sj.canon
print "assigned isoforms = ", obj_sj.assigned_isoform
print "isoform_aberrants = ", obj_sj.isoform_aberrants

tuple_exons = obj_sj.spliced_elems( isoform_id, False )
print "exons ligated by SJ = ", tuple_exons
hash_sj_exon_info = obj_sj.spliced_elems_position_range( isoform_id, False, True )
print "show info about exons ligated by SJ: "
print "range of exons = ", hash_sj_exon_info['str_sj_pos']
print "prev_elem = ", hash_sj_exon_info['prev_elem']
print "next_elem = ", hash_sj_exon_info['next_elem']
print "prev_exon_canon = ", hash_sj_exon_info['prev_exon_canon']
print "next_exon_canon = ", hash_sj_exon_info['next_exon_canon']



##TEST::
# isoform_id = 'NR_024540'
hash_pos = {'chrom': 'chr1', 'pos_oi': 14829}
obj_iso = Isoform( isoform_id, hash_pos )


print "UCSC gene = ", obj_iso.get_ucsc_gene_name( isoform_id )
print "Ensembl transcript ID = ", obj_iso.get_ensembl_gene_name( isoform_id )
ensembl_isoforms = obj_iso.get_ensembl_isoforms( isoform_id )
print "Ensembl gene ID = ", ensembl_isoforms

for i, each_isoform in enumerate( ensembl_isoforms ):
    print "Ensembl gene ", i, " | gene name (ENSG) = ", each_isoform.name2, " | isoform (ENST) = ", each_isoform.name, " | start = ", each_isoform.start, " | end = ", each_isoform.end

#test MultiIsoform & Isoform - make sure the correct isoform
# obj_mi = MultiIsoform( chrom, start, end )
# print "sj_pos = ", sj_pos
# for k,v in obj_mi.hash_isoforms.iteritems():        #k = isoform ID, v = Isoform instance
#     print "isoform = ", k
#     for k2,v2 in v.hashExonList.iteritems():
#         print "k2 = ", k2, " & v2 = ", v2



#check the exon positions of the gene of interest
# gene_name = 'WDR1'
# gene = g.refGene.filter_by( name2 = gene_name ).all()
# for i, each_isoform in enumerate( gene ):
#     print i, " - ", each_isoform.name2, ": ", each_isoform.name, " & "
#     for i2, exon in enumerate( each_isoform.exons ):
#         intron = None if ( len( each_isoform.introns ) <= i2 ) else each_isoform.introns[i2]
#         print each_isoform.chrom, " | exon ", i2, ": ", exon, " & intron = ", intron



print "------------ TDD Completed: 161103_Identify_SJ.py ------------"