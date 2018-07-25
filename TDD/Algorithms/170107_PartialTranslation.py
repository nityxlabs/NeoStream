#/usr/bin/python
import os
import sys

import numpy as np

from cruzdb import Genome
import HTSeq

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv5 import Isoform, SpliceJunction, IsoformSJ, SimpleTranscribe, SimpleTranslate, TranscribeTranscript, TranslateTranscript

#Constants: file paths
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv5"
DIR_DATA = DIR_CURR + "/TestData"
DIR_RESULTS = DIR_CURR + "/TestResults"
DIR_GENOME = DIR_PROJ + "/ArchiveData/hg19.fa"

print "------------ TDD: 170107_PartialTranslation.py ------------"
"""
Algorithm: this algorithm is meant to test the idea of identifying early stop codon, truncated protein, mutated amino acids based on the fact I know the position & reading frame - this is to test the idea and so far it works!!!
"""

g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )

gene_name = 'BRAF'
obj_gene = Isoform.obj_cruzdb.refGene.filter_by( name2 = gene_name ).first()

print "obj_gene.introns: "
print obj_gene.introns


# iso_sj = IsoformSJ( obj_gene.name, [], -2, None, True, 0 )

# print "SJs for ", gene_name
# for i, each_sj in enumerate( iso_sj.list_sj ):
#     print i, ": ", each_sj


#EXPERIMENT: retrieve the genomic position & the reading frame
strand = -1
str_a = "chr7:140534651-140534661*chr7:140534661-140534672"

obj_tt = SimpleTranslate( str_a, strand, DIR_GENOME )

print "show transcript sequence: "
print obj_tt.obj_seq
print obj_tt.obj_seq.reverse_complement()

print "show protein seq: "
print obj_tt.translate_seq()


# hash_transcript_notes = {140534657: 'some_note'}
# obj_tt2 = TranscribeTranscript( str_a, strand, DIR_GENOME, {} )
# print "show hash_seq: "
# print obj_tt2.hash_seq_notes

#EXPERIMENT: retrieve where the frameshift occurred, find previous nucleotides, and then translate
strand = -1
str_a = "chr7:140534651-140534661*chr7:140534661-140534672"
str_sj_canon = "chr7:140481493-140482820*chr7:140482957-140487347"
str_sj_aberrant = "chr7:140481493-140482820*chr7:140482958-140487347"
frameshift_pos = "chr7:140482820-140482958"


#EXPERIMENT: build transcript around a transcript
sample_prevalence = 0
control_prevalence = 0
bool_intronic = True
# sj_canon_1 = { "sj_id": "SJ_1", "chrom": "chr7", "start": 140454033, "end": 140476711, "strand": '-', "gene_sym": 'BRAF' }
sj_canon_1 = { "sj_id": "SJ_1", "chrom": "chr7", "start": 140550011, "end": 140624365, "strand": '-', "gene_sym": 'BRAF' }
obj_sj = SpliceJunction( sj_canon_1['sj_id'], sj_canon_1['chrom'], sj_canon_1['start'], sj_canon_1['end'], sj_canon_1['strand'], 12345, sj_canon_1['gene_sym'], sample_prevalence, control_prevalence, bool_intronic )
sj_isoform_id = obj_sj.assigned_isoform[0]

list_sj = [obj_sj]
iso_sj = IsoformSJ( sj_isoform_id, list_sj, -10, None, True, 0 )
transcript_sj = iso_sj.build_transcript_around_sj( obj_sj )

print "gene_sym = ", obj_sj.gene_sym_oi, " | isoform ID = ", sj_isoform_id
for i, each_sj in enumerate( transcript_sj ):
    tuple_spliced_exons = each_sj.spliced_elems( sj_isoform_id )
    print i, ": ", each_sj, " & connected_exons = exon ", tuple_spliced_exons[0].exonNum, " to exon ", tuple_spliced_exons[1].exonNum


