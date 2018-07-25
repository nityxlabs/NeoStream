#/usr/bin/python
import sys

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv5 import SpliceJunction, Isoform, MultiIsoform, IsoformSJ, TranscribeTranscript, TranslateTranscript

#Constants - directories
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_GENOME = DIR_PROJ + '/ArchiveData/hg19.fa'

def retrieve_aa_around_sj( obj_tt, pos_aberrant, bool_before ):
    """
    Args:
        obj_tt = instance of TranscribeTranscript instance
        pos_aberrant = the 5' end of the aberrant SJ event (for + strand, the lower genomic position, for - strand, the higher genomic position) 
    Function: retrieves information about the amino acid sequence before the 5' end of the aberrant SJ event
    """
    #retrieve the amino acid sequence before mutation
    num_aa = 14      #number of amino acids to retrieve before aberrant SJ event
    hash_range_aa_around_sj = obj_tt.find_aa_pos_surrounding( pos_aberrant, num_aa, bool_before )

    ##TEST::
    print "hash_range_aa_around_sj = ", hash_range_aa_around_sj

    if not hash_range_aa_around_sj:
        return {'genome_range': '-', 'nuc_seq': '-', 'prot_seq': '-', 'len': '-'}

    #retrieve nucleotide sequence for aberrant SJ, NMD-sensitive region, and NMD-irrelevant region
    nuc_seq_aa_around_sj = obj_tt.retrieve_nuc_seq( hash_range_aa_around_sj['i_genome_start'], hash_range_aa_around_sj['i_genome_end'], True )

    #sequence for amino acids before aberrant SJ
    protein_aa_around_sj = Seq( nuc_seq_aa_around_sj ).translate( to_stop = False ) if nuc_seq_aa_around_sj else '-'
    len_protein_aa_around_sj = len( protein_aa_around_sj )

    return {'genome_range': hash_range_aa_around_sj['str_genome_range'],'nuc_seq': nuc_seq_aa_around_sj, 'prot_seq': protein_aa_around_sj, 'len': len( protein_aa_around_sj )}

print "------------ TDD: 170411_AASeq_Mutation.py ------------"
"""
Algorithm: this will test retrieve amino acids before and after a mutation position
Functions to be tested:
    -TranscribeTranscript
        -get_mutated_codon()
"""
#mutation 1 - minus strand
# snv_genome_pos = "chr8:41573229-41573229"
# gene_sym = 'ANK1'
# #NOTE: for mutation, assume it is on the + strand, regardless if the gene is on the + or - strand -> I take care of the strand difference later
# snv_strand = 1      #column 'variant genotype' reports mutation on + strand only, regardless of what strand the mutation or gene is on
# base_orig = 'G'
# base_mut = 'A'

#mutation 2 - plus strand
snv_genome_pos = "chr2:189872304-189872304"
gene_sym = 'COL3A1'
#NOTE: for mutation, assume it is on the + strand, regardless if the gene is on the + or - strand -> I take care of the strand difference later
snv_strand = 1      #column 'variant genotype' reports mutation on + strand only, regardless of what strand the mutation or gene is on
base_orig = 'C'
base_mut = 'T'


#mutation 3
# snv_genome_pos = "chr7:151704932-151704932"
# gene_sym = "GALNTL5"
# snv_strand = 1
# base_orig = 'G'
# base_mut = 'A'

g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )

hash_snv_pos = Isoform.split_genome_pos( snv_genome_pos )
obj_mi = MultiIsoform( hash_snv_pos['chrom'], hash_snv_pos['start'], hash_snv_pos['end'], gene_sym )
hash_pos = { 'chrom': hash_snv_pos['chrom'], 'pos_oi': hash_snv_pos['start'] }      #use this for IsoformSJ - as an isoform may contain many "versions" (one isoform actually still has multiple isoforms funny enough), find the isoform closest to this position
for i3, (k3,v3) in enumerate( obj_mi.hash_isoforms.iteritems() ):     #k3 = isoform ID, v3 = Isoform Instance

    if i3 > 0:
        break

    #create the canonical transcript
    iso_sj = IsoformSJ( k3, [], -10, hash_pos, True )
    canon_transcript = iso_sj.create_canon_transcript( False )
    obj_tt = TranslateTranscript( canon_transcript, iso_sj, DIR_GENOME, {} )

    #retrieve current mutation
    # varied_codon_info = obj_tt.get_mutated_codon( base_orig, base_mut, snv_genome_pos, snv_strand ) 

    # num_aa = 1
    # # nuc_seq_aa_before_sj = obj_tt.find_aa_pos_surrounding( hash_snv_pos['start'], num_aa, True )
    # nuc_seq_aa_after_sj = obj_tt.find_aa_pos_surrounding( hash_snv_pos['start'], num_aa, False )

    # print "snv position = ", hash_snv_pos
    # # print "varied_codon_info = ", varied_codon_info
    # # print "nuc_seq_aa_before_sj = ", nuc_seq_aa_before_sj
    # print "nuc_seq_aa_after_sj = ", nuc_seq_aa_after_sj

    #create mutation in transcript
    try:
        i_genome_pos = obj_tt.arr_genome_pos.index( hash_snv_pos['start'] )
    except ValueError:
        print "Could not find position in gene = ", hash_snv_pos['start']
        continue

    #retrieve mutated codon & resulting amino acid
    # snv_strand = 1 if row['strand'] == '+' else -1        #column 'HGVSc' report the mutation that is strand-corrected (i.e. finds base complement if on other strand)
    varied_codon_info = obj_tt.get_mutated_codon( base_orig, base_mut, snv_genome_pos, snv_strand )     #hash that contains both the original & mutated codon

    genome_range_mut = varied_codon_info['str_genome_range']
    # nuc_seq_aa_mut = obj_tt.arr_nuc_seq[i_genome_pos-1:i_genome_pos+1][::-1] if iso_sj.strand < 0 else obj_tt.arr_nuc_seq[i_genome_pos-1:i_genome_pos+1]
    codon_seq_aa_orig = varied_codon_info['codon_orig']
    protein_aa_orig = varied_codon_info['aa_orig']
    codon_seq_aa_mut = varied_codon_info['codon_mut']
    protein_aa_mut = varied_codon_info['aa_mut']
    len_protein_aa_mut = len( protein_aa_mut )

    ##TEST:: see the mutated amino acid
    print i3, ": mutated codon pos = ", genome_range_mut, " & codon_orig = ", codon_seq_aa_orig, " codon_mut = ", codon_seq_aa_mut, " & AA change = ", protein_aa_orig, " to ", protein_aa_mut
    print i3, "-----------------------\n"


    #retrieve AA sequence before the mutation (towards the 5' end of the gene -> for + strand this is the lower genomic position, for - strand this is the higher genomic position)
    hash_aa_before_mut = retrieve_aa_around_sj( obj_tt, hash_snv_pos['start'], True )
    genome_range_before_mut = hash_aa_before_mut['genome_range']
    nuc_seq_aa_before_mut = hash_aa_before_mut['nuc_seq']
    protein_aa_before_mut = hash_aa_before_mut['prot_seq']
    len_protein_aa_before_mut = hash_aa_before_mut['len']

    ##TEST::
    print "AA BEFORE MUT:"
    print "BEFORE genome pos = ", genome_range_before_mut
    print "BEFORE Prot Seq = ", protein_aa_before_mut, " & len = ", len_protein_aa_before_mut

    #retrieve AA sequence before the mutation (towards the 3' end of the gene -> for + strand this is the higher genomic position, for - strand this is the lower genomic position)
    hash_aa_after_mut = retrieve_aa_around_sj( obj_tt, hash_snv_pos['start'], False )
    genome_range_after_mut = hash_aa_after_mut['genome_range']
    nuc_seq_aa_after_mut = hash_aa_after_mut['nuc_seq']
    protein_aa_after_mut = hash_aa_after_mut['prot_seq']
    len_protein_aa_after_mut = hash_aa_after_mut['len']

    ##TEST::
    print "AA AFTER MUT:"
    print "AFTER genome pos = ", genome_range_after_mut
    print "AFTER Prot Seq = ", protein_aa_after_mut, " & len = ", len_protein_aa_after_mut

print "------------ TDD Completed: 170411_AASeq_Mutation.py ------------"