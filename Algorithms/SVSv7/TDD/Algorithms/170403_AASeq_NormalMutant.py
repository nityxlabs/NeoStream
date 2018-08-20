#/usr/bin/python
import sys

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv5 import SpliceJunction, Isoform, IsoformSJ, TranscribeTranscript, TranslateTranscript

#Constants - directories
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_GENOME = DIR_PROJ + '/ArchiveData/hg19.fa'

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


def get_index_tp_pux( obj_tt ):         #tp_pux = Truncated Protein Penultimate Exon
    """
    Args:
        -obj_tt = instance of class TranscribeTranscript/TranslateTranscript
    Function: retrieves the position where truncated protein boundary occurs. Returns the index of that boundary position for array obj_tt.arr_genome_pos()
    """

    #get the position 55nt from the 3' end of the penultimate exon
    boundary_tp_pux = 55        #boundary_tp_pux = boundary Truncated Protein Penultimate Exon
    if obj_tt.iso_sj.strand < 0:

        ##TEST::
        print "MAIN: get_index_tp_pux show exons: "
        for i, exons in enumerate( obj_tt.list_exons ):
            print "exon #", i, " - ", exons


        pu_exon_prime_3 = obj_tt.list_exons[1].exonPos.location.start + 1       #for - strand gene, get penultimate exon need to add +1 because of 0-based genomic coordinates
        i_tp_pux = obj_tt.arr_genome_pos.index( pu_exon_prime_3 ) + boundary_tp_pux

        ##TEST::
        # print "GET_I_TP_PUX minus: i_prime3 = ", obj_tt.arr_genome_pos.index( pu_exon_prime_3 ), " & i_tp_pux = ", i_tp_pux, " & len( obj_tt.arr_genome_pos ) = ", len( obj_tt.arr_genome_pos )
        # for i, a in enumerate( obj_tt.arr_genome_pos ):
        #     print "i = ", i, " & a = ", a, " & pu_exon_prime_3 = ", pu_exon_prime_3, " & index of 3' = ", obj_tt.arr_genome_pos.index( pu_exon_prime_3 )
        # for i2, each_sj in enumerate( obj_tt.transcript_sj ):
        #     print "i2 = ", i2, " & each_sj = ", str( each_sj ), " & gene = ", obj_tt.iso_sj.gene_sym

    else:
        pu_exon_prime_3 = obj_tt.list_exons[-2].exonPos.location.end        #for + strand genes, get penultimate exon
        i_tp_pux = obj_tt.arr_genome_pos.index( pu_exon_prime_3 ) - boundary_tp_pux

    """
    if i_tp_pux is out of range, then return None as the index will not be found
        -i_tp_pux < 0 should apply to + strand genes
        -i_tp_pux >= len( obj_tt.arr_genome_pos ) should apply to - strand genes )
    """
    if i_tp_pux < 0 or i_tp_pux >= len( obj_tt.arr_genome_pos ):
        return None

    return i_tp_pux

def find_nearest_rf0( obj_tt, pos_oi ):
    """
    Args:
        -obj_tt = instance of class TranscribeTranscript/TranslateTranscript
        -pos_oi = integer that is the position to find the nearest position that has a reading frame of 0
    Function: finds the nearest 0 reading frame based on the position 'pos_oi' and the strand sign 'obj_tt.iso_sj.strand'
    """
    #retrieve the position of the initial nucleotide in the codon closet the boundary for the TP (truncated protein) boundary
    pos_nearest_start = obj_tt.find_codon_beginning_prev( pos_oi )
    if not pos_nearest_start:
        return (None, None)

    i_nearest_start = obj_tt.arr_genome_pos.index( pos_nearest_start )        #this is the 0-rf closes to the truncate protein position of the penultimate exon

    return ( pos_nearest_start, i_nearest_start )

##check class TranslateTranscript from TranscribeTranslate_V4
def get_nmd_sensitive_region( obj_tt, pos_aberrant ):
    """
    Args:
        -obj_tt = instance of class TranscribeTranscript/TranslateTranscript
        -pos_aberrant = integer that is the position that contains the aberrant position of interest. This will usually be the start of aberrant SJ (lower position for + genes & higher position for the - genes)
    Function: retrieve the region of the gene that, if it contains an early stop codon, will lead to degradation of the transcript. Returns the nucleotide sequence of this region
    """
    #retrieve the position of the initial nucleotide in the codon
    direction = 1 if obj_tt.iso_sj.strand < 0 else -1       #if - strand, then find 0-rf nucleotide at higher position, else if + strand, then find 0-rf nucleotide at lower position
    pos_aberrant_start = obj_tt.find_codon_beginning_prev( pos_aberrant )

    if not pos_aberrant_start:
        return None

    i_aberrant_start = obj_tt.arr_genome_pos.index( pos_aberrant_start )
    # ( pos_aberrant_start, i_aberrant_start ) = find_nearest_rf0( obj_tt, pos_aberrant )

    #get the position 55nt from the 3' end of the penultimate exon
    i_tp_pux = get_index_tp_pux( obj_tt )
    if not i_tp_pux:
        return None

    ##TEST:: get start of codon previous to SJ
    diff = 4
    print "GET_NMD: pos_aberrant_start = ", pos_aberrant_start, " & i_aberrant_start = ", i_aberrant_start, " & i_tp_pux = ", i_tp_pux
    print "GET_NMD: genome_pos_start = ", obj_tt.arr_genome_pos[i_aberrant_start - diff : i_aberrant_start], " & 55 bp upstream = ", obj_tt.arr_genome_pos[i_tp_pux : i_tp_pux + diff]
    print "GET_NMD: nucleotide start = ", obj_tt.arr_nuc_seq[i_aberrant_start - diff : i_aberrant_start ], " & nucleotide 55 bp upstream = ", obj_tt.arr_nuc_seq[i_tp_pux : i_tp_pux + diff]
    print "GET_NMD: start reading frame = ", obj_tt.arr_rf[i_aberrant_start - diff: i_aberrant_start ], " & nucleotide 55 bp upstream reading frame = ", obj_tt.arr_rf[i_tp_pux : i_tp_pux + diff]

    ##TEST::
    # print "NMD_SENS: strand = ", obj_tt.iso_sj.strand, " | i_aberrant_start = ", i_aberrant_start, " | aberr_genome_pos = ", obj_tt.arr_genome_pos[i_aberrant_start]
    # print "NMD_SENS: strand = ", obj_tt.iso_sj.strand, " | i_tp_pux = ", i_tp_pux, " | tp_genome_pos = ", obj_tt.arr_genome_pos[i_tp_pux]

    if obj_tt.iso_sj.strand < 0:
        nmd_sensitive_seq = ''.join( obj_tt.arr_nuc_seq[i_tp_pux : i_aberrant_start + 1] )      #need to add +1 to include last base
        nmd_sensitive_seq = nmd_sensitive_seq[::-1]
    else:
        nmd_sensitive_seq = ''.join( obj_tt.arr_nuc_seq[i_aberrant_start : i_tp_pux + 1] )

    ##TEST:: print "NMD_SENS: seq = ", nmd_sensitive_seq

    return nmd_sensitive_seq


def get_nmd_irrelevant_region( obj_tt, pos_aberrant ):
    """
    Args:
        -obj_tt = instance of TranscribeTranslate_V4 class
        -pos_aberrant = the 3' end of the aberrant splicing event (for + strand: genomically higher position, for - strand: genomically lower position)
    Function: retireves the region of the gene that, if it contains an early stop codon, will escape NMD. However, if no stop codon, then will lead to NSD. Returns the nucleotide sequence of this region  
    """
    #get the position 55nt from the 3' end of the penultimate exon
    i_tp_pux = get_index_tp_pux( obj_tt )       #tp_pux = Truncated Protein Penultimate Exon
    if not i_tp_pux:
        return None
    pos_tp_pux = obj_tt.arr_genome_pos[ i_tp_pux ]          #tp_pux = Truncated Protein Penultimate Exon 

    #retrieve the position that should translation in the "NMD-irrelevant region", depending on which position is closer to the end of the gene - tp_pux or the end of the aberrant SJ
    if obj_tt.iso_sj.strand < 0:        #for minus strand genes
        pos_translate_start = pos_tp_pux if pos_tp_pux < pos_aberrant else pos_aberrant
    else:           #for plus strand genes
        pos_translate_start = pos_tp_pux if pos_tp_pux > pos_aberrant else pos_aberrant

    #retrieve the position of the initial nucleotide in the codon closet the boundary for the 55nt Truncated Protein Boundary
    # direction = 1 if obj_tt.iso_sj.strand < 0 else -1       #if - strand, then find 0-rf nucleotide at higher position, else if + strand, then find 0-rf nucleotide at lower position
    # pos_tp_start = obj_tt.find_codon_beginning_prev( pos_tp_pux )
    # i_tp_start = obj_tt.arr_genome_pos.index( pos_tp_start )        #this is the 0-rf closes to the truncate protein position of the penultimate exon
    ( pos_tp_start, i_tp_start ) = find_nearest_rf0( obj_tt, pos_translate_start )
    ##TEST:: get start of codon previous to SJ
    # diff = 4
    # print "GET_NMD_IRRELEVANT: pos_tp_start = ", pos_tp_start, " & i_tp_start = ", i_tp_start
    # print "GET_NMD_IRRELEVANT: before 55 bp upstream = ", obj_tt.arr_genome_pos[i_tp_start - diff : i_tp_start], " & after 55 bp upstream = ", obj_tt.arr_genome_pos[i_tp_start : i_tp_start + diff]
    # print "GET_NMD_IRRELEVANT: before 55 bp nucleotide = ", obj_tt.arr_nuc_seq[i_tp_start - diff : i_tp_start], " & after 55 bp nucleotide = ", obj_tt.arr_nuc_seq[i_tp_start : i_tp_start + diff]
    # print "GET_NMD_IRRELEVANT: before 55 bp reading frame = ", obj_tt.arr_rf[i_tp_start - diff : i_tp_start], " & after 55 bp reading frame = ", obj_tt.arr_rf[i_tp_start : i_tp_start + diff]

    ##TEST:: show the starting codon before the truncated protein boundary (55 bp upstream of 3' end of penultimate exon)

    if not i_tp_start:      #if position is None
        return None

    if obj_tt.iso_sj.strand < 0:
        nmd_irrelevant_seq = ''.join( obj_tt.arr_nuc_seq[0 : i_tp_start + 1] )
        nmd_irrelevant_seq = nmd_irrelevant_seq[::-1]

        ##TEST::
        print "NMD Irrelevant Seq MINUS STRAND: chrom = ", obj_tt.iso_sj.chrom, "& start = ", obj_tt.arr_genome_pos[0], " & end = ", obj_tt.arr_genome_pos[i_tp_start + 1]

    else:
        nmd_irrelevant_seq = ''.join( obj_tt.arr_nuc_seq[i_tp_start:] )     #retrieve from tp_start until the last nucleotide

        ##TEST::
        print "NMD Irrelevant Seq PLUS STRAND: chrom = ", obj_tt.iso_sj.chrom, "& start = ", obj_tt.arr_genome_pos[i_tp_start], " & end = ", obj_tt.arr_genome_pos[-1]

    ##TEST:: print "NMD_IRREL: strand = ", obj_tt.iso_sj.strand, " | i_tp_pux = ", i_tp_pux, " | i_tp_start = ", i_tp_start, " seq = ", nmd_irrelevant_seq

    return nmd_irrelevant_seq

print "------------ TDD: 170403_AASeq_NormalMutant.py ------------"
"""
Algorithm: this tests
    -this algorithm will determine the mutated amino acid sequence based on aberrant splicing.
"""

g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )

#show the end positions for each SJ
get_iso_id = 'NM_001008844'
get_isoform = g.refGene.filter_by( name = get_iso_id ).first()
print "for gene = ", get_isoform.name2, " & isoform = ", get_isoform.name
for i, intron in get_isoform.introns:
    print "intron ", i, " - ", intron


#Test SJ 2
# sj_id = 'JUNC00090748'
# sj_range = 'chr12:56121123-56122061'
# strand = -1
# read_count = 13
# gene_sym = 'CD63'
# isoform_id = 'NM_001257400'
# sample_prevalence = 0
# control_prevalence = 0

#Test SJ 3 - Having issues with this SJ & TranscribeTranslate_v4
# sj_id = 'JUNC00063919'
# sj_range = 'chr6:7585628-7585783'
# strand = 1
# read_count = 125
# gene_sym = 'DSP'
# isoform_id = 'NM_001008844' 
# sample_prevalence = 0
# control_prevalence = 0

##TEST SJ 4 - Check to see if I can extract amino acid sequences before aberrant SJ
sj_id = 'JUNC00124397'
sj_range = 'chr12:49666480-49666608'
strand = 1
read_count = 11
gene_sym = 'TUBA1C'
isoform_id = 'NM_001303115' 
sample_prevalence = 1
control_prevalence = 0

##TEST SJ 5 - same as TEST SJ 4 but for minus strand
# sj_id = 'JUNC00040454'
# sj_range = 'chr3:48668573-48668685'
# strand = -1
# read_count = 44
# gene_sym = 'SLC26A6'
# isoform_id = 'NM_001281733' 
# sample_prevalence = 0
# control_prevalence = 0

hash_sj_info = {'sj_id': sj_id,
'sj_range': sj_range,
'strand': strand,
'read_count': read_count,
'gene_name': gene_sym,
'isoform_id': isoform_id,
'prevalence_all': sample_prevalence,
'prevalence_control': control_prevalence
}


obj_sj = create_obj_sj( hash_sj_info )

#create transcript
list_sj = [obj_sj]
hash_pos = None
simulant_sj = False
group_sj = 0            #do not perform any grouping of SJs
iso_sj = IsoformSJ( isoform_id, list_sj, -10, hash_pos, simulant_sj, group_sj )
# print "reconstruct transcript from single SJ - version 1..."
# transcript_1 = iso_sj.build_transcript_around_sj( obj_sj )
# list_transcripts_ssj = iso_sj.reconstruct_transcript_single_sj( obj_sj )
# canon_transcript = iso_sj.create_canon_transcript()

##TEST:: see each SJ that makes up transcript - this is for function iso_sj.reconstruct_transcript_single_sj()
# print "Real transcript:"
# for i, each_transcript in enumerate( list_transcripts_ssj ):
#     print "transcript ", i, " - ", len( each_transcript )
#     for i2, each_sj in enumerate( each_transcript ):
#         print "real ", i2, " - ", each_sj

print "reconstruct transcript from single SJ - version 2..."
transcript_ssj = iso_sj.reconstruct_transcript_single_sj_v2( obj_sj )


bool_possible = TranscribeTranscript.is_obj_possible( transcript_ssj, iso_sj )
if bool_possible:
    obj_tt = TranslateTranscript( transcript_ssj, iso_sj, DIR_GENOME, {} )

    ##TEST:: show list of exons connected by SJs
    # print "MAIN: get_index_tp_pux show exons: "
    # for i, exons in enumerate( obj_tt.list_exons ):
    #     print "exon #", i, " - ", exons

    ##TEST:: show all indices (genome position, nucleotide seq)
    # obj_tt.test_display_all_indices()

    #retrieve the protein sequences to determine if NMD will occur or not
    pos_aberrant = obj_sj.end + 1 if iso_sj.strand < 0 else obj_sj.start        #get the starting point of the SJ based on the strand sign, need to add +1 for - strand because of 0-based genomic coordinates

    # #nmd_sensitive_seq = obj_tt.get_nmd_sensitive_region( pos_aberrant )
    nmd_sensitive_seq = get_nmd_sensitive_region( obj_tt, pos_aberrant )
    i_tp_pux = get_index_tp_pux( obj_tt )

    
    print "nmd_sensitive_seq = ", nmd_sensitive_seq

    # nmd_irrelevant_seq = obj_tt.get_nmd_irrelevant_region( pos_aberrant )
    nmd_irrelevant_seq = get_nmd_irrelevant_region( obj_tt, pos_aberrant )
    print "nmd_irrelevant_seq = ", nmd_irrelevant_seq

    
    #retrieve the amino acid sequence before mutation
    num_aa = 8
    hash_range_aa_before_mut = obj_tt.find_aa_pos_before( pos_aberrant, num_aa )
  
    nuc_seq_aa_before_mut = obj_tt.retrieve_nuc_seq( hash_range_aa_before_mut['i_genome_start'], hash_range_aa_before_mut['i_genome_end'], True )


    # rev_comp = True if obj_tt.iso_sj.strand < 0 else False
    # adjust_0_base = False
    # nuc_seq_aa_before_mut = ''.join( obj_tt.arr_genome_pos[ hash_range_aa_before_mut['i_genome_pos_start'] : hash_range_aa_before_mut['i_genome_end'] ] )

    #translate sequence
    protein_aa_before_mut = Seq( nuc_seq_aa_before_mut ).translate( to_stop = False ) if nuc_seq_aa_before_mut else '-'
    len_protein_aa_before_mut = len( protein_aa_before_mut )

    protein_nmd_sen = Seq( nmd_sensitive_seq ).translate( to_stop = False ) if nmd_sensitive_seq else '-'
    len_nmd_sen = '-' if not '*' in protein_nmd_sen else len( str( protein_nmd_sen ).split( '*' )[0] )
    protein_nmd_irrel = Seq( nmd_irrelevant_seq ).translate( to_stop = False ) if nmd_irrelevant_seq else '-'
    len_nmd_irrel = '-' if not '*' in protein_nmd_irrel else len( str( protein_nmd_irrel ).split( '*' )[0] )


    #compare to AA sequence of transcript
    str_genome_pos_aa_before_mut = obj_tt.iso_sj.chrom + ":" + str(hash_range_aa_before_mut['genome_start']) + "-" + str(hash_range_aa_before_mut['genome_end'])
    print "position of AA before - ", str_genome_pos_aa_before_mut
    print "nuc_seq_aa_before_mut = ", nuc_seq_aa_before_mut
    print "protein_aa_before_mut = ", protein_aa_before_mut, " & len = ", len( protein_aa_before_mut )
    print "hash for AA before SJ = ", hash_range_aa_before_mut
    print "************\n"
    print "protein_nmd_sen = ", protein_nmd_sen
    print "~~~~~~~~~~~~\n"

    ##TEST:: check the genomic range for NMD sequence
    pos_aberrant = obj_sj.end + 1 if iso_sj.strand < 0 else obj_sj.start        #get the starting point of the SJ based on the strand sign, need to add +1 for - strand because of 0-based genomic coordinates
    if obj_tt.arr_genome_pos[i_tp_pux] > pos_aberrant:
        print "FOR NMD_SENSITIVE: YES, SJ position is GREATER THAN 55 boundary"
        str_nmd_sen_pos = obj_tt.iso_sj.chrom + ":" + str( obj_tt.arr_genome_pos[i_tp_pux] ) + "-" + str( pos_aberrant )
        print "str_nmd_sen_pos = ", str_nmd_sen_pos
    else:
        print "FOR NMD_SENSITIVE: YES, SJ position is GREATER THAN 55 boundary"
        str_nmd_sen_pos = obj_tt.iso_sj.chrom + ":" + str( pos_aberrant ) + "-" + str( obj_tt.arr_genome_pos[i_tp_pux] )
        print obj_tt.iso_sj.chrom + ":" + str( obj_tt.arr_genome_pos[i_tp_pux] ) + "-" + str( pos_aberrant )



    ##TEST:: see nucleotide sequence for amino acid sequence
    print "nuc_seq_aa_before_mut = ", nuc_seq_aa_before_mut
    obj_tt.test_display_indices_range( 48668712, -30 )


print "------------ TDD: 170403_AASeq_NormalMutant.py ------------"