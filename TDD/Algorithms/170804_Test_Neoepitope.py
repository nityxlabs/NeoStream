#/usr/bin/python
import sys
import time

from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv6 import Isoform, IsoformSJ, GenomicVariant, TranscribeTranscript, TranslateTranscript
from mokhaPy import mokhaPy

#Constants - directories
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv6"
DIR_DATA = DIR_CURR + "/TestData"
DIR_RESULTS = DIR_CURR + "/TestResults"

DIR_GENOME = DIR_PROJ + '/ArchiveData/hg19.fa'      #directory for samtool-indexed genome

def find_containing_isoform( check_pos, genomic_range, all_isoforms ):
    """
    Finds the isoform that contains the position 'check_pos'
    """
    list_containing_indices = []        #array of the index of isoforms that contain isoforms that position 'check_pos'
    for i, each_iso in enumerate( all_isoforms ):
        isoform_id = each_iso.name
        #skip 
        if 'NR_' in isoform_id:
            continue

        obj_tt = create_tt_instance( isoform_id, genomic_range )

        if check_pos in obj_tt.arr_genome_pos:
            i_pos = obj_tt.arr_genome_pos.index( check_pos )
            # print i, ": i_pos = ", i_pos, " in isoform ", isoform_id, " & strand = ", each_iso.strand, " & obj_tt.strand = ", obj_tt.iso_sj.strand, " contains position ", check_pos
            list_containing_indices.append( i )
        else:
            print "Position ", check_pos, " not found in isoform ", isoform_id
            # if not check_pos in obj_tt.arr_genome_pos:
                # print "----- the position ", check_pos, " is not in isoform ", isoform_id, " -----"


        # try:
        #     #Testing new function
        #     i_pos = obj_tt.arr_genome_pos.index( check_pos )
        #     print i, ": i_pos = ", i_pos, " in isoform ", isoform_id
        # except ValueError:
        #     print "Position not found in isoform ", isoform_id
        #     # if not check_pos in obj_tt.arr_genome_pos:
        #     #     print "----- the position ", check_pos, " is not in isoform ", isoform_id, " -----"

    return list_containing_indices



def create_tt_instance( isoform_id, genomic_range ):
    """
    Create instance of class TranslateTranscript
    """

    #retrieve position info
    split_genome_pos = Isoform.split_genome_pos( genomic_range )
    hash_pos = { 'chrom': split_genome_pos['chrom'], 'pos_oi': split_genome_pos['start'] }      #used to find the closest position for an isoform

    #create instance of IsoformSJ
    bool_simulant_sj = False
    group_sj = 0        #this means splicing events will NOT be grouped into 5' or 3' competitive splicing
    iso_sj = IsoformSJ( isoform_id, [], -10, hash_pos, bool_simulant_sj, group_sj )

    #create instance of TranslateTranscript
    canon_transcript = iso_sj.create_canon_transcript()
    #reconstruct canonical transcript
    obj_tt = TranslateTranscript( canon_transcript, iso_sj, DIR_GENOME, {} )

    return obj_tt

"""
Functions: Algorithms for making sure my algorithm is correct
"""
def find_nearest_rf0( obj_tt, pos_oi, bool_before = True ):
    """
    Args:
        -obj_tt = instance of class TranscribeTranscript/TranslateTranscript
        -pos_oi = integer that is the position to find the nearest position that has a reading frame of 0
        -bool_before = boolean where if
            -True = retrieves nucleotides & AA before position 'given_pos' (by 'before' I mean closer to the 5' end of the gene - this takes strand sign into consideration)
            -False = retrieves nucleotides & AA after position 'given_pos' (by 'after' I mean closer to the 3' end of the gene - this takes strand sign into consideration)
    Function: finds the nearest 0 reading frame based on the position 'pos_oi' and the strand sign 'obj_tt.iso_sj.strand'
    """
    #retrieve the position of the initial nucleotide in the codon closet the boundary for the TP (truncated protein) boundary
    pos_nearest_start = obj_tt.find_codon_beginning_nearest( pos_oi, bool_before )
    if not pos_nearest_start:
        return (None, None)

    i_nearest_start = obj_tt.arr_genome_pos.index( pos_nearest_start )        #this is the 0-rf closes to the truncate protein position of the penultimate exon

    return ( pos_nearest_start, i_nearest_start )

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

def get_nmd_sensitive_region( obj_tt, pos_aberrant ):
    """
    Args:
        -obj_tt = instance of class TranscribeTranscript/TranslateTranscript
        -pos_aberrant = integer that is the position that contains the aberrant position of interest. This will usually be the start of aberrant SJ (lower position for + genes & higher position for the - genes)
    Function: retrieve the region of the gene that, if it contains an early stop codon, will lead to degradation of the transcript. Returns the nucleotide sequence of this region
    """
    #retrieve the position of the initial nucleotide in the codon
    direction = 1 if obj_tt.iso_sj.strand < 0 else -1       #if minus strand, then find 0-rf nucleotide at higher position, else if + strand, then find 0-rf nucleotide at lower position
    bool_before = True
    pos_aberrant_start = obj_tt.find_codon_beginning_nearest( pos_aberrant, bool_before )

    if not pos_aberrant_start:
        return [None, None, None]

    i_aberrant_start = obj_tt.arr_genome_pos.index( pos_aberrant_start )
    # ( pos_aberrant_start, i_aberrant_start ) = find_nearest_rf0( obj_tt, pos_aberrant, True )

    #get the position 55nt from the 3' end of the penultimate exon
    i_tp_pux = get_index_tp_pux( obj_tt )
    if not i_tp_pux:
        return [None, None, None]
    #FROM FEEDBACK: do not overlap between NMD-sensitive region & NMD-irrelevant region, therefore NMD-sensitive region should have the codon that contains the 55nt penultimate boundary
    #FROM FEEDBACK: retrieve the the codon that contains the 55nt penultimate position
    rf_tp_pux = obj_tt.arr_rf[ i_tp_pux ]
    diff_tp_pux_codon = 2 - rf_tp_pux
    i_codon_before_tp_pux = i_tp_pux - diff_tp_pux_codon if obj_tt.iso_sj.strand < 0 else i_tp_pux + diff_tp_pux_codon
    

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
        nmd_sensitive_seq = ''.join( obj_tt.arr_nuc_seq[i_codon_before_tp_pux : i_aberrant_start + 1] )      #need to add +1 to include last base
        nmd_sensitive_seq = nmd_sensitive_seq[::-1]
        nmd_sensitive_pos = obj_tt.iso_sj.chrom + ':' + str(obj_tt.arr_genome_pos[i_codon_before_tp_pux]) + '-' + str(obj_tt.arr_genome_pos[i_aberrant_start])
    else:
        nmd_sensitive_seq = ''.join( obj_tt.arr_nuc_seq[i_aberrant_start : i_codon_before_tp_pux + 1] )
        nmd_sensitive_pos = obj_tt.iso_sj.chrom + ':' + str(obj_tt.arr_genome_pos[i_aberrant_start]) + '-' + str(obj_tt.arr_genome_pos[i_codon_before_tp_pux])

    ##TEST:: print "NMD_SENS: seq = ", nmd_sensitive_seq

    str_pos_tp_pux = obj_tt.iso_sj.chrom + ':' + str( obj_tt.arr_genome_pos[i_tp_pux] )
    return [nmd_sensitive_seq, nmd_sensitive_pos, str_pos_tp_pux ]


def get_nmd_irrelevant_region( obj_tt, pos_aberrant ):
    """
    Args:
        -obj_tt = instance of TranscribeTranslate_V4 class
        -pos_aberrant = the 3' end of the aberrant splicing event (for + strand: genomically higher position, for - strand: genomically lower position)
    Function: retireves the region of the gene that, if it contains an early stop codon, will escape NMD. However, if no stop codon, then will lead to NSD. Returns the nucleotide sequence of this region  
    """
    #get the position 55nt from the 3' end of the penultimate exon
    i_tp_pux = get_index_tp_pux( obj_tt )       #tp_pux = Truncated Protein Penultimate Exon


    ##TEST::
    print "get_nmd_irrelevant_region: i_tp_pux = ", i_tp_pux

    if not i_tp_pux:
        return [None, None]
    pos_tp_pux = obj_tt.arr_genome_pos[ i_tp_pux ]          #tp_pux = Truncated Protein Penultimate Exon 

    #retrieve the position that should translation in the "NMD-irrelevant region", depending on which position is closer to the end of the gene - tp_pux or the end of the aberrant SJ
    if obj_tt.iso_sj.strand < 0:        #for minus strand genes
        pos_translate_start = pos_tp_pux if pos_tp_pux < pos_aberrant else pos_aberrant
    else:           #for plus strand genes
        pos_translate_start = pos_tp_pux if pos_tp_pux > pos_aberrant else pos_aberrant

    #retrieve the next codon
    if pos_translate_start == pos_tp_pux:
        rf_tp_pux = obj_tt.arr_rf[ i_tp_pux ]
        diff_after_tp_pux = 3 - rf_tp_pux
        i_codon_after_tp_pux = i_tp_pux - diff_after_tp_pux if obj_tt.iso_sj.strand < 0 else i_tp_pux + diff_after_tp_pux
        pos_translate_start = obj_tt.arr_genome_pos[i_codon_after_tp_pux]

    #retrieve the position of the initial nucleotide in the codon closet the boundary for the 55nt Truncated Protein Boundary
    # direction = 1 if obj_tt.iso_sj.strand < 0 else -1       #if - strand, then find 0-rf nucleotide at higher position, else if + strand, then find 0-rf nucleotide at lower position
    # pos_tp_start = obj_tt.find_codon_beginning_prev( pos_tp_pux )
    # i_tp_start = obj_tt.arr_genome_pos.index( pos_tp_start )        #this is the 0-rf closes to the truncate protein position of the penultimate exon
    ( pos_tp_start, i_tp_start ) = find_nearest_rf0( obj_tt, pos_translate_start, True )
    ##TEST:: get start of codon previous to SJ
    # diff = 4
    # print "GET_NMD_IRRELEVANT: pos_tp_start = ", pos_tp_start, " & i_tp_start = ", i_tp_start
    # print "GET_NMD_IRRELEVANT: before 55 bp upstream = ", obj_tt.arr_genome_pos[i_tp_start - diff : i_tp_start], " & after 55 bp upstream = ", obj_tt.arr_genome_pos[i_tp_start : i_tp_start + diff]
    # print "GET_NMD_IRRELEVANT: before 55 bp nucleotide = ", obj_tt.arr_nuc_seq[i_tp_start - diff : i_tp_start], " & after 55 bp nucleotide = ", obj_tt.arr_nuc_seq[i_tp_start : i_tp_start + diff]
    # print "GET_NMD_IRRELEVANT: before 55 bp reading frame = ", obj_tt.arr_rf[i_tp_start - diff : i_tp_start], " & after 55 bp reading frame = ", obj_tt.arr_rf[i_tp_start : i_tp_start + diff]

    ##TEST:: show the starting codon before the truncated protein boundary (55 bp upstream of 3' end of penultimate exon)

    ##TEST::
    print "get_nmd_irrelevant_region: i_tp_start = ", i_tp_start

    if not i_tp_start:      #if position is None
        return [None, None]

    if obj_tt.iso_sj.strand < 0:
        nmd_irrelevant_seq = ''.join( obj_tt.arr_nuc_seq[0 : i_tp_start + 1] )
        nmd_irrelevant_seq = nmd_irrelevant_seq[::-1]
        nmd_irrelevant_pos = obj_tt.iso_sj.chrom + ':' + str(obj_tt.arr_genome_pos[0]) + '-' + str(obj_tt.arr_genome_pos[i_tp_start])
    else:
        nmd_irrelevant_seq = ''.join( obj_tt.arr_nuc_seq[i_tp_start:] )     #retrieve from tp_start until the last nucleotide
        nmd_irrelevant_pos = obj_tt.iso_sj.chrom + ':' + str(obj_tt.arr_genome_pos[i_tp_start]) + '-' + str(obj_tt.arr_genome_pos[-1])

    ##TEST:: print "NMD_IRREL: strand = ", obj_tt.iso_sj.strand, " | i_tp_pux = ", i_tp_pux, " | i_tp_start = ", i_tp_start, " seq = ", nmd_irrelevant_seq

    return [nmd_irrelevant_seq, nmd_irrelevant_pos]



print "------------ TDD: 170804_Test_Neoepitope.py ------------"

start_time = time.time()

g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )



##----- EXPERIMENT - Determine if I can retrieve surrounding nucleotides (& resulting AA).
#Aberrant Splicing: Plus-strand gene
# gene_sym = "KIF20A"
# isoform_id = "NM_005733"
# sj_pos = "chr5:137522852-137522947"

#Aberrant Splicing: Minus-strand gene
# gene_sym = "ZBTB17"
# isoform_id = "NM_001324137"
# sj_pos = "chr1:16271231-16271575"

# hash_sj_pos = Isoform.split_genome_pos( sj_pos )
# obj_tt = create_tt_instance( isoform_id, sj_pos )

# pos_oi = hash_sj_pos['start']
# index_or_pos = 2
# bool_find_nearest = False
# pos_codon_prev = obj_tt.find_pos_x_codons_before( pos_oi, 1, index_or_pos, bool_find_nearest )
# pos_codon_next = obj_tt.find_pos_x_codons_after( pos_oi, 1, index_or_pos, bool_find_nearest )

# print "pos_codon_prev = ", pos_codon_prev
# print "pos_codon_next = ", pos_codon_next

# [b_genomic_range, b_nuc_seq, b_aa_seq] = obj_tt.retrieve_nuc_aa_flanking( pos_codon_prev, 4, True )
# [a_genomic_range, a_nuc_seq, a_aa_seq] = obj_tt.retrieve_nuc_aa_flanking( pos_codon_next, 5, False )

# str_pos = hash_sj_pos['chrom'] + ":" + str(pos_oi)
# print "pos_oi = ", str_pos

# print "Position BEFORE:" 
# print "genomic_range = ", b_genomic_range
# print "nuc_seq = ", b_nuc_seq
# print "aa_seq = ", b_aa_seq

# print "Position AFTER:"
# print "genomic_range = ", a_genomic_range
# print "nuc_seq = ", a_nuc_seq
# print "aa_seq = ", a_aa_seq



##------- UNCOMMENT LATER -------

##----- EXPERIMENT - Be able to retrieve neoepitopes

# #example - nonsynonymous mutation - correct answer: AA (Original/mutated) =  N/K; Codon (original/mutated) = aaT/aaG
# gene_sym = 'RIMKLB'
# genomic_range = 'chr12:8902471-8902471'
# strand = 1
# base_orig = 'T'
# base_alt = 'G'
# change_type = 0         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

# #example - nonsynonymous mutation - correct answer: AA (Original/mutated) =  E/Q; Codon (original/mutated) = Gag/Cag
# gene_sym = 'SLC6A13'
# genomic_range = 'chr12:369101-369101'
# strand = 1
# base_orig = 'C'
# base_alt = 'G'
# change_type = 0         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - DNP (2 nonsynonymous mutation) - correct answer: AA (Original/mutated) =  S/F; Codon (original/mutated) = tCC/tTT
gene_sym = 'CHRNA4'
genomic_range = 'chr20:61981956-61981957'
strand = 1
base_orig = 'GG'
base_alt = 'AA'
change_type = 0         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#I THINK I CAN DELETE THIS
# hash_sj_pos = Isoform.split_genome_pos( sj_pos )
# obj_tt = create_tt_instance( isoform_id, sj_pos )


#retrieve the isoforms that contain position
all_isoforms = Isoform.obj_cruzdb.refGene.filter_by( name2 = gene_sym ).all()
split_genome_pos = Isoform.split_genome_pos( genomic_range )
isoform_indices = find_containing_isoform( split_genome_pos['start'], genomic_range, all_isoforms )

#go through each isoform
for each_iso_i in isoform_indices:
    #retrieve isoform
    isoform_id = all_isoforms[each_iso_i].name
    #create instance TranslateTranscript
    obj_tt = create_tt_instance( isoform_id, genomic_range )

    print each_iso_i, " - ", all_isoforms[each_iso_i].name, ": obj_tt.iso_sj.strand = ", obj_tt.iso_sj.strand, " & strand = ", strand


    #add mutation to obj_tt.list_alterations
    # obj_tt.add_alteration( self, nuc_change, genomic_range, strand, change_type )
    obj_tt.add_alteration( base_alt, genomic_range, strand, change_type )

    # print "obj_tt.list_alterations = ", obj_tt.list_alterations

    [hash_gen_orig, hash_gen_alt] = obj_tt.alterations_to_neoepitope( 9 )

    ##TEST::
    # list_pos = sorted( hash_gen_orig.keys(), reverse = True ) if obj_tt.iso_sj.strand < 0 else sorted( hash_pos_nuc.keys() )
    list_pos = sorted( hash_gen_orig.keys() )
    for k in list_pos:
        print "isoform = ", obj_tt.iso_sj.isoform_id, " | k = ", k, " -> original = ", hash_gen_orig[k]
        print "isoform = ", obj_tt.iso_sj.isoform_id, " | k = ", k, " -> mutation = ", hash_gen_alt[k]
        print ">>>>>>>>>>>\n"

    [hash_window_orig, hash_window_mut] = obj_tt.create_neoepitope_window( hash_gen_orig, hash_gen_alt, 9 )

    for k2,v2 in hash_window_orig.iteritems():
        print "isoform = ", obj_tt.iso_sj.isoform_id, " | k = ", k2, " -> original = ", v2
        print "isoform = ", obj_tt.iso_sj.isoform_id, " | k = ", k2, " -> mutation = ", hash_window_mut[k2]
        print "~~~~~~~~~~~~~~\n"

    print "############################\n\n"


# print "Basic gene information: "
# print "Info about gene = ", obj_tt
# pu_exon = obj_tt.list_exons[1] if obj_tt.iso_sj.strand < 0 else obj_tt.list_exons[::-1][1]
# last_exon = obj_tt.list_exons[0] if obj_tt.iso_sj.strand < 0 else obj_tt.list_exons[::-1][0]
# print "Penultimate exon = ", pu_exon
# print "Last exon = ", last_exon 

# #retrieve the protein sequences to determine if NMD will occur or not
# pos_aberrant_start = hash_sj_pos['end'] + 1 if obj_tt.iso_sj.strand < 0 else hash_sj_pos['start']        #get the starting point (5' end of SJ) of the SJ based on the strand sign, need to add +1 for - strand because of 0-based genomic coordinates\
# pos_aberrant_end = hash_sj_pos['end'] if obj_tt.iso_sj.strand < 0 else hash_sj_pos['start'] + 1        #get the end point (3' end of SJ) of the SJ based on the strand sign, need to add +1 for - strand because of 0-based genomic coordinates


# [sens_genomic_range, sens_nuc_seq, sens_aa_seq] = obj_tt.retrieve_nmd_sensitive_nuc_aa( pos_aberrant_start )
# [irrel_genomic_range, irrel_nuc_seq, irrel_aa_seq] = obj_tt.retrieve_nmd_irrelevant_nuc_aa( pos_aberrant_end )

# print "sens_genomic_range = ", sens_genomic_range
# print "len nuc = ", len(sens_nuc_seq) if sens_nuc_seq else "None"
# print "sens_aa_seq = ", sens_aa_seq, " & len AA = ", len(sens_aa_seq) if sens_aa_seq else "None"

# # print "irrel_genomic_range = ", irrel_genomic_range
# # print "irrel_aa_seq = ", irrel_aa_seq


# #retrieve nucleotide sequence for aberrant SJ, NMD-sensitive region, and NMD-irrelevant region
# [nmd_sensitive_seq, nmd_sensitive_pos, str_pos_tp_pux] = get_nmd_sensitive_region( obj_tt, pos_aberrant_start )
# [nmd_irrelevant_seq, nmd_irrelevant_pos] = get_nmd_irrelevant_region( obj_tt, pos_aberrant_end )

# print "THESE SHOULD BE THE RIGHT ANSWER:"
# print "nmd_sensitive_pos = ", nmd_sensitive_pos
# print "str_pos_tp_pux = ", str_pos_tp_pux
# print "len( NMD sensitive nuc ) = ", len( nmd_sensitive_seq )

elapse_time = time.time() - start_time
mokhaPy.timeElapse_convertToHMS( elapse_time )
print "------------ TDD: 170804_Test_Neoepitope.py ------------"