#/usr/bin/python
import sys

from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv7 import Isoform, IsoformSJ, GenomicVariant, TranscribeTranscript, TranslateTranscript, EnsemblVEP
from mokhaPy import mokhaPy

#Constants - directories
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv7"
DIR_DATA = DIR_CURR + "/TestData"
DIR_RESULTS = DIR_CURR + "/TestResults"

DIR_GENOME = DIR_PROJ + '/ArchiveData/hg19.fa'      #directory for samtool-indexed genome

def create_obj_tt( isoform_id, genome_pos ):
    """
    Creates an instance of TranslateTranscript
    Args:
        -isoform_id = string that is the isoform, usually in the form of an Ensembl ID (e.g. ENST000..)
        -row = from pandas Dataframe, a row from the file contain mutation position
    """
    ##TEST:: print "MAIN: start of cott: isoform_id = ", isoform_id

    db_type = 2     #this means the database is Ensembl
    hash_gp = Isoform.split_genome_pos( genome_pos )        #hash_gp = hash genome pos
    hash_pos = {'chrom': 'chr' + str(hash_gp['chrom']), 'pos_oi': hash_gp['start'] }
    iso_sj = IsoformSJ( db_type, isoform_id, [], -10, hash_pos, False, 0, True )
    canon_transcript = iso_sj.create_canon_transcript( False, False )

    obj_tt = TranslateTranscript( canon_transcript, iso_sj, DIR_GENOME, {} )

    return obj_tt


def retrieve_nmd_nsd_info( obj_tt, pos_start, pos_end, complete_rf ):
    """
    Retrieves the NMD-sensitive & NMD-irrelevant regions based on the 3' position of the genomic alteration 'pos_oi'
    Args:
        -obj_tt = instance of TranslateTranscript
        -pos_start, pos_end = integers that are genomic positions, where pos_start <= pos_end regardless of strand sign
        -complete_rf = boolean where
            -True = will retrieve the a complete reading frame, meaning the 5' position starts with rf = 0 & 3' position ends with rf = 2
            -False = will retrieve from genome_pos (or the nearest position BEFORE this position) to X codons away
    """
    #for reference, this is the column header
    # header+= "NMD_SENS_pos\tNMD_SENS_nuc_seq\tNMD_SENS_aa\t"
    # header+= "NMD_IRREL_pos\tNMD_IRREL_nuc_seq\tNMD_IRREL_aa\t"
    # header+= "bool_NMD\tbool_NSD\t"


    #retrieve the 3' of the alteration
    pos_oi = pos_start if obj_tt.iso_sj.strand < 0 else pos_end       #the reason for the -1 & +1 is so it is the position "after" (closer to the 5' position) of the genomic alteration
    
    # # [nmd_sens_genomic_range, nuc_sens, aa_sens, bool_nmd] = obj_tt.retrieve_nmd_sensitive_nuc_aa( pos_oi )
    # [nmd_irrel_genomic_range, nuc_irrel, aa_irrel, bool_nsd] = obj_tt.retrieve_nmd_irrelevant_nuc_aa( pos_oi )

    # #retrieve the 3' of the alteration
    # pos_oi = row['Start_Position'] if obj_tt.iso_sj.strand < 0 else row['End_Position']
    
    [nmd_sens_genomic_range, nuc_sens, aa_sens, bool_nmd] = obj_tt.retrieve_nmd_sensitive_nuc_aa( pos_oi, complete_rf )
    [nmd_irrel_genomic_range, nuc_irrel, aa_irrel, bool_nsd] = obj_tt.retrieve_nmd_irrelevant_nuc_aa( pos_oi, complete_rf )


    line = ""
    line+= "nmd_sens_genomic_range = " + str(nmd_sens_genomic_range) + '\n'
    # line+= "nuc_sens = " + str(nuc_sens) + '\n'
    line+= "aa_sens = " + str(aa_sens)+ '\n'

    line+= "nmd_irrel_genomic_range = " + str(nmd_irrel_genomic_range) + '\n'
    # line+= "nuc_irrel = " + str(nuc_irrel) + '\n'
    line+= "aa_irrel = " + str(aa_irrel) + '\n'

    line+= "bool_nmd = " + str(bool_nmd) + '\n'
    line+= "bool_nsd = " + str(bool_nsd) + '\n'

    print "NMD & NSD Information:\n", line
    # return line



"""
Functions: retrieve nucleotide positions around 
"""

def retrieve_seq_before_after_alt( obj_tt, row ):
    """
    retrieves the position, nucleotide, and amino acids before & after genomic alteration
    Args:
        -obj_tt = instance of TranslateTranscript
        -row = row from original file containing positions for mutations
        -pos_start & pos_end = genomic positions where start <= end regardless of strand sign 
    """
    #for reference - the columns in the file:
    # header+= "before_alt_pos\tbefore_alt_nuc\tbefore_alt_aa\tbefore_aa_len\t"
    # header+= "after_alt_pos\tafter_alt_nuc\tafter_alt_aa\tafter_aa_len\t"
    # header+= "gene_strand\t"

    [pos_codon_start, pos_codon_end] = obj_tt.find_complete_rf( pos_end, pos_start ) if obj_tt.iso_sj.strand < 0 else obj_tt.find_complete_rf( pos_start, pos_end )
    curr_codon_pos = obj_tt.iso_sj.chrom + ':' + str( pos_codon_start ) + '-' + str( pos_codon_start )

    [ hash_seq_before, hash_seq_after ] = find_surrounding_seq( obj_tt, pos_start, pos_end )

    ##TEST::
    # if 'missense' in hash_gen['most_severe_consequence']:
    #     print "at row ", i, " - there is a missense mutation --> ", hash_gen
    # if hash_seq_before['genomic_range']:
    #     print "row ", i, " - hash_seq_before = ", hash_seq_before
    # if hash_seq_after['genomic_range']:
    #     print "row ", i, " - hash_seq_after = ", hash_seq_after

    # #record position, nucleotides, & amino acids before genomic alteration
    # line = ""
    # line += curr_codon_pos + '\t'
    # line += str(hash_seq_before['genomic_range']) + '\t' if hash_seq_before['genomic_range'] else '-\t'
    # line += str(hash_seq_before['nuc_seq']) + '\t' if hash_seq_before['nuc_seq'] else '-\t'
    # line += str(hash_seq_before['aa_seq']) + '\t' if hash_seq_before['aa_seq'] else '-\t'
    # line += '-\t' if not hash_seq_before['aa_seq'] else str(len(hash_seq_before['aa_seq'])) + '\t'

    # #record position, nucleotides, & amino acids after genomic alteration
    # line += str(hash_seq_after['genomic_range']) + '\t' if hash_seq_after['genomic_range'] else '-\t'
    # line += str(hash_seq_after['nuc_seq']) + '\t' if hash_seq_after['nuc_seq'] else '-\t'
    # line += str(hash_seq_after['aa_seq']) + '\t' if hash_seq_after['aa_seq'] else '-\t'
    # line += '-\t' if not hash_seq_after['aa_seq'] else str(len(hash_seq_after['aa_seq'])) + '\t'

    # line += str( obj_tt.iso_sj.strand ) + '\t'

    # return line

    return [ hash_seq_before, hash_seq_after ]

def find_surrounding_seq( obj_tt, alt_start, alt_end, complete_rf ):
    """
    finds the position, nucleotide sequence, & amino acid sequence surrounding genomic alteration
    Args:
        -obj_tt = instance of the TranslateTranscript class
        -alt_start & alt_end = integers that are the genomic start & end positions for the alteration, note that alt_start <= alt_end in terms of genomic position, regardless of strand sign.
        -complete_rf = boolean where
            -True = will retrieve the a complete reading frame, meaning the 5' position starts with rf = 0 & 3' position ends with rf = 2
            -False = will retrieve from genome_pos (or the nearest position BEFORE this position) to X codons away
    """
    num_codons = 14
    bool_find_nearest = False

    #need to correct for position based on strand sign
    pos_five_prime = alt_end if obj_tt.iso_sj.strand < 0 else alt_start         #this is the position "start"
    pos_three_prime = alt_start if obj_tt.iso_sj.strand < 0 else alt_end        #this is the position "end"


    ##TEST::
    print "\t\t!find_surroudning_seq - pos_five_prime = ", pos_five_prime
    print "\t\t!find_surroudning_seq - pos_three_prime = ", pos_three_prime

    #determine if the flanking codons should be retrieved - do this if complete_rf = True (meaning I want the complete reading frame)
    get_codon_before = complete_rf
    get_codon_after = complete_rf


    #retrieve the position, nucleotide sequence, & aa sequence before the genomic alteration    
    hash_seq_before = obj_tt.get_seq_before_pos( pos_five_prime, num_codons, bool_find_nearest, get_codon_before, complete_rf )
    #retrieve the position, nucleotide sequence, & aa sequence after the genomic alteration
    hash_seq_after = obj_tt.get_seq_after_pos( pos_three_prime, num_codons, bool_find_nearest, get_codon_after, complete_rf )


    # #NOTE: May want to delete this section because I think I have it covered with "complete_rf" in def get_seq_before_pos() & get_seq_after_pos()
    # #NOTE: For frameshifting events (e.g. insertions, deletions), need to find the position immediately before and after
    # index_or_pos = 2
    # [curr_pos_five, pos_codon_before] = self.find_pos_x_codons_before( pos_five_prime, num_codons, index_or_pos, bool_find_nearest )
    # [curr_pos_three, pos_codon_after] = self.find_pos_x_codons_after( pos_three_prime, num_codons, index_or_pos, bool_find_nearest )

    ##I THINK I CAN DELETE ALL THIS
    # #need to correct for position based on strand sign
    # pos_five_prime = alt_end if obj_tt.iso_sj.strand < 0 else alt_start         #this is the position "start"
    # pos_three_prime = alt_start if obj_tt.iso_sj.strand < 0 else alt_end        #this is the position "end"

    # #retrieve the positions before & after the genomic alteration
    # num_codons = 14
    # bool_find_nearest = False
    # [pos_x_codon_before, pos_adj_before] = obj_tt.find_x_codons_before( pos_start, num_codons, bool_find_nearest )
    # [pos_x_codon_after, pos_adj_after] = obj_tt.find_x_codons_after( pos_end, num_codons, bool_find_nearest  )

    # if not hash_range_aa_around_sj:
    #     return {'genome_range': '-', 'nuc_seq': '-', 'prot_seq': '-', 'len': '-'}

    # #retrieve nucleotide sequence for aberrant SJ, NMD-sensitive region, and NMD-irrelevant region
    # nuc_seq_aa_around_sj = obj_tt.retrieve_nuc_seq( hash_range_aa_around_sj['i_genome_start'], hash_range_aa_around_sj['i_genome_end'], True )

    # return [ hash_seq_before, hash_seq_after ]

    print "gene_info: "
    print obj_tt.iso_sj

    print "hash_seq_before:"
    for k,v in hash_seq_before.iteritems():
        print "k = ", k, " & v = ", v

    print "hash_seq_after:"
    for k,v in hash_seq_after.iteritems():
        print "k = ", k, " & v = ", v



def get_flanking_positions( obj_tt, genomic_range, flank_type ):
    """
    Find the position flanking the genomic alteration
    -if mutation, then complete_rf = True for codons flanking codon affected
    -if insertion, then complete_rf = False, need to retrieve codons immediately before & after genomic alteration
    -if deletion, then complete_rf = False, need to retrieve codons immediately before & after genomic alteration
    Args:
        -obj_tt
        -genomic_range = string that is genomic range of genomic alteration ()
        -alt_type = integer that is the "flank_type", either to stay on the position in 'genomic_range' or to get the position before & after the genomic alteration (-1 & +1 for before & after genomic alteration, respectively)
            -0 = stay on genomic position of genomic alteration
                -should do this for SNVs and events where genomic alteration does not alter the reading frame.
            -1 = retrieve the position before (-1) & after (+1)
                -should do this for frameshifting events as the reading frame has been altered
    Returns:
        returns a string 
    """
    hash_pos = Isoform.split_genome_pos( genomic_range )
    if flank_type == 0:
        return hash_pos
    else:
        return {'chrom': hash_pos['chrom'], 'start': hash_pos['start']-1, 'end': hash_pos['end']+1 }


print "------------ TDD: 170907_Test_NMD_NSD.py ------------"

g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )

complete_rf = True          #use this for events that preserve reading frame (e.g. mutations, in-frame indels)
flank_type = 0          #flank_type: basically where does the surrounding nucleotide positions should start with respect to genomic alteration. 0 = should be for frame-preserving events, 1 = for frameshifting events

# complete_rf = False         #use this with a frameshifting event (e.g. indel) because a new reading frame is set with this frameshift event
# flank_type = 1          #flank_type: basically where does the surrounding nucleotide positions should start with respect to genomic alteration. 0 = should be for frame-preserving events, 1 = for frameshifting events

#NOTE: for complete_rf, use True = for events that preserve reading frame (e.g. mutations, in-frame indels), & 


#simulate single point mutation
isoform_id = "ENST00000502297"
genome_pos = "6:17765077-17765077"

# #simulate deletion
# #PROTOCOL: determine gene strand -> retrieve the 3' end after the genomic alteration -> construct the NMD-sensitive & NMD-irrelevant region -> retrieve the nucleotide on both sides of genomic alteration
# isoform_id = "ENST00000164024"
# genome_pos = "3:48698480-48698481"      #for deletion, use position 48698481 if + gene or 48698480 if - gene
# nuc_change = "GC/T"     #original_nucleotide/new_nucleotide -> I don't need this, this is more a reference

# #simulate insertion
# #PROTOCOL: determine gene strand -> retrieve the position after the genomic alteration, in this case ON the 3' most position (rightmost for + genes & leftmost for - genes) -> construct the NMD-sensitive & NMD-irrelevant region -> retrieve the nucleotide on both sides of genomic alteration
# isoform_id = "ENST00000417314"
# genome_pos = "3:67059913-67059914"      #for insertion, use position 67059914 if + gene or 67059913 if - gene
# nuc_change = "TA/AATTTAATGTCTTTAA"

# #simulate insertion
# #PROTOCOL: determine gene strand -> retrieve the position after the genomic alteration, in this case ON the 3' most position (rightmost for + genes & leftmost for - genes) -> construct the NMD-sensitive & NMD-irrelevant region -> retrieve the nucleotide on both sides of genomic alteration
# isoform_id = "ENST00000360128"
# genome_pos = "8:33346656-33346657"      #for insertion, use position 67059914 if + gene or 67059913 if - gene
# nuc_change = "CA/AGAGGAAACTTGTTCCTTTG"


obj_tt = create_obj_tt( isoform_id, genome_pos )
#retrieve the positions for retrieving the surrounding AA
hash_pos = get_flanking_positions( obj_tt, genome_pos, flank_type )

print "\t\t!!hash_pos = ", hash_pos, " & flank_type = ", flank_type


# [ hash_seq_before, hash_seq_after ] = find_surrounding_seq( obj_tt, hash_pos['start'], hash_pos['end'], complete_rf )
find_surrounding_seq( obj_tt, hash_pos['start'], hash_pos['end'], complete_rf )

#evaluate if susceptible to NMD or not
retrieve_nmd_nsd_info( obj_tt, hash_pos['start'], hash_pos['end'], complete_rf )


#DELETE SECTION - MAY WANT TO DO THIS BECAUSE I HAVE "find_surrounding_seq"
# retrieve_nmd_nsd_info( obj_tt, genome_pos )
# #also can I retrieve the nucleotides around this position
# [ hash_seq_before, hash_seq_after ] = retrieve_seq_before_after_alt( obj_tt, row )



#EXPERIMENT - Retrieve effect of 
#simulate single point mutation
isoform_id = "ENST00000502297"
genome_pos = "6:17765077-17765077"

print "------------ TDD Completed: 170907_Test_NMD_NSD.py ------------"