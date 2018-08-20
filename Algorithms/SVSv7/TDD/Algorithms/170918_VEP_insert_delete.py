#/usr/bin/python
import sys
import re
import copy

from Bio.Seq import Seq

from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv7 import Isoform, IsoformSJ, TranscribeTranscript, TranslateTranscript, EnsemblVEP
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

def determine_alt_rf( nuc_change ):
    """
    Determines if a genomic alteration causes a frameshift or not
    Args:
        nuc_change = string, 'allele_string' column evaluated from Ensembl VEP, format: nucleotide_original/nucleotide_mutated (e.g. A/G, T/-, GT/C, GGTGC/-, C/TTAT)
    Returns: an array where
        -[0] = count_change, the nucleotide difference before & after alteration (considered to determine if mutation, insertion, or deletion)
        -[1] = in_frame = boolean where True = alteration DOES NOT perturb gene's reading frame, False = alteartion does perturb gene's reading frame
        -[2] = change_type, the type of alteration it is (mutation, deletion, insertion)
    """
    if not '/' in nuc_change:
        return None
    list_nuc = [str( re.sub('[^ACTG]', '', x) ) for x in nuc_change.split('/')]      #make sure to remove any character that are not nucleotides

    #calculate the change, where count_change = 0 means no change in the number of nucleotides, count_change > 0 means deletion (less nucleotides after alteration), & count_change < 0 means insertion (more nucleotides after alteration)
    count_change = len(list_nuc[0]) - len(list_nuc[1])
    in_frame = False if count_change % 3 else True

    #retrieve the alteration type
    if count_change == 0:        #mutation (same # of nucleotides altered as original)
        change_type = 0
    elif count_change < 0:      #insertion (more nucleotides after alteration)
        change_type = 1
    elif count_change > 0:      #deletion (less nucleotides after alteration)
        change_type = 2

    ##TEST::
    print "D_ALT 1: nuc_change = ", nuc_change, " & list_nuc = ", list_nuc
    print "D_ALT 2: count_change = ", count_change, " & in_frame? - ", count_change % 3, " & ", in_frame

    return [count_change, in_frame, change_type]

def retrieve_affected_aa( obj_tt, obj_vep, genomic_range, nuc_change, aa_change ):
    """
    retrieves the nucleotide sequence 
    Args:
        -obj_tt = instance of TranslateTranscript class, USE: retrieve nucleotide sequence
        -obj_vep = instance of VEP class, USE: for now, only need it for the strand sign where the genomic alteration takes place, which is usually "+" strand anyways
        -genomic_range = string that is the genomic range (not "chr" needed) (format: #:start-end, e.g. 3:)
        -nuc_change = string, 'allele_string' column evaluated from Ensembl VEP, format: nucleotide_original/nucleotide_mutated (e.g. A/G, T/-, GT/C, GGTGC/-, C/TTAT)
        -aa_change = string, 'amino_acids' column evaluated from Ensembl VEP, format: aa_original/aa_mutated
    """
    [count_change, in_frame, change_type] = determine_alt_rf( nuc_change )

    ##TEST::
    print "R_A_AA 0: nuc_change = ", nuc_change, " & aa_change = ", aa_change
    print "R_A_AA 1: count_change = ", count_change, " & in_frame = ", in_frame

    if in_frame:
        list_aa = [ str(re.sub('[JOUX]', '', x)) for x in aa_change.split('/')]     #JOUX are not letters designated to 

        print "R_A_AA InFrame: list_aa = ", list_aa

        # return list_aa[1]
    else:
        #find the rf = 0 previous of the genomic alteration
        #
        if count_change == 0:        #mutation (same # of nucleotides altered as original)
            change_type = 0
        elif count_change < 0:      #insertion (more nucleotides after alteration)
            change_type = 1
        elif count_change > 0:      #deletion (less nucleotides after alteration)
            change_type = 2

        complete_rf = True
        find_closest = False
        hash_pos_nuc = obj_tt.retrieve_genomic_subseq( genomic_range, complete_rf, find_closest )
        # copy_hash_1 = hash_pos_nuc.copy()
        # copy_hash_2 = dict( hash_pos_nuc )
        copy_hash_3 = copy.deepcopy( hash_pos_nuc )     ##This is the only one that does not reference the original "hash_pos_nuc", therefore does not change the original hash

        #need to extract the mutated AA and evaluate the change
        list_nuc = [str( re.sub('[^ACTG]', '', x) ) for x in nuc_change.split('/')]      #make sure to remove any character that are not nucleotides, array where 0 = original nucleotide, 1 = mutated nucleotide
        hash_pos_mut = obj_tt.create_alteration( copy_hash_3, list_nuc[1], genomic_range, obj_vep.strand, change_type )

        ##TEST::
        hash_changetype = {0: "mutation", 1: "insertion", 2: "deletion", 3: "aberrant splicing"}
        print "R_A_AA Frameshift: change_type = ", hash_changetype[change_type]
        print "R_A_AA Frameshift: list_nuc = ", list_nuc
        print "R_A_AA Frameshift: hash_pos_nuc = ", hash_pos_nuc
        # print "R_A_AA COPY: copy_hash_1 = ", copy_hash_1
        # print "R_A_AA COPY: copy_hash_2 = ", copy_hash_2
        print "R_A_AA COPY: copy_hash_3 = ", copy_hash_3
        print "R_A_AA Frameshift: hash_pos_mut = ", hash_pos_mut

        return hash_pos_mut

def retrieve_aa_after_frameshift( obj_tt, hash_pos_nuc ):
    """
    Retrieve the nucleotides & AA affected due to a frameshift
    Args:
        -obj_tt = instance of TranslateTranscript class, USE: retrieve nucleotide sequence
        -hash_pos_nuc = retrieved from def retrieve_genomic_subseq(), where k = genomic position, & v = hash where k2 = 'nuc' & 'rf' & v2 = the nucleotide & the reading frame at that position, respectively
    PROTOCOL:
        -retrieve the codon affected by the frameshift event (insertion, deletion) --> should also consider aberrant splicing
            -use retrieve_affected_aa(), as this uses obj_tt.create_alteration()
            -this should have a data structure along the lines of hash_pos_nuc (key = genomic position, value = hash where k2 = 'nuc' & 'rf' & v2 = nucleotide & reading frame, respectively)
        -make the changes to the position & combine all the nucletide sequences
        -Should I calculate the # of nucleotides missing to make this a complete reading frame (e.g. if delete 1 nucleotide, then need to replace that one. If add X nucleotides so it needs 2 more to make a complete reading frame, then retrieve 2 more nucleotides)
    """
    nuc_seq_mod_codon = ''.join( [ v['nuc'] for k,v in hash_pos_nuc.iteritems() if v['nuc'] != '-' ] )
    #retrieve the nucleotides needed to make 'nuc_seq_mod_codon' have a complete reading frame
    remainder_codon_num = len( nuc_seq_mod_codon ) % 3
    pos_alt_end = min( hash_pos_nuc.keys() ) if obj_tt.iso_sj.strand < 0 else max( hash_pos_nuc.keys() )        #retrieve the 3' end
    

    if not pos_alt_end in obj_tt.arr_genome_pos:
        direction = -1 if obj_tt.iso_sj.strand < 0 else 1     #find a position closer to the 3' end of the gene
        pos_alt_end = obj_tt.find_nearest_pos( pos_alt_end, direction )
        if not pos_alt_end:
            return None
    i_remainder_start = obj_tt.arr_genome_pos.index( pos_alt_end )

    #retrieve the nucleotides to complete 'mod_codon'
    add_nuc_seq = ""
    i_end_nuc_pos = i_remainder_start
    for i in range(0, remainder_codon_num):
        curr_i = i_remainder_start + (i * obj_tt.iso_sj.strand)
        # if curr_i < 0 or curr_i >= len( obj_tt.arr_nuc_seq ):
        if curr_i >= len( obj_tt.arr_nuc_seq ):
            break
        #retrieve the end position & the nucleotide sequence
        i_end_nuc_pos = curr_i
        add_nuc_seq += obj_tt.arr_nuc_seq[curr_i]

    #retrieve the range of the modified codon
    pos_alt_start = max( hash_pos_nuc.keys() ) if obj_tt.iso_sj.strand < 0 else min( hash_pos_nuc.keys() )        #retrieve the 3' end
    pos_mod_codon_end = obj_tt.arr_genome_pos[i_end_nuc_pos]
    if pos_alt_start < pos_mod_codon_end:
        range_mod_codon = obj_tt.iso_sj.chrom + ":" + str( pos_alt_start ) + "-" + str( pos_mod_codon_end )
    else:
        range_mod_codon = obj_tt.iso_sj.chrom + ":" + str( pos_mod_codon_end ) + "-" + str( pos_alt_start )

    ##TEST:: see the modified sequence
    print "add_nuc_seq = ", add_nuc_seq
    print "nuc_seq_mod_codon = ", nuc_seq_mod_codon

    #retrieve sequences for modified codon
    nuc_seq_mod_codon += add_nuc_seq
    aa_seq_mod_codon = str( Seq(nuc_seq_mod_codon).translate( to_stop = False ) )

    #retrieve the sequence after the modified codon
    i_after = i_remainder_start + ((remainder_codon_num + 1) * obj_tt.iso_sj.strand)
    # if i_after < 0 or i_after >= len( obj_tt.arr_nuc_seq ):
    if i_after >= len( obj_tt.arr_nuc_seq ):
        range_after = None
        nuc_seq_after = None
        aa_seq_after = None
    else:
        if obj_tt.iso_sj.strand < 0:
            range_after = obj_tt.iso_sj.chrom + ":" + str(obj_tt.arr_genome_pos[0]) + "-" + str(obj_tt.arr_genome_pos[i_after])
        else:
            range_after = obj_tt.iso_sj.chrom + ":" + str(obj_tt.arr_genome_pos[i_after]) + "-" + str(obj_tt.arr_genome_pos[-1])
        nuc_seq_after = ''.join( obj_tt.arr_nuc_seq[i_after:] )
        aa_seq_after = str( Seq(nuc_seq_after).translate( to_stop = False ) )
    
    return {
    'range_mod_codon': range_mod_codon,
    'nuc_seq_mod_codon': nuc_seq_mod_codon,
    'aa_seq_mod_codon': aa_seq_mod_codon,
    'range_after': range_after,
    'nuc_seq_after': nuc_seq_after,
    'aa_seq_after': aa_seq_after,
    }


def retrieve_aa_after_frameshift_v2( obj_tt, obj_vep, genomic_range, nuc_change ):
    """
    Similar to def retrieve_aa_after_frameshift(), but retrieves the AA sequence for both the original & mutated sequence
    Args:
        -obj_tt = instance of TranslateTranscript class, USE: retrieve nucleotide sequence
        -obj_vep = instance of VEP class, USE: for now, only need it for the strand sign where the genomic alteration takes place, which is usually "+" strand anyways
        -genomic_range = string that is the genomic range (not "chr" needed) (format: #:start-end, e.g. 3:)
        -nuc_change = string, 'allele_string' column evaluated from Ensembl VEP, format: nucleotide_original/nucleotide_mutated (e.g. A/G, T/-, GT/C, GGTGC/-, C/TTAT)
    PROTOCOL:
        -find the rf = 0 before the genomic alteration
        -add the nucleotide into the position of interest using obj_tt.create_alteration
        -then retrieve from that position until the end of the gene
    """
    [count_change, in_frame, change_type] = determine_alt_rf( nuc_change )
    
    #find the rf = 0 previous of the genomic alteration
    complete_rf = True
    find_closest = False
    hash_pos_nuc = obj_tt.retrieve_genomic_subseq( genomic_range, complete_rf, find_closest )
    hash_pos_mut = copy.deepcopy( hash_pos_nuc )     ##This is the only one that does not reference the original "hash_pos_nuc", therefore does not change the original hash

    #need to extract the mutated AA and evaluate the change
    list_nuc = [str( re.sub('[^ACTG]', '', x) ) for x in nuc_change.split('/')]      #make sure to remove any character that are not nucleotides, array where 0 = original nucleotide, 1 = mutated nucleotide
    # hash_pos_mut = obj_tt.create_alteration( hash_pos_mut, list_nuc[1], genomic_range, obj_vep.strand, change_type )
    # nuc_seq_mod_codon = ''.join( [ v['nuc'] for k,v in hash_pos_mut.iteritems() if v['nuc'] != '-' ] )


    [range_seq, alt_nuc_seq, alt_aa_seq] = obj_tt.create_alteration_nuc_aa( hash_pos_nuc, list_nuc[1], genomic_range, obj_vep.strand, change_type )
    
    #retrieve the range of the modified codon
    # pos_codon_end = min( hash_pos_nuc.keys() ) if obj_tt.iso_sj.strand < 0 else max( hash_pos_nuc.keys() )  
    # if obj_tt.iso_sj.strand < 0 :
    #     i_pos_after = obj_tt.arr_genome_pos.index( pos_codon_end ) - 1
    #     mut_nuc_seq = obj_tt.arr_nuc_seq[0:i_pos_after + 1][::-1] 
    # else:
    #     i_pos_after = obj_tt.arr_genome_pos.index( pos_codon_end ) + 1
    #     mut_nuc_seq = obj_tt.arr_nuc_seq[i_pos_orig:]
    # mut_nuc_seq = ''.join( mut_nuc_seq )

    ##TEST:: see the modified sequence
    print "F_V2: nuc_seq_mod_codon = ", nuc_seq_mod_codon, " & mut_nuc_seq = ", mut_nuc_seq

    #combine the modified codon & the amino acids after
    # mut_nuc_seq = nuc_seq_mod_codon + mut_nuc_seq
    # mut_aa_seq = str( Seq(mut_nuc_seq).translate( to_stop = False ) )
    
    #retrieve original sequence
    # pos_codon_start = max( hash_pos_nuc.keys() ) if obj_tt.iso_sj.strand < 0 else min( hash_pos_nuc.keys() )   
    # i_pos_orig = obj_tt.arr_genome_pos.index( pos_codon_start )
    # orig_nuc_seq = obj_tt.arr_nuc_seq[0:i_pos_orig + 1][::-1] if obj_tt.iso_sj.strand < 0 else obj_tt.arr_nuc_seq[i_pos_orig:]
    # orig_nuc_seq = ''.join( orig_nuc_seq )
    # orig_aa_seq = str( Seq(orig_nuc_seq).translate( to_stop = False ) )

    ##TEST::
    # mut_nuc_seq = orig_nuc_seq + mut_nuc_seq
    # mut_aa_seq = str( Seq(mut_nuc_seq).translate( to_stop = False ) )

    ##TEST:: see the original sequence
    print "orig_nuc_seq = ", orig_nuc_seq
    print "pos_codon_start = ", pos_codon_start, " & i_pos_orig = ", i_pos_orig, " & hash_pos_mut = ", hash_pos_mut
    minus_val = 10
    print "SAMPLE NUCS = ", obj_tt.arr_genome_pos[i_pos_orig - minus_val:i_pos_orig][::-1], " || ", obj_tt.arr_nuc_seq[i_pos_orig - minus_val:i_pos_orig], " || ", obj_tt.arr_rf[i_pos_orig - minus_val:i_pos_orig]

    # if obj_tt.iso_sj.strand < 0:
    #     last_pos = obj_tt.arr_genome_pos[0]
    #     mutated_genomic_range = obj_tt.iso_sj.chrom + ":" + str(last_pos) + "-" + str(obj_tt.arr_genome_pos[i_pos_orig])
    # else:
    #     last_pos = obj_tt.arr_genome_pos[-1]
    #     mutated_genomic_range = obj_tt.iso_sj.chrom + ":" + str(obj_tt.arr_genome_pos[i_pos_orig]) + "-" + str(last_pos)

    return {
    "mutated_genomic_range": range_seq,
    "mut_nuc_seq": mut_nuc_seq,
    "mut_aa_seq": mut_aa_seq,
    "orig_nuc_seq": orig_nuc_seq,
    "orig_aa_seq": orig_aa_seq
    }

def retrieve_aa_after_frameshift_v2_B( obj_tt, obj_vep, genomic_range, nuc_change ):
    """
    Similar to def retrieve_aa_after_frameshift(), but retrieves the AA sequence for both the original & mutated sequence
    Args:
        -obj_tt = instance of TranslateTranscript class, USE: retrieve nucleotide sequence
        -obj_vep = instance of VEP class, USE: for now, only need it for the strand sign where the genomic alteration takes place, which is usually "+" strand anyways
        -genomic_range = string that is the genomic range (not "chr" needed) (format: #:start-end, e.g. 3:)
        -nuc_change = string, 'allele_string' column evaluated from Ensembl VEP, format: nucleotide_original/nucleotide_mutated (e.g. A/G, T/-, GT/C, GGTGC/-, C/TTAT)
    PROTOCOL:
        -find the rf = 0 before the genomic alteration
        -add the nucleotide into the position of interest using obj_tt.create_alteration
        -then retrieve from that position until the end of the gene
    """
    # [count_change, in_frame, change_type] = determine_alt_rf( nuc_change )
    
    #find the rf = 0 previous of the genomic alteration
    complete_rf = True
    find_closest = False
    hash_pos_nuc = obj_tt.retrieve_genomic_subseq( genomic_range, complete_rf, find_closest )
    hash_pos_mut = copy.deepcopy( hash_pos_nuc )     ##This is the only one that does not reference the original "hash_pos_nuc", therefore does not change the original hash
    [range_alt, alt_nuc_seq, alt_aa_seq] = obj_tt.retrieve_alteration_nuc_aa( genomic_range, nuc_change, obj_vep.strand, complete_rf, find_closest )


    #need to extract the mutated AA and evaluate the change
    list_nuc = [str( re.sub('[^ACTG]', '', x) ) for x in nuc_change.split('/')]      #make sure to remove any character that are not nucleotides, array where 0 = original nucleotide, 1 = mutated nucleotide
    [range_alt, alt_nuc_seq, alt_aa_seq] = obj_tt.create_alteration_nuc_aa( hash_pos_nuc, list_nuc[1], genomic_range, obj_vep.strand, change_type )
    
    #retrieve the range of the modified codon
    # pos_codon_end = min( hash_pos_nuc.keys() ) if obj_tt.iso_sj.strand < 0 else max( hash_pos_nuc.keys() )  
    # if obj_tt.iso_sj.strand < 0 :
    #     i_pos_after = obj_tt.arr_genome_pos.index( pos_codon_end ) - 1
    #     mut_nuc_seq = obj_tt.arr_nuc_seq[0:i_pos_after + 1][::-1] 
    # else:
    #     i_pos_after = obj_tt.arr_genome_pos.index( pos_codon_end ) + 1
    #     mut_nuc_seq = obj_tt.arr_nuc_seq[i_pos_orig:]
    # mut_nuc_seq = ''.join( mut_nuc_seq )

    [orig_range_seq ,orig_nuc_seq, orig_aa_seq] = obj_tt.retrieve_subseq_nuc_aa( genomic_range, complete_rf, find_closest )

    ##TEST:: see the modified sequence
    print "F_V2: nuc_seq_mod_codon = ", nuc_seq_mod_codon, " & mut_nuc_seq = ", mut_nuc_seq

    #combine the modified codon & the amino acids after
    # mut_nuc_seq = nuc_seq_mod_codon + mut_nuc_seq
    # mut_aa_seq = str( Seq(mut_nuc_seq).translate( to_stop = False ) )
    
    #retrieve original sequence
    # pos_codon_start = max( hash_pos_nuc.keys() ) if obj_tt.iso_sj.strand < 0 else min( hash_pos_nuc.keys() )   
    # i_pos_orig = obj_tt.arr_genome_pos.index( pos_codon_start )
    # orig_nuc_seq = obj_tt.arr_nuc_seq[0:i_pos_orig + 1][::-1] if obj_tt.iso_sj.strand < 0 else obj_tt.arr_nuc_seq[i_pos_orig:]
    # orig_nuc_seq = ''.join( orig_nuc_seq )
    # orig_aa_seq = str( Seq(orig_nuc_seq).translate( to_stop = False ) )

    ##TEST::
    mut_nuc_seq = orig_nuc_seq + mut_nuc_seq
    mut_aa_seq = str( Seq(mut_nuc_seq).translate( to_stop = False ) )

    ##TEST:: see the original sequence
    print "orig_nuc_seq = ", orig_nuc_seq
    print "pos_codon_start = ", pos_codon_start, " & i_pos_orig = ", i_pos_orig, " & hash_pos_mut = ", hash_pos_mut
    minus_val = 10
    print "SAMPLE NUCS = ", obj_tt.arr_genome_pos[i_pos_orig - minus_val:i_pos_orig][::-1], " || ", obj_tt.arr_nuc_seq[i_pos_orig - minus_val:i_pos_orig], " || ", obj_tt.arr_rf[i_pos_orig - minus_val:i_pos_orig]

    # if obj_tt.iso_sj.strand < 0:
    #     last_pos = obj_tt.arr_genome_pos[0]
    #     mutated_genomic_range = obj_tt.iso_sj.chrom + ":" + str(last_pos) + "-" + str(obj_tt.arr_genome_pos[i_pos_orig])
    # else:
    #     last_pos = obj_tt.arr_genome_pos[-1]
    #     mutated_genomic_range = obj_tt.iso_sj.chrom + ":" + str(obj_tt.arr_genome_pos[i_pos_orig]) + "-" + str(last_pos)

    return {
    "mutated_genomic_range": orig_range_seq,
    "alt_nuc_seq": alt_nuc_seq,
    "alt_aa_seq": alt_aa_seq,
    "orig_nuc_seq": orig_nuc_seq,
    "orig_aa_seq": orig_aa_seq
    }

def retrieve_aa_inframe( obj_tt, obj_vep, genomic_range, nuc_change, aa_len ):
    """
    Retrieve the nucleotides & AA affected due to a frameshift
    Args:
        -obj_tt = instance of TranslateTranscript class, USE: retrieve nucleotide sequence
        -obj_vep = instance of VEP class, USE: for now, only need it for the strand sign where the genomic alteration takes place, which is usually "+" strand anyways
        -genomic_range = string that is the genomic range (not "chr" needed) (format: #:start-end, e.g. 3:)
        -nuc_change = string, 'allele_string' column evaluated from Ensembl VEP, format: nucleotide_original/nucleotide_mutated (e.g. A/G, T/-, GT/C, GGTGC/-, C/TTAT)
        -aa_len = the length of the flanking amino acids to retrieve before & after the genomic alteration
    PROTOCOL:
        -retrieve the codon affected by the frameshift event (insertion, deletion) --> should also consider aberrant splicing
            -use retrieve_affected_aa(), as this uses obj_tt.create_alteration()
            -this should have a data structure along the lines of hash_pos_nuc (key = genomic position, value = hash where k2 = 'nuc' & 'rf' & v2 = nucleotide & reading frame, respectively)
        -make the changes to the position & combine all the nucletide sequences
        -Should I calculate the # of nucleotides missing to make this a complete reading frame (e.g. if delete 1 nucleotide, then need to replace that one. If add X nucleotides so it needs 2 more to make a complete reading frame, then retrieve 2 more nucleotides)
    """
    [count_change, in_frame, change_type] = determine_alt_rf( nuc_change )
    
    #find the rf = 0 previous of the genomic alteration
    complete_rf = True
    find_closest = False
    hash_pos_nuc = obj_tt.retrieve_genomic_subseq( genomic_range, complete_rf, find_closest )
    key_pos = sorted( hash_pos_nuc.keys(), reverse = True ) if obj_tt.iso_sj.strand < 0 else sorted( hash_pos_nuc.keys(), reverse = False )

    hash_pos_mut = copy.deepcopy( hash_pos_nuc )     ##This is the only one that does not reference the original "hash_pos_nuc", therefore does not change the original hash

    #need to extract the mutated AA and evaluate the change
    list_nuc = [str( re.sub('[^ACTG]', '', x) ) for x in nuc_change.split('/')]      #make sure to remove any character that are not nucleotides, array where 0 = original nucleotide, 1 = mutated nucleotide
    hash_pos_mut = obj_tt.create_alteration( hash_pos_mut, list_nuc[1], genomic_range, obj_vep.strand, change_type )

    #retrieve the modified nucleotide sequence
    nuc_seq_mod_codon = ''.join( [ hash_pos_mut[k]['nuc'] for k in key_pos if hash_pos_mut[k]['nuc'] != '-' ] )
    #range for each genomic position
    range_mod_codon = obj_tt.iso_sj.chrom + ":" + str( min(hash_pos_nuc.keys()) ) + "-" + str( max(hash_pos_nuc.keys()) )

    #retrieve the range of the modified codon
    pos_codon_end = min( hash_pos_nuc.keys() ) if obj_tt.iso_sj.strand < 0 else max( hash_pos_nuc.keys() )   
    i_pos_after = obj_tt.arr_genome_pos.index( pos_codon_end ) - 1 if obj_tt.iso_sj.strand < 0 else obj_tt.arr_genome_pos.index( pos_codon_end ) + 1
    if obj_tt.iso_sj.strand < 0:
        i_pos_after_end = i_pos_after - (3 * aa_len)
        nuc_seq_after = ''.join( obj_tt.arr_nuc_seq[i_pos_after_end:i_pos_after + 1][::-1] )
        range_after = obj_tt.iso_sj.chrom + ":" + str( obj_tt.arr_genome_pos[i_pos_after_end] ) + "-" + str( obj_tt.arr_genome_pos[i_pos_after] )
    else:
        i_pos_after_end = i_pos_after + (3 * aa_len)
        nuc_seq_after = ''.join( obj_tt.arr_nuc_seq[i_pos_after:i_pos_after_end + 1] )
        range_after = obj_tt.iso_sj.chrom + ":" + str( obj_tt.arr_genome_pos[i_pos_after] ) + "-" + str( obj_tt.arr_genome_pos[i_pos_after_end] )

    #translate sequence
    aa_seq_mod_codon = str( Seq(nuc_seq_mod_codon).translate( to_stop = False ) )
    aa_seq_after = str( Seq(nuc_seq_after).translate( to_stop = False ) )

    nuc_seq_orig = ''.join( [ hash_pos_nuc[k]['nuc'] for k in key_pos if hash_pos_nuc[k]['nuc'] != '-' ] )
    aa_seq_orig = str( Seq(nuc_seq_orig).translate( to_stop = False ) )

    
    return {
    'range_mod_codon': range_mod_codon,
    'nuc_seq_orig': nuc_seq_orig,
    'aa_seq_orig': aa_seq_orig,
    'nuc_seq_mod_codon': nuc_seq_mod_codon,
    'aa_seq_mod_codon': aa_seq_mod_codon,
    'range_after': range_after,
    'nuc_seq_after': nuc_seq_after,
    'aa_seq_after': aa_seq_after,
    }


print "------------ TDD: 170918_VEP_insert_delete.py ------------"

g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )


##DELETIONS
# #simulate single deletion - RESULT: codon =  aAt/at & amino acids =  N/X
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765077-17765077"
# alt = "-"

# #simulate single deletion - RESULT: codon =  Aat/at & amino acids =  N/X
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765076-17765076"
# alt = "-"

# #simulate 2-base deletions - RESULT: codon =  AAt/t  & amino acids =  N/X
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765077-17765078"
# alt = "-"

# #simulate 3-base deletions - RESULT: codon =  taCAAt/tat  & amino acids =  YN/Y
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765077-17765079"
# alt = "-"

##MUTATION
# # simulate single point mutation - RESULT: codon =  aAt/aCt  & amino acids =  N/T
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765077-17765077"
# alt = "G"     #this is just a SNV, no insertion

##INSERTIONS
# #simulate single insertion - RESULT: codon =  aAt/aCCt  & amino acids =  N/TX
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765077-17765077"
# alt = "GG"      #one base for replace nucleotide at position, and +1 for single insertion

# #simulate double insertion - RESULT: codon =  aAt/aCCCt  & amino acids =  N/TX
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765077-17765077"
# alt = "GGG"     #one base for replace nucleotide at position, and +2 for double insertion

# #simulate triple insertion - RESULT: codon =  aAt/aCCCCt  & amino acids =  N/TP 
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765077-17765077"
# alt = "GGGG"     #one base for replace nucleotide at position, and +3 for double insertion

obj_vep = EnsemblVEP( genomic_range, alt )
obj_vep.test_display_vep()


##EXPERIMENT - Test retrieval of amino acids of affected regions


obj_tt = create_obj_tt( isoform_id, genomic_range )

list_info = obj_vep.get_all_alt_info()
for i, (hash_gen, hash_transcript) in enumerate( list_info ):

    ##TEST::
    print i, ": Codon change = ", hash_transcript.values()[i]['codons']


    nuc_change = hash_gen['allele_string']
    aa_change = hash_transcript.values()[i]['amino_acids']
    codon_change = hash_transcript.values()[i]['codons']
    [count_change, in_frame, change_type] = TranscribeTranscript.determine_alt_rf( nuc_change )
    

    ##version 2:
    if not in_frame:
        hash_new_seq_v2 = obj_tt.retrieve_alt_frameshift( genomic_range, nuc_change , obj_vep.strand )
        print "Frameshifted event: "
        for k,v in hash_new_seq_v2.iteritems():
            print "FRAMESHIFT: k = ", k, " - ", v
        print "FRAMESHIFT: Range -> range_full_subseq = ", hash_new_seq_v2['range_full_subseq'] , " || range_mut = ", hash_new_seq_v2['range_mut'], " || range_codon = ", hash_new_seq_v2['range_codon']
        print "FRAMESHIFT: Codon - ", hash_new_seq_v2['orig_codon_nuc_seq'] , "/", hash_new_seq_v2['mut_codon_nuc_seq']
        print "FRAMESHIFT: AA - ", hash_new_seq_v2['orig_codon_aa'] , "/", hash_new_seq_v2['mut_codon_aa']
        print "FRAMESHIFT: Orig Full AA - ", hash_new_seq_v2['orig_aa']
        print "FRAMESHIFT: Mut Full AA - ", hash_new_seq_v2['mut_full_aa']
    else:
        num_codons = 14
        hash_new_seq_v2 = obj_tt.retrieve_alt_inframe( genomic_range, nuc_change, obj_vep.strand, num_codons )
        for k,v in hash_new_seq_v2.iteritems():
            print "INFRAME: k = ", k, " - ", v
        print "IN-FRAME: Range -> range_codon ", hash_new_seq_v2['range_codon'] , " || after_range = ", hash_new_seq_v2['after_range']
        print "IN-FRAME: Codon - ", hash_new_seq_v2['orig_nuc'] , "/", hash_new_seq_v2['mut_nuc']
        print "IN-FRAME: AA - ", hash_new_seq_v2['orig_aa'] , "/", hash_new_seq_v2['mut_aa']


    ##version 1:
    # hash_pos_nuc = retrieve_affected_aa( obj_tt, obj_vep, genomic_range, nuc_change, aa_change )
    # if not in_frame:
    #     # hash_new_seq = retrieve_aa_after_frameshift( obj_tt, hash_pos_nuc )
    #     hash_new_seq_v2 = retrieve_aa_after_frameshift_v2( obj_tt, obj_vep, genomic_range, nuc_change )

    #     print "ORIG GENE INFO: ", obj_tt
    #     print "ORIG gen alt range = ", genomic_range

    #     # print i, ": hash_new_seq - "
    #     # print "range_mod_codon = ", hash_new_seq["range_mod_codon"]
    #     # print "aa_seq_mod_codon = ", hash_new_seq["aa_seq_mod_codon"]
    #     # print "range_after = ", hash_new_seq["range_after"]
    #     # print "aa_seq_after = ", hash_new_seq["aa_seq_after"]
    #     print "--------"
    #     print i, ": hash_new_seq_v2: "
    #     print "mutated_genomic_range = ", hash_new_seq_v2["mutated_genomic_range"]
    #     print "mut_aa_seq = ", hash_new_seq_v2["mut_aa_seq"]
    #     print "orig_aa_seq = ", hash_new_seq_v2["orig_aa_seq"]

    #     # print "obj_tt.arr_genome_pos = ", obj_tt.arr_genome_pos
    #     # print "obj_tt.arr_rf = ", obj_tt.arr_rf
    # else:
    #     print "THIS ALT IS IN_FRAME: count_change = ", count_change, " & in_frame = ", in_frame
    #     aa_len = 16
    #     hash_new_seq_v2 = retrieve_aa_inframe( obj_tt, obj_vep, genomic_range, nuc_change, aa_len )

    #     for k,v in hash_new_seq_v2.iteritems():
    #         print "INFRAME: k = ", k, " - ", v
    #     print "nuc_change = ", hash_new_seq_v2['nuc_seq_orig'], "/", hash_new_seq_v2['nuc_seq_mod_codon'], " || ", codon_change
    #     print "aa_change = ", hash_new_seq_v2['aa_seq_orig'], "/", hash_new_seq_v2['aa_seq_mod_codon'], " || ", aa_change



print "------------ TDD Completed: 170918_VEP_insert_delete.py ------------"