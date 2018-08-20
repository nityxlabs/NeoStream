#/usr/bin/python

import re
import subprocess
import copy

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Isoform import Isoform

# from ExonConnect import ExonConnect
# from ExonConnectMap import ExonConnectMap

delim_exons = '*'       #separates exons ligated by splice junctions
delim_sj = '>>'         #separates splice junctions
class TranscribeTranscript():

    def __init__( self, str_ec, strand_sign, path_genomeidx, hash_transcript_notes ):
        """
        Args:
            str_ec (string exons connected) = string that denotes the exons ligated by splice junctions, in the format  'chrA1:startA1-endA1*chrA2:startA2-endA2>>chrA2:startA2-endA2*chrA3:startA3-endA3>>chrA3:startA3-endA3*chrA4:startA4-endA4>>' and so on
            strand_sign = integer where 1 if '+' or -1 if '-'. Can use def retrieve_exon_from_sj() to create this string from an array of SpliceJunction instances
            path_genomeidx = string that is the path to the samtools-indexed genome
            hash_transcript_notes = hash where key & value are the following:
                -key = "sj", value = array of tuples, where each tuple is the start & end 
            hash_transcript_notes = hash where k = numerical position, v = array of annotations
        NOTE:
            for str_ec, do not need to reverse exon order for minus genes because .reverse_complement takes care of reversing the order the nucleotide sequence
        """
        self.str_ec = str( str_ec )     #if do not convert to string, it will store it as unicode format
        self.arr_pos = TranscribeTranscript.convert_exon_str_to_arr( self.str_ec )

        ##TEST::
        print "TT: arr_pos"
        print self.arr_pos


        self.strand = strand_sign
        self.path_genomeidx = path_genomeidx

        #record events associated with transcript
        self.hash_transcript_notes = hash_transcript_notes
        self.transcript_notes_record_sj()
        
        # self.obj_seq = TranscribeTranscript.compile_dnaseq( self.arr_pos, self.strand, self.path_genomeidx )
        self.hash_seq_notes = self.compile_dna_seq_v2()       #hash_seq_notes -> k = numerical position, v = hash where k2 = category (nucleotide, annotation) & v = corresponding value ('nuc' is the nucleotide base & 'annots' is an array of changes or things to consider)

    def transcript_notes_record_sj( self ):
        """
        Function: records position of SJ in hash 'self.hash_transcript_notes'
        """
        arr_sj = self.find_sj_pos()     #array of hashes that contain information about position 

        for each_sj in arr_sj:
            if not each_sj['start'] in self.hash_transcript_notes:
                self.hash_transcript_notes[ each_sj['start'] ] = []     #will record a set of annotations
            self.hash_transcript_notes[ each_sj['start'] ].append( 'sj_low' )      #sj_low means it is located at a lower genomic position

            if not each_sj['end'] in self.hash_transcript_notes:
                self.hash_transcript_notes[ each_sj['end']] = []     #will record a set of annotations
            self.hash_transcript_notes[ each_sj['end'] ].append( 'sj_high' )       #sj_high means it is located at a higher genomic position


    def find_sj_pos( self ):
        """
        Function: determines the positions of the splice junction based on array "arr_ec"
        """
        arr_sj = []     #will record all SJ positions in hash format 
        for i in range( 0, len( self.arr_pos) - 1 ):
            exon_start = Isoform.split_genome_pos( self.arr_pos[i] )
            exon_end = Isoform.split_genome_pos( self.arr_pos[i+1] )
            hash_sj = {'chrom': exon_start['chrom'], 'start': exon_start['end'], 'end': exon_end['start'] }
            arr_sj.append( hash_sj )

        return arr_sj



    @staticmethod
    def retrieve_exon_from_sj( list_sj, isoform_id ):
        """
        Args:
            list_sj = array of SpliceJunction instances that constitutes an mRNA transcript. This assumes the order the SJs are from least to greatest
        Function: retrieve the exon positions from the array of SpliceJunction instances
        """
        #order from least to greatest - JUST IN CASE, DO NOT REMOVE THIS because this helps sort SpliceJunctions from least to greatest
        # list_sj = sort( list_sj, key = lambda j: j.start )

        #retrieve the exons ligated by the SJs and then return that string
        record_exons_transcript = ""
        for i, each_sj in enumerate( list_sj ):
            
            hash_elems = each_sj.spliced_elems_position_range( isoform_id, False, True )
            if not hash_elems:      #if this is empty, then this means an element (e.g. exon, intron) was not found - skip as this will not be able to be transcribed or translated
                return None
            
            #if not the first element to be added, then add delimiter that separates exons ligated by single SpliceJunction
            if i > 0:
                record_exons_transcript += delim_sj

            record_exons_transcript += hash_elems['prev_elem'].str_genomic_pos( True )
            record_exons_transcript += delim_exons
            record_exons_transcript += hash_elems['next_elem'].str_genomic_pos( True )

        return record_exons_transcript

    @staticmethod
    def retrieve_sj_stat( list_sj, hash_defects = {} ):
        """
        Args:
            list_sj = array of SpliceJunction instances that constitutes an mRNA transcript. This assumes the order the SJs are from least to greatest
            hash_defects = hash where key = integer that is the genomic position & value = information about the defect of interest (e.g. mutation, insertion, deletion). Make sure 
                -mutation: format 'm:B:strand_sign' (e.g. 'm:C:+', 'm:A:-')
                -insertion: format 'in:bases:strand_sign' (e.g. 'in:ACG:+', 'in:GGTC:-')
                -deletion: format 'd:-:strand_sign' (e.g. 'd:-:+', 'd:-:-')
        Function: returns a hash where key = SJ # (could be random, but different from other keys) & value = hash where k2:v2 (sj_range : tuple that is the start & end), (canon : 'True' or 'False'). This function will generate a hash that will contain information for the SJ
        """
        hash_transcript_notes = {}      #k = number that is genomic position, v = array where each element contains
        for each_sj in list_sj:
            if not each_sj.start in hash_transcript_notes:
                hash_transcript_notes[ each_sj.start ] = []
            #record if the SJ is canonical or not
            canon_label_1 = 'sj_after_ab' if not each_sj.canon else 'sj_after'  
            hash_transcript_notes[ each_sj.start ].append( canon_label_1 )

            if not each_sj.end in hash_transcript_notes:
                hash_transcript_notes[ each_sj.end ] = []
            #record if the SJ is canonical or not
            canon_label_2 = 'sj_before_ab' if not each_sj.canon else 'sj_before'  
            hash_transcript_notes[ each_sj.end ].append( canon_label_2 )

        #incorporate each defect into the hash of genomic positions as well
        for k,v in hash_defects.iteritems():        #k = integer that is the genomic position, v = defect of interest
            if not k in hash_transcript_notes:
                hash_transcript_notes[k] = []
            hash_transcript_notes[k].append( v )

        return hash_transcript_notes


    @staticmethod
    def convert_exon_str_to_arr( str_sj_exons ):
        """
        Args:
            str_sj_exons = string that denotes the exons ligated by splice junctions, in the format 'chrA1:startA1-endA1*chrA2:startA2-endA2>>chrA2:startA2-endA2*chrA3:startA3-endA3>>chrA3:startA3-endA3*chrA4:startA4-endA4>>' and so on. DO NOT END with '>>' 
            str_ec = string of ExonConnect objects in the format 'chrA1:startA1-endA1*chrA2:startA2-endA2>>chrA2:startA2-endA2*chrA3:startA3-endA3>>chrA3:startA3-endA3*chrA4:startA4-endA4>>' and so on. DO NOT END with '>>' 
                -Assumes that the positions are sorted from least to greatest
        Function: converts str_pos to an array of genomic ranges, where each element is the genomic range of an exon.
        """
        #split string into array of exons ligated by exon
        arr_sj_exons = str_sj_exons.split( delim_sj )

        #get the first exon, as later in the loop only the second exon will be retrieved
        arr_exons = []
        arr_exons.append( arr_sj_exons[0].split(delim_exons)[0] )

        #go through each exon, and retrieve the latter exon
        for each_sj_exon in arr_sj_exons:
            arr_exons.append( each_sj_exon.split(delim_exons)[1] )

        return arr_exons

    @staticmethod
    def sequence_dna(str_pos, path_genomeidx):
        """ 
        Args:
            str_pos = string that is the position of interest, in the format 'chrom:posStart-posEnd'
            path_genomeidx = string that is the path to the samtools faidx genome file
        Function:
            retrieves the nucleotide sequence from the samtools index file in the path path_genome
        """
        #retrieve the genomic sequence
        exon_seq = subprocess.check_output( ["samtools", "faidx", path_genomeidx, str_pos] )

        #manipulate output to retrieve nucleotide sequence
        exon_seq = re.sub( r"\r\n?", "\n", exon_seq )     #replace any different carriage returns with the unix carriage return
        arr_exon_seq = exon_seq.split("\n")
        del arr_exon_seq[0]            #remove the first element, as this contains ">chr#:posStart-posEnd"
        exon_seq = "".join(arr_exon_seq)

        ##TEST:: print "str_pos = ", str_pos, ": ", Seq(exon_seq).reverse_complement().upper()

        return exon_seq

    ##DELETE - CAN USE Isoform.split_genome_pos
    @staticmethod
    def adjust_position( str_pos ):
        """
        Args:
            str_pos = string that is a genomic position, in the format 'chrom:posStart-posEnd'
        Function: adjusts the position of the string position 'str_pos' by adding +1 to starting position. This is needed when retrieving the nucleotide sequence
        Assumption: 
            -in UCSC genome browser, the starting position of each exon (lower nucleotide position of exon) is not referring to the first nucleotide of the exon but instead the last nucleotide of the intron. Therefore adding +1 gives it the next position 
        """
        hash_exon_pos = Isoform.split_genome_pos( str_pos )
        hash_exon_pos['start'] += 1

        return hash_exon_pos['chrom'] + ':' + str( hash_exon_pos['start'] ) + '-' + str( hash_exon_pos['end'] )


    @classmethod
    def compile_dnaseq( cls_obj, arr_pos, strand_sign, path_genomeidx ):
        """ 
        Args:
            arr_pos = array of positions where each element has a position in the format 'chrom:posStart-posEnd'
            strand_sign = integer where 1 if '+' or -1 if '-'
            path_genomeidx = string that is the path to the samtools faidx genome file
        Function:
            will compile the genomic sequences ( DNA ) from an array contain multiple positions 'arrPos', and returns a Bio.Seq.Seq object
        """
        seq_transcript = ""
        for i, pos in enumerate( arr_pos ):
            adjust_exon_pos = TranslateTranscript.adjust_position( pos )
            seq_transcript += cls_obj.sequence_dna( adjust_exon_pos, path_genomeidx )

        #create BioSeq sequence object
        if strand_sign < 0:     #if strand sign is negative, then get reverse complement
            # seq_transcript = seq_transcript[::-1]     #reverse the string of nucleotides (assuming the exons entered are from right to left aka higher to lower)
            return Seq( seq_transcript.upper(), IUPAC.unambiguous_dna ).reverse_complement()      ##QUES: what is IUPAC.unambiguous_dna?? WHAT DOES IT DO?
        else:       #else strand is positive
            return Seq( seq_transcript.upper(), IUPAC.unambiguous_dna )

    @classmethod
    def demarcate_dna_seq( cls_obj, exon_range, nuc_seq, hash_transcript_notes, hash_seq_notes ):
        """
        Args:
            -exon_range = hash where key:value, 'chrom', 'start', & 'end'
            -nuc_seq = string that is the retrieved nucleotide sequence (from def compile_dna_seq_v2() ). If strand_sign = '-', then the nuc_seq is the reverse complement of the '+' where the genomic position of the first nucleotide > genomic position of the last nucleotide, so will need to flip it here.
            -strand_sign = string that is either '+' or '-', denoting whether the gene is on the plus or minus side of the gene
            -hash_transcript_notes = hash that keeps track of events (e.g. aberrant splicing, mutations, insertions, deletions, etc.) for each position. k = integer position, v = array of events that occur in this position.
            -hash_seq_notes = hash that keeps track of nucleotides at a position and associated annotations, where k = numerical genomic position, v = hash with the following key:value pair - nucleotide:character that is nucleotide, annot:array where each element
        Function: this algorithm will mark specific events (e.g. splicing, aberrant splicing, mutation, etc.) that accompany the nucleotide sequence at a specific position.
        """
        
        #I DON'T NEED THIS TO TAKE CARE OF THE STRAND SIGN - if negative strand sign, then reverse the string so the numerical genomic position is from least to greatest. 
        # if strand_sign == '-':
        #     nuc_seq = nuc_seq[::-1]

        #retrieve features that are contained in this range
        key_pos = [x for x in hash_transcript_notes.keys() if exon_range['start'] <= x <= exon_range['end'] ]
        #traverse along the nucleotide sequence
        for i, each_nuc in enumerate( list(nuc_seq) ):
            key_pos = exon_range['start'] + i
            hash_seq_notes[key_pos] = { 'nuc': each_nuc, 'annots': [] }
            if not key_pos in hash_transcript_notes:
                continue

            for each_note in hash_transcript_notes[key_pos]:
                hash_seq_notes[key_pos]['annots'].append( each_note )

        return hash_seq_notes

    def compile_dna_seq_v2( self ):
        """
        Args:
            -arr_pos = array of exon positions, where each element has format 'chrom:start-end'
            -strand_sign = string that is either '+' or '-', denoting whether the gene is on the plus or minus side of the gene
            -hash_transcript_notes = hash that keeps track of events (e.g. aberrant splicing, mutations, insertions, deletions, etc.) for each position. k = integer that is genomic position, v = array of events that occur in this position.
        Function: retrieves the nucleotide sequences in a position and returns a hash the contains the genomic position, nucleotide, and the associated events (e.g. splicing, aberrant splicing, mutation, etc.) that occur at this position
        """
        hash_seq_notes = {}     #key = integer that is the genomic position, v = hash with following key:value pairs - nucleotide:character that is nucleotide, annot:array of defects associated
        #retrieve the nucleotide sequence from a given exon
        for each_pos in self.arr_pos:
            # adjust_exon_pos = TranslateTranscript.adjust_position( pos )
            exon_range = Isoform.split_genome_pos( each_pos )       #returns hash 
            exon_seq = TranscribeTranscript.sequence_dna( each_pos, self.path_genomeidx )
            #create BioSeq sequence object
            if self.strand < 0:     #if strand sign is negative, then get reverse complement
                ##METHOD 1 - retrieve complementary sequence
                exon_seq = Seq( exon_seq.upper(), IUPAC.unambiguous_dna ).complement()      ##QUES: what is IUPAC.unambiguous_dna?? WHAT DOES IT DO?
                ##METHOD 2 - retrieve complementary sequence:
                # exon_seq = Seq( exon_seq.upper(), IUPAC.unambiguous_dna ).reverse_complement()      ##QUES: what is IUPAC.unambiguous_dna?? WHAT DOES IT DO?
                # exon_seq = exon_seq[::-1]       #reverse the order of the nucleotides so the genomic position refers to the correct nucleotide sequence

            hash_seq_notes = TranscribeTranscript.demarcate_dna_seq( exon_range, exon_seq, self.hash_transcript_notes, hash_seq_notes )

        return hash_seq_notes

    # @classmethod
    # def compile_dna_seq_v2( cls_obj, arr_pos, strand_sign, path_genomeidx, hash_transcript_notes ):
    #     """
    #     Args:
    #         -arr_pos = array of exon positions, where each element has format 'chrom:start-end'
    #         -strand_sign = string that is either '+' or '-', denoting whether the gene is on the plus or minus side of the gene
    #         -hash_transcript_notes = hash that keeps track of events (e.g. aberrant splicing, mutations, insertions, deletions, etc.) for each position. k = integer that is genomic position, v = array of events that occur in this position.
    #     Function: retrieves the nucleotide sequences in a position and returns a hash the contains the genomic position, nucleotide, and the associated events (e.g. splicing, aberrant splicing, mutation, etc.) that occur at this position
    #     """
    #     hash_seq_notes = {}     #key = integer that is the genomic position, v = hash with following key:value pairs - nucleotide:character that is nucleotide, annot:array of defects associated
    #     #retrieve the nucleotide sequence from a given exon
    #     for each_pos in arr_pos:
    #         # adjust_exon_pos = TranslateTranscript.adjust_position( pos )
    #         exon_range = Isoform.split_genome_pos( each_pos )
    #         exon_seq = TranscribeTranscript.sequence_dna( each_pos, path_genomeidx )
    #         #create BioSeq sequence object
    #         if strand < 0:     #if strand sign is negative, then get reverse complement
    #             exon_seq = Seq( exon_seq.upper(), IUPAC.unambiguous_dna ).reverse_complement()      ##QUES: what is IUPAC.unambiguous_dna?? WHAT DOES IT DO?

    #         hash_seq_notes = cls_obj.demarcate_dna_seq( exon_range, exon_seq, hash_transcript_notes, hash_seq_notes )

    #     return hash_seq_notes

    def compile_dna_seq_only( self ):
        """
        Function: returns only the nucleotide sequence record in 'hash_seq_notes'
        """
        return ''.join( [v['nucleotide'] for k,v in self.hash_seq_notes.iteritems()] )


    ##NOTE: May need to adjust this for strand sign (self.strand)
    @classmethod
    def apply_seq_changes( cls_obj, hash_seq_notes ):
        """
        Args:
            orig_transcript = hash where k = numerical position, v = hash where k2 = category (nucleotide, annotation) & v = corresponding value (nucleotide is a character & 'annots' is an array of changes or things to consider)
            strand = string that is either '+' or '-'. This will determine how to treat the changes
        Function: applies the changes recorded in the hash 
        """
        hash_transcript = copy.deepcopy( hash_seq_notes )      #use copy.deepcopy() to copy hashes, especially hashes nested within hashes
        #go through each element of the hash, and make the changes recorded
        for k,v in hash_transcript.iteritems():     #k = numerical position, v = hash where k2 = category (nucleotide, annotation) & v = corresponding value (nucleotide is a character & 'annots' is an array of changes or things to consider)
            for each_change in v['annots']:
                arr_change = each_change.split(':')
                if arr_change[0].lower() == 'm':        #for mutations
                    v['nuc'] = arr_change[1]
                elif arr_change[0].lower() == 'in':     #for insertions
                    v['nuc'] = v['nuc'] + arr_change[1]
                elif arr_change[0].lower() == 'd':      #for deletions
                    v['nuc'] = arr_change[1]         #make sure it is a non-word character, such as '-'

        return hash_transcript


    @classmethod
    def display_seq_changes_nucleotide( self ):
        """
        Args:
            hash_seq_notes = hash where k = numerical position, v = hash where k2 = category (nucleotide, annotation) & v = corresponding value (nucleotide is a character & 'annots' is an array of changes or things to consider)
            strand = string that is either '+' or '-'. This will determine how to treat the changes
        Function: this will mark up the nucleotide sequence and the amino acid sequence (if annot_aa is True)
        """
        #make all changes to the nucleotide
        hash_transcript = copy.deepcopy( self.hash_seq_notes )

        str_mod_nuc = ""
        for k,v in hash_transcript.iteritems():     #k = numerical position, v = hash where k2 = category (nucleotide, annotation) & v = corresponding value (nucleotide is a character & 'annots' is an array of changes or things to consider)
            #if no annotations, then no changes that need to be recorded
            if not v['annots']:
                str_mod_nuc += v['nucleotide']
                continue

            #make all necessary changes to nucleotide base
            sj_stat = None      #None = no SJ, 1 = before nucleotide position, 2 = after nucleotide position, 3 = aberrant before nucleotide position, 4 = aberrant after nucleotide position
            record_changes = ""
            for each_change in v['annots']:
                arr_change = each_change.split(':')
                if arr_change[0].lower() == 'sj_before':
                    sj_stat = 1 
                elif arr_change[0].lower() == 'sj_after':
                    sj_stat = 2
                elif arr_change[0].lower() == 'sj_before_ab':
                    sj_stat = 3
                elif arr_change[0].lower() == 'sj_after_ab':
                    sj_stat = 4
                elif arr_change[0].lower() == 'm':        #for mutations
                    # v['nucleotide'] = arr_change[1]
                    record_changes = "~mut(" + v['nucleotide'] + "->" + arr_change[1] + ")" if not record_changes else ":mut(" + v['nucleotide'] + "->" + arr_change[1] + ")"
                elif arr_change[0].lower() == 'in':     #for insertions
                    # v['nucleotide'] = v['nucleotide'] + arr_change[1]
                    record_changes = "~add(" + arr_change[1] + ")" if not record_changes else ":add(" + arr_change[1] + ")"
                elif arr_change[0].lower() == 'd':      #for deletions
                    # v['nucleotide'] = arr_change[1]         #make sure it is a non-word character, such as '-'
                    record_changes = "~del(" + v['nucleotide'] + ")" if not record_changes else ":del(" + v['nucleotide'] + ")"

            #record version of string if 'record_changes' contains information
            record_nuc = v['nucleotide'] if not record_changes else record_changes + "[" + v['nucleotide'] + "]~"

            #check the splicing status
            if not sj_stat:
                str_mod_nuc += record_nuc
            elif sj_stat == 1:      #this mean SJ is before nucleotide
                str_mod_nuc += "|" + record_nuc
            elif sj_stat == 2:      #this mean SJ is after nucleotide
                str_mod_nuc += record_nuc + "|"
            elif sj_stat == 3:      #this mean aberrant SJ is before nucleotide
                str_mod_nuc += "//" + record_nuc
            elif sj_stat == 4:      #this mean aberrant SJ is after nucleotide
                str_mod_nuc += record_nuc + "//"

        return str_mod_nuc

    @staticmethod
    def search_nuc_defects( mod_transcript, i_start, i_end, strand ):
        """
        Args:
            -mod_transcript = hash where k = numerical position, v = hash where k2 = category (nucleotide, annotation) & v = corresponding value (nucleotide is a character & 'annots' is an array of changes or things to consider)
            strand = string that is either '+' or '-'. This will determine how to treat the changes
            -i_start & i_end = integer that is the starting & ending index for mod_transcript, respectively. These will be used to retrieve the genomic position at the index location and access the annotations in those positions 
            -strand = integer where 1 means '+' gene & -1 means '-' gene
        Function: returns the defects associated with the range of positions with nucleotide positions in hash 'mod_transcript'. NOTE that this function will be usually be used in conjunction with def display_seq_changes_aa()
        """
        list_annot = []
        i_end_2 = i_end + 1 if strand > 0 else i_end - 1        #need to adjust for strand sign
        for i in range( i_start, i_end_2, strand ):     #strand here is used as the step
            genome_pos = mod_transcript.keys()[i]       #retrieve the genome position of interest
            for each_change in mod_transcript[genome_pos]['annots']:
                arr_change = each_change.split(':')
                if 'sj_' in arr_change[0].lower():
                    list_annot.append( arr_change[0] )
                elif arr_change[0].lower() == 'm':        #for mutations
                    list_annot.append( 'mut' )
                elif arr_change[0].lower() == 'in':     #for insertions
                    list_annot.append( 'insert' )
                elif arr_change[0].lower() == 'd':      #for deletions
                    list_annot.append( 'delete' )

        return ':'.join( list_annot )


    @classmethod
    def display_seq_changes_aa( cls_obj, hash_seq_notes, strand ):
        """
        Args:
            orig_transcript = hash where k = numerical position, v = hash where k2 = category (nucleotide, annotation) & v = corresponding value (nucleotide is a character & 'annots' is an array of changes or things to consider)
            strand = string that is either '+' or '-'. This will determine how to treat the changes
            -p_start = integer that is the absolute genomic starting position 
            -strand = strand = integer where 1 means '+' gene & -1 means '-' gene
                -if '-', will need reverse complement of nucleotide strand
        Function: annotates the amino acid
        """
        mod_transcript = cls_obj.apply_seq_changes( hash_seq_notes )
        #translate the nucleotide seqeuence
        transcript_seq = ''.join( [v['nucleotide'] for k,v in mod_transcript.iteritems()] ).replace( '-', '' )
        if strand < 0:      #if on minus strand, then need to flip since the transcript will be translated from right to left, not left to right
            transcript_seq = transcript_seq[::-1]
        aa_seq = Seq( str( transcript_seq ) ).translate( to_stop = True )

        annot_aa = ""
        for i, each_aa in enumerate( aa_seq ):
            if strand < 0:
                aa_nuc_start = ( len( mod_transcript.keys() ) - i - 1 ) * 3
                aa_nuc_end = aa_nuc_start - 3 if ( aa_nuc_start - 3 ) >= 0 else 0
            else:
                aa_nuc_start = i * 3
                aa_nuc_end = aa_nuc_start + 3 if (aa_nuc_start + 3) < len( mod_transcript.keys() ) else len( mod_transcript.keys() )

            str_annots = TranscribeTranscript.search_nuc_defects( mod_transcript, aa_nuc_start, aa_nuc_end, strand )
            if not str_annots:
                annot_aa += each_aa
            else:
                annot_aa += str_annots + "[" + each_aa + "]"

        return annot_aa


    @staticmethod
    def search_nuc_defects_v3( mod_transcript, list_pos ):
        """
        Args:
            -mod_transcript = hash where k = numerical position, v = hash where k2 = category (nucleotide, annotation) & v = corresponding value (nucleotide is a character & 'annot' is an array of changes or things to consider)
            strand = string that is either '+' or '-'. This will determine how to treat the changes
            -list_pos = array of positions of interest. These values should correspond to keys in hash 'mod_transcript'
        Function: Same as search_nuc_defects() but uses relative positions for pos_start & pos_end (e.g. )
        """
        list_annot = []
        for each_pos in list_pos:
            for each_change in mod_transcript[each_pos]['annot']:
                arr_change = each_change.split(':')
                if 'sj_' in arr_change[0].lower():
                    list_annot.append( arr_change[0] )
                elif arr_change[0].lower() == 'm':        #for mutations
                    list_annot.append( 'mut' )
                elif arr_change[0].lower() == 'in':     #for insertions
                    list_annot.append( 'insert' )
                elif arr_change[0].lower() == 'd':      #for deletions
                    list_annot.append( 'delete' )

        return ':'.join( list_annot )

    @classmethod
    def display_seq_changes_aa_V2( cls_obj, orig_transcript, strand = None ):
        """
        Args:
            orig_transcript = hash where k = numerical position, v = hash where k2 = category (nucleotide, annotation) & v = corresponding value (nucleotide is a character & 'annots' is an array of changes or things to consider)
            strand = string that is either '+' or '-'. This will determine how to treat the changes
            -strand = integer that is 1 if '+' strand or -1 if '-' strand
        Function: this function will translate and annotate amino acid sequence. This is similar to def display_seq_changes_aa(), but doesn't need to define starting position
        """
        mod_transcript = cls_obj.apply_seq_changes( orig_transcript )

        #retrieve all keys from hash "orig_transcript", and sort based on strand sign (+ strand: ascending order, - strand: descending order)
        key_pos = orig_transcript.keys() if strand > 0 else sorted( orig_transcript.keys(), reverse = True )

        #go through each codon, translate, and record annotations for each amino acid
        i = 0
        annot_aa = ""       #this will record the annotated amino acid sequence
        while i < len(key_pos):

            #retrieve the nucleotide sequences that make up the codon, and all the genomic events that occur within this range
            list_pos = []       #records keys from mod_transcript (these are the genomic positions)
            codon_seq = []      #records nucleotides that make up the a codon, which should be 3 nucleotides long
            while len( codon_seq ) < 3 and i < len(key_pos):
                list_pos.append( key_pos[i] )
                if mod_transcript[ key_pos[i] ]['nucleotide'] in ['A', 'T', 'C', 'G']:
                    codon_seq.append( mod_transcript[ key_pos[i] ]['nucleotide'] )

                #increment counter - this will keep track of which position to start next 
                i = i + 1

            #retrieve the amino acid
            str_codon_seq = ''.join( codon_seq )
            if len( codon_seq ) == 3:
                obj_dna = simpleDNA( str_codon_seq )
                single_aa = ''.join( obj_dna.convert_aa() )
            else:
                single_aa = '-'

            #retrieve annotations associated amino acid (e.g. mutation, insertion, deletion, SJ)
            str_annots = search_nuc_defects_v3( mod_transcript, list_pos )
            if not str_annots:
                annot_aa += single_aa
            else:
                annot_aa += str_annots + "[" + single_aa + "]"
            

        return annot_aa

    
    @staticmethod
    def convert_complementary_bases( list_bases ):
        """
        Args:
            list_bases = array where each element is a nucleotide base. Make sure all bases are upper case
        Function: converts bases to complementary bases
        """
        hash_complementary_base = {'A':'T', 'T':'A', 'G':'C', 'C':'G', '-':'-'}     #NOTE: '-' is there just in case
        return [hash_complementary_base[x] for x in list_bases]

    @staticmethod
    def check_defect_sign( hash_defects, strand ):
        """
        Args:
            hash_defects = hash where k = integer that is genomic position, v = array of defects
            strand = integer where 1 is "+ strand" and -1 is "- strand"
        Function: check defects with respect to strand sign, and adjust for strand sign
        """
        str_strand = '+' if strand > 0 else '-'     #string for comparing which strand the defect is on

        hash_adj_defects = {}       #this hash will record the new, adjusted defects based on strand sign
        for k,v in hash_defects.iteritems():     #k = integer that is genomic position, v = array of defects
            list_adj_defects = []       #this will record all defects, those that are adjusted for strand sign as well as unmodified defects
            for each_defect in v:       #each_defect = string in the format type:change:strand_sign
                d_info = each_defect.split(':')     #0 = defect type (e.g. mutation, insertion, deletion), 1 = change to make, 2 = strand sign
                if d_info[2] != strand:
                    #change bases to complementary bases IF mutation, insertion
                    comp_strand = ''.join( convert_complementary_bases(list(d_info[1])) ) if d_info[0] != 'd' else d_info[1]
                    adj_defect = d_info[0] + ':' + comp_strand + ':' + strand
                    list_adj_defects.append( adj_defect )
                else:
                    list_adj_defects.append( each_defect )

            #record the array of defects to the new hash
            hash_adj_defects[k] = list_adj_defects

        return hash_adj_defects


    @staticmethod
    def sequence_transcribe():
        """ transcribes DNA sequence into RNA """
        return Seq.transcribe()



class TranslateTranscript( TranscribeTranscript ):
    #Protocol for translating transcript
    # *Perhaps should make this a class that records all possible translated protein, maybe in "TranslateTranscript"
    # -STEP: compile the entire nucleotide sequence
    # -STEP: find potential start codon sites
    #     -find exon & the relative position within the exon
    #     -also get the absolute position so i can use this to translate the protein
    # -STEP: from the potential start codon sites, retrieve the nucleotide sequence & translate it - split based on the first occurrence of a stop codon
    # -STEP: quantify the number of amino acids per possible protein translated
    # -QUES: should I have a way to find the start & end position of the protein based on exons? I think so, but how will I define the end exon? 
    #     -CONJ 1 (NOPE): perhaps should go through each exon and translate the proteins - NO WON'T WORK because exons are not discrete units for protein domains
    #     -CONJ 2: based on the start position & the number amino acids, back calculate to exon position - could this be error prone?
    #         -could subtract the length of each exon until subtraction leads to negative number, which means the last amino acid lands within the exon
    #NOTE: stated that "coding regions of transcript tend to be longer than by chance", which means what is the random chance of bumping into a stop codon? It is 3/64, which means usually 21.333 amino acids will have to be translated before a stop codon is encountered

    def __init__( self, str_pos, strand_sign, path_genomeidx, hash_transcript_notes ):
        TranscribeTranscript.__init__( self, str_pos, strand_sign, path_genomeidx, hash_transcript_notes )

    def translate_transcript( self ):
        """
        Args:
            arrPos = array of positions where each element has a position in the format 'chrom:posStart-posEnd'
        Function:
            translates nucleotide sequence into protein sequence, and returns an array of the longest amino acid sequences
        """

        #find all potential start positions
        pos_startcodons = self.find_start_codons()

        if not pos_startcodons:
            return []

        ##QUES: DO I NEED THIS??
        pos_startcodons.append(-1)      #for finding last position

        #find the longest protein sequence
        final_protein = Seq( str( self.obj_seq[pos_startcodons[0]:-1] ) ).translate( to_stop = True )
        all_proteins = []       #records all proteins that are 

        ##TEST:: show protein sequence at first start codon found
        # print "FIRST: start = ", pos_startcodons[0], " & protein length = ", len(final_protein)
        # print "FIRST: nuc = ", self.obj_seq[pos_startcodons[0] : -1]
        # print "FIRST: protein = ", final_protein

        for i in range( 1, len(pos_startcodons) ):
            
            #.translate(to_stop) will only show protein sequence up to the stop signal
            protein_seq = Seq( str( self.obj_seq[pos_startcodons[i]:-1] ) ).translate( to_stop = True )


            ##TEST:: show protein sequences subsequent start codons found
            # print "start = ", pos_startcodons[i], " & protein length = ", len(protein_seq)
            # print "nuc = ", self.obj_seq[pos_startcodons[i] : -1]
            # print "protein = ", protein_seq

            #record longest protein sequence
            if len( protein_seq ) > len( final_protein ):
                final_protein = protein_seq
                all_proteins = []       #reset all recorded proteins
            #else if match in length, then just append to array
            elif len( protein_seq ) == len( final_protein ):
                all_proteins.append( protein_seq )

        all_proteins.insert( 0, final_protein )

        return all_proteins


    def find_start_codons( self ):
        """
        Args:
            obj_seq = Seq object from Biopython. This object contains the nucleotide sequences
        Function: finds all potential positions with a start codon
        """
        #pick start codon
        start_codon = 'ATG' if 'T' in self.obj_seq.upper() else 'AUG'

        #return all positions where start codon occurs
        pos = 0
        start_pos = []
        while self.obj_seq.upper().find( start_codon, pos ) > -1:
            next_pos = self.obj_seq.upper().find( start_codon, pos )
            start_pos.append( next_pos )

            pos = next_pos + 1

        return start_pos
