#/usr/bin/python

import re
import subprocess

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Isoform import Isoform

delim_exons = '*'       #separates exons ligated by splice junctions
delim_sj = '>>'         #separates splice junctions
class SimpleTranscribe( object ):

    def __init__( self, str_ec, strand_sign, path_genomeidx ):
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
        self.arr_pos = SimpleTranscribe.convert_exon_str_to_arr( self.str_ec )

        self.strand = strand_sign
        self.path_genomeidx = path_genomeidx

        self.obj_seq = self.compile_dnaseq()


    @staticmethod
    def convert_exon_str_to_arr( str_sj_exons ):
        """
        Args:
            str_sj_exons = string that denotes the exons ligated by splice junctions, in the format 'chrA1:startA1-endA1*chrA2:startA2-endA2>>chrA2:startA2-endA2*chrA3:startA3-endA3>>chrA3:startA3-endA3*chrA4:startA4-endA4>>' and so on. DO NOT END with '>>' 
            str_ec = string of ExonConnect objects in the format 'chrA1:startA1-endA1*chrA2:startA2-endA2>>chrA2:startA2-endA2*chrA3:startA3-endA3>>chrA3:startA3-endA3*chrA4:startA4-endA4>>' and so on. DO NOT END with '>>' 
        Function: converts str_pos to an array of genomic ranges, where each element is the genomic range of an exon.
        """
        #split string into array of exons ligated by exon
        arr_sj_exons = str_sj_exons.split( delim_sj )

        #get the first exon, as later in the loop only the second exon will be retrieved
        arr_exons = []
        arr_exons.append( arr_sj_exons[0].split(delim_exons)[0] )

        ##TEST:: print "arr_sj_exons = ", arr_sj_exons

        #go through each exon, and retrieve the latter exon
        for each_sj_exon in arr_sj_exons:
            arr_exons.append( each_sj_exon.split(delim_exons)[1] )

        return arr_exons

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

    @staticmethod
    def sequence_dna( str_pos, path_genomeidx ):
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

    #I'M NOT SURE IF I NEED THIS FUNCTION BECAUSE OF def compile_dnaseq()
    # @classmethod
    # def retrieve_nuc_seq( cls_obj, pos, path_genomeidx, reverse_comp ):
    #     """
    #     Args:
    #         pos = string that is the genomic position, in the format chrom:start-end
    #         path_genomeidx = string that is the path to the samtools-indexed genome
    #     Function: will retrieve nucleotide sequence
    #     """
    #     adjust_exon_pos = cls_obj.adjust_position( pos )
    #     seq_transcript += cls_obj.sequence_dna( adjust_exon_pos, path_genomeidx )

    #     #create BioSeq sequence object - note, can take the reverse complement using Seq( seq_transcript.upper(), IUPAC.unambiguous_dna ).reverse_complement()
    #     return Seq( seq_transcript.upper(), IUPAC.unambiguous_dna )


    def compile_dnaseq( self ):
        """ 
        Function:
            will compile the genomic sequences ( DNA ) from an array contain multiple positions 'arrPos', and returns a Bio.Seq.Seq object
        """
        seq_transcript = ""
        for i, pos in enumerate( self.arr_pos ):
            adjust_exon_pos = SimpleTranscribe.adjust_position( pos )
            seq_transcript += SimpleTranscribe.sequence_dna( adjust_exon_pos, self.path_genomeidx )
            seq_transcript_TEST = SimpleTranscribe.sequence_dna( pos, self.path_genomeidx )

            ##TEST::
            # print "DNASeq adjusted part ", i, " - ", adjust_exon_pos, " -> Seq = ", seq_transcript 
            # print "DNASeq normal part ", i, " - ", pos, " -> Seq = ", seq_transcript_TEST

        #create BioSeq sequence object
        if self.strand < 0:     #if strand sign is negative, then get reverse complement
            # seq_transcript = seq_transcript[::-1]     #reverse the string of nucleotides (assuming the exons entered are from right to left aka higher to lower)
            return Seq( seq_transcript.upper(), IUPAC.unambiguous_dna ).reverse_complement()      ##QUES: what is IUPAC.unambiguous_dna?? WHAT DOES IT DO?
        else:       #else strand is positive
            return Seq( seq_transcript.upper(), IUPAC.unambiguous_dna )

    @classmethod
    def compile_dnaseq_array( cls_obj, arr_pos, path_genomeidx, strand_sign ):
        """
        Args:
            -arr_pos = array of genomic positions (format = chrom:start-end)
            -path_genomeidx = string that is the path to the samtools-indexed genome
            -strand_sign = integer where 1 if '+' or -1 if '-'. Can use def retrieve_exon_from_sj() to create this string from an array of SpliceJunction instances
        Function: Similar to def compile_dnaseq(), but will return an array of DNA sequences instead of string.
            -the DNA sequence for each position saved in array "arr_pos", and returns the nucleotide sequence for each nucleotide range. Note that if "strand_sign" is minus, it will return the reverse complement sequence and will reverse the array, meaning:
               -first nucleotide range will be the last element in the returned array
               -last nucleotide range will be the first element in the returned array
        """
        arr_nuc = []        #records nucleotide sequence
        for i, pos in enumerate( arr_pos ):
            adjust_exon_pos = cls_obj.adjust_position( pos )
            seq_dna = cls_obj.sequence_dna( adjust_exon_pos, path_genomeidx )
            # obj_seq = Seq( seq_dna.upper(), IUPAC.unambiguous_dna ).complement() if strand_sign < 0 else Seq( seq_dna.upper(), IUPAC.unambiguous_dna )
            obj_seq = Seq( seq_dna.upper(), IUPAC.unambiguous_dna ).reverse_complement() if strand_sign < 0 else Seq( seq_dna.upper(), IUPAC.unambiguous_dna )

            #append nucleotide sequence 
            arr_nuc.append( obj_seq )

        return arr_nuc[::-1] if strand_sign < 0 else arr_nuc



class SimpleTranslate( SimpleTranscribe ):

    def __init__( self, str_ec, strand_sign, path_genomeidx ):
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
        super( SimpleTranslate, self ).__init__( str_ec, strand_sign, path_genomeidx )


    def translate_seq( self ):
        """
        Function:
            translates nucleotide sequence into protein sequence, and returns an array of the longest amino acid sequences
        """
        return Seq( str(self.obj_seq) ).translate( to_stop = False )

