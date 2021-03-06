#/usr/bin/python

import re
import subprocess

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from ExonConnect import ExonConnect
from ExonConnectMap import ExonConnectMap

class TranscribeTranscript():

    def __init__(self, str_ec, strand_sign, path_genomeidx):
        """
        Args:
            str_ec = string of ExonConnect objects in the format 'chrA1:startA1-endA1*chrA2:startA2-endA2>>chrA2:startA2-endA2*chrA3:startA3-endA3>>chrA3:startA3-endA3*chrA4:startA4-endA4>>' and so on
            strand_sign = integer where 1 if '+' or -1 if '-'
            path_genomeidx = string that is the path to the samtools-indexed genome
        NOTE:
            for str_ec, do not need to reverse exon order for minus genes because .reverse_complement takes care of reversing the order the nucleotide sequence
        """
        self.str_ec = str(str_ec)     #if do not convert to string, it will store it as unicode format
        self.arr_pos = TranscribeTranscript.create_arr_pos(self.str_ec)
        self.strand = strand_sign
        self.path_genomeidx = path_genomeidx
        self.obj_seq = TranscribeTranscript.compile_dnaseq(self.arr_pos, self.strand, self.path_genomeidx)

    @staticmethod
    def create_arr_pos(str_ec):
        """
        Args:
            str_ec = string of ExonConnect objects in the format 'chrA1:startA1-endA1*chrA2:startA2-endA2>>chrA2:startA2-endA2*chrA3:startA3-endA3>>chrA3:startA3-endA3*chrA4:startA4-endA4>>' and so on. DO NOT END with '>>' 
        Function: converts str_pos to an array of positions
        """
        #split string of ExonConnect strings into individual strings
        arr_ec = str_ec.split(ExonConnectMap.delimMultiEC)

        #get the first exon from the first ExonConnect, as later in the loop only the second exon will be retrieved
        arr_pos = []
        arr_pos.append(arr_ec[0].split(ExonConnect.delimEC)[0])
        
        #go through each ExonConnect, and retrieve the latter exon
        for each_ec in arr_ec:
            arr_pos.append(each_ec.split(ExonConnect.delimEC)[1])

        return arr_pos



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
        exon_seq = subprocess.check_output(["samtools", "faidx", path_genomeidx, str_pos])

        #manipulate output to retrieve nucleotide sequence
        exon_seq = re.sub(r"\r\n?", "\n", exon_seq)     #replace any different carriage returns with the unix carriage return
        arr_exon_seq = exon_seq.split("\n")
        del arr_exon_seq[0]            #remove the first element, as this contains ">chr#:posStart-posEnd"
        exon_seq = "".join(arr_exon_seq)

        ##TEST:: print "str_pos = ", str_pos, ": ", Seq(exon_seq).reverse_complement().upper()

        return exon_seq

    @staticmethod
    def adjust_position(str_pos):
        """
        Args:
            str_pos = string that is a genomic position, in the format 'chrom:posStart-posEnd'
        Function: adjusts the position of the string position 'str_pos' by adding +1 to starting position. This is needed when retrieving the nucleotide sequence
        Assumption: 
            -in UCSC genome browser, the starting position of each exon (lower nucleotide position of exon) is not referring to the first nucleotide of the exon but instead the last nucleotide of the intron. Therefore adding +1 gives it the next position 
        """
        #split by ':' and then by '-', and add +1 to starting position
        chrom = str_pos.split(':')[0]
        pos_start = str( int( str_pos.split(':')[1].split('-')[0] ) + 1 )
        pos_end = str_pos.split(':')[1].split('-')[1]

        return chrom + ':' + pos_start + '-' + pos_end

    @classmethod
    def compile_dnaseq(clsObj, arr_pos, strand_sign, path_genomeidx):
        """ 
        Args:
            arr_pos = array of positions where each element has a position in the format 'chrom:posStart-posEnd'
            strand_sign = integer where 1 if '+' or -1 if '-'
            path_genomeidx = string that is the path to the samtools faidx genome file
        Function:
            will compile the genomic sequences ( DNA ) from an array contain multiple positions 'arrPos', and returns a Bio.Seq.Seq object
        """
        seqTranscript = ""
        for i, pos in enumerate(arr_pos):
            pos_adjust = clsObj.adjust_position(pos)
            seqTranscript += clsObj.sequence_dna(pos_adjust, path_genomeidx)

        #create BioSeq sequence object
        if strand_sign < 0:     #if strand sign is negative, then get reverse complement
            # seqTranscript = seqTranscript[::-1]     #reverse the string of nucleotides (assuming the exons entered are from right to left aka higher to lower)
            return Seq(seqTranscript.upper(), IUPAC.unambiguous_dna).reverse_complement()      ##QUES: what is IUPAC.unambiguous_dna?? WHAT DOES IT DO?
        else:       #else strand is positive
            return Seq(seqTranscript.upper(), IUPAC.unambiguous_dna)


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

    def __init__(self, str_pos, strand_sign, path_genomeidx):
        TranscribeTranscript.__init__(self, str_pos, strand_sign, path_genomeidx)

    def translate_transcript(self):
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
        final_protein = Seq(str(self.obj_seq[pos_startcodons[0] : -1])).translate(to_stop = True)
        all_proteins = []       #records all proteins that are 

        ##TEST:: show protein sequence at first start codon found
        # print "FIRST: start = ", pos_startcodons[0], " & protein length = ", len(final_protein)
        # print "FIRST: nuc = ", self.obj_seq[pos_startcodons[0] : -1]
        # print "FIRST: protein = ", final_protein

        for i in range( 1, len(pos_startcodons) ):
            
            #.translate(to_stop) will only show protein sequence up to the stop signal
            protein_seq = Seq(str(self.obj_seq[pos_startcodons[i] : -1])).translate(to_stop = True)


            ##TEST:: show protein sequences subsequent start codons found
            # print "start = ", pos_startcodons[i], " & protein length = ", len(protein_seq)
            # print "nuc = ", self.obj_seq[pos_startcodons[i] : -1]
            # print "protein = ", protein_seq

            #record longest protein sequence
            if len(protein_seq) > len(final_protein):
                final_protein = protein_seq
                all_proteins = []       #reset all recorded proteins
            #else if match in length, then just append to array
            elif len(protein_seq) == len(final_protein):
                all_proteins.append(protein_seq)

        all_proteins.insert(0, final_protein)

        return all_proteins


    def find_start_codons(self):
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
        while self.obj_seq.upper().find(start_codon, pos) > -1:
            next_pos = self.obj_seq.upper().find(start_codon, pos)
            start_pos.append(next_pos)

            pos = next_pos + 1

        return start_pos
