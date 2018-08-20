#/usr/bin/python

import re
import subprocess
import copy

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Isoform import Isoform
from Exon import Exon

# from ExonConnect import ExonConnect
# from ExonConnectMap import ExonConnectMap

delim_exons = '*'       #separates exons ligated by splice junctions
delim_sj = '>>'         #separates splice junctions
class TranscribeTranscript():

    def __init__( self, transcript_sj, iso_sj, path_genomeidx, hash_transcript_notes ):
        """
        Args:
            -transcript_sj = an array of SpliceJunction instances that constitute a transcript. The SpliceJunctions should be sorted from lowest numerical position to highest numerical position.
            -iso_sj = Instance of the IsoformSJ class
            -path_genomeidx = string that is the path to the samtools-indexed genome
            -hash_transcript_notes = hash where key & value are the following:
                -key = "sj", value = array of tuples, where each tuple is the start & end 
            -hash_transcript_notes = hash where k = numerical position, v = array of annotations
        NOTE:
            for str_ec, do not need to reverse exon order for minus genes because .reverse_complement takes care of reversing the order the nucleotide sequence
        """
        self.transcript_sj = sorted( transcript_sj, key = lambda s: s.start, reverse = False )      #this is a list of splicing positions
        self.iso_sj = iso_sj        #Instance of the IsoformSJ class

        #self.list_exons = list of exons ordered by genomic position (least to greatest)
        if self.transcript_sj:      #if splice junctions are recorded for isoform
            self.list_exons = TranscribeTranscript.create_exons_from_sj( self.transcript_sj, self.iso_sj )  
        else:
            self.list_exons = TranscribeTranscript.create_exons_no_sj( self.iso_sj )
        self.path_genomeidx = path_genomeidx

        #record events associated with transcript
        self.hash_transcript_notes = hash_transcript_notes
        # self.transcript_notes_record_sj()           #not sure if I need this, could just use "self.transcript_sj"

        # self.obj_seq = TranscribeTranscript.compile_dnaseq( self.arr_pos, self.iso_sj.strand, self.path_genomeidx, True, True )
        # self.hash_seq_notes = self.compile_dna_seq_v2()       #hash_seq_notes -> k = numerical position, v = hash where k2 = category (nucleotide, annotation) & v = corresponding value ('nuc' is the nucleotide base & 'annots' is an array of changes or things to consider)
        """
        self.arr_genome_pos = an array where each element is a genomic position. 
        self.arr_nuc_seq = an array where each element is a nucleotide. The reason for this is so each nucleotide will have a reading frame associated with it (based on self.arr_rf)
        self.arr_rf = an array where each element is a number 0, 1, or 2. Each is a reading frame that corresponds to each a specific position
        """
        self.arr_genome_pos, self.arr_nuc_seq, self.arr_rf = self.create_array_indices()

    def test_display_all_indices( self ):
        """
        Function: just a test function that shows all indices created by def create_array_indices() - genomic position, nucleotide sequence, & reading frame (self.arr_genome_pos, self.arr_nuc_seq, & self.arr_rf, respectively)
        """
        print "show indices of genomic position, nucleotide sequence, & reading frame"
        for i in range( 0, len(self.arr_genome_pos) ):
            print i, ": pos = ", self.arr_genome_pos[i], " | base = ", self.arr_nuc_seq[i], " | frame = ", self.arr_rf[i]

    def test_display_indices_range( self, genome_start, num_after ):
        """
        Args:
            genome_start = integes that are genomic positions
            num_after = integer that is:
                -if positive, the number of indices after genome_start
                -if negative, the number of indices before genome_start
        Function: similar to TranscribeTranscript.test_display_all_indices(), but takes into account a starting genomic position & ending genomic position
        """
        #retrieve the indices corresponding to position
        i_start = self.arr_genome_pos.index( genome_start )
        index_range = range( i_start, i_start + num_after + 1 ) if num_after > 0 else range( i_start + num_after, i_start + 1 )     #added +1 as to consider the last position as well
        for i in index_range:
            print i, ": pos = ", self.arr_genome_pos[i], " | base = ", self.arr_nuc_seq[i], " | frame = ", self.arr_rf[i], " & strand sign = ", self.iso_sj.strand 

    def retrieve_coding_region( self ):
        """
        NOTE: THIS HAS NOT BEEN TESTED YET!!
        Function: retrieve the nucleotides that have a reading frame associated with it (entails these nucleotides will translate into protein sequence)
        """
        #find first occurrence
        i_first_occur = next(i for i, v in enumerate( self.arr_rf ) if v != -1)

        #find last occurrence
        reverse_arr_rf = self.arr_rf[::-1] 
        i_last_occur_A = next(i for i, v in enumerate( reverse_arr_rf ) if v != -1)
        i_last_occur_B = len( self.arr_rf ) - i_last_occur_A - 1     #it's -1 because of iterating through index 0 in i_last_occur_A  

        return {'i_coding_start': i_first_occur, 'i_coding_end': i_last_occur_B}



    @staticmethod
    def is_obj_possible( transcript_sj, iso_sj ):
        """
        Args:
            -transcript_sj = array of SpliceJunction instances, sorted by start position of each SJ from least to greatest
            -iso_sj = IsoformSJ instance, will be used for isoform_id & list of canonical exons in iso_sj.hashExonList
        Function: see if it is possible to make an instance of this class (TranscribeTranscript)
            -Step: checks to see if the isoform from IsoformSJ instance "iso_sj" is present across SJs (e.g. for the aberrant SJ, perhaps it is assigned to a different isoform)
        """
        #check if isoform present in any of the SpliceJunction instances in transcript_sj
        multi_list = [x.hash_isoforms.keys() for x in transcript_sj]
        flat_list = [item for sublist in multi_list for item in sublist]
        if not iso_sj.isoform_id in flat_list:
            return False

        return True


    @staticmethod
    def compile_dnaseq_list( list_pos, strand_sign, path_genomeidx ):
        """
        Args:
            list_pos = array of genomic positions, where each element has to be in the format "chrom:start-end"
            strand_sign = integer where 1 if '+' or -1 if '-'
            path_genomeidx = string that is the path to the genome.fa file where the nucleotide sequences will be extracted (e.g. hg19.fa)
        Function: retrieves all DNA sequences from the list of positions in 'list_pos'. Returns an array where each element is a string of nucleotide sequences (based on the range for each element in array 'list_pos'). Note for plus strand the sequence from the plus strand, and the 
        """
        arr_elem_seq = []
        for each_elem in list_pos:
            rev_comp = False        #even if rev_comp is False, compile_dnaseq() will still get the complement nucleotide (not the reverse complement)
            adjust_0_base = False
            nucleotide_seq = TranscribeTranscript.compile_dnaseq( [each_elem], strand_sign, path_genomeidx, rev_comp, adjust_0_base )
            arr_elem_seq.append( nucleotide_seq )

        return arr_elem_seq



    def create_array_indices( self ):
        """
        Function: Creates arrays that will help match different types of information (e.g. genomic position to reading frame to nucleotide). Can use these arrays to find corresponding information between arrays. Creates array indices for the following:
            -genomic position = an array where each element is a numerical genomic position
            -nucleotide = an array where each element is nucleotide character
            -reading frame = an array where each element is a reading frame (0, 1, or 2)
            -exon # -> this is self.list_exons
            -amino acid sequence?
        """

        ##TEST::
        # for i, each_exon in enumerate(self.list_exons):
        #     print "TTV4 TEST ", i, " - exon ", each_exon.exonNum, " -> ", str( each_exon )

        arr_genome_pos = []     #records the genomic position
        arr_nuc_seq = []
        for each_exon in self.list_exons:
            exon_start = each_exon.exonPos.location.start + 1                   #added +1 because 0-based position
            exon_end = each_exon.exonPos.location.end
            arr_genome_pos+= [i for i in range( exon_start, exon_end + 1 )]     #added +1 because range excludes the last position

            #retrieve the nucleotides assigned to this exon
            str_pos = each_exon.chrom + ':' + str( exon_start ) + '-' + str( exon_end )     #need to add +1 else it doesn't retrieve the last nucleotide
            nucleotide_seq = TranscribeTranscript.compile_dnaseq( [str_pos], self.iso_sj.strand, self.path_genomeidx, False, False )
            arr_nuc_seq+= list( nucleotide_seq )


        #calculate reading frame for each position
        arr_rf = [-1 for i in range(0, len(arr_genome_pos) )]
        #retrieve the start and end cds positions
        try:
            start_i = arr_genome_pos.index( self.iso_sj.coding_boundary[0] + 1 )        #added +1 because 0-based position
        except (ValueError, TypeError) as e:
            start_i = 0
        #try to retrieve the end position
        try:
            end_i = arr_genome_pos.index( self.iso_sj.coding_boundary[1] )
        except (ValueError, TypeError) as e:
            end_i = len( arr_genome_pos ) - 1
        if self.iso_sj.strand < 0:
            for i in range( end_i, start_i-1, -1 ):
                arr_rf[i] = (end_i - i) % 3
        else:
            for i in range( start_i, end_i + 1 ):
                arr_rf[i] = (i - start_i) % 3


        ##TEST:: print "TT_V4 Check Lens = ", len( arr_genome_pos), " = ", len( arr_nuc_seq ), " = ", len( arr_rf )
        return [arr_genome_pos, arr_nuc_seq, arr_rf]

    def retrieve_nuc_seq( self, i_start, i_end, str_form = False ):
        """
        Args:
            i_start, i_end = integers that are indices for self.arr_nuc_seq, where i_start < i_end
        Function: retrieves the nucleotide sequence from self.arr_nuc_seq based on i_start & i_end. Not similar to function TranscribeTrascript.sequence_dna() since this uses self.arr_nuc_seq - returns array of nucleotide sequences
        NOTE: if I want to include the last index "i_end", need to add +1, therefore I did it in the code in this function
        """
        try:
            nuc_seq = self.arr_nuc_seq[ i_start: i_end + 1 ]        #add +1 to include the last index
            if self.iso_sj.strand < 0:
                nuc_seq = nuc_seq[::-1]     #if on minus strand
        except IndexError:
            print "Error in retrieve_nuc_seq: out of index range"
            return None

        return nuc_seq if not str_form else ''.join( nuc_seq )


    def find_codon_beginning_nearest( self, genome_pos, bool_before = True ):
        """
        Args:
            -genome_pos = integer position that will be used to find the reading frame for this given position
            -bool_before = boolean where if
                -True = retrieves nucleotides & AA before position 'given_pos' (by 'before' I mean closer to the 5' end of the gene - this takes strand sign into consideration)
                -False = retrieves nucleotides & AA after position 'given_pos' (by 'after' I mean closer to the 3' end of the gene - this takes strand sign into consideration)
        Function: finds the nucleotide position with reading frame = 0 that comes before genome_pos (or nearest position to genome_pos if genome_pos not in self.arr_genome_pos) -> means lower nucleotide position for + strand genes and higher nucleotide position for - strand genes 
            -when looking for nearest position, looks for position based on strand sign
                -nearest position 
        """
        ##TEST::
        # for i, a in enumerate( self.arr_genome_pos ):
        #     print "TTV4 FCBP i = ", i, " & a = ", a, " & genome_pos = ", genome_pos
        # for i2, each_sj in enumerate( self.transcript_sj ):
        #     print "TTV4. SJ i2 = ", i2, " & each_sj = ", str( each_sj ), " & gene = ", self.iso_sj.gene_sym
        # for i3, each_exon in enumerate( self.list_exons ):
        #     print "TTV4. Exon = ", i3, " & each_exon = ", str( each_exon )

        #if genome_pos is not in self.arr_genome_pos, find the nearest genomic position to genome_pos
        if not genome_pos in self.arr_genome_pos:
            if bool_before:     #look for position before 'genome_pos' (by 'before', I mean closer to the 5' end - this depends on the strand sign)
                nearest_pos = min( [x for x in self.arr_genome_pos if x > genome_pos]) if self.iso_sj.strand < 0 else max( [x for x in self.arr_genome_pos if x < genome_pos])
            else:       #look for position after 'genome_pos' (by 'after', I mean closer to the 3' end - this depends on the strand sign)
                nearest_pos = min( [x for x in self.arr_genome_pos if x < genome_pos]) if self.iso_sj.strand < 0 else max( [x for x in self.arr_genome_pos if x > genome_pos])
            #record new genome_pos
            genome_pos = nearest_pos
        
        #retrieve the reading frame based on teh position
        i_pos = self.arr_genome_pos.index( genome_pos )
        get_rf = self.arr_rf[ i_pos ]
            

        ##TEST:: 
        # print "TTV4 - FCB: genome_pos = ", genome_pos, " && self.arr_genome_pos[i_pos] = ", self.arr_genome_pos[i_pos], " && i_pos = ", i_pos, " && get_rf = ", get_rf

        if get_rf < 0:      #this means there is no reading frame associated with genome_pos, therefore no 0 reading frame near genome_pos
            return None
        elif get_rf == 0:
            return genome_pos
        else:
            if self.iso_sj.strand < 0:       #search for reading frame 0 in lower position
                ##TEST::
                # print "TTV4 - FCB (-1): genome_pos = ", self.arr_genome_pos[ i_pos + get_rf ], " && get_rf = ", self.arr_rf[ i_pos + get_rf ]

                return self.arr_genome_pos[ i_pos + get_rf ] if bool_before else self.arr_genome_pos[ i_pos - get_rf ]
            else:                   #search for reading frame 0 in higher position
                ##TEST:: print "TTV4 - FCB (+1): genome_pos = ", self.arr_genome_pos[ i_pos - get_rf ], " && get_rf = ", self.arr_rf[ i_pos - get_rf ]

                return self.arr_genome_pos[ i_pos - get_rf ] if bool_before else self.arr_genome_pos[ i_pos + get_rf ]

    def find_codon_rf_nearest( self, genome_pos, rf_oi = 2, bool_before = True ):
        """
        Args:
            -genome_pos = integer position that will be used to find the reading frame for this given position
            -rf_oi = the reading frame of interest. It should be either 0, 1, or 2, but since find_codon_beginning_prev() finds position for reading frame 0, should just use either 1 or 2
            -bool_before = boolean where if
                -True = retrieves nucleotides & AA before position 'given_pos' (by 'before' I mean closer to the 5' end of the gene - this takes strand sign into consideration)
                -False = retrieves nucleotides & AA after position 'given_pos' (by 'after' I mean closer to the 3' end of the gene - this takes strand sign into consideration)
        Function: similar to TranscribeTranscript.find_codon_beginning_prev(), but instead of just looking for the start of the codon (i.e. where reading frame of nucleotide is 0), it will look for the previous genomic position with the reading frame 'rf_oi' (for + strand = position will be <= genome_pos, for - strand = position will be >= genome_pos)
        """
        #determine if the current position has the reading frame    
        if genome_pos in self.arr_genome_pos:
            i_pos = self.arr_genome_pos.index( genome_pos )
            get_rf = self.arr_rf[ i_pos ]
            if get_rf == rf_oi:
                return genome_pos


        #calculate the distance to the next reading frame 
        dist_rf_oi_prev = {0: 0, 1: 2, 2: 1}     #the distance to reach the previous reading frame, where key = the desired reading frame (rf_oi) and the value is the distance to reach that previous reading frame
        dist_rf_oi_next = {0: 0, 1: 1, 2: 2}     #use if bool_before = False -> the distance to reach the next reading frame, where key = the desired reading frame (rf_oi) and the value is the distance to reach that previous reading frame

        #find position of nearest nucleotide with reading frame = 0 (start of codon)
        genome_pos_rf_zero = self.find_codon_beginning_nearest( genome_pos, bool_before )
        if not genome_pos_rf_zero:
            return None

        i_genome_pos_rf_zero = self.arr_genome_pos.index( genome_pos_rf_zero )

        #calculate the position that contains that reading frame (for previous exon - use dist_rf_oi_prev; for next exon - use dist_rf_oi_next)
        if self.iso_sj.strand < 0:      #find the previous or next codon for minus strand
            i_new_rf = i_genome_pos_rf_zero + dist_rf_oi_prev[rf_oi] if bool_before else i_genome_pos_rf_zero - dist_rf_oi_next[rf_oi] 
        else:           #find the previous or next codon for plus strand
            i_new_rf = i_genome_pos_rf_zero - dist_rf_oi_prev[rf_oi] if bool_before else i_genome_pos_rf_zero + dist_rf_oi_next[rf_oi] 

        if i_new_rf < 0 or i_new_rf > ( len(self.arr_genome_pos) - 1 ):
            return None
        else:
            return self.arr_genome_pos[ i_new_rf ]

    ##MAY DELETE THIS FUNCTION AS I HAVE def find_codon_beginning_nearest()
    # def find_codon_beginning_prev( self, genome_pos ):
    #     """
    #     Args:
    #         -genome_pos = integer position that will be used to find the reading frame for this given position
    #     Function: finds the nucleotide position with reading frame = 0 that comes before genome_pos (or nearest position to genome_pos if genome_pos not in self.arr_genome_pos) -> means lower nucleotide position for + strand genes and higher nucleotide position for - strand genes 
    #         -when looking for nearest position, looks for position based on strand sign
    #             -nearest position 
    #     """
    #     ##TEST::
    #     # for i, a in enumerate( self.arr_genome_pos ):
    #     #     print "TTV4 FCBP i = ", i, " & a = ", a, " & genome_pos = ", genome_pos
    #     # for i2, each_sj in enumerate( self.transcript_sj ):
    #     #     print "TTV4. SJ i2 = ", i2, " & each_sj = ", str( each_sj ), " & gene = ", self.iso_sj.gene_sym
    #     # for i3, each_exon in enumerate( self.list_exons ):
    #     #     print "TTV4. Exon = ", i3, " & each_exon = ", str( each_exon )

    #     #if genome_pos is not in self.arr_genome_pos, find the nearest genomic position to genome_pos
    #     if not genome_pos in self.arr_genome_pos:
    #         nearest_pos = min( [x for x in self.arr_genome_pos if x > genome_pos]) if self.iso_sj.strand < 0 else max( [x for x in self.arr_genome_pos if x < genome_pos])
    #         genome_pos = nearest_pos
        
    #     i_pos = self.arr_genome_pos.index( genome_pos )
    #     get_rf = self.arr_rf[ i_pos ]
            

    #     ##TEST:: 
    #     # print "TTV4 - FCB: genome_pos = ", genome_pos, " && self.arr_genome_pos[i_pos] = ", self.arr_genome_pos[i_pos], " && i_pos = ", i_pos, " && get_rf = ", get_rf

    #     if get_rf < 0:      #this means there is no reading frame associated with genome_pos, therefore no 0 reading frame near genome_pos
    #         return None
    #     elif get_rf == 0:
    #         return genome_pos
    #     else:
    #         if self.iso_sj.strand < 0:       #search for reading frame 0 in lower position
    #             ##TEST:: print "TTV4 - FCB (-1): genome_pos = ", self.arr_genome_pos[ i_pos + get_rf ], " && get_rf = ", self.arr_rf[ i_pos + get_rf ]

    #             return self.arr_genome_pos[ i_pos + get_rf ]
    #         else:                   #search for reading frame 0 in higher position
    #             ##TEST:: print "TTV4 - FCB (+1): genome_pos = ", self.arr_genome_pos[ i_pos - get_rf ], " && get_rf = ", self.arr_rf[ i_pos - get_rf ]

    #             return self.arr_genome_pos[ i_pos - get_rf ]

    ##MAY DELETE THIS FUNCTION AS I HAVE def find_codon_rf_nearest()
    # def find_codon_rf_prev( self, genome_pos, rf_oi = 2 ):
    #     """
    #     Args:
    #         -genome_pos = integer position that will be used to find the reading frame for this given position
    #         -rf_oi = the reading frame of interest. It should be either 0, 1, or 2, but since find_codon_beginning_prev() finds position for reading frame 0, should just use either 1 or 2
    #     Function: similar to TranscribeTranscript.find_codon_beginning_prev(), but instead of just looking for the start of the codon (i.e. where reading frame of nucleotide is 0), it will look for the previous genomic position with the reading frame 'rf_oi' (for + strand = position will be <= genome_pos, for - strand = position will be >= genome_pos)
    #     """
    #     #determine if the current position has the reading frame    
    #     if genome_pos in self.arr_genome_pos:
    #         i_pos = self.arr_genome_pos.index( genome_pos )
    #         get_rf = self.arr_rf[ i_pos ]
    #         if get_rf == rf_oi:
    #             return genome_pos


    #     #calculate the distance to the next reading frame 
    #     dist_rf_oi = {0: 0, 1: 2, 2: 1}     #the distance to reach the previous reading frame, where key = the desired reading frame (rf_oi) and the value is the distance to reach that previous reading frame

    #     #find position of nearest nucleotide with reading frame = 0 (start of codon)
    #     genome_pos_rf_zero = self.find_codon_beginning_prev( genome_pos )
    #     if not genome_pos_rf_zero:
    #         return None

    #     i_genome_pos_rf_zero = self.arr_genome_pos.index( genome_pos_rf_zero )

    #     #calculate the position that contains that reading frame IN THE PREVIOUS CODON
    #     i_new_rf = i_genome_pos_rf_zero + dist_rf_oi[rf_oi] if self.iso_sj.strand < 0 else i_genome_pos_rf_zero - dist_rf_oi[rf_oi] 
    #     if i_new_rf < 0 or i_new_rf > ( len(self.arr_genome_pos) - 1 ):
    #         return None
    #     else:
    #         return self.arr_genome_pos[ i_new_rf ]


    def find_containing_codon( self, given_pos ):
        """
        PROBLEM: This function is nearly identical to def find_aa_pos_before() (in this script), need to find a way to consolidate these functions
        Args:
            -given_pos = integer that is the genomic position
            -num_aa = integer that is the number of amino acids that needs to be retrieved 
            -bool_before = boolean where if
                -True = retrieves nucleotides & AA before position 'given_pos' (this takes strand sign into consideration)
                -False = retrieves nucleotides & AA after position 'given_pos' (this takes strand sign into consideration)
        Function: finds x nucleotide positions before that correspond to the number of amino acids before the position.
            -This is useful for identifying the amino acid sequence before a mutation or a position where the amino acid sequence changes
        """
        #I do not want to consider the amino acid at position "given_pos" as this usually is the position of mutation/5' end of aberrant splicing
        # given_pos = given_pos + 3 if self.iso_sj.strand < 0 else given_pos - 3
        try:
            i_given_pos = self.arr_genome_pos.index( given_pos )
        except (ValueError, TypeError) as e:
            print "Error in find_containing_codon(): could not find given position ", given_pos
            return {}

    
        #retrieve the start & end position for the amino acid sequence, and sort the positions based on lower & higher positions
        bool_before = True      #set to True because I want to find the start of the current codon
        pos_start = self.find_codon_beginning_nearest( given_pos, bool_before )
        i_pos_start = self.arr_genome_pos.index( pos_start )
        i_pos_end = i_pos_start - 2 if self.iso_sj.strand < 0 else i_pos_start + 2
        try:
            pos_end = self.arr_genome_pos[ i_pos_end ]
        except IndexError:
            print "Error in find_containing_codon(): cannot retrieve codon because end position at index does not exist = ", i_pos_end
            pos_end = None

        if not pos_start or not pos_end:
            return {}

        ##TEST::
        # print "\n<<<<<<<<<<<"
        # print "TT_V4.find_containing_codon() PART 1 - given_pos = ", given_pos, " & pos_start = ", pos_start, " & pos_end = ", pos_end
        # print "TT_V4.find_containing_codon() PART 2 - RF for given_pos = ", self.arr_rf[ i_given_pos ], " & RF for pos_start ", self.arr_rf[ i_pos_start ], " & RF for end pos = ", self.arr_rf[ i_pos_end ]
        # print ">>>>>>>>>>>\n"


        if pos_start <= pos_end:        #this usually applies to + strand genes, where start position < end position
            genome_pos_lower = pos_start
            genome_pos_higher = pos_end
            # i_genome_pos_lower = self.arr_genome_pos.index( genome_pos_lower )
            # i_genome_pos_higher = self.arr_genome_pos.index( genome_pos_higher )
            str_genome_range = self.iso_sj.chrom + ':' + str( pos_start ) + '-' + str( pos_end )
        else:
            genome_pos_lower = pos_end
            genome_pos_higher = pos_start
            # i_genome_pos_lower = self.arr_genome_pos.index( genome_pos_lower )
            # i_genome_pos_higher = self.arr_genome_pos.index( genome_pos_higher )
            str_genome_range = self.iso_sj.chrom + ':' + str( pos_end ) + '-' + str( pos_start )
        i_genome_pos_lower = self.arr_genome_pos.index( genome_pos_lower )
        i_genome_pos_higher = self.arr_genome_pos.index( genome_pos_higher )

        return {'genome_start': genome_pos_lower, 'genome_end': genome_pos_higher, 'str_genome_range': str_genome_range, 'i_genome_start': i_genome_pos_lower, 'i_genome_end': i_genome_pos_higher }

    
    def retrieve_containing_codon( self, given_pos, strand_sign = None ):
        """
        returns the nucleotide sequence for the containing codon

        Args:
            -given_pos = integer that is the genomic position
            -strand_sign = integer that is the strand sign of interest (+1 = plus strand, -1 = negative strand). If None, then will take on the strand sign of self.iso_sj.strand
                -if iso_sj.strand & strand_sign do not match, will need to find the reverse_complement of the codon sequence
        Output: returns the nucleotide sequence for the codon of interest
        """
        if not strand_sign:
            strand_sign = self.iso_sj.strand

        hash_codon = self.find_containing_codon( given_pos )
        #if 'given_pos' is not found in def find_containing_codon(), then this means hash_codon_pos is empty
        if not hash_codon:
            return None

        if self.iso_sj.strand < 0:
            str_nuc_codon = ''.join( self.arr_nuc_seq[ hash_codon['i_genome_start']:hash_codon['i_genome_end'] + 1 ][::-1] )
        else:
            str_nuc_codon = ''.join( self.arr_nuc_seq[ hash_codon['i_genome_start']:hash_codon['i_genome_end'] + 1 ] )

        #if the isoform strand sign does not match the parameter strand_sign, then need to take the reverse complement
        if strand_sign != self.iso_sj.strand:
            str_nuc_codon = str( Seq(str_nuc_codon).reverse_complement().upper() )

        return str_nuc_codon


    def find_aa_pos_surrounding( self, given_pos, num_aa, bool_before = True ):
        """
        PROBLEM: AS OF NOW IT DOESN'T GET THE EXACT NUMBER OF AMINO ACIDS "num_aa", it is either the nucleotide sequence is equal to or less than num_aa (usually 1 less)
        PROBLEM: This function is nearly identical to def find_containing_codon() (in this script), need to find a way to consolidate these functions
        Args:
            -given_pos = integer that is the genomic position
            -num_aa = integer that is the number of amino acids that needs to be retrieved 
            -bool_before = boolean where if
                -True = retrieves nucleotides & AA before position 'given_pos' (by 'before' I mean closer to the 5' end of the gene - this takes strand sign into consideration)
                -False = retrieves nucleotides & AA after position 'given_pos' (by 'after' I mean closer to the 3' end of the gene - this takes strand sign into consideration)
        Function: finds x nucleotide positions before that correspond to the number of amino acids before the position.
            -This is useful for identifying the amino acid sequence before a mutation or a position where the amino acid sequence changes
        """
        #I do not want to consider the amino acid at position "given_pos" as this usually is the position of mutation/5' end of aberrant splicing
        # given_pos = given_pos + 3 if self.iso_sj.strand < 0 else given_pos - 3
        try:
            i_given_pos = self.arr_genome_pos.index( given_pos )
        except (ValueError, TypeError) as e:
            print "Error in find_aa_pos_before: could not find given position ", given_pos
            return {}

        #calculate the number of nucleotides needed amino acid bases 
        num_nuc_before = num_aa * 3       #do -1 because will include the amino acid of codon in 'given_pos'
        # if bool_before:
        #     i_aa_start_pos = i_given_pos + num_nuc_before if self.iso_sj.strand < 0 else i_given_pos - num_nuc_before 
        # else:
        #     i_aa_start_pos = i_given_pos - num_nuc_before if self.iso_sj.strand < 0 else i_given_pos + num_nuc_before

        #NOTE: I use '- 3' because I want to step out of the position of the current codon
        if self.iso_sj.strand < 0:
            if bool_before:     #retrieve position closer to 5' end
                # bool_before_2_start = True
                # bool_before_2_end = False
                i_aa_start_pos = i_given_pos + num_nuc_before
                i_aa_end_pos = i_given_pos + 3      #need to '+ 3' because I want to step out of the position of the current codon
            else:
                # bool_before_2_start = True
                # bool_before_2_end = False
                i_aa_start_pos = i_given_pos - 3
                i_aa_end_pos = i_given_pos - num_nuc_before
        else:           #look at plus strands
            if bool_before:
                # bool_before_2_start = True
                # bool_before_2_end = False
                i_aa_start_pos = i_given_pos - num_nuc_before
                i_aa_end_pos = i_given_pos - 3      #need to '- 3' because I want to step out of the position of the current codon
            else:
                # bool_before_2_start = True
                # bool_before_2_end = False
                i_aa_start_pos = i_given_pos + 3
                i_aa_end_pos = i_given_pos + num_nuc_before
        
        try:
            aa_start_pos = self.arr_genome_pos[ i_aa_start_pos ]
            aa_end_pos = self.arr_genome_pos[ i_aa_end_pos ]
        except IndexError:
            print "Error in TranscribeTranscript.find_aa_pos_surrounding: genomic positions for AA are out of range -> index of start = ", i_aa_start_pos, " & index of end = ", i_aa_end_pos
            return {}

        #retrieve the start & end position for the amino acid sequence, and sort the positions based on lower & higher positions

        ##TEST::
        print "TT_V4.find_aa_pos_surrounding(): aa_start_pos = ", aa_start_pos, " & aa_end_pos = ", aa_end_pos

        bool_before_2_start = True          #for aa_start_pos, retrieve position closer to 5' end
        bool_before_2_end = False           #for aa_end_pos, retrieve position closer to 3' end
        pos_start = self.find_codon_beginning_nearest( aa_start_pos, bool_before_2_start )
        pos_end = self.find_codon_rf_nearest( aa_end_pos, 2, bool_before_2_end )

        if not pos_start or not pos_end:
            return {}

        ##TEST::
        # print "\n<<<<<<<<<<<"
        # print "TT_V4.find_aa_pos_before() PART 1 - given_pos = ", given_pos, " & aa_start_pos = ", aa_start_pos, " &  aa_end_pos = ", aa_end_pos, " & pos_start = ", pos_start, " & pos_end = ", pos_end
        # i_pos_end = self.arr_genome_pos.index( pos_end )
        # print "TT_V4.find_aa_pos_before() PART 2 - RF for given_pos = ", self.arr_rf[ i_given_pos ], " & RF for aa_start_pos ", self.arr_rf[ i_aa_start_pos ], " & RF for end pos = ", self.arr_rf[ i_pos_end ]
        # print ">>>>>>>>>>>\n"


        if pos_start <= pos_end:        #this usually applies to + strand genes, where start position < end position
            genome_pos_lower = pos_start
            genome_pos_higher = pos_end
            str_genome_range = self.iso_sj.chrom + ':' + str( pos_start ) + '-' + str( pos_end )
        else:
            genome_pos_lower = pos_end
            genome_pos_higher = pos_start
            str_genome_range = self.iso_sj.chrom + ':' + str( pos_end ) + '-' + str( pos_start )
        i_genome_pos_lower = self.arr_genome_pos.index( genome_pos_lower )
        i_genome_pos_higher = self.arr_genome_pos.index( genome_pos_higher )

        return {'genome_start': genome_pos_lower, 'genome_end': genome_pos_higher, 'str_genome_range': str_genome_range, 'i_genome_start': i_genome_pos_lower, 'i_genome_end': i_genome_pos_higher }


    ##MAY DELETE THIS FUNCTION BECAUSE OF def find_aa_pos_surrounding()
    # def find_aa_pos_before( self, given_pos, num_aa, bool_before = True ):
    #     """
    #     PROBLEM: AS OF NOW IT DOESN'T GET THE EXACT NUMBER OF AMINO ACIDS "num_aa", it is either the nucleotide sequence is equal to or less than num_aa (usually 1 less)
    #     PROBLEM: This function is nearly identical to def find_containing_codon() (in this script), need to find a way to consolidate these functions
    #     Args:
    #         -given_pos = integer that is the genomic position
    #         -num_aa = integer that is the number of amino acids that needs to be retrieved 
    #         -bool_before = boolean where if
    #             -True = retrieves nucleotides & AA before position 'given_pos' (by 'before' I mean closer to the 5' end of the gene - this takes strand sign into consideration)
    #             -False = retrieves nucleotides & AA after position 'given_pos' (by 'after' I mean closer to the 3' end of the gene - this takes strand sign into consideration)
    #     Function: finds x nucleotide positions before that correspond to the number of amino acids before the position.
    #         -This is useful for identifying the amino acid sequence before a mutation or a position where the amino acid sequence changes
    #     """
    #     #I do not want to consider the amino acid at position "given_pos" as this usually is the position of mutation/5' end of aberrant splicing
    #     # given_pos = given_pos + 3 if self.iso_sj.strand < 0 else given_pos - 3
    #     try:
    #         i_given_pos = self.arr_genome_pos.index( given_pos )
    #     except (ValueError, TypeError) as e:
    #         print "Error in find_aa_pos_before: could not find given position ", given_pos
    #         return {}

    #     #calculate the number of nucleotides needed amino acid bases 
    #     num_nuc_before = num_aa * 3       #do -1 because will include the amino acid of codon in 'given_pos'
    #     if bool_before:
    #         i_aa_start_pos = i_given_pos + num_nuc_before if self.iso_sj.strand < 0 else i_given_pos - num_nuc_before
    #     else:
    #         i_aa_start_pos = i_given_pos - num_nuc_before if self.iso_sj.strand < 0 else i_given_pos + num_nuc_before
        
    #     try:
    #         aa_start_pos = self.arr_genome_pos[ i_aa_start_pos ]
    #     except IndexError:
    #         return {}

    #     #retrieve the start & end position for the amino acid sequence, and sort the positions based on lower & higher positions
    #     pos_start = self.find_codon_beginning_prev( aa_start_pos )
    #     pos_end = self.find_codon_rf_prev( given_pos, 2 )

    #     if not pos_start or not pos_end:
    #         return {}

    #     ##TEST::
    #     print "\n<<<<<<<<<<<"
    #     print "TT_V4.find_aa_pos_before() PART 1 - given_pos = ", given_pos, " & aa_start_pos = ", aa_start_pos, " & pos_start = ", pos_start, " & pos_end = ", pos_end
    #     i_pos_end = self.arr_genome_pos.index( pos_end )
    #     print "TT_V4.find_aa_pos_before() PART 2 - RF for given_pos = ", self.arr_rf[ i_given_pos ], " & RF for aa_start_pos ", self.arr_rf[ i_aa_start_pos ], " & RF for end pos = ", self.arr_rf[ i_pos_end ]
    #     print ">>>>>>>>>>>\n"


    #     if pos_start <= pos_end:        #this usually applies to + strand genes, where start position < end position
    #         genome_pos_lower = pos_start
    #         genome_pos_higher = pos_end
    #         # i_genome_pos_lower = self.arr_genome_pos.index( genome_pos_lower )
    #         # i_genome_pos_higher = self.arr_genome_pos.index( genome_pos_higher )
    #         str_genome_range = self.iso_sj.chrom + ':' + str( pos_start ) + '-' + str( pos_end )
    #     else:
    #         genome_pos_lower = pos_end
    #         genome_pos_higher = pos_start
    #         # i_genome_pos_lower = self.arr_genome_pos.index( genome_pos_lower )
    #         # i_genome_pos_higher = self.arr_genome_pos.index( genome_pos_higher )
    #         str_genome_range = self.iso_sj.chrom + ':' + str( pos_end ) + '-' + str( pos_start )
    #     i_genome_pos_lower = self.arr_genome_pos.index( genome_pos_lower )
    #     i_genome_pos_higher = self.arr_genome_pos.index( genome_pos_higher )

    #     return {'genome_start': genome_pos_lower, 'genome_end': genome_pos_higher, 'str_genome_range': str_genome_range, 'i_genome_start': i_genome_pos_lower, 'i_genome_end': i_genome_pos_higher }


    def transcript_notes_record_sj( self ):
        """
        Function: records position of SJ in hash 'self.hash_transcript_notes'
        """
        for each_sj in self.transcript_sj:
            if not each_sj.start in self.hash_transcript_notes:
                self.hash_transcript_notes[ each_sj.start ] = []     #will record a set of annotations
            self.hash_transcript_notes[ each_sj.start ].append( 'sj_low' )      #sj_low means it is located at a lower genomic position

            if not each_sj.end in self.hash_transcript_notes:
                self.hash_transcript_notes[ each_sj.end ] = []     #will record a set of annotations
            self.hash_transcript_notes[ each_sj.end ].append( 'sj_high' )       #sj_high means it is located at a higher genomic position

    ##NOTE: This maybe obsolete now because I do not use an array of positions, but instead an array of SpliceJunction instances that constitute a transcript (i.e. replaced self.arr_pos with self.transcript_sj)
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

            #record the previous and next exon position ranges based on 
            record_exons_transcript += hash_elems['prev_elem'].str_genomic_pos( True )
            record_exons_transcript += delim_exons
            record_exons_transcript += hash_elems['next_elem'].str_genomic_pos( True )

        return record_exons_transcript


    #MAY DELETE BECAUSE THIS DOESN'T WORK CORRECTLY
    # @staticmethod
    # def retrieve_exon_from_sj_v2( list_sj, isoform_id ):
    #     """
    #     Args:
    #         list_sj = array of SpliceJunction instances that constitutes an mRNA transcript. This assumes the order the SJs are from least to greatest
    #         isoform_id = string that is isoform ID
    #     Function: similar to def retrieve_exon_from_sj() where retrieves the exon positions from the array of SpliceJunction instances, but returns a list of Exon instances
    #     """
    #     #record all canonical exons based on canonical SJs
    #     hash_exons = {}     #k = exon number, v = Exon instance
    #     canon_sj = [x for x in list_sj if x.canon]
    #     for each_sj in canon_sj:
    #         tuple_exons = each_sj.spliced_elems( isoform_id, True )
    #         hash_exons[ tuple_exons[0].exonNum ] = tuple_exons[0]
    #         hash_exons[ tuple_exons[1].exonNum ] = tuple_exons[1]

    #     #replace canonical exons with modified exons based on aberrant splicing
    #     # aberrant_sj = [x for x in list_sj if not x.canon]
    #     aberrant_sj = [x for x in list_sj if not x in canon_sj]     #prefer this way so I get all possible exons 
    #     for each_sj in aberrant_sj:
    #         tuple_exons = each_sj.spliced_elems_position_range_v2( isoform_id, False, True )

    #         ##TEST::
    #         print "TTV4 - refsjv2: start = ", tuple_exons['exon_start'].exonNum, " -> ", tuple_exons['exon_start']
    #         print "TTV4 - refsjv2: end = ", tuple_exons['exon_end'].exonNum, " -> ", tuple_exons['exon_end']

    #         hash_exons[ tuple_exons['exon_start'].exonNum ] = tuple_exons['exon_start']
    #         hash_exons[ tuple_exons['exon_end'].exonNum ] = tuple_exons['exon_end']

    #     list_exons = sorted( [v for k,v in hash_exons.iteritems() ], key = lambda x: x.exonPos.location.start )
    #     return list_exons

    @staticmethod
    def create_exons_from_sj( transcript_sj, iso_sj ):
        """
        Args:
            -transcript_sj = array of SpliceJunction instances, sorted by start position of each SJ from least to greatest
            -iso_sj = IsoformSJ instance, will be used for isoform_id & list of canonical exons in iso_sj.hashExonList
        Function: create exons from 
        """
        dist_intergenic = 100         #use this if the SJ is before the start of the gene or after the end of the gene. This will tend to fall into the intergenic region 
        transcript_sj = sorted( transcript_sj, key = lambda s: s.start, reverse = False )

        #get start & end position
        obj_iso = transcript_sj[0].hash_isoforms[iso_sj.isoform_id]        #get instance of Isoform object
        if iso_sj.strand < 0:       #if minus-strand gene
            start_exon = obj_iso.get_exon_num( obj_iso.last_exon_num )
            end_exon = obj_iso.get_exon_num( 1 )            #use 1 if 1-based exon numbering (first exon has exon number 1), and 0-based exon numbering (first exon has exon number 0)
        else:           #else if plus-strand gene
            start_exon = obj_iso.get_exon_num( 1 )          #use 1 if 1-based exon numbering (first exon has exon number 1), and 0-based exon numbering (first exon has exon number 0)
            end_exon = obj_iso.get_exon_num( obj_iso.last_exon_num )

        ##TEST:: show start & end exons
        # print "TTV4 CEFSJ 2 - start_exon = ", start_exon, " & end_exon = ", end_exon

        list_exons = []     #list of exons ordered by genomic position (least to greatest)
        for i in range( 0, len(transcript_sj) + 1 ):

            if i == 0:      #if the first SJ in the transcript
                pos_start = start_exon.exonPos.location.start 
                pos_end = transcript_sj[i].start
                exon_num = start_exon.exonNum

                #if the starting position of the SJ is before the start of the exon
                if pos_end <= pos_start:
                    pos_start = iso_sj.boundary[0] if pos_end > iso_sj.boundary[0] else (pos_end - dist_intergenic)
            elif i == ( len(transcript_sj) ):       #if the last SJ in the transcript
                pos_start = transcript_sj[i - 1].end
                pos_end = end_exon.exonPos.location.end
                exon_num = end_exon.exonNum

                #if the ending position of the SJ is after the end of the last exon
                if pos_start >= pos_end:
                    pos_end = iso_sj.boundary[1] if pos_start < iso_sj.boundary[1] else (pos_start + dist_intergenic)
            else:
                pos_start = transcript_sj[i - 1].end 
                pos_end = transcript_sj[i].start
                exon_num = None

            #if the start & end position are the same, then there is no exon as the ends of the splice junction have overlapped
            if pos_start == pos_end:
                continue

            exon_info = {
            "pos_start": int( pos_start ),
            "pos_end": int( pos_end ),
            "strand_sign": iso_sj.strand,
            "chrom": iso_sj.chrom,
            "isoform_id": iso_sj.isoform_id,
            "exon_frame": None,
            "exon_num": exon_num,
            "canonical": None,
            "feat_type": 'exon' }

            #record the Exon instance, and see if it is canonical or not
            obj_exon = Exon(exon_info)
            canon_exon_match = [v for k,v in iso_sj.hashExonList.iteritems() if v == obj_exon ]
            if canon_exon_match:
                list_exons.append( canon_exon_match[0] )
            else:
                list_exons.append( obj_exon )

            """
            Things to consider in the future:
                -calculate exon number
                -get reading frame
            """

        return list_exons

    @staticmethod
    def create_exons_no_sj( iso_sj ):
        """
        Args:
            -iso_sj = instance of IsoformSJ class. USE = retrieve exons from hashExonList
        Function: returns a list of exons for IsoformSJ's hashExonList. Use this when there are no splice junctions for gene (i.e. self.transcript_sj is empty)
        """
        return [v for k,v in iso_sj.hashExonList.iteritems()]

        

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
            -UCSC genome browser uses 0-based coordinates instead of 1-based coordinates, the starting position of each exon (lower nucleotide position of exon) is not referring to the first nucleotide of the exon but instead the last nucleotide of the intron. Therefore adding +1 gives it the next position 
        """
        hash_exon_pos = Isoform.split_genome_pos( str_pos )
        hash_exon_pos['start'] += 1

        return hash_exon_pos['chrom'] + ':' + str( hash_exon_pos['start'] ) + '-' + str( hash_exon_pos['end'] )


    @classmethod
    def compile_dnaseq( cls_obj, arr_pos, strand_sign, path_genomeidx, rev_comp = False, adjust_0_base = False ):
        """ 
        Args:
            -arr_pos = array of positions where each element has a position in the format 'chrom:pos_start-pos_end'
            -strand_sign = integer where 1 if '+' or -1 if '-'
            -path_genomeidx = string that is the path to the samtools faidx genome file
            -rev_comp = means "reverse_complement". This should only be implemented if strand_sign = '-'. Boolean where
                -True = returns the reverse complement of the nucleotide sequence if strand_sign is '-'. Use this if I want to retrieve the transcribed sequence from 5' to 3' direction
                -False = returns the complement of the nucleotide sequence. Use this if I want to keep retrieve the complement sequence and have the genomic positions preserved.
            -adjust_0_base = boolean where:
                -True = will add +1 to the starting position of each exon. Use this for 0-based position
                -False = will keep the original position of each exon
        Function:
            will compile the genomic sequences ( DNA ) from an array contain multiple positions 'arrPos', and returns a Bio.Seq.Seq object
        """
        seq_transcript = ""
        for i, pos in enumerate( arr_pos ):
            exon_pos = pos if not adjust_0_base else TranslateTranscript.adjust_position( pos )
            seq_transcript += cls_obj.sequence_dna( exon_pos, path_genomeidx )

        #create BioSeq sequence object
        if strand_sign < 0:     #if strand sign is negative, then get reverse complement
            # seq_transcript = seq_transcript[::-1]     #reverse the string of nucleotides (assuming the exons entered are from right to left aka higher to lower)
            return Seq( seq_transcript.upper(), IUPAC.unambiguous_dna ).reverse_complement() if rev_comp else Seq( seq_transcript.upper(), IUPAC.unambiguous_dna ).complement()      ##QUES: what is IUPAC.unambiguous_dna?? WHAT DOES IT DO?
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
            -transcript_sj = an array of SpliceJunction instances that constitute a transcript. The SpliceJunctions should be sorted from lowest numerical position to highest numerical position.
            -strand_sign = string that is either '+' or '-', denoting whether the gene is on the plus or minus side of the gene
            -hash_transcript_notes = hash that keeps track of events (e.g. aberrant splicing, mutations, insertions, deletions, etc.) for each position. k = integer that is genomic position, v = array of events that occur in this position.
        Function: retrieves the nucleotide sequences in a position and returns a hash the contains the genomic position, nucleotide, and the associated events (e.g. splicing, aberrant splicing, mutation, etc.) that occur at this position
        """
        hash_seq_notes = {}     #key = integer that is the genomic position, v = hash with following key:value pairs - nucleotide:character that is nucleotide, annot:array of defects associated
        #retrieve the nucleotide sequence from a given exon
        for each_exon in self.list_exons:
            exon_seq = TranscribeTranscript.sequence_dna( each_exon.str_genomic_pos(True), self.path_genomeidx )

            if self.iso_sj.strand < 0:     #if strand sign is negative, then get reverse complement
                ##METHOD 1 - retrieve complementary sequence
                exon_seq = Seq( exon_seq.upper(), IUPAC.unambiguous_dna ).complement()      ##QUES: what is IUPAC.unambiguous_dna?? WHAT DOES IT DO?
                ##METHOD 2 - retrieve complementary sequence:
                # exon_seq = Seq( exon_seq.upper(), IUPAC.unambiguous_dna ).reverse_complement()      ##QUES: what is IUPAC.unambiguous_dna?? WHAT DOES IT DO?
                # exon_seq = exon_seq[::-1]       #reverse the order of the nucleotides so the genomic position refers to the correct nucleotide sequence

            hash_seq_notes = TranscribeTranscript.demarcate_dna_seq( each_exon.str_genomic_pos(False), exon_seq, self.hash_transcript_notes, hash_seq_notes )

        return hash_seq_notes


    def compile_dna_seq_only( self ):
        """
        Function: returns only the nucleotide sequence record in 'hash_seq_notes'
        """
        return ''.join( [v['nuc'] for k,v in self.hash_seq_notes.iteritems()] )


    ##NOTE: May need to adjust this for strand sign (self.iso_sj.strand)
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

    """
    Functions: make single base mutations
    """
    def get_mutation_BOS( self, base_orig, base_mut, snv_strand ):        #get_mutation_BOS = get mutation based on strand
        """
        Args:
            base_orig = character that is the original nucleotide base
            base_mut = character that is the mutated nucleotide base
            snv_strand = integer that is the strand sign -> -1 for minus strand & +1 for plus strand
        Function: determines the correct nucleotide sequence based on strand sign
        """
        if self.iso_sj.strand != snv_strand:
            base_orig = str( Seq(base_orig, IUPAC.unambiguous_dna).complement() )
            base_mut = str( Seq(base_mut, IUPAC.unambiguous_dna).complement() )

        return [base_orig, base_mut]

    def make_mutation( self, genome_pos, base_change ):
        """
        Args:
            genome_pos = integer that is the genomic position of interest
            base_change = character that is the change to be made to the nucleotide sequence. Assumes that this is the correct nucleotide based on the strand sign, so make sure to use TranscribeTranscript.get_mutation_BOS() before using this function.
        Function: create a mutation, changing based at positon 'genome_pos' to character 'base_change'
        """
        try:
            i_genome_pos = self.arr_genome_pos.index( genome_pos )
            self.arr_nuc_seq[ i_genome_pos ] = base_change
            return True
        except ValueError:
            print "Making mutation - could not make mutation at position ", genome_pos
            return False


    def get_mutated_codon( self, base_orig, base_mut, snv_pos, snv_strand, bool_before = True ):
        """
        Args:
            -self = instance of TranscribeTranslate
            -base_orig = character that is the original nucleotide at position 'snv_pos' (assumes plus strand, so need to find nucleotide complement if 'snv_strand' is minus)
            -base_mut = character that is the mutated nucleotide at position 'snv_pos' (assumes plus strand, so need to find nucleotide complement if 'snv_strand' is minus)
            -snv_pos = string in the format 'chr:start-end' (e.g. chr9:4-25), this is the position of the SNV
            -bool_before = boolean where if
                -True = retrieves nucleotides & AA before position 'snv_pos' (by 'before' I mean closer to the 5' end of the gene - this takes strand sign into consideration)
                -False = retrieves nucleotides & AA after position 'snv_pos' (by 'after' I mean closer to the 3' end of the gene - this takes strand sign into consideration)
        Function: returns information about the original and mutated codon. Includes information from def find_aa_pos_before() (index position, genomic positions)
        """
        #get original codon
        hash_genome_pos = Isoform.split_genome_pos( snv_pos )
        pos_codon_start = self.find_codon_beginning_nearest( hash_genome_pos['start'], bool_before )
        codon_pos = self.find_containing_codon( pos_codon_start )
        if not codon_pos:       #this means no codon was found associated with position 'pos_codon_start'
            hash_null = {'codon_orig': '-', 'codon_mut': '-', 'aa_orig': '-', 'aa_mut': '-'}
            hash_null.update( {'genome_start': '-', 'genome_end': '-', 'str_genome_range': '-', 'i_genome_start': '-', 'i_genome_end': '-'} )     #this is a hash that would have been returned by def find_containing_codon()
            return hash_null

        #get original codon
        codon_orig = ''.join( self.arr_nuc_seq[ codon_pos['i_genome_start'] : codon_pos['i_genome_end'] + 1 ] )

        #get mutated codon
        [base_orig, base_mut] = self.get_mutation_BOS( base_orig, base_mut, snv_strand )
        stat_mutation = self.make_mutation( hash_genome_pos['start'], base_mut )
        codon_mut = ''.join( self.arr_nuc_seq[ codon_pos['i_genome_start'] : codon_pos['i_genome_end'] + 1 ] )

        #reverse codons if on minus strand
        if self.iso_sj.strand < 0:
            codon_orig = codon_orig[::-1] 
            codon_mut = codon_mut[::-1] 

        ##TEST::
        # print "TT_V4.get_mutated_codon(): base_orig = ", base_orig, " & base_mut = ", base_mut
        # print "TT_V4.get_mutated_codon(): codon_orig = ", codon_orig, " & codon_mut = ", codon_mut

        #get resulting amino acid
        aa_orig = Seq( codon_orig ).translate( to_stop = False )
        aa_mut = Seq( codon_mut ).translate( to_stop = False )

        #return information about codon
        hash_codon_info = {'codon_orig': codon_orig, 'codon_mut': codon_mut, 'aa_orig': aa_orig, 'aa_mut': aa_mut}
        hash_codon_info.update( codon_pos )     #hash structure = {'genome_start': genome_pos_lower, 'genome_end': genome_pos_higher, 'str_genome_range': str_genome_range, 'i_genome_start': i_genome_pos_lower, 'i_genome_end': i_genome_pos_higher }
        return hash_codon_info


    """
    Functions: create annotations to nucleotide & amino acid sequence based on genomic alterations & aberrant SJ present in annotations
    """
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

    def __init__( self, transcript_sj, iso_sj, path_genomeidx, hash_transcript_notes ):
        TranscribeTranscript.__init__( self, transcript_sj, iso_sj, path_genomeidx, hash_transcript_notes )

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

    # def get_nmd_sensitive_region( self, pos_aberrant ):
    #     """
    #     Args:
    #         -self = instance of class TranscribeTranscript/TranslateTranscript
    #         -pos_aberrant = integer that is the position that contains the aberrant position of interest. This will usually be the start of aberrant SJ (lower position for + genes & higher position for the - genes)
    #     Function: retrieve the region of the gene that, if it contains an early stop codon, will lead to degradation of the transcript. Returns the nucleotide sequence of this region
    #     """
    #     #retrieve the position of the initial nucleotide in the codon
    #     direction = 1 if self.iso_sj.strand < 0 else -1       #if - strand, then find 0-rf nucleotide at higher position, else if + strand, then find 0-rf nucleotide at lower position
    #     pos_aberrant_start = self.find_codon_beginning_prev( pos_aberrant )

    #     if not pos_aberrant_start:
    #         return None

    #     i_aberrant_start = self.arr_genome_pos.index( pos_aberrant_start )
    #     # ( pos_aberrant_start, i_aberrant_start ) = find_nearest_rf0( self, pos_aberrant )

    #     #get the position 55nt from the 3' end of the penultimate exon
    #     i_tp_pux = get_index_tp_pux( self )
    #     if not i_tp_pux:
    #         return None

    #     ##TEST:: get start of codon previous to SJ
    #     diff = 4
    #     print "GET_NMD: pos_aberrant_start = ", pos_aberrant_start, " & i_aberrant_start = ", i_aberrant_start, " & i_tp_pux = ", i_tp_pux
    #     print "GET_NMD: genome_pos_start = ", self.arr_genome_pos[i_aberrant_start - diff : i_aberrant_start], " & 55 bp upstream = ", self.arr_genome_pos[i_tp_pux : i_tp_pux + diff]
    #     print "GET_NMD: nucleotide start = ", self.arr_nuc_seq[i_aberrant_start - diff : i_aberrant_start ], " & nucleotide 55 bp upstream = ", self.arr_nuc_seq[i_tp_pux : i_tp_pux + diff]
    #     print "GET_NMD: start reading frame = ", self.arr_rf[i_aberrant_start - diff: i_aberrant_start ], " & nucleotide 55 bp upstream reading frame = ", self.arr_rf[i_tp_pux : i_tp_pux + diff]

    #     ##TEST::
    #     # print "NMD_SENS: strand = ", self.iso_sj.strand, " | i_aberrant_start = ", i_aberrant_start, " | aberr_genome_pos = ", self.arr_genome_pos[i_aberrant_start]
    #     # print "NMD_SENS: strand = ", self.iso_sj.strand, " | i_tp_pux = ", i_tp_pux, " | tp_genome_pos = ", self.arr_genome_pos[i_tp_pux]

    #     if self.iso_sj.strand < 0:
    #         nmd_sensitive_seq = ''.join( self.arr_nuc_seq[i_tp_pux : i_aberrant_start + 1] )      #need to add +1 to include last base
    #         nmd_sensitive_seq = nmd_sensitive_seq[::-1]
    #     else:
    #         nmd_sensitive_seq = ''.join( self.arr_nuc_seq[i_aberrant_start : i_tp_pux + 1] )

    #     ##TEST:: print "NMD_SENS: seq = ", nmd_sensitive_seq

    #     return nmd_sensitive_seq


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
        else:
            nmd_irrelevant_seq = ''.join( obj_tt.arr_nuc_seq[i_tp_start:] )     #retrieve from tp_start until the last nucleotide

        ##TEST:: print "NMD_IRREL: strand = ", obj_tt.iso_sj.strand, " | i_tp_pux = ", i_tp_pux, " | i_tp_start = ", i_tp_start, " seq = ", nmd_irrelevant_seq

        return nmd_irrelevant_seq
