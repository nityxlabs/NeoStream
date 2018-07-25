#/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from cruzdb import Genome

from IsoformSJ import IsoformSJ
from TranscribeTranslate_V5 import TranscribeTranscript, TranslateTranscript
from NeoepitopeMHC import NeoepitopeMHC
from SimpleNeoepitopeV2 import SimpleNeoepitopeAllV2, SimpleNeoepitopeIsoformV2

class SimpleNeoepitopeIsoformSJ():

    def __init__( self, db_type, isoform_id, obj_sj_aberr, path_genomeidx ):
        """
        NOTE: For this class, it is only considering/incorporating single aberrant splicing when reconstructing a transcript isoform
        Args:
            -iso_sj_aberr = IsoformSJ instance that contains the aberrant SJ of interest
            -path_genomeidx = string that is the path to the samtools-indexed genome
        """
        self.db_type = db_type
        self.isoform_id = isoform_id
        self.obj_sj_aberr = obj_sj_aberr
        self.is_coding = self.obj_sj_aberr.hash_isoforms[isoform_id].is_coding         #is this isoform protein coding? True if it is protein-coding, otherwise False
        self.path_genomeidx = path_genomeidx

        #create IsoformSJ & TranscribeTranslate instances for both the aberrant transcript isoform & its canonical transcript isoform counterpart
        [self.iso_sj_aberr, self.obj_tt_aberr] = self.create_instances_aberrant()
        [self.iso_sj_canon, self.obj_tt_canon] = self.create_instances_canonical()


        #DELETE THIS
        # list_sj = [obj_sj_aberr]
        # sj_thres_for_isosj = -10
        # hash_pos = None
        # simulant_sj = False
        # group_sj = 0            #do not perform any grouping
        # strictly_isoform = True         #only interested in isoform_id, not just closest position, right? I'm pretty sure this is right
        # iso_sj_aberr = IsoformSJ( db_type, isoform_id, list_sj, sj_thres_for_isosj, hash_pos, simulant_sj, group_sj, strictly_isoform )


        ##DO I DELETE THIS?? - I'm pretty sure I do
        # obj_tt_aberr = TranslateTranscript( list_transcripts_ssj, iso_sj_aberr, path_genomeidx, {} )

        # canon_transcript = iso_sj_aberr.create_canon_transcript()
        # iso_sj_canon = IsoformSJ( db_type, isoform_id, canon_transcript, sj_thres_for_isosj, hash_pos, simulant_sj, group_sj, strictly_isoform )
        # obj_tt_canon = TranslateTranscript( canon_transcript, iso_sj_canon, path_genomeidx, {} )

    @staticmethod
    def obj_possible_IsoformSJ( isoform_id, obj_sj ):
        """
        Function: determines if it possible to make IsoformSJ instance. It will not be possible if the ends of the SJ instance 'obj_sj' are both outside the range of the isoform 'isoform_id'
        """
        #db_type = this will select the genomic database - used for class SpliceJunction & IsoformSJ (1 is for RefSeq, 2 is Ensembl, 3 is UCSC, & 4 is GENCODE)
        # db_type = 1     #for RefSeq
        db_type = 4     #for GENCODE

        str_sj_pos = obj_sj.str_genomic_pos()
        obj_possible_for_isosj = IsoformSJ.is_obj_possible( isoform_id, str_sj_pos, db_type )
        #if not possible to create the object, then return None
        if not obj_possible_for_isosj:
            return None

        #else return instance of IsoformSj
        list_sj = [obj_sj]
        sj_thres_for_isosj = -10
        hash_pos = obj_possible_for_isosj
        simulant_sj = False
        group_sj = 0            #do not perform any grouping
        strictly_isoform = True         #only interested in isoform_id, not just closest position, right? I'm pretty sure this is right
        
        #check if it is possible to create an IsoformSJ instance
        return IsoformSJ( db_type, isoform_id, list_sj, sj_thres_for_isosj, hash_pos, simulant_sj, group_sj, strictly_isoform )

    @staticmethod
    def is_possible_TranscribeTranscript( iso_sj, obj_sj ):
        """
        Args:
            iso_sj = IsoformSJ instance made from def obj_possible_IsoformSJ()
        Function: checks if it is possible to create a is_possible_TranscribeTranscript_V4 instance using IsoformSJ instance 'iso_sj'
        """
        list_transcripts_ssj = False
        if not iso_sj:      #if iso_sj is None, then not possible to create a TranscribeTranscript_V4  instance
            bool_possible = False
        else:
            # list_transcripts_ssj = iso_sj.reconstruct_transcript_single_sj_v2( obj_sj )        #ssj = single SJ, all the transcripts reconstructed by a single SJ
            list_transcripts_ssj = iso_sj.reconstruct_transcript_single_sj_v3( obj_sj )        #ssj = single SJ, all the transcripts reconstructed by a single SJ, this is the same as reconstruct_transcript_single_sj_v2(), but instead of looking for canonical SJ overlapped by aberrant SJ, just groups them based on 5' competitive splicing (reason I did this is because there is a chance a canonical SJ can overlap the aberrant SJ 'obj_sj', therefore missing the aberrant SJ effect on the transcript)
            if not list_transcripts_ssj:
                bool_possible = False
            else:
                bool_possible = TranscribeTranscript.is_obj_possible( list_transcripts_ssj, iso_sj )

        return {'bool_possible': bool_possible, 'list_transcripts_ssj': list_transcripts_ssj}

    def create_instances_aberrant( self ):
        """
        Creates a IsoformSJ instance for the aberrant transcript isoform for "self.isoform_id" as it includes the aberrant transcript "self.obj_sj_aberr"
        """
        list_sj = [self.obj_sj_aberr]
        sj_thres = -10
        hash_pos = None
        simulant_sj = False
        group_sj = 0            #do not perform any grouping
        strictly_isoform = True         #only interested in isoform_id, not just closest position, right? I'm pretty sure this is right

        iso_sj_aberr = SimpleNeoepitopeIsoformSJ.obj_possible_IsoformSJ( self.isoform_id, self.obj_sj_aberr )
        if not iso_sj_aberr:
            return [None, None]

        hash_info_tt = SimpleNeoepitopeIsoformSJ.is_possible_TranscribeTranscript( iso_sj_aberr, self.obj_sj_aberr )
        if not hash_info_tt['bool_possible']:
            return [None, None]

        #create the TranslateTranscript instance of the aberrant transcript
        obj_tt_aberr = TranslateTranscript( hash_info_tt['list_transcripts_ssj'], iso_sj_aberr, self.path_genomeidx, {} )

        return [iso_sj_aberr, obj_tt_aberr]


    def create_instances_canonical( self ):
        """
        Creates a IsoformSJ instance for the canonical transcript isoform for "self.isoform_id"
        """
        #create canonical form
        sj_thres = -10
        hash_pos = None
        simulant_sj = False
        group_sj = 0
        strictly_isoform = True         #only interested in isoform_id, not just closest position, right? I'm pretty sure this is right

        if not self.iso_sj_aberr:
            return [None, None]

        canon_transcript = self.iso_sj_aberr.create_canon_transcript()
        if not canon_transcript:        #if "canon_transcript", this means isoform does not have any splice junctions - could be because only 1 exon
            return [None, None]
        
        iso_sj_canon = IsoformSJ( self.db_type, self.isoform_id, canon_transcript, sj_thres, hash_pos, simulant_sj, group_sj, strictly_isoform )

        # ##TEST::
        # print "SNIsoSJ.CIC 0: show the canon_transcript - ", canon_transcript, " & len = ", len( canon_transcript )
        # for i, t_sj in enumerate( canon_transcript ):
        #     print "SNIsoSJ.CIC - ", i, " & t_sj = ", t_sj

        #create the TranslateTranscript instance of the canonical transcript
        obj_tt_canon = TranslateTranscript( canon_transcript, iso_sj_canon, self.path_genomeidx, {} )

        return [iso_sj_canon, obj_tt_canon]


    def retrieve_fiveprime_sj_pair( self ):
        """
        Retrieves the 5' most-positioned aberrant SJ and its corresponding 5' most-positioned overlapped, canonical SJ
        """

        #retrieve the canonical SJ overlapped by the aberrant SJ
        aberr_sj = self.obj_tt_aberr.retrieve_fiveprime_aberr_sj()
        canon_overlap_sj = self.obj_tt_canon.retrieve_canonical_counterpart( aberr_sj )

        #these are both SpliceJunction Instances
        return [aberr_sj, canon_overlap_sj]


    ##IS THIS DONE???
    def retrieve_aa_before_sj( self, pos_oi = None, num_aa = 14, bool_aberr = True ):
        """
        This is the same as def retrieve_aberr_canon_after_seq(), but I wanted to back up and consider amino acids before the aberrant SJ position as I will need this window for neoepitope evaluation
        OR SHOULD I JUST DO obj_tt.retrieve_aa_before_sj() --> maybe.....
        Args:
            pos_aberrant = the 5' end of the aberrant SJ event (for + strand, the lower genomic position, for - strand, the higher genomic position) . If not given, then the 5' position of the 5' most-positioned aberrant SJ is used
            num_aa = number of amino acids to retrieve before aberrant SJ event
            bool_aberr = boolean that determines which self.obj_tt_[aberr/canon] to use
                -True = use self.obj_tt_aberr
                -False = use self.obj_tt_canon
        """
        if bool_aberr:
            return self.obj_tt_aberr.retrieve_aa_before_sj( pos_oi, num_aa )
        else:
            return self.obj_tt_canon.retrieve_aa_before_sj( pos_oi, num_aa )

    
    def retrieve_aberr_canon_after_seq( self ):
        """
        retrieves the genomic position, nucleotide sequence, & protein sequence after the 5' most-position aberrant splicing event and its corresponding canonical transcript - NOTE there there should only be a single aberrant SJ incorporated in the transcript isoform
        -NOTE: the strand for the both self.iso_sj_aberr & self.iso_sj_canon should be the same
        """

        [aberr_sj, canon_overlap_sj] = self.retrieve_fiveprime_sj_pair()
        if not aberr_sj or not canon_overlap_sj:
            return [None, None]

        #retrieves the genomic position, nucleotide sequence, and protein sequence for the aberrant transcript
        pos_aberr_oi = aberr_sj.start if self.iso_sj_aberr.strand < 0 else aberr_sj.end + 1      #retrieves SJ end position
        hash_after_aberr = self.obj_tt_aberr.retrieve_aa_after_sj( pos_aberr_oi )

        #retrieves the genomic position, nucleotide sequence, and protein sequence for the aberrant transcript
        #retrievs SJ start position
        #NOTE: the strand for both iso_sj_aberr & iso_sj_canon should be the same, which is why I used "self.iso_sj_aberr" here instead of "self.iso_sj_canon"
        pos_canon_oi = canon_overlap_sj.start if self.iso_sj_aberr.strand < 0 else canon_overlap_sj.end + 1       #retrieves SJ end position
        hash_after_canon = self.obj_tt_canon.retrieve_aa_after_sj( pos_canon_oi )

        return [hash_after_aberr, hash_after_canon]

    """
    Evaluate NMD & NSD
    """
    def get_nmd_sensitive_region( self ):
        """
        Args:
            -self = instance of class TranscribeTranscript/TranslateTranscript
            -pos_aberrant_fiveprime = integer that is the position that contains the aberrant position of interest. This will usually be the start of aberrant SJ (lower position for + genes & higher position for the - genes)
        Function: retrieve the region of the gene that, if it contains an early stop codon, will lead to degradation of the transcript. Returns the nucleotide sequence of this region
        """
        [aberr_sj, canon_overlap_sj] = self.retrieve_fiveprime_sj_pair()
        if not aberr_sj or not canon_overlap_sj:
            return {}

        pos_fiveprime_aberr = aberr_sj.end + 1 if self.iso_sj_aberr.strand < 0 else aberr_sj.start
        pos_fiveprime_canon = canon_overlap_sj.end + 1 if self.iso_sj_canon.strand < 0 else canon_overlap_sj.start

        hash_nmd_sens = {}
        #returns hash element contains an array where [0] = NMD-sensitive nucleotide sequence, [1] = the genomic range of NMD-sensitive region, [2] = the genomic position of the TP boundary (50-55 nt rule)
        hash_nmd_sens['aberrant'] = self.obj_tt_aberr.get_nmd_sensitive_region( pos_fiveprime_aberr )
        hash_nmd_sens['canon'] = self.obj_tt_canon.get_nmd_sensitive_region( pos_fiveprime_canon )

        return hash_nmd_sens

    def get_nmd_irrelevant_region( self ):
        """
        Args:
            -self = instance of TranscribeTranslate_V4 class
            -pos_aberrant_threeprime = the 3' end of the aberrant splicing event (for + strand: genomically higher position, for - strand: genomically lower position)
        Function: retrieves the region of the gene that, if it contains an early stop codon, will escape NMD. However, if no stop codon, then will lead to NSD. Returns the nucleotide sequence of this region  
        """
        [aberr_sj, canon_overlap_sj] = self.retrieve_fiveprime_sj_pair()
        if not aberr_sj or not canon_overlap_sj:
            return {}
        
        pos_threeprime_aberr = aberr_sj.start if self.iso_sj_aberr.strand < 0 else aberr_sj.end + 1
        pos_threeprime_canon = canon_overlap_sj.start if self.iso_sj_canon.strand < 0 else canon_overlap_sj.end + 1

        hash_nmd_irrel = {}
        #returns array where [0] = nucleotide sequence of NMD-irrelevant region & [1] = genomic range of NMD-irrelevant region (genomic position after TP boundary (50-55 nt rule) )
        hash_nmd_irrel['aberrant'] = self.obj_tt_aberr.get_nmd_irrelevant_region( pos_threeprime_aberr )
        hash_nmd_irrel['canon'] = self.obj_tt_canon.get_nmd_irrelevant_region( pos_threeprime_canon )

        return hash_nmd_irrel

    @staticmethod
    def retrieve_corresponding_neoeps( seq_aberr, seq_canon, neoep_len = 9 ):
        """
        Retrieves the neoepitopes generated by the aberrant splicing event and its corresponding canonical counterparts
        Args:
            -seq_aberr & seq_canon = protein sequence for the aberrant transcript & canonical transcript, respectively. NOTE that "seq_aberr" may contain a stop signal "*", meaning that is where translation should stop.
            -neoep_len = integer that is neoepitope length to consider. 9 is the default as it often seen the 9-mer peptides are the peptides presented
        """
        obj_neoep = NeoepitopeMHC( neoep_len )
        seq_aberr_shorter = str( seq_aberr.split('*')[0] )      #up until stop signal
        [hash_neoep_mut, hash_neoep_orig] = obj_neoep.sliding_window_neoepitope_v3( seq_aberr_shorter, seq_canon, 3 )

    def retrieve_mRNA_neoep_subseq_v3( self ):
        """
        compares the nucleotide sequence between "mrna_orig" & "mrna_alt" and returns the subsequence that is different between them, and returns the start & end integer position referring to the position in "subseq_mut"
        Args:
            -mrna_orig & mrna_alt = array of nucleotide sequences, where each element is a nucleotide base
            -aa_len = integer that is the number of codons (aa) to retrieve before nucleotide difference between mrna_orig & mrna_alt
        Returns: an array of 2 elements of the start [0] and end [1] of the position in the string form of "self.obj_tt_aberr.arr_nuc_seq"
        NOTE: to find the true genomic position, need to remember that "i_coding_start_orig" & "i_coding_start_alt" are offsetting the values "diff_start" & "diff_end", therefore i_start = i_coding_start_orig (or i_coding_start_alt) + diff_start & i_end = i_coding_start_alt (or i_coding_start_alt) + diff_end
        """
        #DELETE THIS!
        # mrna_orig = ''.join( self.obj_tt_canon.arr_nuc_seq )
        # mrna_alt = ''.join( self.obj_tt_aberr.arr_nuc_seq )

        #retrieve the starting coding position - only want to extract the coding region of the mRNA for comparison
        i_coding_start_orig = self.obj_tt_canon.retrieve_coding_start_pos()
        i_coding_start_alt = self.obj_tt_aberr.retrieve_coding_start_pos()
        #if there is no start coding position for either, then return None
        if not i_coding_start_orig or not i_coding_start_alt:
            return [None, None]
        #retrieve the full mRNA sequence from the start of the coding region
        mrna_orig = ''.join( self.obj_tt_canon.arr_nuc_seq[i_coding_start_orig:] )
        mrna_alt = ''.join( self.obj_tt_aberr.arr_nuc_seq[i_coding_start_alt:] )
        [ diff_start, diff_end ] = SimpleNeoepitopeIsoformV2.retrieve_mRNA_neoep_subseq_v3( mrna_orig, mrna_alt )
        
        return [ diff_start, diff_end ]


    ##CAUTION: these positions may not be exact - only use this for finding where sequence differences (amino acid) occur between original transcript & aberrant transcript
    def retrieve_mRNA_neoep_subseq_v4( self ):
        """
        NOTE: this is similar to function retrieve_mRNA_neoep_subseq_v3(), but returns the indices & genomic positions associated with the sequences that differ between the original & altered sequence
        compares the nucleotide sequence between "mrna_orig" & "mrna_alt" and returns the subsequence that is different between them, and returns the start & end integer position referring to the position in "subseq_mut"
        Args:
            -mrna_orig & mrna_alt = array of nucleotide sequences, where each element is a nucleotide base
            -aa_len = integer that is the number of codons (aa) to retrieve before nucleotide difference between mrna_orig & mrna_alt
        Returns: an array of 2 elements of the start [0] and end [1] of the position in the string form of "self.obj_tt_aberr.arr_nuc_seq"
        Retunrs: returns a nested hash element, where the outer hash has 2 elements, where each element is also a hash. The outer hash has 'canon' & 'aberr' for the features associated with the canonical & aberrant SJ, respectively. Each hash has the index for the start & end position that corresponds to obj_tt_[canon/aberr].arr_genome_pos, and the genomic range for where the canonical & aberrant transcript differ (in terms of AA seq)
        """
        #DELETE THIS!
        # mrna_orig = ''.join( self.obj_tt_canon.arr_nuc_seq )
        # mrna_alt = ''.join( self.obj_tt_aberr.arr_nuc_seq )

        #retrieve the starting coding position - only want to extract the coding region of the mRNA for comparison
        i_coding_start_orig = self.obj_tt_canon.retrieve_coding_start_pos()
        i_coding_start_alt = self.obj_tt_aberr.retrieve_coding_start_pos()
        #if there is no start coding position for either, then return None
        if not i_coding_start_orig or not i_coding_start_alt:
            return {}
        #retrieve the full mRNA sequence from the start of the coding region
        mrna_orig = ''.join( self.obj_tt_canon.arr_nuc_seq[i_coding_start_orig:] )
        mrna_alt = ''.join( self.obj_tt_aberr.arr_nuc_seq[i_coding_start_alt:] )
        [ diff_start, diff_end ] = SimpleNeoepitopeIsoformV2.retrieve_mRNA_neoep_subseq_v3( mrna_orig, mrna_alt )

        #need to calculate with the offset of the starting position for the coding region
        #canon SJ range
        i_start_canon = i_coding_start_orig + diff_start - 1

        # ##TEST::
        # print "\tSNIsoSJ: i_coding_start_orig = ", i_coding_start_orig, " & genome = ", self.obj_tt_canon.arr_genome_pos[i_coding_start_orig]
        # print "\tSNIsoSJ: i_start_canon = ", i_start_canon, " & genome = ", self.obj_tt_canon.arr_genome_pos[i_start_canon]


        i_end_canon = len(self.obj_tt_canon.arr_genome_pos) - 1 if (i_coding_start_orig + diff_end) >= len(self.obj_tt_canon.arr_genome_pos) else i_coding_start_orig + diff_end
        range_canon = self.iso_sj_canon.chrom + ":" + str(self.obj_tt_canon.arr_genome_pos[i_start_canon]) + "-" + str(self.obj_tt_canon.arr_genome_pos[i_end_canon]) 
        hash_canon = {'i_start': i_start_canon, 'i_end': i_end_canon, 'genome_range': range_canon}
        #aberrant SJ range
        i_start_aberr = i_coding_start_alt + diff_start
        i_end_aberr = len(self.obj_tt_aberr.arr_genome_pos) - 1 if (i_coding_start_alt + diff_end) >= len(self.obj_tt_aberr.arr_genome_pos) else i_coding_start_alt + diff_end

        # ##TEST::
        # print "\tSNIsoSJ: i_start_aberr = ", i_start_aberr, " & i_end_aberr = ", i_end_aberr, " & len(self.obj_tt_aberr.arr_genome_pos) = ", len(self.obj_tt_aberr.arr_genome_pos) 

        range_aberr = self.iso_sj_aberr.chrom + ":" + str(self.obj_tt_aberr.arr_genome_pos[i_start_aberr]) + "-" + str(self.obj_tt_aberr.arr_genome_pos[i_end_aberr])
        hash_aberr = {'i_start': i_start_aberr, 'i_end': i_end_aberr, 'genome_range': range_canon}
        
        return {'canon': hash_canon, 'aberr': hash_aberr}

    def retrieve_orig_alt_neoeps_v3( self, aa_len = 9 ):
        """
        compares the mRNA sequence from the original transcript and aberrant transcript and returns the mRNA sequence that is different between them
        Returns: an array of 2 elements where each is the subsequence "aa_len" positions before the genomic alteration, where [0] = array that is original mRNA sequence (each element is a nucleotide base), [1] = same as [0] but for altered mRNA sequence
        """
        #DELETE THIS - THIS IS WRONG
        # mrna_orig = ''.join( self.obj_tt_canon.arr_nuc_seq )
        # mrna_alt = ''.join( self.obj_tt_aberr.arr_nuc_seq )

        #retrieve the starting coding position - only want to extract the coding region of the mRNA for comparison
        i_coding_start_orig = self.obj_tt_canon.retrieve_coding_start_pos()
        i_coding_start_alt = self.obj_tt_aberr.retrieve_coding_start_pos()
        #if there is no start coding position for either, then return None
        if not i_coding_start_orig or not i_coding_start_alt:
            return [None, None]
        #retrieve the full mRNA sequence from the start of the coding region
        mrna_orig = ''.join( self.obj_tt_canon.arr_nuc_seq[i_coding_start_orig:] )
        mrna_alt = ''.join( self.obj_tt_aberr.arr_nuc_seq[i_coding_start_alt:] )
        # [subseq_orig, subseq_mut] = SimpleNeoepitopeIsoformV2.retrieve_mRNA_neoep_subseq_v2( mrna_orig, mrna_alt, aa_len )
        [subseq_orig, subseq_mut] = SimpleNeoepitopeIsoformV2.retrieve_orig_alt_neoeps_v3( mrna_orig, mrna_alt, aa_len )        #this retrieves the subsequence mRNA that differs between the original & altered mRNA

        subseq_orig = ''.join( subseq_orig )
        subseq_mut = ''.join( subseq_mut )

        return [subseq_orig, subseq_mut]


    ##CAUTION: I've tested this function and so far it does reconstitute the correct AA sequences for both + & - strand genes, but may require more testing
    def retrieve_comparative_neoeps_v2( self, aa_len = 9 ):
        """
        This is basically like the "SimpleNeoepitopeIsoformV2.retrieve_comparative_neoeps_v2()" where this retrieves the different amino acid sequences based on the mRNA of the original transcript and the aberrant SJ transcript
        Args:
            -aa_len = integer that is the number of codons (aa) to retrieve before nucleotide difference between mrna_orig & mrna_alt
        Returns: same output as algorithm SimpleNeoepitopeIsoformV2.retrieve_comparative_neoeps_v2(), where:
            -[0] = original, unmutated AA sequence
            -[1] = altered AA sequence due to mutation
            -[2] = original, unmutated AA sequence but is the same length as the mutated sequence [3] if [3] has a termination codon
            -[3] = altered AA sequence due to mutation but is only up to the termination codon
        """
        #DELETE THIS - THIS IS WRONG
        # mrna_orig = ''.join( self.obj_tt_canon.arr_nuc_seq )
        # mrna_alt = ''.join( self.obj_tt_aberr.arr_nuc_seq )

        #retrieve the starting coding position - only want to extract the coding region of the mRNA for comparison
        i_coding_start_orig = self.obj_tt_canon.retrieve_coding_start_pos()
        i_coding_start_alt = self.obj_tt_aberr.retrieve_coding_start_pos()

        # ##TEST::
        # print "\tSNIsoSJ.RCNv2 | 0: i_coding_start_orig = ", i_coding_start_orig, " &  ", self.obj_tt_canon.arr_genome_pos[i_coding_start_orig]
        # print "\tSNIsoSJ.RCNv2 | 1: i_coding_start_alt = ", i_coding_start_alt, " &  ", self.obj_tt_aberr.arr_genome_pos[i_coding_start_alt], " & ", self.obj_tt_aberr.arr_genome_pos[i_coding_start_alt-10:i_coding_start_alt + 1]
        # TEST_nuc_seq = self.obj_tt_aberr.arr_nuc_seq[i_coding_start_alt-10:i_coding_start_alt + 1]
        # print "\tSNIsoSJ.RCNv2 | 1: TEST_nuc_seq = ", TEST_nuc_seq, " & reverse = ", TEST_nuc_seq[::-1]


        #if there is no start coding position for either, then return None
        if not i_coding_start_orig or not i_coding_start_alt:
            return [None, None, None, None]
        #retrieve the full mRNA sequence from the start of the coding region
        if self.iso_sj_aberr.strand < 0:        #for minus strand genes
            nuc_seq_orig = self.obj_tt_canon.arr_nuc_seq[:i_coding_start_orig + 1]
            nuc_seq_alt = self.obj_tt_aberr.arr_nuc_seq[:i_coding_start_alt + 1]
            mrna_orig = ''.join( nuc_seq_orig[::-1] )
            mrna_alt = ''.join( nuc_seq_alt[::-1] )
        else:       #for plus strand genes
            mrna_orig = ''.join( self.obj_tt_canon.arr_nuc_seq[i_coding_start_orig:] )
            mrna_alt = ''.join( self.obj_tt_aberr.arr_nuc_seq[i_coding_start_alt:] )
        [aa_orig, aa_alt, aa_orig_2, aa_alt_2] = SimpleNeoepitopeIsoformV2.retrieve_comparative_neoeps_v2( mrna_orig, mrna_alt, aa_len )

        return [aa_orig, aa_alt, aa_orig_2, aa_alt_2]
