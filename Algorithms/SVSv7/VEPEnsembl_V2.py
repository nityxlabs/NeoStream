#/usr/bin/python
import re
import requests

from VEPTranscribeTranslate_V5 import VEPTranscribeTranscript, VEPTranslateTranscript

"""
Questions:
-Q: Should I retrieve the gene information associated with alteration (e.g. )
    -CONJ: I think I should create a hash that will record the following: Ensembl ID (gene & trancript), CDS position, reading frame, strand sign for transcript ID - WHAT ELSE?
-Q: what about frameshift-inducing events (e.g. indels) - how do I retrieve the AAs after the frameshift?
    -IDEA 1: perhaps extend the original event to include the # of AAs after - NO, this won't work because I need to specify the specific position where the genomic alteration occurs.
    -IDEA 2: I may need to retrieve the nucleotide sequence after the alteration, treat the nucleotide sequence as if it starts from reading frame 0, & then translate (need to consider strand sign, for minus strand need to get reverse complement sequence unless it is already returned by VEP) 


-Q: which strand should I choose when retrieve the correct codon? -> use the strand opposing the gene strand (e.g. if gene is on minus strand, use plus strand)
-Q: which strand should I choose when retrieve the correct AA? -> Need to use the strand opposing the gene strand (e.g. if gene is on minus strand, use plus strand), I'm assuming because the gene's strand acts as the template & the opposite strand is the "coding strand"
-Q: Does VEP always retrieve the correct codon frame based on the position retrieved? -> YES! It does retrieve the correct codon and the correct nucleotide position that is modified (according to the gene isoform)
"""

class EnsemblVEP():

    def __init__( self, genome_pos, alt, strand, path_genomeidx, build_hg38 = False ):
        """
        Args:
            -genome_pos = string that is genomic position of genomic alteration (format: chrom:start-end, BUT exclude 'chr' (e.g. 9:25-40, NOT chr9:25-40))
            -alt = string that is the genomic string change. Can be either nucleotide string (A,C,G,T) or '-' (meaning a deletion)
            -strand = integer that is either 1 (plus strand) or -1 (minus strand). NOTE: In my experience, only use the "+" strand for retrieve the correct the codon information for the gene
            -path_genomeidx = string that is the path to the samtools-indexed genome

        """
        self.genome_pos = genome_pos
        self.alt = alt
        self.strand = strand
        self.ext = "vep/human/region/" + self.genome_pos + ":" + str( self.strand ) + "/" + self.alt
        #examples for URL for 
        # url_pattern_grch38 = 'http://rest.ensembl.org/vep/Homo_sapiens/'        #pattern for hg38 (GRCh38 Build)
        # url_pattern_grch37 = 'http://grch37.rest.ensembl.org/vep/Homo_sapiens/' #pattern for hg19 (GRCh19 Build)
        if build_hg38:
            # self.server_build = "http://rest.ensembl.org/vep/Homo_sapiens/"
            # self.server_build = "http://rest.ensembl.org/vep/human/"
            self.server_build = "http://rest.ensembl.org/"
        else:
            # self.server_build = "http://grch37.rest.ensembl.org/vep/Homo_sapiens/"
            # self.server_build = "http://grch37.rest.ensembl.org/vep/human/"
            self.server_build = "http://grch37.rest.ensembl.org/"

        #create TranscribeTranslate instance - retrieve associated isoform ID
        r = requests.get( self.server_build + self.ext, headers = {'Content-Type': 'application/json'} )
        
        #retrieve info for making IsoformSJ & TranscribeTranslate instance
        alt_pos = split_genome_pos( self.genome_pos )
        hash_pos = {'chrom': 'chr' + str(alt_pos['chrom']), 'pos_oi': alt_pos['start'] }       #this is for IsoformSJ to find the closest isoform to this position
        list_info = self.get_all_alt_info()
        self.hash_isoform = {}      #k = isoform_id, v = TranslateTranscript instance of that isoform
        for each_iso in list_info:      #each element is hash
            iso_sj = IsoformSJ( each_iso['transcript_id'], 2, [], -10, hash_pos, False )
            canon_transcript = iso_sj.create_canon_transcript( False )
            self.hash_isoform[each_iso['transcript_id']] = TranslateTranscript( canon_transcript, iso_sj, path_genomeidx, {} )



    @staticmethod
    def retrieve_vep_info( r_json ):
        """
        Prints annotations based on genomic alterations 
        Args:
            r_json = r.json() that is retrieve as a request from Ensembl's REST API - VEP 
        """
        print "r_json.keys() = ", r_json.keys()
        #retrieve specific information about genomic alteration
        print "\n"
        print "gene_symbol = ", r_json['transcript_consequences'][0]['gene_symbol']
        print "gene_id = ", r_json['transcript_consequences'][0]['gene_id']
        print "transcript_id = ", r_json['transcript_consequences'][0]['transcript_id']
        print "codons = ", r_json['transcript_consequences'][0]['codons']
        print "amino acid consequence = ", r_json['transcript_consequences'][0]['amino_acids']
        print "consequence_terms = ", r_json['transcript_consequences'][0]['consequence_terms']
        print "cds_start = ", r_json['transcript_consequences'][0]['cds_start']
        print "cds_end = ", r_json['transcript_consequences'][0]['cds_end']
        print "protein_start = ", r_json['transcript_consequences'][0]['protein_start']
        print "protein_end = ", r_json['transcript_consequences'][0]['protein_end']
        print "most severe consequence = ", r_json['most_severe_consequence']
        allele_info = r_json['allele_string']      #allele_string = this is the complementary sequence to the bases on the strand. For example, if the + strand is called, then the allele_string is the complementary strand on minus strand
        print "allele_string = ", allele_info, " | original base = ", allele_info.split('/')[0], " | new base = ", allele_info.split('/')[1]
        #calculate reading frame
        hash_rf = {1: 0, 2: 1, 0: 2}        #calculate reading frame based on (cds_pos % 3), where key = remainder of 3 (cds_pos % 3),value = the reading frame (0 = start of codon, 1 = middle of codon, 2 = end of codon)
        cds_start_pos = r_json['transcript_consequences'][0]['cds_start']
        cds_end_pos = r_json['transcript_consequences'][0]['cds_end']
        print "Is this a way to calculate the reading frame?? - ", cds_start_pos, " & type = ", type( cds_start_pos ), " & reading frame = ", hash_rf[cds_start_pos % 3]
        print "Is this a way to calculate the reading frame?? - ", cds_end_pos, " & type = ", type( cds_end_pos ), " & reading frame = ", hash_rf[cds_end_pos % 3]

    def get_vep_request( self ):
        """
        Retrieves the genomic alteration information as a VEP request
        """
        r = requests.get( self.server_build + self.ext, headers = {'Content-Type': 'application/json'})
        return r.json()


    @staticmethod
    def get_each_alt_info( r_json ):
        """
        retrieves Ensembl's REST API VEP request information for each entry provided.
        """
        hash_info = {}
        # keys_transcript = ['gene_symbol', 'gene_id', 'transcript_id', 'codons', 'amino_acids', 'consequence_terms', 'cds_start', 'cds_end', 'protein_start', 'protein_end', 'cds_start', 'cds_end']
        keys_transcript = r_json['transcript_consequences'][0].keys()
        ##TEST::
        print "EnsemblVEP - GEAI keys = ", r_json['transcript_consequences'][0].keys()
        # keys_transcript = r_json['transcript_consequences'][0].keys()

        for k in keys_transcript:
            hash_info[k] = r_json['transcript_consequences'][0][k]
        #more keys
        hash_info['most_severe_consequence'] = r_json['most_severe_consequence']
        hash_info['allele_string'] = r_json['allele_string']


        #calculate reading frame
        hash_rf = {1: 0, 2: 1, 0: 2}        #calculate reading frame based on (cds_pos % 3), where key = remainder of 3 (cds_pos % 3),value = the reading frame (0 = start of codon, 1 = middle of codon, 2 = end of codon)
        cds_start_pos = r_json['transcript_consequences'][0]['cds_start']
        cds_end_pos = r_json['transcript_consequences'][0]['cds_end']
        hash_info['rf_cds_start'] = hash_rf[cds_start_pos % 3]
        hash_info['rf_cds_end'] = hash_rf[cds_end_pos % 3]

        ##TEST::
        print "Is this a way to calculate the reading frame?? - ", cds_start_pos, " & type = ", type( cds_start_pos ), " & reading frame = ", hash_rf[cds_start_pos % 3]
        print "Is this a way to calculate the reading frame?? - ", cds_end_pos, " & type = ", type( cds_end_pos ), " & reading frame = ", hash_rf[cds_end_pos % 3]

        return hash_info


    def get_all_alt_info( self ):
        """
        retrieves annotation for genomic information based on genomic alteration
        """
        r_json = self.get_vep_request()

        list_info = []
        for i in range( 0, len( r_json ) ):
            list_info.append( self.get_each_alt_info( r_json[i] ) )

        return list_info

    def hgvs_nearest_rf_start( self, pos_oi, bool_before ):
        """
        finds the nearest start of the reading frame (rf = 0) by finding the modulus of the CDS position of 3. This is possible because Ensembl's VEP returns the relative CDS position of the gene isoform. If the current position is already at rf = 0 (pos_oi % 3 == 0), then returns pos_oi.
        Args:
            -pos_oi = integer that should the current CDS position (based on the relative CDS position returned by VEP)
            -bool_before = boolean where:
                -True = looks for position before 'pos_oi' that has reading frame 0 (new_pos < pos_oi)
                -False = looks for position after 'pos_oi' that has reading frame 0 (new_pos > pos_oi)
        Returns: relative CDS integer that -> compatible with use with HGVS
        """
        curr_rf = pos_oi % 3        #calculate current reading frame
        if curr_rf == 0:
            return pos_oi
        else:
            return pos_oi - curr_rf if bool_before else pos_oi + (3 - curr_rf)      #this is new_pos


    def get_aa_before( self, num_aa ):
        """
        retrieve the amino acid sequence before the defect of interest
        Args:
            num_aa = integer that is the number of AAs to retrieve before defect (this is strand dependent)
        -IDEA 1: use the HGVS method - example: /vep/human/hgvs/ENST00000003084:c.1431_1433delTTC?content-type=application/json
            -need to retrieve the Ensembl transcript ID (ENST), or will the genomic ID also work? (ENSG)
            -need the CDS position, from start to end
        -Q: do I need the following: transcript ID, CDS position, gene strand sign? I may have to assign different AAs based on different transcript IDs
        -Q: What information should be returned by this function?
            -CONJ: the transcript ID, the CDS position range, the AA sequence, & the AA sequence length
        """
        #PROTOCOL: calculate the CDS range based on the number of AAs needed -> use the HGVS position and the transcript ID to retrieve the nucleotide sequence -> translate the nucleotide sequence in the AA sequence.
        #Q: for def get_aa_before(), how do I know where the CDS starts? Where the CDS ends? I want to make sure I'm retrieving the nucleotide sequence for AA without going out of bounds for the gene. I can't include introns (or maybe with aberrant splicing...)

        #retrieve following information: CDS position, Ensembl transcript ID
        list_nuc_seq = []       #array where each element will be a hash that will record 
        list_info = self.get_all_alt_info()
        for each_info in list_info:
            #retrieve the start of the neoepitope sequence
            bool_before = False if self.strand < 0 else True
            cds_before_start = each_info['cds_start'] - (num_aa * 3)
            cds_before_start = self.hgvs_nearest_rf_start( cds_before_start, bool_before )
            
            #create request for Ensembl's REST API - VEP
            hgvs_info = each_info['transcript_id'] + ":c." + str( cds_before_start ) + "_" + str( each_info['cds_start'] ) + "X>X?"
            ext = "vep/human/hgvs/" + hgvs_info
            #request information from VEP - this should lead to an error that will allow to be extract the nucleotide sequence that is before the range
            r = requests.get( self.server_build + ext, headers = {'Content-Type': 'application/json'} )


            ##TEST:: see the error
            print "EVEP_GAAB: repr(r.json()) = ", repr(r.json())


            diff = abs( each_info['cds_start'] - cds_before_start )
            regex_nuc = re.search( r'\([ACGT]{' + str(diff) + ',}\)', repr(r.json()) )
            nuc_seq = re.sub( '[\(\)]', '', regex_nuc.group() )     #need to remove the parentheses at the ends

            #record information about transcript_id & nucleotide sequence
            hash_nuc_seq = { 'hgvs_info': ext, 'transcript_id': each_info['transcript_id'], 'nuc_seq': nuc_seq }
            list_nuc_seq.append( hash_nuc_seq )

        return list_nuc_seq


    def get_aa_after( self, num_aa ):
        """
        see def get_aa_before() - very similar but looks for AAs after (looks )
        """
        pass

    """
    Generate Neoepitopes
    """
    def retrieve_neoepitope( self, num_aa ):
        """
        Retrieves the potential neoepitopes based on surrounding position 
        """
        pass

    def sliding_window_neoepitope_v3( neoep_seq_mut, neoep_seq_orig, window_size ):
        """
        Arg:
            neoep_seq = string that the neoepitope sequence (amino acid sequence. USE = extract sublists of length 'window_size'
            window_size = integer that is the "window size" - the length of each sublist, where each sublist is an array 
        Function: sliding window algorithm that slides across an array, retrieving a sublist of length 'window_size'
        NOTE: to calculate the number of neoepitopes generated based on sliding window algorithm. Number of epitopes = (length of amino acid sequence - length of sliding window) + 1. Examples: 4 amino acids, sliding window is 2 -> 4 - 2 + 1 = 3, 6 - 3 + 1 = 4, 17 - 9 + 1 = 9
        """
        #split string into an array of amino acid characters
        list_neoep_seq_mut = list( neoep_seq_mut )
        list_neoep_seq_orig = list( neoep_seq_orig )
        len_list = len( list_neoep_seq_mut )
        slide_count = 0
        hash_window_mut = {}        #For the mutated peptide sequence, k = integer (assigned from slide_count), v = hash that contains neoepitope sequence & frame
        hash_window_orig = {}        #For the original peptide sequence, k = integer (assigned from slide_count), v = hash that contains neoepitope sequence & frame
        while (slide_count + window_size) <= len_list:
            window_frame = str( slide_count ) + ':' + str ( (slide_count + window_size) )

            #for mutated peptide
            hash_window_mut[slide_count] = {}       #k2 = neoep_seq & window_frame, v2 = neoep_seq = neoepitope sequence & window_frame = range of window for sublist (e.g. 0:8, 1:9, 2:10)
            hash_window_mut[slide_count]['neoep_seq'] = ''.join( list_neoep_seq_mut[slide_count:(slide_count + window_size)] )
            hash_window_mut[slide_count]['window_frame'] = window_frame

            #for original peptide
            hash_window_orig[slide_count] = {}       #k2 = neoep_seq & window_frame, v2 = neoep_seq = neoepitope sequence & window_frame = range of window for sublist (e.g. 0:8, 1:9, 2:10)
            hash_window_orig[slide_count]['neoep_seq'] = ''.join( list_neoep_seq_orig[slide_count:(slide_count + window_size)] )
            hash_window_orig[slide_count]['window_frame'] = window_frame

            slide_count += 1

        return [hash_window_mut, hash_window_orig]


        