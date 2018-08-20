#/usr/bin/python
import re
import requests
import collections
import time

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

#Examples of extension - RAF1 (minus strand, multiple isoforms)
# ext_ensembl_1 = "/vep/human/region/3:12641705-12641707:1/AAA"            #cds_pos = 934, codons =  GTG/TTT, AA = V/F, read_frame = 0 --> this is correct
# ext_ensembl_1 = "/vep/human/region/3:12641705-12641707:-1/AAA"            #cds_pos = 934, codons =  CAC/TTT, AA = H/F, read_frame = 0
# ext_ensembl_1 = "/vep/human/region/3:12641706-12641706:1/G"            #cds_pos = 289, codons =  Acg/Ccg, AA = T/P, read_frame = 0

class EnsemblVEP():

    def __init__( self, genome_pos, orig, alt, strand = 1, opt_param = None, build_hg38 = False ):
        """
        Args:
            -genome_pos = string that is genomic position of genomic alteration (format: chrom:start-end, BUT exclude 'chr' (e.g. 9:25-40, NOT chr9:25-40))
            -orig = string that is suppose to be the reference nucleotide base, but is either None or '-', where '-' means it is an insertion. If it is any other alteration (mutation, deletion), then this should be "None"
            -alt = string that is the genomic string change. Can be either nucleotide string (A,C,G,T) or '-' (meaning a deletion)
            -strand = integer that is either 1 (plus strand) or -1 (minus strand). NOTE: In my experience, only use the "+" strand for retrieve the correct the codon information for the gene
            -opt_param = Optional Parameters. This will be appended on the url after "?"
                -Optional parameters list for vep/human/region: https://rest.ensembl.org/documentation/info/vep_region_get
                -Example 1: refseq=1
                -Example 2: refseq=1&hgvs=1
                -Example 3: refseq=1&hgvs=1&variant_class=1
                -"hgvs" is HIGHLY RECOMMENDED since this determines the nucleotide & AA consequence of an alteration
                -"variant_class" is HIGHLY RECOMMENDED since this determines the type of alteration
                -use "refseq" if I need RefSeq ID, else Ensembl ID will be returned
            -build_hg38 = boolean -> True = use hg38 build, False = use hg19 build
        Examples of URL:
            -http://grch37.rest.ensembl.org/vep/human/region/3:12641706-12641706:1/G?refseq=1&hgvs=1
                -this contains optional parameter after "?" -> "refseq=1&hgvs=1"
        """
        if "chr" in genome_pos:
            genome_pos = genome_pos.replace( "chr", "" )

        if orig == '-':
            [genome_pos, correct_vep_call] = EnsemblVEP.insertion_change_genome_pos( genome_pos )
        else:
            correct_vep_call = True

        self.genome_pos = genome_pos
        """
        boolean that keeps track if the VEP request is made correctly. True = VEP request is correct, else False = incorrectly made. False can happen if:
            -for insertions, the start & end position have a numerical difference > 1, this will cause an issue with VEP calling
        """
        self.correct_vep_call = correct_vep_call
        self.alt = alt
        self.strand = strand
        # self.ext = "vep/human/region/" + self.genome_pos + ":" + str( self.strand ) + "/" + self.alt

        self.ext = "vep/human/region/" + self.genome_pos + ":" + str(1) + "/" + self.alt
        if opt_param:
            self.ext += "?" + opt_param
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

        # #create TranscribeTranslate instance - retrieve associated isoform ID
        # r = requests.get( self.server_build + self.ext, headers = {'Content-Type': 'application/json'} )
        

    @classmethod
    def insertion_change_genome_pos( cls_obj, genome_pos ):
        """
        this function adapts the genomic position if it is experiencing an insertion, where the start position is +1 greater than the end position (e.g. chr9:5-4)
        """
        hash_pos = cls_obj.split_genome_pos( genome_pos )
        if abs( hash_pos['end'] - hash_pos['start'] ) > 1:
            correct_vep_call = False
            str_genome_pos = genome_pos
        elif hash_pos['end'] > hash_pos['start']:
            correct_vep_call = True
            str_genome_pos = hash_pos['chrom'] + ':' + str(hash_pos['end']) + '-' + str(hash_pos['start'])
        else:
            correct_vep_call = True
            str_genome_pos = genome_pos

        return [str_genome_pos, correct_vep_call]

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

    ##MAY DELETE
    # def get_vep_request( self ):
    #     """
    #     Retrieves the genomic alteration information as a VEP request
    #     """
    #     r = requests.get( self.server_build + self.ext, headers = {'Content-Type': 'application/json'})
    #     #check with r.ok to see if the request 
    #     return r.json() if r.ok else None

    def get_vep_request( self, time_out = -1 ):
        """
        Retrieves the genomic alteration information as a VEP request
        Args:
            -time_out = the amount of time a request should wait for a response. Usually I've seen timeouts = 0.001, but if timeout < 0, then no timeout parameter is defined
        Returns: returns an array where [0] = integer that is status of request & [1] = the request object if the request is successful, else None
            -[0] values:
                - -1 = request error
                -0 = request not ok (r.ok is False)
                -1 = request is ok
        """
        ##TEST:: print "EnsemblVEP.get_vep_request 1: url = ", self.server_build + self.ext

        try:
            if time_out < 0:
                r = requests.get( self.server_build + self.ext, headers = {'Content-Type': 'application/json'} )
            else:
                r = requests.get( self.server_build + self.ext, headers = {'Content-Type': 'application/json'}, timeout = time_out )
        except requests.exceptions.RequestException as e:
            print "VEP Request Error: ", e
            return [-1, None]
        #check with r.ok to see if the request 
        # return r.json() if r.ok else None
        return [1, r] if r.ok else [0, None]


    ##NEED TO FIX THIS ASAP TO REPLACE OTHER get_each_alt_info
    # @staticmethod
    # def get_each_alt_info( r_json ):
    #     """
    #     retrieves Ensembl's REST API VEP request information for each entry provided.
    #     """
    #     hash_info = EnsemblVEP.flatten_nested_hash( dict(r_json), '', '_' )

    #     # #calculate each reading frame for each 'transcript_consequences'
    #     # hash_rf = {1: 0, 2: 1, 0: 2}        #calculate reading frame based on (cds_pos % 3), where key = remainder of 3 (cds_pos % 3),value = the reading frame (0 = start of codon, 1 = middle of codon, 2 = end of codon)
    #     # for i, tc in enumerate( r_json['transcript_consequences'] ):

    #     #     key_rf_start = 'rf_cds_start_' + str(i)
    #     #     key_rf_end = 'rf_cds_end_' + str(i)

    #     #     if not 'cds_start' in tc:
    #     #         cds_start_pos = -1
    #     #         hash_info[key_rf_start] = -1
    #     #     else:
    #     #         cds_start_pos = tc['cds_start']
    #     #         hash_info[key_rf_start] = hash_rf[cds_start_pos % 3]
    #     #     if not 'cds_end' in tc:
    #     #         cds_end_pos = -1
    #     #         hash_info['rf_cds_end'] = -1
    #     #     else:
    #     #         cds_end_pos = tc['cds_end']
    #     #         hash_info['rf_cds_end'] = hash_rf[cds_end_pos % 3]

    #     return hash_info

    @staticmethod
    def get_each_alt_info( r_json ):
        """
        retrieves Ensembl's REST API VEP request information for each entry provided.
        """
        ##TEST::
        # print "EnsemblVEP - GEAI keys = ", r_json['transcript_consequences'][0].keys()
        # print "EnsemblVEP - GEAI ALL = ", r_json
        # keys_transcript = r_json['transcript_consequences'][0].keys()

        hash_gen = {}       #will record general information about variant, k = category, v = value of category
        #general keys
        # keys_gen = ['seq_region_name', 'start', 'end', 'regulatory_feature_consequences', 'most_severe_consequence', 'allele_string', 'input', 'variant_class']
        keys_gen = r_json.keys()        #may want to remove key "transcript_consequences"
        #remove the transcript_consequences as I will record it separately later
        if 'transcript_consequences' in keys_gen:
            keys_gen.remove( 'transcript_consequences' )
        for k in keys_gen:
            hash_gen[k] = r_json[k] if k in r_json else '-'

        #if transcript_consequences are not present, then 
        if not 'transcript_consequences' in r_json:
            return [hash_gen, {}]

        hash_transcript = {}        #k = isoform ID, v = another hash where k2 = category, v2 = value of category
        # keys_transcript = ['gene_symbol', 'gene_id', 'transcript_id', 'codons', 'amino_acids', 'consequence_terms', 'cds_start', 'cds_end', 'protein_start', 'protein_end', 'cds_start', 'cds_end']

        keys_transcript = [str(x) for x in r_json['transcript_consequences'][0].keys()]
        for i, tc in enumerate( r_json['transcript_consequences'] ):        #tc = transcript consequences, which is a hash that contains information about consequences of genomic alterations
            hash_transcript[ tc['transcript_id'] ] = {}

            ##TEST:: see the keys associated with each transcript
            # print "EVEP - GEAI: transcript keys = ", tc.keys()

            for k,v in tc.iteritems():      #k = name of specific transcript consequence, v = value of that consequence
                hash_transcript[ tc['transcript_id'] ][k] = v

        return [hash_gen, hash_transcript]


    def get_all_alt_info( self, time_out = -1 ):
        """
        retrieves annotation for genomic information based on genomic alteration
        Args:
            -time_out = the amount of time a request should wait for a response. Usually I've seen timeouts = 0.001, but if timeout < 0, then no timeout parameter is defined
        Returns: returns an array where [0] = integer that is status of request (look at def get_vep_request()), where a value = 1 means this is fine else not fine, [1] = array of all info for the variant 
        """
        ##TEST:: use to calculate time for response 
        # start_time = time.time()

        [r_stat, r] = self.get_vep_request( time_out )

        # ##TEST::
        # elapse_time = time.time() - start_time
        # print "EnsemblVEP.GAAI 2: after VEP request - time_out = ", time_out, " | r_stat = ", r_stat, " & time for VEP request = ", elapse_time
        

        if r_stat < 1:      #if below 1, this means there is an error with the request
            return [r_stat, []]

        r_json = r.json()
        list_info = []
        for i in range( 0, len( r_json ) ):
            list_info.append( self.get_each_alt_info( r_json[i] ) )

        return [r_stat, list_info]

    def test_display_vep( self ):
        """
        Displays results requested from Ensembl VEP based on mutation --> used for testing purposes only
        """
        list_info = self.get_all_alt_info()
        for i, (hash_gen, hash_transcript) in enumerate( list_info ):
            print i, ": ~~~~~hash_gen --> ", hash_gen
            print i, ": *****hash_transcript --> ", hash_transcript

            for i_t, (k_t, v_t) in enumerate( hash_transcript.iteritems() ):
            #TEST:: see the 
                print "\t", i, "_", i_t, ": gene_id = ", v_t['gene_id']
                print "\t", i, "_", i_t, ": transcript_id = ", v_t['transcript_id']
                print "\t", i, "_", i_t, ": codon = ", v_t['codons']
                print "\t", i, "_", i_t, ": amino acids = ", v_t['amino_acids']
                print "\t", i, "_", i_t, ": codon = ", v_t['codons'], " & amino acids = ", v_t['amino_acids']
                print "------------"

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

    @staticmethod
    def calc_reading_frame( cds_pos ):
        """
        Calculate the reading frame based on the relative CDS position. If the position is in the CDS of the gene, then VEP will give the relative CDS position, and the reading frame can be calculated by (CDS % 3) - therefore parameter 'cds_pos' needs to be the relative CDS position in order to be relevant.
        Args:
            -cds_pos = the relative CDS position given by VEP
        """
        #calculate each reading frame for each 'transcript_consequences'
        hash_rf = {1: 0, 2: 1, 0: 2}        #calculate reading frame based on (cds_pos % 3), where key = remainder of 3 (cds_pos % 3),value = the reading frame (0 = start of codon, 1 = middle of codon, 2 = end of codon)
        return hash_rf[cds_pos % 3]

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
            # print "EVEP_GAAB: repr(r.json()) = ", repr(r.json())


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


    @staticmethod
    #I did not write this function - Source: https://stackoverflow.com/questions/6027558/flatten-nested-python-dictionaries-compressing-keys
    def flatten_nested_hash( d, parent_key='', sep='_' ):
        """
        Flattens an nested hash into an a 1D hash
        Args:
            -d = "dictionary", which is basically the hash
            -parent_key = a string that may be prepended to a 
        """
        items = []
        for k, v in d.items():
            new_key = parent_key + sep + k if parent_key else k
            if isinstance(v, collections.MutableMapping):
                items.extend(flatten(v, new_key, sep=sep).items())
            else:
                items.append((new_key, v))
        return dict(items)

    @staticmethod
    def split_genome_pos( str_gene_range ):
        """
        Args:
            str_gene_range = string in the format (chrom:start-end)
        Function: splits string 'str_gene_range' (format: chrom:start-end) into individual elements
        """
        num_pos = str_gene_range.split(':')[1].split('-')
        exon_len = int( num_pos[1] ) - int( num_pos[0] )
        return {'chrom': str_gene_range.split(':')[0], 'start': int( num_pos[0] ), 'end': int( num_pos[1] ), 'exon_len': exon_len }
        