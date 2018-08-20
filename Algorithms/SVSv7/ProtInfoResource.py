#/usr/bin/python
import re
import time

import requests
from requests.packages.urllib3.exceptions import InsecureRequestWarning     #use this when "InsecureRequestWarning" is called, which usually happens when requests uses "verify = False". To silence "InsecureRequestWarning", use line: requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

import xml.etree.ElementTree as ET      #use this to convert XML to hash table

import pandas as pd

class ProtInfoResource():

    def __init__( self ):
        """
        Source of Documentation: https://research.bioinformatics.udel.edu/peptidematch/api/v2/
        """
        # self.url_get = "https://research.bioinformatics.udel.edu/peptidematchapi2/match_get"        #for the "GET" method
        self.url = "https://research.bioinformatics.udel.edu/peptidematchapi2/match_post"       #for the "POST" method --> this is my default

        self.taxonids = "9606"          #comma-delimited string that is all the organism IDs to be used (source of taxonomy IDs: http://www.uniprot.org/docs/speclist) --> 9606: this is the taxonomy ID for Humans, Homosapiens
        self.swissprot = "false"        #if "true" = only use SwissProt Database, "false" = use other databses including SwissProt
        self.uniref100 = "false"        #if "true" = only use UniRef100 Database, "false" = use other databses including UniRef100
        self.isoform = "true"          #if "true" = will retrieve the isoform ID with the peptide, "false" = do not report the isoform ID (Q: does isoform = "true" lead to more matches for a peptide?)
        self.leqi = "false"             #Treat Leucine (L) and Isoleucine (I) equivalent ("leqi" means Leucine equal to Isoleucine) --> I think I should always keep this false

    def get_pir_request( self, peptides, time_out = -1 ):
        """
        sends HTTP request to PIR (Protein Information Resource)
        Args:
            -peptides = string, a list of amino acid sequence that is delimited by commas ",". But as of now I think ONLY USE 1 PEPTIDE AT A TIME FOR NOW
        """
        hash_data = {
            "peptides": peptides,
            "taxonids": self.taxonids,
            "swissprot": self.swissprot,
            "uniref100": self.uniref100,
            "isoform": self.isoform,
            "leqi": self.leqi,
        }

        ##TEST:: print "MHC_IEDB.retrieve_epitope_analysis: url = ", url, " & hash_data = ", hash_data
        #use this line to suppress "InsecureRequestWarning" warning
        requests.packages.urllib3.disable_warnings( InsecureRequestWarning )

        try:
            if time_out < 0:
                # r = requests.post( self.url, data = hash_data )
                r = requests.post( self.url, data = hash_data, verify = False )     #using verify = False because I'm receiving error " (Caused by SSLError(SSLError(1, u'[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed"
            else:
                # r = requests.post( self.url, data = hash_data, timeout = time_out )
                r = requests.post( self.url, data = hash_data, timeout = time_out, verify = False )     #using verify = False because I'm receiving error " (Caused by SSLError(SSLError(1, u'[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed"
        except requests.exceptions.RequestException as e:
            print "PIR Request Error with url ", self.url, " & hash_data = ", hash_data, " --> ", e
            return [ -1, None ]

        if not r.ok:
            print "ProtInfoResource.get_pir_request() Error: ", self.url, " & hash_data = ", hash_data, " & r.ok = ", r.ok, " & r.status_code = ", r.status_code
            return [ 0, None ]      #return empty dataframe

        #check with r.ok to see if the request 
        return [1, r] if r.ok else [0, None]


    # ##MAY DELETE THIS BECAUSE OF def get_pep_match_count()
    # def get_pep_match_count_BACKUP( self, peptides, bool_json = False, time_out = -1 ):
    #     """
    #     retrieves information about how frequently a peptide is found across all protein-coding genes
    #     Args:
    #         -peptides = string, a list of amino acid sequence that is delimited by commas ",". But as of now I think ONLY USE 1 PEPTIDE AT A TIME FOR NOW
    #         -bool_json = boolean where:
    #             -True = will use the r.json to retrieve the peptide matches
    #             -False = will use the r.text to retrieve the peptide matches
    #     """
    #     [r_stat, r] = self.get_pir_request( peptides, time_out )

    #     if r_stat < 1:      #if below 1, this means there is an error with the request
    #         return [r_stat, {}]

    #     ##TEST::
    #     # print "PIR.get_pep_match_count: r.text = ", r.text
    #     # print "PIR.get_pep_match_count: r.json = ", r.json()

    #     #need to go through each gene match to which genes contain this peptide sequence as well
    #     """
    #     -matches = number of documents where peptide sequence found, which I'm pretty sure means number of genes containing this peptide sequence
    #     -qtime = amount of time needed to perform query in milliseconds
    #     -status = integer that is the status of the query. I know that '0' means the status is fine
    #     """
    #     if bool_json:
    #         r_json = r.json()
    #         matches = r_json['numberFound']
    #         qtime = r_json['qtime']
    #         status = r_json['status']
    #     else:
    #         obj_et = ET.fromstring( r.text )
    #         matches = -1 if obj_et.find( 'numberFound' ) == None else int( obj_et.find( 'numberFound' ).text )
    #         qtime = -1 if obj_et.find( 'qtime' ) == None else int( obj_et.find( 'qtime' ).text )
    #         status = -1 if obj_et.find( 'status' ) == None else int( obj_et.find( 'status' ).text )

    #     return [r_stat, { 'matches': matches, 'qtime': qtime, 'status': status }]


    def get_pep_match_count( self, peptides, time_out = -1 ):
        """
        retrieves information about how frequently a peptide is found across all protein-coding genes
        Args:
            -peptides = string, a list of amino acid sequence that is delimited by commas ",". But as of now I think ONLY USE 1 PEPTIDE AT A TIME FOR NOW
            -bool_json = boolean where:
                -True = will use the r.json to retrieve the peptide matches
                -False = will use the r.text to retrieve the peptide matches
        """
        [r_stat, r] = self.get_pir_request( peptides, time_out )

        if r_stat < 1:      #if below 1, this means there is an error with the request
            return [r_stat, {}]

        ##TEST::
        # print "PIR.get_pep_match_count: r.text = ", r.text
        # print "PIR.get_pep_match_count: r.json = ", r.json()

        #need to go through each gene match to which genes contain this peptide sequence as well
        """
        -matches = number of documents where peptide sequence found, which I'm pretty sure means number of genes containing this peptide sequence
        -qtime = amount of time needed to perform query in milliseconds
        -status = integer that is the status of the query. I know that '0' means the status is fine
        """
        #check if return string is in XML format
        try:
            obj_et = ET.fromstring( r.text )
            matches = -1 if obj_et.find( 'numberFound' ) == None else int( obj_et.find( 'numberFound' ).text )
            qtime = -1 if obj_et.find( 'qtime' ) == None else int( obj_et.find( 'qtime' ).text )
            status = -1 if obj_et.find( 'status' ) == None else int( obj_et.find( 'status' ).text )
            
            return [r_stat, { 'matches': matches, 'qtime': qtime, 'status': status }]
        except Exception as e:
            print "PIR.get_pep_match_count(): Return text is not in XML format: ", e

        try:
            r_json = r.json()
            matches = r_json['numberFound']
            qtime = r_json['qtime']
            status = r_json['status']

            return [r_stat, { 'matches': matches, 'qtime': qtime, 'status': status }]
        except Exception as e:
            print "PIR.get_pep_match_count(): Return text is not in JSON format: ", e

        return [r_stat, {}]



    ##NOTE: I need to add use the XML extraction method just as sometimes r.json does not work
    def get_found_pep_info( self, peptides, time_out = -1 ):
        """
        similar to def get_pep_match_count(), but also retrievs the genes that contain the peptide seqeunce "peptide"
        Args:
            -peptides = string, a list of amino acid sequence that is delimited by commas ",". But as of now I think ONLY USE 1 PEPTIDE AT A TIME FOR NOW
            -time_out = the amount of time a request should wait for a response. Usually I've seen timeouts = 0.001, but if timeout < 0, then no timeout parameter is defined
        """
        [r_stat, r] = self.get_pir_request( peptides, time_out )

        if r_stat < 1:      #if below 1, this means there is an error with the request
            return [r_stat, {}]

        #need to go through each gene match to which genes contain this peptide sequence as well
        r_json = r.json()
        list_info = []      #will record all genes that contains peptides
        for i, result in enumerate( r_json['results'] ):
            for i2, (k2, v2) in enumerate( result.iteritems() ):        #result has 2 keys associated with it, "queryPeptide" & "proteins", where "queryPeptide" is the peptide submitted for a match & "proteins" contains information of genes containing peptide sequence
                if not re.match( r'prot', k2, re.I ):
                    continue

                #retrieve information about genes containing peptide sequene
                for i3, v3 in enumerate( result['proteins'] ):
                    query_pep = result['queryPeptide']
                    gene_name = v3['name']      #gene name
                    uniprot_ac = v3['ac']       #Uniprot Accession Code
                    peptide_pos_match = v3['matchingPeptide'][0]['matchRange'][0]
                    peptide_match_start = peptide_pos_match['start']
                    peptide_match_end = peptide_pos_match['end']

                    hash_gene = {
                    'query_peptide': query_pep,
                    'gene_name': gene_name, 
                    'uniprot_ac': uniprot_ac,
                    'peptide_match_start': peptide_match_start,
                    'peptide_match_end': peptide_match_end
                    }

                    list_info.append( hash_gene )

        return [r_stat, {'matches': r_json['numberFound'], 'qtime': r_json['qtime'], 'status': r_json['status'], 'list_genes': list_info}]


    def tester_display_gene_matches( self, peptides ):
        """
        This is a tester function to display results from def get_found_pep_info()
        Args:
            -peptides = string, a list of amino acid sequence that is delimited by commas ",". But as of now I think ONLY USE 1 PEPTIDE AT A TIME FOR NOW
        """
        r_stat, hash_gene_matches = self.get_found_pep_info( peptides )

        if r_stat == 1:
            print "number of matches = ", hash_gene_matches['matches']
            print "status of query = ", hash_gene_matches['status']
            print "query time = ", hash_gene_matches['qtime']

            for i, g in enumerate( hash_gene_matches['list_genes'] ):       #g = hash that contains information about gene that contains peptide sequence
                print "gene match ", i
                for k, v in g.iteritems():
                    print "\t", i, " - ", k, ": ", v


    """
    Function: functions to process peptides
    """
    ##DO NOT DELETE THIS FUNCTION - this shows how to return multiple outputs using dataframe.apply -> this will be useful for returning number of matches for peptide & the genes that contain this peptide sequence
    # def process_peptides_endogen_v0( self, df_mhc, bool_json = False ):
    #     """
    #     NOTE: This is the same as def process_peptides_endogen() but this will return multiple outputs from inner function "eval_pep()"
    #     Processes a list of peptides for frequency of peptides as endogenously occurring, and returns the pandas dataframe with the appended columns with frequency of the peptide as endogenous and MAYBE the genes that contain each peptide
    #     Args:
    #         -df_mhc = a pandas dataframe coming from class MHC_IEDB.retrieve_epitope_analysis(), where this dataframe contains a column titled 'peptide' that contains peptide sequence of interest
    #     """
    #     #METHOD 1: assign multiple values to dataframe by returning "row" from apply function
    #     def eval_pep( row ):
    #         """
    #         evaluate affinity, where, as a rough guideline, peptides with IC50 values <50 nM are considered high affinity, <500 nM intermediate affinity and <5000 nM low affinity --> --> Source: http://tools.iedb.org/mhci/help/
    #         """
    #         col_pep = "peptide"
    #         [r_stat, hash_matches] = self.get_pep_match_count( row[col_pep], bool_json )

    #         retry_counter = 0
    #         num_retry = 4
    #         if r_stat < 1:
    #             while r_stat < 1 and retry_counter < num_retry:
    #                 retry_counter += 1
    #                 print "Failed in process_peptides_endogen(): retry ", retry_counter, " of ", num_retry
    #                 [r_stat, hash_matches] = self.get_pep_match_count( row[col_pep], bool_json )
    #             if r_stat < 1:      #if still erroring, then just print error
    #                 print "Error with retrieving matches for peptide sequence: ", pep_seq, " & retry = ", retry_counter

    #         if r_stat == 1 and hash_matches['status'] == 0:
    #             row['pep_matches'] = hash_matches['matches']
    #             row['pep_match_qstat'] = hash_matches['status']
    #         else:
    #             row['pep_matches'] = 'error'
    #             row['pep_match_qstat'] = 'error'

    #         return row

    #     df_mhc = df_mhc.apply( eval_pep, axis = 1 )

    #     return df_mhc

    ##DO NOT DELETE THIS FUNCTION - this shows how to only send peptide sequence & just retrieve the result, so this is thes simplest version of my algorithm
    # def process_peptides_endogen_v1( self, df_mhc ):
    #     """
    #     Processes a list of peptides for frequency of peptides as endogenously occurring, and returns the pandas dataframe with the appended columns with frequency of the peptide as endogenous and MAYBE the genes that contain each peptide
    #     Args:
    #         -df_mhc = a pandas dataframe coming from class MHC_IEDB.retrieve_epitope_analysis(), where this dataframe contains a column titled 'peptide' that contains peptide sequence of interest
    #     """
    #     #Method: only send peptide sequence & retrieve only the number of matches
    #     def eval_pep( pep_seq ):
    #         """
    #         evaluate affinity, where, as a rough guideline, peptides with IC50 values <50 nM are considered high affinity, <500 nM intermediate affinity and <5000 nM low affinity --> Source: http://tools.iedb.org/mhci/help/
    #         """
    #         [r_stat, hash_matches] = self.get_pep_match_count( pep_seq, bool_json )

    #         retry_counter = 0
    #         num_retry = 4
    #         if r_stat < 1:
    #             while r_stat < 1 and retry_counter < num_retry:
    #                 retry_counter += 1
    #                 print "Failed in process_peptides_endogen(): retry ", retry_counter, " of ", num_retry
    #                 [r_stat, hash_matches] = self.get_pep_match_count( pep_seq, bool_json )
    #             if r_stat < 1:      #if still erroring, then just print error
    #                 print "Error with retrieving matches for peptide sequence: ", pep_seq, " & retry = ", retry_counter

    #         return hash_matches['matches'] if r_stat == 1 and hash_matches['status'] == 0 else 'error'
    #         ##TEST:: just make sure the peptide sequences match
    #         # return str(hash_matches['matches']) + "_" + pep_seq if r_stat == 1 and hash_matches['status'] == 0 else 'error' + "_" + pep_seq

    #     #I think I can remove "pep_match_qstat", I was just using this to make sure PIR result didn't return an error
    #     df_mhc['pep_matches'] = df_mhc['peptide'].apply( eval_pep )       #axis = used to select either 'index' (row) or 'column' -> 0 or 'index': apply function to each column, 1 or 'columns': apply function to each row

    #     return df_mhc


    ##NOTE: This function does not need the "bool_json" parameter because def get_pep_match_count() does not need it
    def process_peptides_endogen( self, df_mhc, bool_json = False ):
        """
        Processes a list of peptides for frequency of peptides as endogenously occurring, and returns the pandas dataframe with the appended columns with frequency of the peptide as endogenous and MAYBE the genes that contain each peptide
        Args:
            -df_mhc = a pandas dataframe coming from class MHC_IEDB.retrieve_epitope_analysis(), where this dataframe contains a column titled 'peptide' that contains peptide sequence of interest
        """
        #METHOD 1: assign multiple values to dataframe by returning "row" from apply function
        def eval_pep( row ):
            """
            evaluate affinity, where, as a rough guideline, peptides with IC50 values <50 nM are considered high affinity, <500 nM intermediate affinity and <5000 nM low affinity --> --> Source: http://tools.iedb.org/mhci/help/
            """
            sec_break = 10       #take this many seconds to take a break if there is an error with a request
            col_pep = "peptide"
            [r_stat, hash_matches] = self.get_pep_match_count( row[col_pep] )

            ##TEST::
            # print "start ", row['start'], " & row[col_pep] = ", row[col_pep]
            # sec_break = 2
            # print "PIR.process_peptides_endogen() - take ", sec_break, " second break"
            # time.sleep( sec_break )

            retry_counter = 0
            num_retry = 4
            if r_stat < 1:
                while r_stat < 1 and retry_counter < num_retry:
                    #Take a break first, maybe overloading with
                    retry_counter += 1
                    print "Failed in process_peptides_endogen(): retry ", retry_counter, " of ", num_retry, ". Take a ", sec_break, " second break"
                    time.sleep( sec_break )     #NOTE: THIS BREAK IS ESSENTIAL, ELSE the program will relentlessly send requests to API, therefore need to take a break
                    [r_stat, hash_matches] = self.get_pep_match_count( row[col_pep] )
                if r_stat < 1:      #if still erroring, then just print error
                    print "FINAL ERROR IN process_peptides_endogen(): with retrieving matches for peptide sequence: ", row[col_pep], " & retry = ", retry_counter, " of ", num_retry

            #record information for peptide (number of matches & status of query)
            if r_stat < 1 or not hash_matches:
                row['pep_matches'] = 'error'
                row['pep_match_qstat'] = 'error'
            elif r_stat == 1 and hash_matches['status'] == 0:
                row['pep_matches'] = hash_matches['matches']
                row['pep_match_qstat'] = hash_matches['status']
            else:
                row['pep_matches'] = 'error'
                row['pep_match_qstat'] = 'error'

            #record the peptide sequence just to make sure the correct peptide sequence is returned
            row['pep_doublecheck_seq'] = row[col_pep]

            return row
        
        df_mhc = df_mhc.apply( eval_pep, axis = 1 )     #axis = used to select either 'index' (row) or 'column' -> 0 or 'index': apply function to each column, 1 or 'columns': apply function to each row

        return df_mhc

    def evaluate_pep_endogfreq( self, pep_seq, num_retry = 4 ):
        """
        NOTE: this function is very similar to nested function "def eval_pep()" in "def process_peptides_endogen()" in this case
        evaluate affinity, where, as a rough guideline, peptides with IC50 values <50 nM are considered high affinity, <500 nM intermediate affinity and <5000 nM low affinity --> --> Source: http://tools.iedb.org/mhci/help/
        Args:
            -pep_seq = string that is the peptide sequence
            -num_retry = integer that is the number of retries to attempt before quitting
        """
        sec_break = 10       #take this many seconds to take a break if there is an error with a request
        [r_stat, hash_matches] = self.get_pep_match_count( pep_seq )

        ##TEST::
        # print "start ", row['start'], " & pep_seq = ", pep_seq
        # sec_break = 2
        # print "PIR.process_peptides_endogen() - take ", sec_break, " second break"
        # time.sleep( sec_break )

        retry_counter = 0
        if r_stat < 1:
            while r_stat < 1 and retry_counter < num_retry:
                #Take a break first, maybe overloading with
                retry_counter += 1
                print "Failed in process_peptides_endogen(): retry ", retry_counter, " of ", num_retry, ". Take a ", sec_break, " second break"
                time.sleep( sec_break )     #NOTE: THIS BREAK IS ESSENTIAL, ELSE the program will relentlessly send requests to API, therefore need to take a break
                [r_stat, hash_matches] = self.get_pep_match_count( pep_seq )
            if r_stat < 1:      #if still erroring, then just print error
                print "FINAL ERROR IN process_peptides_endogen(): with retrieving matches for peptide sequence: ", pep_seq, " & retry = ", retry_counter, " of ", num_retry

        #record information for peptide (number of matches & status of query)
        hash_pep_info = {}
        if r_stat < 1 or not hash_matches:
            hash_pep_info['pep_matches'] = 'error'
            hash_pep_info['pep_match_qstat'] = 'error'
        elif r_stat == 1 and hash_matches['status'] == 0:
            hash_pep_info['pep_matches'] = hash_matches['matches']
            hash_pep_info['pep_match_qstat'] = hash_matches['status']
        else:
            hash_pep_info['pep_matches'] = 'error'
            hash_pep_info['pep_match_qstat'] = 'error'

        #record the peptide sequence just to make sure the correct peptide sequence is returned
        hash_pep_info['pep_doublecheck_seq'] = pep_seq

        return hash_pep_info




    # ##----------- JUST FOR REFERENCE (FROM EnsemblVEP.py) - DO NOT USE THIS FOR THIS CLASS -----------##


    # def get_all_alt_info( self, peptides, time_out = -1 ):
    #     """
    #     retrieves annotation for genomic information based on genomic alteration
    #     Args:
    #         -time_out = the amount of time a request should wait for a response. Usually I've seen timeouts = 0.001, but if timeout < 0, then no timeout parameter is defined
    #     Returns: returns an array where [0] = integer that is status of request (look at def get_vep_request()), where a value = 1 means this is fine else not fine, [1] = array of all info for the variant 
    #     """
    #     ##TEST:: use to calculate time for response 
    #     # start_time = time.time()

    #     [r_stat, r] = self.get_vep_request( time_out )

    #     # ##TEST::
    #     # elapse_time = time.time() - start_time
    #     # print "EnsemblVEP.GAAI 2: after VEP request - time_out = ", time_out, " | r_stat = ", r_stat, " & time for VEP request = ", elapse_time

    #     if r_stat < 1:      #if below 1, this means there is an error with the request
    #         return [r_stat, []]

    #     r_json = r.json()
    #     list_info = []
    #     for i in range( 0, len( r_json ) ):
    #         list_info.append( self.get_each_alt_info( r_json[i] ) )

    #     return [r_stat, list_info]
    

    # @staticmethod
    # def get_each_alt_info( r_json ):
    #     """
    #     retrieves Ensembl's REST API VEP request information for each entry provided.
    #     """
    #     ##TEST::
    #     # print "EnsemblVEP - GEAI keys = ", r_json['transcript_consequences'][0].keys()
    #     # print "EnsemblVEP - GEAI ALL = ", r_json
    #     # keys_transcript = r_json['transcript_consequences'][0].keys()

    #     hash_gen = {}       #will record general information about variant, k = category, v = value of category
    #     #general keys
    #     # keys_gen = ['seq_region_name', 'start', 'end', 'regulatory_feature_consequences', 'most_severe_consequence', 'allele_string', 'input', 'variant_class']
    #     keys_gen = r_json.keys()        #may want to remove key "transcript_consequences"
    #     #remove the transcript_consequences as I will record it separately later
    #     if 'transcript_consequences' in keys_gen:
    #         keys_gen.remove( 'transcript_consequences' )
    #     for k in keys_gen:
    #         hash_gen[k] = r_json[k] if k in r_json else '-'

    #     #if transcript_consequences are not present, then 
    #     if not 'transcript_consequences' in r_json:
    #         return [hash_gen, {}]

    #     hash_transcript = {}        #k = isoform ID, v = another hash where k2 = category, v2 = value of category
    #     # keys_transcript = ['gene_symbol', 'gene_id', 'transcript_id', 'codons', 'amino_acids', 'consequence_terms', 'cds_start', 'cds_end', 'protein_start', 'protein_end', 'cds_start', 'cds_end']

    #     keys_transcript = [str(x) for x in r_json['transcript_consequences'][0].keys()]
    #     for i, tc in enumerate( r_json['transcript_consequences'] ):        #tc = transcript consequences, which is a hash that contains information about consequences of genomic alterations
    #         hash_transcript[ tc['transcript_id'] ] = {}

    #         ##TEST:: see the keys associated with each transcript
    #         # print "EVEP - GEAI: transcript keys = ", tc.keys()

    #         for k,v in tc.iteritems():      #k = name of specific transcript consequence, v = value of that consequence
    #             hash_transcript[ tc['transcript_id'] ][k] = v

    #     return [hash_gen, hash_transcript]