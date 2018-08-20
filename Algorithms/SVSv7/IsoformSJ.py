#/usr/bin/python

import numpy as np
import tabix

from Isoform import Isoform
from MultiIsoform import MultiIsoform
from SpliceJunction import SpliceJunction
from SJPrevalence import SJPrevalence

#CONSTANTS: column in SJ file
#NOTE: Find true left position = COLSJ_LEFTMOST + COLSJ_BLOCKLEN[0], where COLSJ_BLOCKLEN is 2 numbers separated by ",", where [0] = length of first block (see UCSC Genome Browser), [1] = length of second block (see UCSC Genome Browser)
#NOTE: Method 1 - Find true right position = COLSJ_RIGHTMOST - COLSJ_BLOCKLEN[1], where COLSJ_BLOCKLEN is 2 numbers separated by ",", where [0] = first number, [1] = second number
#NOTE: Method 2 - Find true right position = COLSJ_LEFTMOST + COLSJ_BLOCKSTART[1], where COLSJ_BLOCKSTART is 2 numbers separated by ",", where [0] = start position of first block, [1] = start position of second block
COLSJ_CHROM = 0
COLSJ_LEFTMOST = 1
COLSJ_RIGHTMOST = 2
COLSJ_SJID = 3
COLSJ_READCOUNT = 4
COLSJ_SS = 5      #COLSJ_SS  =  column splice junction strand sign
COLSJ_NUMBLOCKS = 9       #COLSJ_NUMBLOCKS  =  this indicates the number of blocks in COLSJ_BLOCKLEN & COLSJ_BLOCKSTART. This is purely a construct for .bed files
COLSJ_BLOCKLEN = 10       #COLSJ_BLOCKLEN  =  this is the block length in UCSC genome browser, NOTE: this contains 2 numbers separated by ",", where the first number is the length for the first block, and the second number is the length for the second block
COLSJ_BLOCKSTART = 11     #COLSJ_BLOCKSTART  =  this is the block start in UCSC genome browser, NOTE: this contains 2 numbers separated by ","

ROWLEN_MAX = COLSJ_BLOCKSTART

exon_base = 1           #if 0, then exons are 0-based, if 1, then exons are 1-based

class IsoformSJ( Isoform ):
    tabix_prevalence = None         #will record the file that contains the sample prevalence of splice junctions (specifically how frequently an aberrant SJ occurs for each sample)


    def __init__( self, db_type, isoform_id, list_sj, sj_thres = -2, hash_pos = None, simulant_sj = False, group_sj = 0, strictly_isoform = False ):
        """
        Args:
            -isoform_id = string that is the isoform ID to recreate the transcript
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
            -list_sj = array of SpliceJunction instances that will be used to reconstruct transcripts. This can be an empty array if need be. Furthermore, if "simulant_sj" = True, it can add canonical SJs associated with isoform.
            -hash_pos = used for Isoform class. hash that will be used to determine the closest position of interest for the Isoform class. It will have 2 keys:
                -"chrom" = string that is the chromosome of interest (format: chr# (like chr9, chr11))
                -"pos_oi" = integer that suggest which isoform to select by selecting the isoform version that is closest to this position 
            -sj_thres = read count threshold that will be used to filter "noisy, background" splice junctions from "signal" splice junctions. The reason it is set to -2 is because sometimes 'list_sj' can contain SJ instances that have a read count less than 0 (e.g. simulant SJ).
            -simulant_sj = boolean
                -True = add simulant SJ to list of known SJs (idea is perhaps sequencing missed these SJs) - will help reconstruct transcripts
                -False = do not add simulant SJs
            -group_sj = integer that will choose one of the following options:
                -0 = no grouping, so do nothing
                -5 = group SJs by competing 5' splice sites using IsoformSJ's def group_splice_donors_bracket_acceptors()
                -3 = group SJs by competing 3' splice sites (not sure when I'd need to use this)
            -strictly_isoform = boolean where
                -True = will strictly retrieve isoform_id - ignores hash_pos
                -False = will retrieve isoform_id closest to hash_pos, assuming hash_pos has been assigned or list_sj has SpliceJunction instances associated with it
        """
        ##TEST:: print "ISOFORM SJ: hash_pos = ", hash_pos

        if strictly_isoform:
            hash_pos = None
        #if hash_pos for Isoform class is not assigned, then retrieve the first SJ 
        elif not hash_pos and list_sj:
            hash_pos = {'chrom': list_sj[0].chrom, 'pos_oi': list_sj[0].start}
        super( IsoformSJ, self ).__init__( db_type, isoform_id, hash_pos )

        self.sj_thres = sj_thres
        self.list_sj = []       #records all SpliceJunction instances that pass the read count threshold
        for each_sj in list_sj:
            if each_sj.read_count >= sj_thres:
                self.list_sj.append( each_sj )

        #if simulant_sj = True, then add simulant SJ to list of SJs
        # if simulant_sj:     #if I want a threshold of a certain number of real SJs before adding simulant SJs
        #     #check if enough support for isoform_id
        #     if self.is_isoform_present( 0.5, 0 ):       #param1 (threshold_missing) = allowable percentage of SJs missing to still be deemed as present, param2 (start_end) = should start exon/end exon/both be present?
        #         #if enough support, create simulant SJ for missing canonical SJ, else cease with transcript reconstruction
        #         self.list_sj += self.sj_canon_absent( True )          #append list of simulant splice junctions to list of record SJs

        #else if I just want to create simulant SJs
        if simulant_sj:     
            self.list_sj += self.sj_canon_absent( True )          #append list of simulant splice junctions to list of record SJs

        #group competing 5' competing SJs
        self.competing_sj_five_prime = {}
        self.sj_bracket_five_prime = {}
        if group_sj == 5:
            # self.group_splice_donors_bracket_acceptors( bool_add_new_range = True, bool_canon = False, bool_canon_isoform = False )
            group_sj_fpc = self.group_splice_donors_bracket_acceptors( True, False, False )     #group_sj_fpc = Group SpliceJunctions Five Prime Competing
            self.competing_sj_five_prime = group_sj_fpc['hash_brackets']        #k = tuple that is the genomic range, where [0] = previous 3' splice site, [1] = next 3' splice site, v = array of SpliceJunctions contained within the range
            self.sj_bracket_five_prime = group_sj_fpc['hash_sj_to_bracket']     #k = SpliceJunction instance, v = tuple that is the range that contains this SpliceJunction instance
        
    @classmethod
    def is_obj_possible( cls_obj, isoform_id, pos_range, db_type = 1 ):
        """
        Args:
            -isoform_id = string that is the isoform ID
            -pos_range = string in the format 'chrom:start-end'. This will usually be the end positions of a splicing event
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        Function: checks if it is possible to create an IsoformSJ instance using the string position 'pos_range'. Returns a hash of the 'chrom' & 'pos_oi' that should be used for Isoform class (hash_pos parameter)
        """
        return Isoform.is_obj_possible( isoform_id, pos_range, db_type )


    def json_allsj_info( self, arr_sj ):
        """
        Args:
            arr_sj = array of SpliceJunction objects that will be converted to string
        Function: returns a hash that contains information about all Splice Junctions - gene name, isoform id, genomic position, exons & introns, etc. 
        """
        hash_all_sj = {}
        for i, sj in enumerate( arr_sj ):
            aberrations = sj.isoform_aberrants if sj.isoform_aberrants else ''

            hash_sj = {
            'chrom': sj.chrom,
            'start': sj.start,
            'end': sj.end,
            'strand': sj.sj_strand,
            'readCount': sj.read_count,
            'aberrant': aberrations
            }

            key_sj = 'sj' + str(i)
            hash_all_sj[key_sj] = hash_sj

        return hash_all_sj


    @staticmethod
    #NOTE: new in SVSv6: This is considering "leftmost" & "rightmost"
    def get_sj_tabix_file( db_type, path_sample_tabix, pos_range, gene_sym, isoform_id = None, bool_record_canon_isoforms = False, sj_thres = 0, bool_intronic = False ):
        """
        Args:
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
            -path_sample_tabix = string that is the path to a tabix-indexed file of the sample's junctions.bed file
            -pos_range = string that is format 'chrom:start-end' - will be used by tabix.querys() to find any elements within range
            -sj_thres = read count threshold that will be used to filter "noisy, background" splice junctions from "signal" splice junctions
            -isoform_id = string that is the isoform ID. NEEDS TO MATCH FORMAT FOR "db_type" (e.g. if db_type = 1, then isoform ID needs to be in RefSeq format)
                -NOTE: Reason for "isoform_id" is for SpliceJunction class. Assigning isoform_id creates only Isoform instance for each SJ -> this can save a lot of time and memory when generating each SpliceJunction instance (one Isoform instance for MultiIsoform as opposed to 70 isoforms. I think gene TTN has this issue)
            -bool_record_canon_isoforms = boolean where
                -True = will record all isoforms where this SJ is canonical. Use this when I want to find constitutive exons or see how prevalent an exon is across isoforms
                -False = will only record isoform in 'isoform_id', assuming it is given.
        Function: retrieves splice junctions within file that are within a genomic range 'pos_range' - use this if I want to retrieve a list of SJs from 
        """
        tabix_file = tabix.open( path_sample_tabix )
        try:
            get_sj = tabix_file.querys( pos_range )
        except:
            print "Error in IsoformSJ.sj_in_range: Tabix Error with range ", pos_range
            return []


        #retrieve all aberrant splice junctions
        hash_aberrant_prevalence = IsoformSJ.aberrant_prevalence( pos_range )       #key = string (format: 'chrom:start-end') & value = tuple where [0] = sample prevalence & [1] = control prevalence

        ##TEST:: used to count the number of SJs in the Tabix file
        # tabix_counter = 0
        # get_sj_test = tabix_file.querys( pos_range )
        # for i, x in enumerate( get_sj_test ):
        #     tabix_counter = i

        list_sj = []       #record all splice junctions found within range
        for i, each_sj in enumerate( get_sj ):
            #filter SJ that do not pass the threshold
            if int( each_sj[COLSJ_READCOUNT] ) < sj_thres:
                continue

            #retrieve SJs true end positions -> find if the SJ is prevalent -> record SJ
            true_pos = SpliceJunction.calc_sjpos( each_sj )

            #see if SJ is prevalent in other samples
            key = each_sj[COLSJ_CHROM] + ':' + str( true_pos[0] ) + '-' + str( true_pos[1] )
            if key in hash_aberrant_prevalence:
                sample_prevalence = hash_aberrant_prevalence[key][0]
                control_prevalence = hash_aberrant_prevalence[key][1]
            else:
                sample_prevalence = 0
                control_prevalence = 0

            #create splice junction object & add to list of SJs
            #NOTE: new in SVSv6: This is considering "leftmost" & "rightmost"
            #NOTE: new in SVSv7: Considering "isoform_id" & "bool_record_canon_isoforms"
            obj_sj = SpliceJunction( db_type, each_sj[COLSJ_SJID], each_sj[COLSJ_CHROM], true_pos[0], true_pos[1], each_sj[COLSJ_LEFTMOST], each_sj[COLSJ_RIGHTMOST], each_sj[COLSJ_SS], int( each_sj[COLSJ_READCOUNT] ), gene_sym, isoform_id, bool_record_canon_isoforms, sample_prevalence, control_prevalence, bool_intronic )
            list_sj.append( obj_sj )

            ##TEST:: print "IsoformSJ.sj_in_range: ", i, " of ", tabix_counter, " & start = ", obj_sj.start, " & end = ", obj_sj.end, " & canon = ", obj_sj.canon, " & isoform = ", obj_sj.assigned_isoform

        return list_sj

    @classmethod
    def set_tabix_prevalence( cls_obj, path_tabix ):
        """ Function: sets path of tabix file that contains prevalence numbers """
        cls_obj.tabix_prevalence = tabix.open( path_tabix )

    @staticmethod
    def get_prevalence( arr_tabix_sj ):
        """
        Args:
            arr_tabix_sj = row that is split into an array by '\t'
                -row from the tabix-indexed file that has information about each aberrant SJ & the number of samples that also have this # of aberrant SJs. This row is already split into an array by '\t'
                -PATH TO FILE:
        Function: retrieve the prevalence of the SJ in all samples
        """
        #get the actual end points 
        true_pos = SpliceJunction.calc_sjpos( arr_tabix_sj )
        #NOTE: increase range by 10 (- 5 & + 5), so 
        sj_range = str( arr_tabix_sj[COLSJ_CHROM] ) + ':' + str( true_pos[0] ) + '-' + str( true_pos[1] )
        prevalence = IsoformSJ.aberrant_prevalence( sj_range )     #retrieve hash where k = genomic range & value = returns tuple where [0] = sample prevalence & [1] = control prevalence

        if not prevalence or not sj_range in prevalence:
            return ( -1, -1, None )
        elif sj_range in prevalence:
            return prevalence[sj_range]

    @classmethod
    def aberrant_prevalence( cls_obj, genomic_range ):
        """
        Args:
            filter_sj = matrix of splice junctions, where each array contains an array of columns
            genomic_range = string in format 'chrom:start-end', used to extract positions from tabix file
        Function:
            formats filter_sj as to be able to find splice junction information easily & information about how many other samples contain sample splice junction & controls
        Output: 
            hash_aberrant_prevalence = hash that contains common aberrant splice junctions (sj) that should be filtered, where key = position range (format: chrom:start-end), value = tuple where [0] = # of samples that contains sj & [1] = # of controls that contains sj
        """
        #columns for tabix file 'tabix_filter'
        col_chrom = 0
        col_start = 1
        col_end = 2
        col_gene_name = 3
        col_list_samples = 4
        col_samplecount = 5
        col_controlcount = 6
        col_readcounts = 7

        try:
            tabix_sj = cls_obj.tabix_prevalence.querys( genomic_range )
        except:
            return {}

        #save information into hash
        hash_aberrant_prevalence = {}
        for arr_tabix_sj in tabix_sj:
            #key = splice junction position (chrom:start-end), value = tuple where [0] = # of total samples with SJ, [1] = # of controls that contains sj
            key = arr_tabix_sj[0] + ':' + arr_tabix_sj[1] + '-' + arr_tabix_sj[2]

            #need to record only the version of the SJ that contains the highest prevalence (I noticed that the tabix file will contain duplicates of an entry, but the prevalence count will be different -> need to make sure I record the instance that contains the highest prevalence count)
            if not key in hash_aberrant_prevalence:
                hash_aberrant_prevalence[key] = ( int(arr_tabix_sj[5]), int(arr_tabix_sj[6]), arr_tabix_sj[4] )      #[0] = # of total samples with SJ, [1] = # of controls that contains sj, [2] = all samples that contain this SJ
            else:
                saved_total_sample_prevalence = hash_aberrant_prevalence[key][0]
                if saved_total_sample_prevalence < int(arr_tabix_sj[5]):        #if saved prevalence is lower than current prevalence, then save
                    hash_aberrant_prevalence[key] = ( int(arr_tabix_sj[5]), int(arr_tabix_sj[6]), arr_tabix_sj[4] )      #[0] = # of total samples with SJ, [1] = # of controls that contains sj, [2] = all samples that contain this SJ

        return hash_aberrant_prevalence


    def get_aberrant_sj( self ):
        """ Function: retrieves aberrant splice junctions for a particular isoform """
        list_aberrant_sj = []
        for sj in self.list_sj:
            if not sj.canon and self.isoform_id in sj.assigned_isoform:
                list_aberrant_sj.append( sj )

        return list_aberrant_sj

    ##DELETE?? - I already have a line to solve for this - self.get_exon( position )
    def pos_to_exon( self, position ):
        """ Function: retrieve exons that contain a position of interest """
        # self.get_exon( position )
        pass

    def pos_to_sj( self, position ):
        """ Function: retrieves the Splice Junction objects that contain the exon """
        #retrieve exon, and retrieve SJ that contain exon
        exon = self.get_exon( position )
        return self.exon_to_sj( exon )

    def exon_to_sj( self, exon, direction = 0 ):
        """
        Args:
            exon = instance of Exon that will be used to see if it is present in splice junction
            isoform_id = string for isoform to look into
            direction = integer that will check if exon is in a specific position, i.e. start or end
                -0 = is if exon is present, position doesn't matter
                -1 = checks if exon is the start exon
                -2 = checks if exon is the end exon
        Function: retrieves the Splice Junction objects that contain the exon
        """
        #if exon is None, then no splice junctions will be connected to it
        if exon == None:
            return []

        sj_with_exon = []       #record all Splice Junction objects that contain exon
        for sj in self.list_sj:
            if sj.has_exon( exon, self.isoform_id, direction ):
                sj_with_exon.append( sj )

        return sj_with_exon

    def exon_modified_to_sj( self, exon, direction = 0 ):
        """
        Args:
            exon = instance of Exon that will be used to see if it is present in splice junction
            list_sj = array of SpliceJunction Objects that will be used to see if they contain 'exon'
            isoform_id = string for isoform to look into
            direction = integer that will check if exon is in a specific position, i.e. start or end
                -0 = if exon is present, position doesn't matter
                -1 = checks if exon is the start exon
                -2 = checks if exon is the end exon
        Function: similar to IsoformSJ.exon_to_sj_v2(), however will also consider SJs that may splice into introns, effectively "stretching an exon" (as the exon will now include the intronic region)
        """
        #if exon is None, then no splice junctions will be connected to it
        if exon == None:
            return []

        sj_with_exon = []       #record all Splice Junction objects that contain exon
        for sj in self.list_sj:
            if sj.has_exon_modified( exon, self.isoform_id, direction ):
                sj_with_exon.append( sj )

        return sj_with_exon

    def exon_to_sj_v2( self, exon, list_sj, direction = 0 ):
        """
        Args:
            exon = instance of Exon that will be used to see if it is present in splice junction
            list_sj = array of SpliceJunction Objects that will be used to see if they contain 'exon'
            isoform_id = string for isoform to look into
            direction = integer that will check if exon is in a specific position, i.e. start or end
                -0 = if exon is present, position doesn't matter
                -1 = checks if exon is the start exon
                -2 = checks if exon is the end exon
        Function: retrieves the Splice Junction objects that contain the exon. Difference between this function & exon_to_sj is that it pulls from SpliceJunction instances from "list_sj"
        """
        #if exon is None, then no splice junctions will be connected to it
        if exon == None:
            return []

        sj_with_exon = []       #record all Splice Junction objects that contain exon
        for sj in list_sj:
            if sj.has_exon( exon, self.isoform_id, direction ):
                sj_with_exon.append( sj )

        return sj_with_exon

    def exon_modified_to_sj_v2( self, exon, list_sj, direction = 0 ):
        """
        Args:
            exon = instance of Exon that will be used to see if it is present in splice junction
            list_sj = array of SpliceJunction Objects that will be used to see if they contain 'exon'
            isoform_id = string for isoform to look into
            direction = integer that will check if exon is in a specific position, i.e. start or end
                -0 = if exon is present, position doesn't matter
                -1 = checks if exon is the start exon
                -2 = checks if exon is the end exon
        Function: similar to IsoformSJ.exon_to_sj_v2(), however will also consider SJs that may splice into introns, effectively "stretching an exon" (as the exon will now include the intronic region)
        """
        #if exon is None, then no splice junctions will be connected to it
        if exon == None:
            return []

        sj_with_exon = []       #record all Splice Junction objects that contain exon
        for sj in list_sj:
            if sj.has_exon_modified( exon, self.isoform_id, direction ):
                sj_with_exon.append( sj )

        return sj_with_exon

    def get_sj_within_range( self, pos_range, inclusive_type = 3, canon_mode = 0 ):
        """
        Args:
            -pos_range = string that is genomic position (format = chr:start-end)
            -inclusive_type = integer that conveys the type of inclusion to consider
                -1 = only get SJs where the start position is within the boundary 'pos_range'
                -2 = only get SJs where the end position is within the boundary 'pos_range'
                -3 = only get SJs where the start and end position is within the boundary 'pos_range'
            -canon_mode = integer taht conveys the type of canonical 
                -0 = will return all SJs, canonical & non-canonical, within the range
                -1 = will only return canonical SJs within the range 'pos_range'
                -2 = will only return canonical SJs within the range 'pos_range' and the same isoform ID is assigned
        Function: retrieves all SJs in the list 'self.list_sj' that within the range, and returns a list of SpliceJunction instances within that range
        """
        hash_genomic_pos = Isoform.split_genome_pos( pos_range )
        #deter
        if inclusive_type == 1:
            contained_sj = [ x for x in self.list_sj if hash_genomic_pos['start'] <= x.start <= hash_genomic_pos['end'] ]
        elif inclusive_type == 2:
            contained_sj = [ x for x in self.list_sj if hash_genomic_pos['start'] <= x.end <= hash_genomic_pos['end'] ]
        else:
            contained_sj = [ x for x in self.list_sj if x.start >= hash_genomic_pos['start'] and x.end <= hash_genomic_pos['end'] ]

        #determine if I should keep canonical SJs or not
        if canon_mode == 1:
            contained_sj = [x for x in contained_sj if x.canon]
        elif canon_mode == 2:
            contained_sj = [x for x in contained_sj if x.canon and self.isoform_id in x.assigned_isoform]

        return contained_sj

    """
    Functions: Isoform Transcript Reconstruction - check if enough support for isoform to be constructed
    """

    def is_start_end_present( self, start_end = 0 ):
        """
        Args:
            start_end = integer to check if start/end is present:
                -0: neither have to be present
                -1: start must be present
                -2: end must be present
                -3: either start/end should be present
                -4: both should be present
        Function: checks if start & end exon are present in recorded SJ
        """
        if start_end == 0:
            return True
        
        #retrieve start and end exon
        start = self.get_exon_num( exon_base )      #get start exon
        end = self.get_exon_num( self.last_exon_num )       #get end exon
        if start_end == 1:
            sj_start = self.exon_to_sj( start )
            return ( len( sj_start ) > 0 )
        if start_end == 2:
            sj_end = self.exon_to_sj( end )
            return ( len( sj_end ) > 0 )
        if start_end == 3:
            sj_start = self.exon_to_sj( start )
            sj_end = self.exon_to_sj( end )
            # return True if len( sj_start ) > 0 or len( sj_end ) > 0 else False
            return ( len( sj_start ) > 0 or len( sj_end ) > 0 )
        if start_end == 4:
            sj_start = self.exon_to_sj( start )
            sj_end = self.exon_to_sj( end )
            return ( len( sj_start ) > 0 and len( sj_end ) > 0 )

    def sj_canon_present( self ):
        """ Function: finds canonical splice junctions present from all recorded splice junctions """
        canon_sj = [x for x in self.list_sj if x.canon and self.isoform_id in x.assigned_isoform]
        canon_present = []      #records list of canonical splice junctions missing. Array contains tuples where [0] = start position & [1] = end position
        for sj in canon_sj:
            if not sj in canon_present:
                canon_present.append( sj )

        return canon_present

    @staticmethod
    def create_sj_list( db_type, list_tuple_pos, chrom, strand, gene_sym, isoform_id = None ):
        """
        Args:
            -list_tuple_pos = array of tuples where each tuple contains 2 elements, [0] = start position for SJ [1] = end position for SJ, where end position > start position.
                -NOTE: one source to get SJs is from Isoform's def get_donor_acceptor_sites()
            -chrom = string in the format 'chr#' (e.g. chr9, chr13)
            -strand = string that is either '+' or '-'
            -gene_sym = string that is the name of the gene of interest
            -isoform_id = string that is the isoform_id. Make sure this isoform is associated with gene_sym else it will not be recorded in MultiIsoform class.
                -If isoform_id = None -> MultiIsoform will record all possible isoforms for gene_sym
                -If isoform_id is defined -> MultiIsoform will record for a single isoform (this will be faster)
        Function: creates a transcript where elements in array 'list_tuple_pos' are used. Creates a list of simulant SJ for the positions recorded in 'list_tuple_pos' -> makes them SpliceJunction instances
        """
        #retrieve all donor-acceptor sites (canonical end points for SJ) -> append 
        record_sj = []
        for i, sj in enumerate( list_tuple_pos ):
            id = 'SIMSJ_' + str( i )
            #NOTE: new in SVSv6: This is considering "leftmost" & "rightmost"
            record_sj.append( SpliceJunction( db_type, id, chrom, sj[0], sj[1], None, None, strand, -1, gene_sym, isoform_id ) )

        return record_sj

    def create_canon_transcript( self, bool_strand = True, bool_record_canon_isoforms = False ):
        """
        Args:
            -bool_strand = Used for self.get_donor_acceptor_sites(), where True will make a transcript with respect to strand sign & False will make the transcript from lowest to highest genomic position, regardless of strand sign.
                -NOTE: for TranscribeTranslate_V4, the 'transcript_sj' parameter (parameter that accepts the reconstructed transcript based on an array of SpliceJunction instances) the array (containing) SpliceJunction instances) should be sorted from lower to higher genomic position.
            -bool_record_canon_isoforms = boolean where
                -True = will record all isoforms where this SJ is canonical
                -False = will only record isoform in 'isoform_id', assuming it is given.
        Function: create canonical version of the transcript
        Returns: returns an array of SpliceJunction instances that compose a canonical transcript
        """

        ##TEST:: print "IsoSJ_create_canon_transcript: ", self.isoform_id


        #retrieve all donor-acceptor sites (canonical end points for SJ) -> append 
        da_sites = self.get_donor_acceptor_sites( bool_strand )
        transcript_sj = []

        ##TEST:: print "IsoformSJ.create_canon_transcript: len( da_sites ) = ", len( da_sites )

        for i, sj in enumerate( da_sites ):
            id = 'SIMSJ_' + str( i )
            convert_strand = '-' if self.strand < 0 else '+'
            #NOTE: I should assign 'self.isoform_id' therefore decreases time & memory use to generate each simulant SJ, especially for genes with many isoforms
            #NOTE: new in SVSv6: This is considering "leftmost" & "rightmost"

            ##TEST:: print "Iso_SJ.CCT: simulant id = ", id

            transcript_sj.append( SpliceJunction( self.db_type, id, self.chrom, sj[0], sj[1], None, None, convert_strand, -1, self.gene_sym, self.isoform_id, bool_record_canon_isoforms ) )

        return transcript_sj

    def create_canon_transcript_real_sj( self ):
        """ 
        Function: same as IsoformSJ's def create_canon_transcript(), but uses the splice junctions recorded in self.list_sj to reconstruct the transcript
        NOTE: need to make sure simulant SJs are made for missing canonical SJs
        """
        canon_only = True       #boolean variable = True means only reconstruct canonical transcript
        t_start = self.build_start( canon_only )
        t_full = self.build_transcript( t_start, [], 0, 0, canon_only )
        return t_full[0] if t_full else None


    def check_canon_transcript( self, sj_transcript ):
        """
        Arg:
            sj_transcript = array where each element is a SpliceJunction instance. Usually this array is created from IsoformSJ's def reconstruct_transcript()
            Function: determine if variable 'sj_transcript' matches the canonical form of a transcript for this specific isoform. Returns True if 'sj_transcript' is the canonical form else return False.
        """
        #create canonical form of transcript
        vct = self.create_canon_transcript()        #vct = virtual canonical transcript

        #if not the same length, then "sj_transcript" is not canonical
        if len( vct ) != len( sj_transcript ):
            return False

        #check each SJ to make sure they are the same
        for i in range( 0, len( vct ) ):
            if not vct[i] == sj_transcript[i]:
                return False

        return True


    def sj_canon_absent( self, create_sj = True ):
        """ 
        Args:
            create_sj = boolean
                -True = returns an array of SpliceJunction objects. If I want to be able to use missing canonical SJ (e.g. reconstruct transcript), then use this
                -False = returns an array of integers that are the positions of missing SJs. If I just want to see how many canonical SJs are missing, I'll use this.
        Function: finds canonical splice junctions missing from all recorded splice junctions
        NOTE: this can be used to retrieve simulant SJs missing from recorded from sample
        NOTE: I used -1 as the read count for simulant SJ because it is an imaginary splice junction. But I now use +1 because I need to calculate the relative abundance of each transcript
        """
        #retrieve all canonical splice junctions that are present in sample
        present_canon = [x for x in self.list_sj if x.canon and self.isoform_id in x.assigned_isoform]
        canon_absent = []      #records list of canonical splice junctions missing. Array contains tuples where [0] = start position & [1] = end position
        
        #find any canonical SJs that are missing
        da_sites = self.get_donor_acceptor_sites()
        for i, da_site in enumerate( da_sites ):
            #see if any of the SJs recorded are canonical
            canon_sj = [ x for x in self.list_sj if x.start == da_site[0] and x.end == da_site[1] ]
            
            #if the array of canonical SJ is empty, then need to make simulant SJ
            if not canon_sj:
                #if parameter 'create_sj' is True, then create SJ object, else just enter the string 'id'
                id = 'SIMSJ_' + str( i )
                #NOTE: new in SVSv6: This is considering "leftmost" & "rightmost"
                elem = SpliceJunction( self.db_type, id, self.chrom, da_site[0], da_site[1], None, None, self.strand, -1, self.gene_sym, self.isoform_id ) if create_sj else id
                canon_absent.append( elem )

        return canon_absent


    # def sj_canon_absent( self ):
    #     """ 
    #     Function: finds canonical splice junctions missing from all recorded splice junctions
    #     NOTE: this can be used to retrieve simulant SJs missing from recorded from sample
    #     """
    #     #retrieve all canonical splice junctions that are present in sample
    #     present_canon = [x for x in self.list_sj if x.canon and self.isoform_id in x.assigned_isoform]
    #     canon_absent = []      #records list of canonical splice junctions missing. Array contains tuples where [0] = start position & [1] = end position

    #     #check which canonical splice junctions are not present in sample
    #     for sj in self.canon_transcript:
    #         if not sj in present_canon and not sj in canon_absent:
    #             canon_absent.append( sj )

    #     return canon_absent


    def is_isoform_present( self, threshold_missing = 0, start_end = 0 ):
        """ 
        Args:
            threshold_missing = percentage value (e.g. ) of missing splice junctions (i.e. ExonConnect objects) to still consider isoform
            start_end = integer to check if start/end is present:
                -0: neither have to be present
                -1: start must be present
                -2: end must be present
                -3: either start/end should be present
                -4: both should be present
        Function: see which isoforms are present depending on canonical splice junctions present
        """
        threshold = np.floor( ( self.last_exon_num + 1 ) * threshold_missing )       #same as len( self.hashExonList) * threshold_missing

        #check if start or end
        bool_start_end = self.is_start_end_present( start_end )

        #retrieve all canonical SJ - list of tuples where [0] = splice start & [1] = splice end
        num_canon_absent = self.sj_canon_absent( False )

        ##TEST:: print "ISO_SJ - num_canon_absent = ", num_canon_absent

        ##TEST:: count the number of SJs (aka donor-acceptor sites)
        # da_sites = self.get_donor_acceptor_sites()
        # print "IsoformSJ.isoform_present: canon_absent = ", len( num_canon_absent ), " & threshold = ", threshold, " & # of SJ = ", len( self.hashExonList ), " & last_exon_num = ", self.last_exon_num, " & len( da_sites ) = ", len( da_sites )


        #if the number of missing SJs is less than the threshold, then there is enough support
        #NOTE: # of SJ = # of exons - 1. Here, self.last_exon = # of exons - 1
        if threshold > len( num_canon_absent ):     #means less SJs missing than threshold allows
            return True
        else:           #else not enough support
            return False

        ##INSTEAD RETURN EXONS MISSING  - THESE WILL BE USED TO GENERATE SIMULANT SJ

    # def check_isoform_present(self, isoform_id, threshold_missing = 0, start_end = 0):
    #     """
    #     Args:
    #         isoform_id = string that is the key to hash self.hashCanonTranscripts
    #         threshold_missing = percentage of missing splice junctions (i.e. ExonConnect objects) to still consider isoform
    #         start_end = integer to check if start/end is present:
    #             -0: neither have to be present
    #             -1: start must be present
    #             -2: end must be present
    #             -3: either start/end should be present
    #             -4: both should be present
    #     Function: checks if isoform is present based on the ExonConnect objects present (i.e. splice junctions) & the allowable amount of splice junctions missing
    #     """
    #     #find length of transcript - used 
    #     threshold = np.floor( len( self.hashCanonTranscripts[isoform_id] ) * threshold_missing )

    #     ##TEST:: print "ECM_CIP - 0: in check_isoform_present 0: threshold = ", threshold, " & isoform = ", isoform_id

    #     if start_end > 0:

    #         ##TEST::
    #         print "iECM_CIP - 1: & isoform = ", isoform_id

    #         if not self.check_isoform_start_end(isoform_id, start_end):
    #             ##TEST::
    #             print "iECM_CIP - 2 - start_end is not present: & isoform = ", isoform_id, "\n\n"
    #             return False


    #     ##TEST:: print "in check_isoform_present 1"

    #     #see if each ExonConnect is present
    #     count_missing = 0
    #     present = True          #boolean that is True if isoform is sufficiently present (i.e. count_missing <= threshold), else False meaning not enough splice junctions to support isoform's presence
    #     for i, each_ec in enumerate(self.hashCanonTranscripts[isoform_id]):       #each_ec = each ExonConnect object in canonical
    #         if not each_ec in self.list_sj:
    #             count_missing += 1

    #         #if number of missing splice junction exceeds threshold, break loop & return false
    #         if count_missing > threshold:
    #             present = False
    #             break


    #     ##TEST:: 
    #     print "ECM.check_isoform_present: ID = ", isoform_id, " & length = ", len(self.hashCanonTranscripts[isoform_id]) , " & threshold = ", threshold, " & missing = ", count_missing



    #     return present

    """
    Functions: Determine relative abundance of isoforms 
    """
    def filter_canon_sj( self, hash_sj, bool_isoform ):
        """
        Args:
            hash_sj = hash where key = integer or tuple, value = array of SJs
            bool_isoform = boolean that determines if canon isoform should be extracted
                -True = only preserve canonical SJs for this specific isoform ID
                -False = only preserve canonical SJs, regardless of isoform ID
        Function: this function will filter all SJs based on if they are canonical and, depending on 'bool_isoform'
        """
        filter_sj = {}      #k = the original key from 'hash_sj', v = array of SJs that are canonical
        for k,v in hash_sj.iteritems():     #k = integer or tuple that is the position
            if bool_isoform:
                get_canon_sj = [x for x in v if x.canon and self.isoform_id in v.assigned_isoform]
            else:
                get_canon_sj = [x for x in v if x.canon]

            #if the array recording canonical SJs
            if get_canon_sj:
                filter_sj[k] = get_canon_sj

        return filter_sj

    def group_splice_donors_exact_pos( self, bool_donor = True, bool_canon = False, bool_canon_isoform = False ):
        """
        Args:
            bool_donor = boolean (donor meaning 'splice donor' aka 5' splice site)
                -True: will group SJs that have the same 5' splice site (splice donor)
                -False: will group SJs that have the same 3' splice site (splice acceptor)
            bool_canon = boolean
                -True: will only record SJs that are canonical
                -False: will retrieve all SJs, regardless if they are canonical or not
            bool_canon_isoform = boolean, this depends on the parameter 'bool_canon'
                -True: will retrieve SJs that are canonical and assigned to isoform 'self.isoform_id'
                -False: will retrieve any SJs, doesn't have to be canonical or assigned to siform 'self.isoform_id'
        Function: group Splice Junction instances that have the same end, based on boolean variable 'bool_donor'. Returns a hash where the key = sj position (5' splice position if bool_donor True else 3' splice position) & v = array of SpliceJunction positions
        """
        if bool_donor:      #record for the same 5' splice site (splice donor)
            hash_group_sj = {x.start:[] for x in self.list_sj} if self.strand > 0 else {x.end:[] for x in self.list_sj}      #k = numerical position of SJ, v = array that contains SpliceJunction instances with same position record in the key of this hash
            bool_sj_start = True if self.strand > 0 else False       #determine which SJ position to get, where True means get sj start position, and False means get sj end position
        else:       #else record for the same 3' splice site (splice acceptor)
            hash_group_sj = {x.end:[] for x in self.list_sj} if self.strand > 0 else {x.start:[] for x in self.list_sj}      #k = numerical position of SJ, v = array that contains SpliceJunction instances with same position record in the key of this hash
            bool_sj_start = False if self.strand > 0 else True       #determine which SJ position to get, where True means get sj start position, and False means get sj end position

        #record each SpliceJunction object depend on position required
        for each_sj in self.list_sj:
            sj_pos = each_sj.start if bool_sj_start else each_sj.end
            hash_group_sj[ sj_pos ].append( each_sj )

        #if bool_canon is True, then retrieve all canonical SJs
        if bool_canon:
            hash_group_sj = self.filter_canon_sj( hash_group_sj, bool_canon_isoform )

        return hash_group_sj

    #NOTE: I use group_splice_donors_bracket_acceptors() over group_splice_donors_bracket_acceptors_v2() because this function not only groups SJ in bracket but also returns hash to determine which bracket an SJ is contained in..
    def group_splice_donors_bracket_acceptors( self, bool_add_new_range = True, bool_canon = False, bool_canon_isoform = False ):
        """
        Args:
            bool_add_new_range = boolean
                -True = will add new range for a SpliceJunction if it is not found in the brackets
                -False = will NOT add new range for a SpliceJunction if it is not found in the brackets
            bool_canon = boolean
                -True: will only record SJs that are canonical
                -False: will retrieve all SJs, regardless if they are canonical or not
            bool_canon_isoform = boolean, this depends on the parameter 'bool_canon'
                -True: will retrieve SJs that are canonical and assigned to isoform 'self.isoform_id'
                -False: will retrieve any SJs, doesn't have to be canonical or assigned to siform 'self.isoform_id'
        Function: groups 5' splice sites (splice donors), constitutive & alternative, together. This will help in calculating relative abundance of each isoform
        Protocol: create brackets between the 3' splice sites -> retrieve the 5' splice site for SpliceJunction instances (plus strand genes: lower position, minus strand gene: higher position)
        """

        #create hash of brackets, where each bracket is a range of 3' splice sites
        list_AS = self.get_acceptor_sites_only( False )     #list_AS = list of Acceptor Sites
        #get the start & end exons, add to the list of 3' splice sites, and then sort from least to greatest
        # start_exon = self.get_exon_num( exon_base )
        # end_exon = self.get_exon_num( self.last_exon_num )
        # pos_start = start_exon.exonPos.location.start if self.strand > 0 else start_exon.exonPos.location.end
        # pos_end = end_exon.exonPos.location.end if self.strand > 0 else end_exon.exonPos.location.start
        pos_start = self.boundary[0]        #get starting position of isoform
        pos_end = self.boundary[1]          #get ending position of isoform
        list_AS.append( pos_start )
        list_AS.append( pos_end )
        list_AS.sort( reverse = False )


        hash_sj_to_bracket = {}     #k = SpliceJunction instance, v = tuple that is the range that contains the SpliceJunction
        hash_brackets = { (list_AS[i], list_AS[i+1] ):[] for i in range(0, len(list_AS)-1) }        #hash where key = tuple that is range between 3' splice sites, v = array of SpliceJunction instances whose start position (5' splice position) is within the 3' splice bracket

        for each_sj in self.list_sj:
            #get the 5' splice site, and then find which bracket range the position lands in and record it
            sj_start = each_sj.start if self.strand > 0 else each_sj.end
            bool_pos_assign = False         #checks if SpliceJunction object is added, if not then will create a new key and value
            for k,v in hash_brackets.iteritems():       #k = tuple of position ranges, where [0] = 3' splice site of previous exon & [1] = 3' splice site of next exon, v = array of SpliceJunction instances that are between the range of interest
                if k[0] < sj_start < k[1]:        #made them both equal so they are the same for both + & - strand genes
                    v.append( each_sj )
                    hash_sj_to_bracket[each_sj] = (k[0], k[1])
                    bool_pos_assign = True
                    break
            
            #if "each_sj" is not added, then add element with new key
            if not bool_pos_assign and bool_add_new_range:
                new_range = ( each_sj.start, each_sj.start )
                hash_brackets[ new_range ] = [each_sj]
                hash_sj_to_bracket[each_sj] = new_range

        #if bool_canon is True, then retrieve all canonical SJs
        if bool_canon:
            hash_brackets = self.filter_canon_sj( hash_brackets, bool_canon_isoform )

        """returns a hash where:
        -"hash_brackets" is a hash where key = bracket range (3' end of previous exon to 3' end of next exon, therefore value = any SJ falling between this region is considered a fiveprime splicing competitor)
        -"hash_sj_to_bracket" is a where key = specific splicing event, value = the corresponding bracket range that contains the specific splicing event
        """
        return {'hash_brackets':hash_brackets, 'hash_sj_to_bracket':hash_sj_to_bracket}


    ##BM1
    def group_splice_donors_bracket_acceptors_listsj( self, list_sj, bool_add_new_range = True, bool_canon = False, bool_canon_isoform = False ):
        """
        Args:
            list_sj = an array of SpliceJunction instances that will be evaluated to determine which SJs are 5' splicing competitors
                -this is useful if there is an aberrant SJ in the list, therefore can see which canonical splicing events are competitors to the aberrant SJ event
            bool_add_new_range = boolean
                -True = will add new range for a SpliceJunction if it is not found in the brackets
                -False = will NOT add new range for a SpliceJunction if it is not found in the brackets
            bool_canon = boolean
                -True: will only record SJs that are canonical
                -False: will retrieve all SJs, regardless if they are canonical or not
            bool_canon_isoform = boolean, this depends on the parameter 'bool_canon'
                -True: will retrieve SJs that are canonical and assigned to isoform 'self.isoform_id'
                -False: will retrieve any SJs, doesn't have to be canonical or assigned to siform 'self.isoform_id'
        Function: groups 5' splice sites (splice donors), constitutive & alternative, together. This will help in calculating relative abundance of each isoform - this is the same as def group_splice_donors_bracket_acceptors() but with this function I can upload a list of SJs constituting a transcript
        Protocol: create brackets between the 3' splice sites -> retrieve the 5' splice site for SpliceJunction instances (plus strand genes: lower position, minus strand gene: higher position)
        """

        #create hash of brackets, where each bracket is a range of 3' splice sites
        list_AS = self.get_acceptor_sites_only( False )     #list_AS = list of Acceptor Sites
        #get the start & end exons, add to the list of 3' splice sites, and then sort from least to greatest
        # start_exon = self.get_exon_num( exon_base )
        # end_exon = self.get_exon_num( self.last_exon_num )
        # pos_start = start_exon.exonPos.location.start if self.strand > 0 else start_exon.exonPos.location.end
        # pos_end = end_exon.exonPos.location.end if self.strand > 0 else end_exon.exonPos.location.start
        pos_start = self.boundary[0]        #get starting position of isoform
        pos_end = self.boundary[1]          #get ending position of isoform
        list_AS.append( pos_start )
        list_AS.append( pos_end )
        list_AS.sort( reverse = False )


        hash_sj_to_bracket = {}     #k = SpliceJunction instance, v = tuple that is the range that contains the SpliceJunction
        hash_brackets = { (list_AS[i], list_AS[i+1] ):[] for i in range(0, len(list_AS)-1) }        #hash where key = tuple that is range between 3' splice sites, v = array of SpliceJunction instances whose start position (5' splice position) is within the 3' splice bracket

        for each_sj in list_sj:
            #get the 5' splice site, and then find which bracket range the position lands in and record it
            sj_start = each_sj.start if self.strand > 0 else each_sj.end
            bool_pos_assign = False         #checks if SpliceJunction object is added, if not then will create a new key and value
            for k,v in hash_brackets.iteritems():       #k = tuple of position ranges, where [0] = 3' splice site of previous exon & [1] = 3' splice site of next exon, v = array of SpliceJunction instances that are between the range of interest
                if k[0] < sj_start < k[1]:        #made them both equal so they are the same for both + & - strand genes
                    v.append( each_sj )
                    hash_sj_to_bracket[each_sj] = (k[0], k[1])
                    bool_pos_assign = True
                    break
            
            #if "each_sj" is not added, then add element with new key
            if not bool_pos_assign and bool_add_new_range:
                new_range = ( each_sj.start, each_sj.start )
                hash_brackets[ new_range ] = [each_sj]
                hash_sj_to_bracket[each_sj] = new_range


        #if bool_canon is True, then retrieve all canonical SJs
        if bool_canon:
            hash_brackets = self.filter_canon_sj( hash_brackets, bool_canon_isoform )

        """returns a hash where:
        -"hash_brackets" is a hash where key = bracket range (3' end of previous exon to 3' end of next exon, therefore value =any SJ falling between this region is considered a fiveprime splicing competitor)
        -"hash_sj_to_bracket" is a where key = SpliceJunction instance that is specific splicing event, value = the corresponding bracket range that contains the specific splicing event
        """
        return {'hash_brackets': hash_brackets, 'hash_sj_to_bracket': hash_sj_to_bracket}

    ##MAYBE DELETE - This is done
    # def group_splice_donors_bracket_acceptors_v3( self, bool_add_new_range = True, bool_canon = False, bool_canon_isoform = False ):
    #     """
    #     Function: retrieves
    #     ASSUME: assumes that self.
    #     """
    #     #create hash of brackets, where each bracket is a range of 3' splice sites
    #     list_AS = self.get_acceptor_sites_only( False )     #list_AS = list of Acceptor Sites
    #     #get the start & end exons, add to the list of 3' splice sites, and then sort from least to greatest
    #     # start_exon = self.get_exon_num( exon_base )
    #     # end_exon = self.get_exon_num( self.last_exon_num )
    #     # pos_start = start_exon.exonPos.location.start if self.strand > 0 else start_exon.exonPos.location.end
    #     # pos_end = end_exon.exonPos.location.end if self.strand > 0 else end_exon.exonPos.location.start
    #     pos_start = self.boundary[0]        #get starting position of isoform
    #     pos_end = self.boundary[1]          #get ending position of isoform
    #     list_AS.append( pos_start )
    #     list_AS.append( pos_end )
    #     list_AS.sort( reverse = False )

    #     #find the range that contains pos_oi
    #     if pos_oi:
    #         hash_brackets = {}      #key = tuple that is the range from the previous 3' splice site to the next 3' splice site, value = array that will contain SpliceJunction instances that within the 3' splice site 
    #         for i in range( 0, len(list_AS)-1 ):
    #             if pos_oi in range(list_AS[i], list_AS[i+1]):
    #                 hash_brackets[ (list_AS[i], list_AS[i+1]) ] = []


    #     pass


    def group_splice_donors_exact_pos_v2( self, bool_donor = True, pos_oi = None, bool_canon = False, bool_canon_isoform = False ):
        """
        Args:
            bool_donor = boolean (donor meaning 'splice donor' aka 5' splice site)
                -True: will group SJs that have the same 5' splice site (splice donor)
                -False: will group SJs that have the same 3' splice site (splice acceptor)
            pos_oi = integer that, if not None, will only find SpliceJunction instances that map to this position only
            bool_canon = boolean
                -True: will only record SJs that are canonical
                -False: will retrieve all SJs, regardless if they are canonical or not
            bool_canon_isoform = boolean, this depends on the parameter 'bool_canon'
                -True: will retrieve SJs that are canonical and assigned to isoform 'self.isoform_id'
                -False: will retrieve any SJs, doesn't have to be canonical or assigned to siform 'self.isoform_id
        Function: group Splice Junction instances that have the same end, based on boolean variable 'bool_donor'. Returns a hash where the key = sj position (5' splice position if bool_donor True else 3' splice position) & v = array of SpliceJunction positions
        """
        if pos_oi:      #this means to only look at one specific position instead of all positions in self.list_sj
            hash_group_sj = { pos_oi: [] }
            if bool_donor:      #if looking for 5' splice site aka splice donor position
                bool_sj_start = True if self.strand > 0 else False       #determine which SJ position to get, where True means get sj start position, and False means get sj end position
            else:               #if looking for 3' splice site aka splice acceptor position
                bool_sj_start = False if self.strand > 0 else True       #determine which SJ position to get, where True means get sj start position, and False means get sj end position
        elif bool_donor:      #record for the same 5' splice site (splice donor)
            hash_group_sj = {x.start:[] for x in self.list_sj} if self.strand > 0 else {x.end:[] for x in self.list_sj}      #k = numerical position of SJ, v = array that contains SpliceJunction instances with same position record in the key of this hash
            bool_sj_start = True if self.strand > 0 else False       #determine which SJ position to get, where True means get sj start position, and False means get sj end position
        else:       #else record for the same 3' splice site (splice acceptor)
            hash_group_sj = {x.end:[] for x in self.list_sj} if self.strand > 0 else {x.start:[] for x in self.list_sj}      #k = numerical position of SJ, v = array that contains SpliceJunction instances with same position record in the key of this hash
            bool_sj_start = False if self.strand > 0 else True       #determine which SJ position to get, where True means get sj start position, and False means get sj end position

        #record each SpliceJunction object depend on position required
        for each_sj in self.list_sj:
            sj_pos = each_sj.start if bool_sj_start else each_sj.end
            
            if sj_pos in hash_group_sj.keys():
                hash_group_sj[ sj_pos ].append( each_sj )

        #if bool_canon is True, then retrieve all canonical SJs
        if bool_canon:
            hash_group_sj = self.filter_canon_sj( hash_group_sj, bool_canon_isoform )

        return hash_group_sj


    def group_splice_donors_bracket_acceptors_v2( self, bool_add_new_range = True, pos_oi = None, bool_canon = False, bool_canon_isoform = False ):
        """
        Args:
            bool_add_new_range = boolean
                -True = will add new range for a SpliceJunction if it is not found in the brackets
                -False = will NOT add new range for a SpliceJunction if it is not found in the brackets
            poi_oi = integer that, if defined, will only record a range that contains this position. Usually if this is set, then bool_add_new_range SHOULD be set to False (as to not add new ranges)
            bool_canon = boolean
                -True: will only record SJs that are canonical
                -False: will retrieve all SJs, regardless if they are canonical or not
            bool_canon_isoform = boolean, this depends on the parameter 'bool_canon'
                -True: will retrieve SJs that are canonical and assigned to isoform 'self.isoform_id'
                -False: will retrieve any SJs, doesn't have to be canonical or assigned to siform 'self.isoform_id
        Function: groups 5' splice sites (splice donors), constitutive & alternative, together. This will help in calculating relative abundance of each isoform
        Protocol: create brackets between the 3' splice sites -> retrieve the 5' splice site for SpliceJunction instances (plus strand genes: lower position, minus strand gene: higher position)
        """

        #create hash of brackets, where each bracket is a range of 3' splice sites
        list_AS = self.get_acceptor_sites_only( False )     #list_AS = list of Acceptor Sites
        #get the start & end exons, add to the list of 3' splice sites, and then sort from least to greatest
        # start_exon = self.get_exon_num( exon_base )
        # end_exon = self.get_exon_num( self.last_exon_num )
        # pos_start = start_exon.exonPos.location.start if self.strand > 0 else start_exon.exonPos.location.end
        # pos_end = end_exon.exonPos.location.end if self.strand > 0 else end_exon.exonPos.location.start
        pos_start = self.boundary[0]        #get starting position of isoform
        pos_end = self.boundary[1]          #get ending position of isoform
        list_AS.append( pos_start )
        list_AS.append( pos_end )
        list_AS.sort( reverse = False )     #sort from least to greatest

        #find the range that contains pos_oi
        if pos_oi:
            hash_brackets = {}      #key = tuple that is the range from the previous 3' splice site to the next 3' splice site, value = array that will contain SpliceJunction instances that within the 3' splice site 
            for i in range( 0, len(list_AS)-1 ):
                if pos_oi in range(list_AS[i], list_AS[i+1]):
                    hash_brackets[ (list_AS[i], list_AS[i+1]) ] = []
            
            # hash_brackets = { (list_AS[i], list_AS[i+1]):[] for i in range(0, len(list_AS)) if pos_oi in range(list_AS[i], list_AS[i+1]) }
        else:
            hash_brackets = { (list_AS[i], list_AS[i+1]):[] for i in range(0, len(list_AS)-1 ) }

        if not hash_brackets:
            return {}

        #go through each SpliceJunction instance, and place them in the correct 3' splice site bracket
        for each_sj in self.list_sj:
            #get the 5' splice site, and then find which bracket range the position lands in and record it
            sj_start = each_sj.start if self.strand > 0 else each_sj.end
            bool_pos_assign = False         #checks if SpliceJunction object is added, if not then will create a new key and value
            for k,v in hash_brackets.iteritems():       #k = tuple of position ranges, where [0] = 3' splice site of previous exon & [1] = 3' splice site of next exon, v = 
                if k[0] <= sj_start < k[1]:
                    v.append( each_sj )
                    bool_pos_assign = True
                    break
            
            #if "each_sj" is not added, then add element with new key
            if not bool_pos_assign and bool_add_new_range:
                new_range = ( each_sj.start, each_sj.start )
                hash_brackets[ new_range ] = [each_sj]

        #if bool_canon is True, then retrieve all canonical SJs
        if bool_canon:
            hash_brackets = self.filter_canon_sj( hash_brackets, bool_canon_isoform )

        return hash_brackets

    """
    Functions: Isoform Transcript Reconstruction - check if enough support for isoform to be constructed
    """
    def transcript_exon_str( self, transcript_sj, display_stat = 1, allow_none = True ):
        """
        Args:
            transcript_sj: array where each element is a SpliceJunction instance. This array is all connected SJs that constitutes an mRNA transcript. This value can be retreive from function def reconstruct_transcript()
            canon = boolean:
                -if True, only look for exons where splicing occurs at the exon ends (aka donor-acceptor sites) - therefore "exonLeft" for starting position, & "exonRight" for ending position
                -if False, then look for all exons that are spliced together ('exonLeft', 'exonRight', 'withinExon')
            display_stat = integer with the following rules
                -1 = show element number (i.e. exon number)
                -2 = show connected positions
                -3 = show information about the Splice Junction
                -4 = show splice junction IDs
            allow_none = boolean (this is for the SpliceJunction def spliced_elems_position_range():
                -if True, then will still output string even if either of the exons from self.spliced_elems() returns None
                -if False, then will only output string only if both exons are not None, else this function will return None
        Function: returns a string of the exon positions in the the transcript
        """
        ##INTERESTING - when a function returns a hash, I can immediately call on one of the elements in the hash after the function declaration (i.e. x.spliced_elems_position_range( self.isoform_id, False, allow_none )['str_sj_pos'])
        if display_stat == 2:       #display genomic positions
            list_pos = [ x.spliced_elems_position_range( self.isoform_id, False, allow_none )['str_sj_pos'] for x in transcript_sj ]
        elif display_stat == 3:
            list_pos = [ str(x) for x in transcript_sj ]
        elif display_stat == 4:
            list_pos = [ x.sj_id for x in transcript_sj ]
        else:               #display element numbers (e.g. exon 9)
            list_exons = [ x.spliced_elems( self.isoform_id, False ) for x in transcript_sj ]
            list_pos = []
            for exon_pair in list_exons:
                if not exon_pair:       #if None, then skip
                    continue
                
                pair_str = exon_pair[0].exonPos.type + '_' + str( exon_pair[0].exonNum ) + '*' if exon_pair[0] else "None*"
                pair_str+= exon_pair[1].exonPos.type + '_' + str( exon_pair[1].exonNum ) if exon_pair[1] else "None*"
                list_pos.append( pair_str )

        if list_pos:
            return ' >> '.join( list_pos )
        else:
            return None

    


    #NOTE: new in SVSv6: This is considering "leftmost" & "rightmost"
    def reconstruct_transcript_single_sj_v2( self, obj_sj, consider_most = False ):
        """
        Args:
            -obj_sj = SpliceJunction instance, usually an aberrant splice junction. USE = will replace canonical SJ with this SJ
            -consider_most = boolean used for obj_sj.compare_overlap_sj_v2
                -True = Will consider self.leftmost & self.rightmost as endpoints to see if 'obj_sj' is contained. These are the lowest & highest genomic position boundaries where the reads that support the splicing event end. 
                -False = will considered self.start & self.end as the endpoints to see if 'obj_sj' is contained. These are the genomic start & end positions of the gapped alignments.
        Function: similar to reconstruct_transcript_single_sj(), but will use a different method for reconstructing the transcript.
        """
        #reconstruct canonical transcript -> determine which canonical transcripts overlap
        #create canonical form of transcript
        ct = self.create_canon_transcript( False )        #ct = canonical transcript

        #find affected canonical SJs by 'obj_sj'
        # overlapped_sj = [ sj for sj in ct if obj_sj.compare_overlap_sj_v2( sj, consider_most ) > 0 or sj.compare_overlap_sj_v2( obj_sj, consider_most ) ]

        overlapped_sj = [ sj for sj in ct if obj_sj.compare_overlap_sj_v2( sj, consider_most ) > 0 ]

        for each_sj in overlapped_sj:
            ct.remove( each_sj )
        #add the SJ of interest
        ct.append( obj_sj )

        #sort the transcript from least to greatest genomic position based on start positoin
        transcript_single_sj = sorted( ct, key = lambda s: s.start )

        ##TEST:: show each SJ instance in transcript_single_sj
        # for i_test, test_each_sj in enumerate( transcript_single_sj ):
        #     str_iso_id = [( ' | '.join(x.hash_isoforms.keys()) ) for x in transcript_single_sj]
        #     print "ISO_SJ - RTSSJV2 0: ", i_test, " - ", test_each_sj, " && isoform IDs = ", str_iso_id


        return transcript_single_sj

    ##BM1
    def reconstruct_transcript_single_sj_v3( self, obj_sj, consider_most = False ):
        """
        Args:
            -obj_sj = SpliceJunction instance, usually an aberrant splice junction. USE = will replace canonical SJ with this SJ
            -consider_most = boolean used for obj_sj.compare_overlap_sj_v2
                -True = Will consider self.leftmost & self.rightmost as endpoints to see if 'obj_sj' is contained. These are the lowest & highest genomic position boundaries where the reads that support the splicing event end. 
                -False = will considered self.start & self.end as the endpoints to see if 'obj_sj' is contained. These are the genomic start & end positions of the gapped alignments.
        Function: similar to reconstruct_transcript_single_sj_v2(), but in addition to looking for aberrant SJs that overlap canonical SJs, this looks for canonical SJs that are competitive 5' splicing events to aberrant SJ 'obj_sj' - so removes competitive 5' SJs & overlapped SJs
        """
        #reconstruct canonical transcript -> determine which canonical transcripts overlap
        #create canonical form of transcript
        ct = self.create_canon_transcript( False )        #ct = canonical transcript

        #Find canonical SJs overlapped by 'obj_sj'
        overlapped_sj = [ sj for sj in ct if obj_sj.compare_overlap_sj_v2( sj, consider_most ) > 0 ]
        for each_sj in overlapped_sj:
            ct.remove( each_sj )

        #retrieve all 5' competitive SJs to aberrant SJ 'obj_sj'
        ct.append( obj_sj )     #append the aberrant SJ 'obj_sj' into list of SJ constituting canoncial transcript isoform
        hash_acceptor_brackets = self.group_splice_donors_bracket_acceptors_listsj( ct )        #returns a hash of hashes -> {'hash_brackets': hash_brackets, 'hash_sj_to_bracket': hash_sj_to_bracket}
 
        key_bracket = hash_acceptor_brackets['hash_sj_to_bracket'][obj_sj]      #retrieve the bracket containing the aberrant SJ
        
        ##TEST::
        # print "hash_acceptor_brackets['hash_sj_to_bracket']: ", hash_acceptor_brackets['hash_sj_to_bracket']

        # for k,v in hash_acceptor_brackets['hash_sj_to_bracket'].iteritems():
        #     print "\tIsoSJ.RTSSJV3 - HSJBRAC: k = ", k, " & v = ", v

        # print "IsoSJ.RTSSJV3 - 0 = ", key_bracket
        # print "IsoSJ.RTSSJV3 - 1 key_bracket = ", key_bracket
        # print "IsoSJ.RTSSJV3 - 2: ", hash_acceptor_brackets['hash_brackets'][key_bracket]
       
        for each_sj in hash_acceptor_brackets['hash_brackets'][key_bracket]:
            # ##TEST::
            # print "\tIsoSJ.RTSSJV3 loop: each_sj = ", each_sj
            if each_sj != obj_sj:
                # ##TEST::
                # print "\t\tREMOVE - IsoSJ.RTSSJV3 loop: each_sj = ", each_sj
                ct.remove( each_sj )


        #sort the transcript from least to greatest genomic position based on start positoin
        transcript_single_sj = sorted( ct, key = lambda s: s.start )

        ##TEST:: show each SJ instance in transcript_single_sj
        # for i_test, test_each_sj in enumerate( transcript_single_sj ):
        #     str_iso_id = [( ' | '.join(x.hash_isoforms.keys()) ) for x in transcript_single_sj]
        #     print "ISO_SJ - RTSSJV2 0: ", i_test, " - ", test_each_sj, " && isoform IDs = ", str_iso_id


        return transcript_single_sj

    ##DON'T DELETE THIS FUNCTION YET, BUT IT HAS BEEN SUPERCEDED BY def reconstruct_transcript_single_sj_v2()
    def reconstruct_transcript_single_sj( self, obj_sj, bool_longest_transcript = True ):
        """
        Args:
            -obj_sj = SpliceJunction instance that will be used to build a transcript around 
            -bool_longest_transcript = boolean that
                -True = will return the longest transcript reconstructed (I've noticed sometimes multiple transcripts can be generated), therefore will only return an array of SpliceJunction instances that constitute a transcript
                -False = will return all results, which means this will return an array of array, where each subarray is an array of SpliceJunction instances that constitute a transcript
        Function: this will reconstruct a transcript for each aberrant SJ, only incorporating 1 aberrant SJ per transcript. This will return an array of transcripts, each transcript containing 1 aberrant SJ. Note that each transcript will be an array of SpliceJunction instances
        """
        #Retrieve all canonical SJs
        da_sites = self.get_donor_acceptor_sites( False )        #retrieve all donor-acceptor sites
        convert_ss = '-' if self.strand < 0 else '+'
        list_simulants = IsoformSJ.create_sj_list( self.db_type, da_sites, self.chrom, convert_ss, self.gene_sym )

        list_sj_oi = list_simulants[:]      #make sure to make a duplicate as to not affect original 'list_simulants'
        list_sj_oi.append( obj_sj )
        t_start = self.build_start_v2( list_sj_oi, False )

        ##TEST:: show starting exons
        print "ISO_SJ RTS.SJ 1 - t_start = ", t_start
        for i, each_t in enumerate( t_start ):
            for i2, each_sj in enumerate( each_t ):
                (x_start, x_end) = each_sj.spliced_elems( self.isoform_id )
                print i, ":", i2, " - ", each_sj, "\n-each_sj exons start = ", x_start, "\n-end exon = ", x_end
            print "****** START SJ ******\n"

        #get ready to build transcript
        max_penalty = 1         #this should be set to just 1 as it will only include 1 aberrant SJ
        canon_only = False
        t_full = self.build_transcript_v2( list_sj_oi, t_start, [], max_penalty, 0, canon_only )       #t_full = array of SJ instances that constitute an mRNA transcript

        ##TEST:: show reconstructed transcripts
        print "ISO_SJ RTS.SJ 3:"
        for i, each_t in enumerate( t_full ):
            for i2, each_sj in enumerate( each_t ):
                print i, ":", i2, " - ", each_sj
            print "-------- Next transcript -------\n"

        #only keep the transcripts that contain the aberrant splicing event
        t_keep = [x for x in t_full if obj_sj in x]

        #get longest transcript
        if t_keep and bool_longest_transcript:
            #Method 1 to retrieve longest transcript
            # max_transcript = max( t_keep, key = lambda x: len( x ) )
            #Method 2 to retrieve longest transcript
            max_transcript = max( t_keep, key = len )

            return max_transcript
        else:
            return t_keep

    ##MAY DELETE THIS AS I NOW USE def reconstruct_transcript_single_sj()
    # def find_nearest_elem( self, elem, list_sj, direction ):
    #     """
    #     Args:
    #         -elem = Exon instance that is the starting point for building one side of the transcript ()
    #         -list_sj = list of splice junctions to select from, will be used in def exon_to_sj_v2()
    #         -direction = integer that is the direction to build the transcript
    #             - -1 = build the "lower side" of the transcript (i.e. get all SJs with lower genomic positions - this is the towards that starting part of the transcript for + strand genes & the ending part of the transcript for - strand genes)
    #             - +1 = build the "higher side" of the transcript (i.e. get all SJs with higher genomic positions - this is the towards that ending part of the transcript for + strand genes & the starting part of the transcript for - strand genes)
    #     Function: builds each side of transcript, depending on variable "direction"
    #     """
    #     transcript_sj = []

    #     #if this element is an intron, then need to find the nearest exon based on parameter "direction"
    #     if sj_exons[0].exonPos.type == 'intron':
    #         pos_oi = sj_exons[0].exonPos.location.start if direction == -1 else sj_exons[0].exonPos.location.end
    #         exon_oi = self.find_closest_element( pos_oi, direction, 1 )
    #     else:
    #         exon_oi = sj_exons[0]
 
    #     bool_continue = True
    #     while( bool_continue ):
    #         prev_sj = self.exon_to_sj_v2( exon_oi, simulant_sj, 2 )
    #         if prev_sj:     #if the next SJ is found, then record SJ and use the lower-positioned as the next exon of interest
    #             transcript_sj.append( prev_sj[0] )
    #             exon_oi = prev_sj[0].spliced_elems( self.isoform_id, False )[0]
    #         else:
    #             bool_continue = False

    #     pass


    # def build_transcript_around_sj( self, obj_sj ):
    #     """
    #     Args:
    #         -obj_sj = instance of SpliceJunction class
    #         -isoform_id = string that is the isoform_id. Make sure this isoform is associated with gene_sym else it will not be recorded in MultiIsoform class.
    #             -If isoform_id = None -> MultiIsoform will record all possible isoforms for gene_sym
    #             -If isoform_id is defined -> MultiIsoform will record for a single isoform (this will be faster)
    #     Function: will build a transcript around SpliceJunction Object 'obj_sj'. Will build it around a specific isoform structure 'self.isoform_id'
    #     Test Cases:
    #         -Case 1: what happens when SJ is intronic (i.e. one or both ends land in intron)
    #             -Potential Solution 1: recreate a SpliceJunction instance to record the intronic positin
    #     """
    #     #create all simulant, canonical SJ
    #     da_sites = self.get_donor_acceptor_sites( False )

    #     str_strand = '+' if self.strand > 0 else '-'        #create_sj_list() requires strand be a string '+' or '-'
    #     simulant_sj = IsoformSJ.create_sj_list( self.db_type, da_sites, self.chrom, str_strand, self.gene_sym, self.isoform_id )

    #     #retrieve the next SJ based on the 
    #     transcript_sj = [obj_sj]      #records all SJs that will reconstruct the transcript based 
    #     sj_exons = obj_sj.spliced_elems( self.isoform_id, False )
    #     #build the "lower side" of the transcript (i.e. get all SJs with lower genomic positions - this is the towards that starting part of the transcript for + strand genes & the ending part of the transcript for - strand genes)

    #     ##TEST::print "IsoSJ BTASJ 1 - build LOWER end: exon[0] ", sj_exons[0], " & exon[1] = ", sj_exons[1]

    #     if sj_exons[0].exonPos.type == 'intron':
    #         self.find_closest_element( pos, direction, feat_search = 1 )
    #     exon_oi = sj_exons[0]
    #     bool_continue = True
    #     while( bool_continue ):
    #         prev_sj = self.exon_to_sj_v2( exon_oi, simulant_sj, 2 )
    #         if prev_sj:     #if the next SJ is found, then record SJ and use the lower-positioned as the next exon of interest
    #             transcript_sj.append( prev_sj[0] )
    #             exon_oi = prev_sj[0].spliced_elems( self.isoform_id, False )[0]
    #         else:
    #             bool_continue = False

    #     ##TEST::print "IsoSJ BTASJ 2 - build HIGHER end - ", len( transcript_sj ), " of # da sites = ", len( da_sites )

    #     #build the "higher side" of the transcript (i.e. get all SJs with higher genomic positions - this is the towards that ending part of the transcript for + strand genes & the starting part of the transcript for - strand genes)
    #     exon_oi = sj_exons[1]
    #     bool_continue = True
    #     while( bool_continue ):
    #         next_sj = self.exon_to_sj_v2( exon_oi, simulant_sj, 1 )

    #         ##TEST:: print "ISO_SJ - BUILD TOWARDS HIGHER: exon_oi = ", exon_oi, " & next_sj = ", next_sj,  " & len = ", len( transcript_sj )

    #         if next_sj:     #if the next SJ is found, then record SJ and use the higher-positioned as the next exon of interest
    #             transcript_sj.append( next_sj[0] )
    #             exon_oi = next_sj[0].spliced_elems( self.isoform_id, False )[1]
    #         else:
    #             bool_continue = False

    #     ##TEST::print "IsoSJ BTASJ 3 - SHOULD HAVE FINISHED TRANSCRIPT - ", len( transcript_sj ), " of # da sites = ", len( da_sites )

    #     #sort list of SJs from lowest to highest position
    #     transcript_sj = sorted( transcript_sj, key = lambda s: s.start )

    #     return transcript_sj


    def reconstruct_transcript( self, canon_only = False, max_penalty = 0 ):
        """
        Args:
           max_penalty = integer that is the number of allowable non-canonical SJs in transcript
           canon_only = boolean that, if True, will only record canonical SJ for transcript. Else if False will record canonical & non-canonical SJs.
        Function: reconstruct the transcript for specific isoform, incorporating only a certain number penalties (this is a recursive function) 
        """
        
        """
        Method 1: check if there is enough support for a specific isoform - if so, then generate simulant SJ for missing splice junctions --> reconstruct transcript one SJ at a time, save different paths in an array
            
            -need to consider the 'next exon' for SJ depending on strand sign
            -need to find SJ that splice at the donor-acceptor sites as opposed to within exon (or within intron)
            -how am I finding a penalty? - 1) compare to each canonical isoform, 2) check the next exon number (not sure about this)
        
        Method 2 (PROBABLY NOT): see if enough support for isoform --> if yes, then reconstruct canonical form of isoform, retrieve aberrant SJ assigned to isoform, and then incorporate combinations of penalties depending on max_penalty score

        Method 3: determine if enough canonical support for isoform -> if true then create simulant SJ -> begin with first exon in isoform -> find SJ that contain that exon & record it -> save different paths in an array -> check the penalty score, if doesn't pass then discard
            -Potential Functions: connect_sj(), 
            -ANS to returning data structure: create a list of SJ objects! This will help retrieve information like frameshift, exon skip, penalties (i.e. # of non-canonical SJs)
            -ANS to finding exon & SJ: retrieve exon -> find the SJ associated with exon, record SJ -> find next exon -> find next SJ, record SJ -> etc. until no more SJs found
            -ANS to strand sign: start with correct side depending on strand sign - SJ needs to be able to go either direction
        """

        #TEST:: see if transcript is present
        isoform_present = self.is_isoform_present( 0.2, 3 )
        print "Reconstruct transcript, iso present? - ", isoform_present

        #check if enough support for isoform_id
        if not self.is_isoform_present( 0.5, 0 ):
            return []


        #if enough support, create simulant SJ for missing canonical SJ, else cease with transcript reconstruction
        # self.list_sj += self.sj_canon_absent( True )          #append list of simulant splice junctions to list of record SJs

        ##TEST:: print "IsoformSJ Step 1 - get absent SJ = ", len( self.sj_canon_absent( False ) )

        ##TEST:: 
        # print "*****Rebuild transcript: # of SJ = ", len( self.list_sj ), "*******\n"
        # print "ISJ.reconstruct_transcript: "
        # for i, sj in enumerate( self.list_sj ):
        #     print i, ": ", str( sj )

        #use this set of functions if do not want to incorporate aberrant SJs that splice intronically (into the introns)
        # t_start = self.build_start( canon_only )
        # # use self.build_transcript() if do not want to incorporate aberrant SJs that splice intronically (into the introns)
        # t_full = self.build_transcript( t_start, [], max_penalty, 0, canon_only )       #t_full = array of SJ instances that constitute an mRNA transcript
        
        #use this set of functions if do not want to incorporate aberrant SJs that splice intronically (into the introns)
        t_start = self.build_start_v2( self.list_sj, False )
        #use self.build_transcript_v2() if want to incorporate aberrant SJs that splice intronically (into the introns)
        t_full = self.build_transcript_v2( self.list_sj, t_start, [], max_penalty, 0, canon_only )       #t_full = array of SJ instances that constitute an mRNA transcript

        
        ##TEST:: observe the number of reconstructed transcripts
        # print "--TEST: Reconstruct Transcript-- --> len( t_full ) = ", len( t_full )
        # for i, each_transcript in enumerate( t_full ):
        #     for i2, each_sj in enumerate( each_transcript ):
        #         exons = each_sj.spliced_elems( self.isoform_id )
        #         print "RT  ", i, ":", i2, " : each_sj = ", each_sj, " & exons[0] = ", exons[0], " & exons[1] = ", exons[1]
        #         # print "RT - exonNum: ", i, ": ", str( each_sj ), " & start = ", exons[0].exonNum, " & end = ", exons[1].exonNum

        #just return all transcripts found
        return t_full

        # #if I want to return "valid" transcripts (transcripts with a start & end)
        # #check if valid transcript
        # t_valid = []
        # for t in t_full:        #t = array of SpliceJunctions
        #     if self.transcript_validity( t ):
        #         t_valid.append( t )

        # ##TEST:: print "reconstruct transcript: t_full = ", len( t_full ), " & t_valid = ", len( t_valid )

        # #returns array of arrays, where the inner array has SpliceJunction object that, once all connected, constitute an mRNA transcript
        # return t_valid


    def build_start( self, canon_only = False ):
        """
        Args:
           canon_only = boolean that, if True, will only record canonical SJ for transcript. Else if False will record canonical & non-canonical SJs.
        Function: find the first splice junction in the transcript
        """
        transcript_start = []
        #retrieve SJs that start with exon
        first_exon = self.get_exon_num( exon_base )     #from class Isoform, exons are 1-based, meaning the first exon is exon 1. If exons were 0-based, the the first exon would be exon 0

        #create a new array for each new start splice junction
        direction = 1 if self.strand == 1 else 2        #if + gene, then starting exon will have lowest numerical genomic position, else if - gene, then starting exon will have highest numerical genomic position
        first_sj = self.exon_to_sj( first_exon, direction )


        #if canon_only == True, then only record next splice junctions that are canonical
        if canon_only:
            first_sj = [x for x in first_sj if x.canon and self.isoform_id in x.assigned_isoform]

        for sj in first_sj:
            transcript_start.append( [sj] )

        return transcript_start


    def build_transcript( self, t_temp, t_full, max_penalty, counter, canon_only = False ):
        """
        Args:
           t_temp = array of SpliceJunction Objects that are still being 'built', where each element is a different path. This variable initially contains the start exons, and then are built from there. 
           t_full = array of ExonConnect objects that have been completed
           max_penalty = integer that is the number of allowable non-canonical SJs in transcript
           canon_only = boolean that, if True, will only record canonical SJ for transcript. Else if False will record canonical & non-canonical SJs.
        Function:
            returns all connections in self.hashExonConnect where each ExonConnect object can connect to another via the connectEnd binds to connectStart
        """
        ##TEST:: print counter, ":: IsoformSJ - build_transcript: t_temp = ", len( t_temp ), " & t_full = ", len( t_full )

        t_next = []
        for t in t_temp:        #t = array of Splice Junction objects that constitute a transcript
            #if transcript contains more aberrant splice junctions than max_penalty allows, then discard transcript
            if len( self.calculate_penalty( t ) ) > max_penalty:
                continue

            #find the next exon for the splice junction, find the next splice junction that starts with the next exon, and then append that splice junction
            sj_exons = t[-1].spliced_elems( self.isoform_id )       #returns tuple of Exon objects, where [0] = lower numerical exon & [1] = higher numerical exon
            if self.strand == 1:        #for + gene, find sj that start where the last SJ left off
                next_sj = self.exon_to_sj( sj_exons[1], 1 )
            else:                       #for - gene, look for sj that end where the last SJ began
                next_sj = self.exon_to_sj( sj_exons[0], 2 )

            #if canon_only == True, then only record next splice junctions that are canonical
            if canon_only:
                next_sj = [x for x in next_sj if x.canon and self.isoform_id in x.assigned_isoform]


            ##TEST:
            # for n in next_sj:
            #     print "each sj = ", n.sj_id
            # print "sj_exons: [0] = ", sj_exons[0], " & [1] = ", sj_exons[1]
            # if sj_exons[0] and sj_exons[1]:
            #     print "exon_num: [0] = ", sj_exons[0].exonNum, " & [1] = ", sj_exons[1].exonNum
            # print '<<<<<<>>>>>>\n\n'



            #do potential next splice junctions exists
            if next_sj:
                for sj in next_sj:      #sj = SpliceJunction object
                    if not sj in t:     #if splice junction 'sj' not already present in transcript, then record 
                        t_next.append( t + [sj] )
                    else:       #else if SJ is recorded, then record as completed transcript
                        if len( self.calculate_penalty( t ) ) <= max_penalty:
                            t_full.append( t )
            #else if next_sj exists but is a part of transcript 't' OR if next_sj does not exist
            else:
                if len( self.calculate_penalty( t ) ) <= max_penalty:
                    t_full.append( t )

        #if more unfinished transcripts exist, then continue building
        if t_next:
            self.build_transcript( t_next, t_full, max_penalty, counter + 1, canon_only )

        return t_full


    def build_start_v2( self, list_sj, canon_only = False ):
        """
        Args:
           canon_only = boolean that, if True, will only record canonical SJ for transcript. Else if False will record canonical & non-canonical SJs.
        Function: find the first splice junction in the transcript
        """
        #retrieve all exons from all SpliceJunctions Objects in 'list_sj'
        all_start_exons = []       #records all the "starting exons" ligated by the SpliceJunction (the exon with the lower numerical genomic position, regardless of exon number)
        all_end_exons = []          #records all the "ending exons" ligated by the SpliceJunction (the exon with the higher numerical genomic position, regardless of exon number)
        for each_sj in list_sj:
            #using SpliceJunction.spliced_elems() instead of SpliceJunction.spliced_elems_position_range() is because I need to see known feature elements (i.e. exons and introns)
            ligated_exons = each_sj.spliced_elems( self.isoform_id, False )
            all_start_exons.append( ligated_exons[0] )
            all_end_exons.append( ligated_exons[1] )

        #find the "first exon" - exon where there is no splice junction
        uniq_start_exons = {x3 for x3 in all_start_exons if not x3 in all_end_exons}       #the curly brackets make this a "set" comprehension, removing duplicates
        uniq_start_exons.add( self.get_exon_num(exon_base) )     #from class Isoform, exons are 1-based, meaning the first exon is exon 1. If exons were 0-based, the the first exon would be exon 0

        ##TEST::
        print "ISO_SJ.build_start_v2() - start exons:"
        for i, use in enumerate( uniq_start_exons ):
            print "start exon ", i, " - ", use

        #retrieve & record the starting SJs by recording all SJs that contain these unique start exons
        first_sj = []       #records all SJs that contain these unique start exons
        for start_exon in uniq_start_exons:
            #create a new array for each new start splice junction
            direction = 1 if self.strand == 1 else 2        #if + gene, then starting exon will have lowest numerical genomic position, else if - gene, then starting exon will have highest numerical genomic position
            
            #retrieve each start SJ & record it in first_sj
            get_start_sjs = self.exon_modified_to_sj_v2( start_exon, list_sj, direction )
            for each_start_sj in get_start_sjs:
                first_sj.append( each_start_sj )

        ##TEST::
        print "ISO_SJ.build_start_v2() - start SJ:"
        for i, each_fsj in enumerate( first_sj ):
            print "start SJ ", i, " - ", each_fsj

        #if canon_only == True, then only record next splice junctions that are canonical
        if canon_only:
            first_sj = [x for x in first_sj if x.canon and self.isoform_id in x.assigned_isoform]
        #remove duplicate elements
        first_sj = list( set(first_sj) )

        #create an array for each sj start, as each will be added to in def build_transcript_v2()
        transcript_start = []
        for sj in first_sj:
            transcript_start.append( [sj] )
            ##TEST::
            print "ISO_SJ start: add first sj = ", sj

        return transcript_start


    def build_transcript_v2( self, list_sj, t_temp, t_full, max_penalty, counter, canon_only = False ):
        """
        Args:
           list_sj = array of all SpliceJunction Objects that will be used to reconstruct the transcripot
           t_temp = array of SpliceJunction Objects that are still being 'built', where each element is a different path. This variable initially contains the start exons, and then are built from there. 
           t_full = array of ExonConnect objects that have been completed
           max_penalty = integer that is the number of allowable non-canonical SJs in transcript
           canon_only = boolean that, if True, will only record canonical SJ for transcript. Else if False will record canonical & non-canonical SJs.
        Function:
            returns all connections in self.hashExonConnect where each ExonConnect object can connect to another via the connectEnd binds to connectStart
        """
        ##TEST:: print counter, ":: IsoformSJ - build_transcript: t_temp = ", len( t_temp ), " & t_full = ", len( t_full )

        t_next = []
        for t in t_temp:        #t = array of Splice Junction objects that constitute a transcript

            #if transcript contains more aberrant splice junctions than max_penalty allows, then discard transcript
            if len( self.calculate_penalty( t ) ) > max_penalty:
                continue

            #find the next exon for the splice junction, find the next splice junction that starts with the next exon, and then append that splice junction
            sj_exons = t[-1].spliced_elems_position_range( self.isoform_id )       #returns tuple of Exon objects, where [0] = lower numerical exon & [1] = higher numerical exon. NOTE: I'm using SJ.spliced_elems_position_range() as I need to make sure both elements are exons. An SJ may splice into an intron, in which case it is "elongating" the next exon
            if self.strand == 1:        #for + gene, find sj that start where the last SJ left off
                next_sj = self.exon_modified_to_sj_v2( sj_exons['next_elem'], list_sj, 1 )
            else:                       #for - gene, look for sj that end where the last SJ began
                next_sj = self.exon_modified_to_sj_v2( sj_exons['prev_elem'], list_sj, 2 )

            #if canon_only == True, then only record next splice junctions that are canonical
            if canon_only:
                next_sj = [x for x in next_sj if x.canon and self.isoform_id in x.assigned_isoform]


            ##TEST:
            # for n in next_sj:
            #     print "each sj = ", n.sj_id
            # print "sj_exons: [0] = ", sj_exons[0], " & [1] = ", sj_exons[1]
            # if sj_exons[0] and sj_exons[1]:
            #     print "exon_num: [0] = ", sj_exons[0].exonNum, " & [1] = ", sj_exons[1].exonNum
            # print '<<<<<<>>>>>>\n\n'



            #do potential next splice junctions exists
            if next_sj:
                for sj in next_sj:      #sj = SpliceJunction object
                    if not sj in t:     #if splice junction 'sj' not already present in transcript, then record 
                        t_next.append( t + [sj] )
                    else:       #else if SJ is recorded, then record as completed transcript
                        if len( self.calculate_penalty( t ) ) <= max_penalty:
                            t_full.append( t )
            #else if next_sj exists but is a part of transcript 't' OR if next_sj does not exist
            else:
                if len( self.calculate_penalty( t ) ) <= max_penalty:
                    t_full.append( t )

        #if more unfinished transcripts exist, then continue building
        if t_next:
            self.build_transcript_v2( list_sj, t_next, t_full, max_penalty, counter + 1, canon_only )

        return t_full

    """
    Function: check validity of transcript created
    """

    def calculate_penalty( self, transcript_sj ):
        """
        Args:
            transcript_sj = array of SpliceJunction objects
        Function: returns an array of non-canonical splice junctions
        """
        # return [x for x in transcript_sj if not x in self.canon_transcript]
        return [x for x in transcript_sj if not x.canon and self.isoform_id in x.assigned_isoform]

    def transcript_validity( self, transcript_sj ):
        """
        Args:
            transcript_sj = array of SpliceJunction objects that constitutes an mRNA transcript
        Function: returns True if transcript_sj meets criteria to be a valid transcript, else returns False
            -checks if first & last exon present in transcript
        """
        #check if first & last exon present in transcript
        exon_start = self.get_exon_num( exon_base )
        exon_end = self.get_exon_num( self.last_exon_num )

        if self.strand == 1:
            has_start = transcript_sj[0].has_exon( exon_start, self.isoform_id, 1 )
            has_end = transcript_sj[-1].has_exon( exon_end, self.isoform_id, 2 )
        else:
            has_start = transcript_sj[0].has_exon( exon_start, self.isoform_id, 2 )
            has_end = transcript_sj[-1].has_exon( exon_end, self.isoform_id, 1 )

        ##TEST:: print "TRANSCRIPT VALID: start = ", exon_start.exonNum, " & end = ", exon_end.exonNum, " & has_start = ", has_start, " & has_end = ", has_end

        return True if has_start and has_end else False

            

    # def isoform_transcript_valid_v3(self, isoform_id, transcripts_full):
    #     """ 
    #     Args:
    #         isoform_id = string that is the isoform ID
    #         transcripts_full = array where each element is an array of ExonConnect Objects that constitutes a transcripts 
    #     Function: 
    #         checks if a transcript is a valid, meaning the transcript contains the start & end exon
    #     """
    #     #get start & end exon for isoform 'isoform_id'
    #     exon_start, exon_end = self.ExonList.isoform_exons_start_end(isoform_id)
    #     return [x for x in transcripts_full if x[0].match_start(exon_start) and x[-1].match_end(exon_end) ]

    # def calculate_penalty_tester(self, isoform_id, transcript):
    #     """
    #     Args:
    #         transcript = array of splice junctions (ExonConnect objects)
    #     Function: calculates penalty score for transcript by calculating the number of splice junctions in transcript that is also present in the canonical isoform 
    #     """
    #     #Version 1: calculate number of canonical splice junctions missing in 'transcript'
    #     # isoform_len = len(self.hashCanonTranscripts[isoform_id])
    #     # return isoform_len - len( [x for x in transcript if x in self.hashCanonTranscripts[isoform_id] ] )

    #     #Version 2: calculate number of non-canonical splice junctions in 'transcript'
    #     return len( [x for x in transcript if not x in self.hashCanonTranscripts[isoform_id] ] )

    # def isoform_reconstruct_v3_tester(self, isoform_id, isoform_sj, transcripts_temp, transcripts_full, max_penalty):
    #     """ 
    #     Args:
    #        transcripts_temp = array of ExonConnect objects that are still being 'built', where each element is a different path. This variable initially contains the start exons, and then are built from there. 
    #        transcripts_full = array of ExonConnect objects that have been completed
    #     Function:
    #         returns all connections in self.hashExonConnect where each ExonConnect object can connect to another via the connectEnd binds to connectStart
    #     Protocol:
    #         1: retrieve an array of arrays, where each element is an array of array
    #             -create GeneConnectMap object based on SJ, add SJs to GeneMap, which will add to ExonMap
    #             -ExonMap will contain ExonList I think - this.ExonList.hashExonList
    #         2: CAN I RETRIEVE THE FIRST EXONS OR TSS EXONS - ExonList.firstExon, or ExonList.tssExons, or retrieve objExon based on position
    #         3: go through each array element, and find the next connection using class ExonConnect
    #         4: find all connections & then go through function recursively
    #         5: once a connection is complete, then save to arrayFullPaths
    #         6: from each connection of arrays, recreate the path
    #     """
    #     #go through all_sj to find which connections should be added
    #     transcripts_next = []
    #     for i, each_transcript in enumerate(transcripts_temp):      #each_transcript = ExonConnectMap object
    #         #calculate number of penalties - if it is above the max penalty, then do not pursue
    #         if self.calculate_penalty(isoform_id, each_transcript) > max_penalty:
    #             continue

    #         #retrieve all connections where the start connection is the exon of interest
    #         next_connect = [ y for y in isoform_sj if y.match_start(each_transcript[-1].connectEnd) ]

    #         if len(next_connect):
    #             for each_ec in next_connect:
    #                 transcripts_next.append(each_transcript + [each_ec])
    #         else:
    #             #if transcript passes valid criteria, then add to array of fully reconstructed transcripts
    #             if self.isoform_transcript_valid_v3(isoform_id, [each_transcript]):
    #                 transcripts_full.append(each_transcript)

    #     if transcripts_next:
    #         self.isoform_reconstruct_v3_tester(isoform_id, isoform_sj, transcripts_next, transcripts_full, max_penalty)

    #     return transcripts_full
