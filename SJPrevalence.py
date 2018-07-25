#/usr/bin/python

from SpliceJunction import SpliceJunction

class SJPrevalence():
    def __init__( self ):
        pass

    """
    Function: determine read counts above a read count threshold
    """
    def check_above_thres( sample_names, sample_readcounts, thres ):
        """
        Retrieves the samples with read counts above threshold
        Args:
            sample_names = string that is list of sample names separated by ','. Assumes this order corresponds to order in 'sample_readcounts'
            sample_readcounts = string that is the list of read counts separated ','. Assumes this order corresponds to order in 'sample_names'
            thres = integer that is the read count threshold
        Returns: returns an array of arrays - read counts above threshold, samples with read counts above threshold, & controls with read counts above threshold
        """
        read_counts = [ int(x) for x in str(sample_readcounts).split(',') ]
        samples = str(sample_names).split(',')

        #retrieve indices above threshold value
        i_above_thres = [ i for i,x in enumerate( read_counts ) if x >= thres ]
        hash_sample_reads = {samples[i]:read_counts[i] for i in i_above_thres}
        reads_above_thres = [ read_counts[i] for i in i_above_thres ]
        samples_above_thres = [ samples[i] for i in i_above_thres ]
        ctrl_above_thres = [ c for c in LIST_CTRLS if c in samples_above_thres ]

        return [hash_sample_reads, reads_above_thres, samples_above_thres, ctrl_above_thres]

    """
    Function: methods to determine the dominance of each SJ
    """
    @staticmethod
    def display_list_sj_v2( list_sj ):
        """
        Args:
            -list_sj = array of SpliceJunction instances
        Function: displays all the SJs recorded in 'hash_sj'
        """
        display_found_sj = ""
        for get_sj in list_sj:
            display_found_sj += str( get_sj ) + " //// "

        return display_found_sj

    @staticmethod
    def compare_acceptor_bracket_sj_v3( obj_sj, iso_sj, canon_only = False ):
        """
        compares the aberrant SJ to the all the other SJs within the same 3' splice bracket (this is to find competing 5' splicing events)

        Args:
            -obj_sj = instance of SpliceJunction, where the position is the aberrant SJ
            -obj_iso_sj = instance of IsoformSJ -> this will be used for stat_overlap_sj = 2 or 3
            -canon_only = boolean that will be used in IsoformSJ def group_splice_donors_bracket_acceptors_v2(). These SJs will be used to calculate if SJ 'obj_sj' is dominant or not
                -True = only finds canonical splice events within the 3' splice site bracket (splice acceptor range)
                -False = finds all splice events (canon & non-canon) within the 3' splice site bracket (splice acceptor range)
        Returns:
            returns a hash including information on
                -ratio = calculate ratio between aberrant splicing read count & max overlapped canonical sk read count
                -list_bracket_readcount = array of all read counts support canonical splicing
                -display_canon_sj = string displayed by function def display_canon_sj()
                -bracket_range = string that is the genomic range of where 'obj_sj' is located
                    -this comes from property iso_sj.sj_bracket_five_prime --> k = SpliceJunction instance, v = tuple that is the range that contains the SpliceJunction
                -max_sj_info = information of max canonical SJ, including SJ ID, assigned isoforms
        """
        
        #retrieve the array of SpliceJunction instances in the same bracket
        # hash_sj_in_brackets = iso_sj.competing_sj_five_prime[obj_sj]
        get_bracket_range_oi = iso_sj.sj_bracket_five_prime[obj_sj]     #records the bracket range for the specific SpliceJunction instance. sj_bracket_five_prime --> k = SpliceJunction instance, v = tuple that is the range that contains the SpliceJunction
        list_overlap_sj = iso_sj.competing_sj_five_prime[get_bracket_range_oi]

        #get the range of the interest
        bracket_range = obj_sj.chrom + ':' + '-'.join( [str(x) for x in get_bracket_range_oi] )      #this should be only one range, the is the bracket range that encapsulates the range of interest

        if not list_overlap_sj:
            max_sj_info = str(obj_sj.read_count) + " - " + obj_sj.sj_id + " - " + ' | '.join( obj_sj.assigned_isoform )
            return  {'ratio': obj_sj.read_count, 'list_bracket_readcount': [1], 'display_canon_sj': obj_sj.sj_id, 'bracket_range': bracket_range, 'max_sj_info': max_sj_info}

        #retrieve the set of SJs depending on 'canon_only' (if 'canon_only' = True, then retrieve only canonial SJ, else if 'canon_only' = False, then retrieve all SJ)
        if canon_only:
            list_bracket_sj = [] if not list_overlap_sj else [x for x in list_overlap_sj if x.canon]
            list_bracket_readcount = [x.read_count for x in list_bracket_sj] if list_bracket_sj else [1]
        else:
            list_bracket_sj = [] if not list_overlap_sj else list_overlap_sj
            list_bracket_readcount = [x.read_count for x in list_bracket_sj] if list_bracket_sj else [1]

        
        ##TEST:: print "list_overlap_sj = ", list_overlap_sj

        #retrieve information from the max sj
        max_sj_info = 'nothing'
        if list_bracket_sj:
            get_max_sj = max( list_bracket_sj, key = lambda x: x.read_count )
            max_sj_info = str(get_max_sj.read_count) + " - " + get_max_sj.sj_id + " - " + ' | '.join( get_max_sj.assigned_isoform )


        #calculate the ratio between the aberrant SJ & the canonical SJ
        ratio = float( obj_sj.read_count ) / max( list_bracket_readcount )

        #display all the canonical SJs
        display_canon_sj = display_list_sj_v2( list_bracket_sj )

        # return [ratio, list_bracket_readcount, display_canon_sj, bracket_range]
        return  {'ratio': ratio, 'list_bracket_readcount': list_bracket_readcount, 'display_canon_sj': display_canon_sj, 'bracket_range': bracket_range, 'max_sj_info': max_sj_info}

    # @staticmethod
    # def create_iso_sj_sample( obj_sj, sample_name ):
    #     """
    #     creates an IsoformSJ instance to retrieve splice junctions for a specific sample (for tabix file)

    #     Returns:
    #         returns an instance of IsoformSJ for splice junctions specifically 
    #     """
    #     ##TEST::
    #     # print "obj_sj.assigned_isoform = ", obj_sj.assigned_isoform, " & gene_sym = ", obj_sj.gene_sym_oi, " & obj_sj.hash_isoforms.keys() = ", obj_sj.hash_isoforms.keys()

    #     #retrieve isoform ID
    #     if obj_sj.assigned_isoform:
    #         isoform_id = obj_sj.assigned_isoform[0]
    #     else:
    #         list_isoforms = [ x for x in obj_sj.hash_isoforms.keys() if 'NM' in x ]
    #         isoform_id = list_isoforms[0] if list_isoforms else obj_sj.hash_isoforms.keys()[0]

    #     path_sample_tabix = DIR_TABIX + "/sorted_" + sample_name + "_junctions.bed.gz" 
    #     isoform_info = obj_sj.hash_isoforms[ isoform_id ]
    #     #get range of isoform
    #     iso_chrom = str( isoform_info.chrom ) if not '_' in isoform_info.chrom else str( isoform_info.chrom.split('_')[0] )
    #     iso_start = str( isoform_info.boundary[0] )
    #     iso_end = str( isoform_info.boundary[1] )
    #     pos_range = iso_chrom + ':' + iso_start + '-' + iso_end
    #     #get the list of splicing events (e.g. list of SpliceJunction instances) for IsoformSJ class
    #     list_sj = IsoformSJ.get_sj_tabix_file( path_sample_tabix, pos_range, obj_sj.gene_sym_oi, 0 )
    #     list_sj+= [obj_sj]      #NOTE: there may be a duplicate of this SJ because of Tabix, but that's fine I think

    #     #create IsoformSJ instance
    #     sj_thres = 0
    #     return IsoformSJ( isoform_id, list_sj, sj_thres, None, False, 5 )             #new version

    @staticmethod
    def determine_sj_dominance( gene_sym, chrom, start, end, hash_sample_readcount ):
        """
        Determines the dominance of aberrant splicing by comparing the read count of the aberrant splicing event vs. overlapped canonical splicing events.

        Args:
            row = a hash that has the following information (or this could be a panda DataFrame):
                -'gene_sym': string that is the gene symbol (e.g. BRAF)
                -'chrom': string in the format 'chr#' (e.g chr9)
                -'start': integer that is the start genomic position for the SJ of interest. 'start' < 'end', regardless of strand sign
                -'end': integer that is the end genomic position for the SJ of interest. 'start' < 'end', regardless of strand sign
                -hash_sample_readcount = hash where key = sample name, value = read count associated with that sample
                -'samples_above_thres': an array of sample names that correspond to the order of 'read_count_above_thres'
                -'read_count_above_thres': an array of read counts that correspond to the order of 'samples_above_thres'
        Returns:
            returns an array where [0] = boolean if SJ is aberrant or not (False if canonical, else True if aberrant), and [1] = hash where key = sample name, value = dominance ratio (aberrant SJ read count / max canonical SJ)
        """
        #PROTOCOL: retrieve position of SJ -> determine if the SJ is aberrant -> if no, then skip, else if True then retrieve overlapped SJ for 
        #parameters that don't matter here for SpliceJunction class:
        sj_id = 'sj_id'
        strand = '+'
        isoform_id = None
        sample_prevalence = 0
        control_prevalence = 0
        bool_intronic = False

        
        ##TEST::
        # print "DSJ_DOM - samples = ", row['samples_above_thres'], " & ", row['read_count_above_thres']

        #go through all samples, & calculate the dominant ratio for each sample
        list_samples = row['samples_above_thres']
        list_readcounts = row['read_count_above_thres']
        hash_sample_ratios = {}     #hash where key = sample name, value = dominance ratio (aberrant SJ read count / max canonical SJ)
        ##TEST::
        # print "determine_sj_dom: list_samples = ", list_samples
        for k,v in hash_sample_readcount.iteritems():       #k = sample name, v = read count
            #retrieve sample name & corresponding read count
            sample_name = k
            sample_readcount = v

            #create SpliceJunction instance for sample
            obj_sj = SpliceJunction( sj_id, chrom, start, end, strand, sample_readcount, gene_sym, isoform_id, sample_prevalence, control_prevalence, bool_intronic )

            if obj_sj.canon:
                return [False, {}]       #[0] = aberrant SJ (False if canonical, else True if aberrant )

            #if no isoforms are recorded, then skip as I can reconstruct any transcript without an isoform
            if not obj_sj.hash_isoforms.keys():
                continue

            #calculate overlapping SJ
            iso_sj = create_iso_sj_sample( obj_sj, sample_name )

            # hash_overlap_canon = compare_overlapped_canon_sj_v2( obj_sj, iso_sj )
            hash_bracket_canon = compare_acceptor_bracket_sj_v3( obj_sj, iso_sj, True )     #find all canonical SJs that overlap the same range as 'obj_sj'
            # hash_bracket_all = compare_acceptor_bracket_sj_v3( obj_sj, iso_sj, False )       #find all SJs (canonical + aberrant) that overlap the same range as 'obj_sj'

            hash_sample_ratios[sample_name] = hash_bracket_canon['ratio']

        return [True, hash_sample_ratios]
