#/usr/bin/python
#Alg: SuiteNeoepStatistics.py

import pandas as pd


# #Version 1: for IEDB MHC Only combined with Proteasome + TAP + MHC
# g_aff_neoep_orig = 'aff_cat_x_1'       #affinity of originaneoepitope to MHC --> I'M NOT SURE, NEED TO DOUBLE-CHECK COLUMN
# g_aff_neoep_alt = 'aff_cat_y_2'       #affinity of altered neoepitope to MHC

#Version 2: for IEDB Proteasome + TAP + MHC
g_aff_neoep_orig = 'aff_cat_1'       #affinity of originaneoepitope to MHC --> YES I DOUBLE-CHECK, THIS IS CORRECT
g_aff_neoep_alt = 'aff_cat_2'       #affinity of altered neoepitope to MHC
g_alt_peptide = 'peptide_2'

df_empty = pd.DataFrame( {'A' : []} )       #just use this just in case I need to send an empty dataframe

g_sample_id = 'case_id'     #this is for anti-PD1 files
# g_sample_id = 'case_id'     #this is for Veliparib files

class NeoepStatistics():

    def __init__( self, df_neoep, hash_filters ):
        """
        Args:
            -hash_filters = a hashtable that contains the following keys & values. NOTE: the filters for proteasome cleavage, TAP transport, & endogenous
                -frame_stat = integer with the following status: (NOTE: Default = 0)
                    -0 = will only consider both the alteration type (mutation, indel) (column = 'my_change_type')
                    -1 = will only consider the reading frame (column = 'in_frame')
                    -2 = will consider both the alteration type (mutation, indel) & frame-preservation (column = 'alt_readframe_type')
                -neoep_stat = integer that will perform a specific type of analysis with neoepitopes. Default = 0
                    -0 = will consider ALL genomic alterations that produce any neoepitope
                    -1 = will only consider high-affinity MHC neoepitopes that result from genomic alterations
                    -2 = will only consider high-efficacious neoepitopes that result from genomic alterations
                -check_NMD = boolean that considers Nonsense-mediate Decay (NMD) -> default is False
                    -True = will only count events (genomic alterations) where NMD should not occur ()
                    -False = will consider all events regardless if it is susceptible to NMD or not
                -thres_percent_exp = number (integer or float) that serves as gene expression threshold, but only if the expression level is greater than -1. Default = -1
                -thres_percent_prot = number (integer or float) that serves as the proteasome percentile threshold, but only if the expression level is greater than -1. Default = -1
                -thres_percent_tap = number (integer or float) that serves as the TAP percentile threshold, but only if the expression level is greater than -1. Default = -1
                -thres_endogenous_freq = the limit of the number times the peptide is found as an endogenous peptide -> I will usually use 0 for this to find peptides that are unique. Default = -1
        """
        #the filters used on dataset
        self.mhc_affinity = hash_filters['mhc_affinity']
        self.frame_stat = hash_filters['frame_stat']
        self.neoep_stat = hash_filters['neoep_stat']
        self.check_NMD = hash_filters['check_NMD']
        self.thres_percent_exp = hash_filters['thres_percent_exp']
        self.thres_percent_prot = hash_filters['thres_percent_prot']
        self.thres_percent_tap = hash_filters['thres_percent_tap']
        self.thres_endogenous_freq = hash_filters['thres_endogenous_freq']


        # alt_column = 'alt_readframe_type' if with_rf else 'my_change_type'
        if self.frame_stat == 0:     #only look at genomic alterations (SNV, insertions, deletions)
            alt_column = 'my_change_type'
        elif self.frame_stat == 1:       #only look at in-frame or frameshift events
            alt_column = 'in_frame'
        elif self.frame_stat == 2:       #look at the combination of genomic alteration (SNV, insertions, deletions) + reading frame status (in-frame or frameshift)
            alt_column = 'alt_readframe_type'

        #information about dataset
        self.df_neoep = df_neoep
        self.list_samples = self.df_neoep[g_sample_id].unique()
        self.list_alt_types = self.df_neoep[alt_column].unique()        #all categories of alterations (e.g. SNV, insertion, deletion, SNV_inframe, insertion_frameshift)

    @staticmethod
    def write_df_to_file( df_oi, df_title, file_summary ):
        """
        writes the output of a dataframe to the file "file_summary"
        Args:
            -df_oi = dataframe that will be written to the file "file_summary"
            -df_title = string that is the title describing "df_oi"
            -file_summary = the file that will be written to, will write to the file by using "file_summary.name"
        """
        df_blank = pd.DataFrame()

        df_blank.to_csv( file_summary.name, sep = '\t', mode = 'a', index_label = df_title )
        df_oi.to_csv( file_summary.name, mode = 'a', header = True, index = True, sep = '\t' )
        df_blank.to_csv( file_summary.name, sep = '\t', mode = 'a' )


    @staticmethod
    def create_hash_filter( mhc_affinity, frame_stat, neoep_stat, check_NMD, thres_gene_exp, thres_prot, thres_tap, thres_endogenous_freq ):
        """
        Creates a hash filter for the function quantify_all_alterations_v2()
        Args:
            -mhc_affinity = integer that filters based on MHC affinity, here is what each integer means:
                -0 = consider all neoepitopes, no need to select neoepitopes with specific affinity to MHC allele
                -1 = consider all neoepitopes with at least "LOW" affinity to MHC allele (so consider LOW, MID, & HIGH)
                -2 = consider all neoepitopes with at least "MID" affinity to MHC allele (so consider MID, & HIGH)
                -3 = consider all neoepitopes with at least "HIGH" affinity to MHC allele (so consider HIGH only)
            -frame_stat = integer with the following status: (NOTE: Default = 0)
                -0 = will only consider both the alteration type (mutation, indel) (column = 'my_change_type')
                -1 = will only consider the reading frame (column = 'in_frame')
                -2 = will consider both the alteration type (mutation, indel) & frame-preservation (column = 'alt_readframe_type')
            -neoep_stat = integer that will perform a specific type of analysis with neoepitopes. Default = 0
                -0 = will consider ALL genomic alterations that produce any neoepitope
                -1 = will only consider high/mid-affinity MHC neoepitopes that result from genomic alterations
                -2 = will only consider high-efficacious neoepitopes (high/mid-affinity neoepitopes) that result from genomic alterations
            -check_NMD = boolean that considers Nonsense-mediate Decay (NMD) -> default is False
                -True = will only count events (genomic alterations) where NMD should not occur ()
                -False = will consider all events regardless if it is susceptible to NMD or not
            -thres_gene_exp = number (integer or float) that serves as gene expression threshold, but only if the expression level is greater than -1. Default = -1
            -thres_prot = number (integer or float) that serves as the proteasome percentile threshold, but only if the expression level is greater than -1. Default = -1
            -thres_tap = number (integer or float) that serves as the TAP percentile threshold, but only if the expression level is greater than -1. Default = -1
            -thres_endogenous_freq = the limit of the number times the peptide is found as an endogenous peptide -> I will usually use 0 for this to find peptides that are unique. Default = -1
        """
        hash_filters = {
            "mhc_affinity": mhc_affinity,
            "frame_stat": frame_stat,
            "neoep_stat": neoep_stat,
            "check_NMD": check_NMD,
            "thres_percent_exp": thres_gene_exp,
            "thres_percent_prot": thres_prot,
            "thres_percent_tap": thres_tap,
            "thres_endogenous_freq": thres_endogenous_freq
        }

        return hash_filters

    @staticmethod
    def filter_correspond_labels( hash_filters ):
        """
        Creates a text label based on the filters being used in dictionary "hash_filters" (hash_filters created by def create_hash_filter())
        Args:
            -hash_filters = hash of filters that will be applied to dataframe of neoepitopes - "hash_filters" is generated by "create_hash_filter()"
        """
        list_label = []
        #look at neoepitope affinity to MHC allele
        if ( hash_filters['mhc_affinity'] == 1 ):
            list_label.append( "MHC affinity >= LOW" )
        elif ( hash_filters['mhc_affinity'] == 2 ):
            list_label.append( "MHC affinity >= MID" )
        elif ( hash_filters['mhc_affinity'] == 3 ):
            list_label.append( "MHC affinity = LOW Only" )
        elif ( hash_filters['mhc_affinity'] == 4 ):
            list_label.append( "MHC affinity = MID Only" )
        elif ( hash_filters['mhc_affinity'] == 5 ):
            list_label.append( "MHC affinity = HIGH ONLY" )


        #look at frame_stat
        if ( hash_filters['frame_stat'] == 0 ):
            list_label.append( "frame_stat = alteration only" )
        elif ( hash_filters['frame_stat'] == 1 ):
            list_label.append( "frame_stat = read-frame only" )
        elif ( hash_filters['frame_stat'] == 2 ):
            list_label.append( "frame_stat = alteration + read-frame" )

        #look at neoep_stat
        if ( hash_filters['neoep_stat'] == 0 ):
            list_label.append( "neoep_stat = all neoeps" )
        elif ( hash_filters['neoep_stat'] == 1 ):
            list_label.append( "neoep_stat = high-affinity neoeps" )
        elif ( hash_filters['neoep_stat'] == 2 ):
            list_label.append( "neoep_stat = high-efficacy neoeps" )

        #look at check_NMD
        list_label.append( "Filter NMD = " + str( hash_filters['check_NMD'] ) )
        list_label.append( "Gene Exp >= " + str( hash_filters['thres_percent_exp'] ) )
        list_label.append( "ProtCleav Thres >= " + str( hash_filters['thres_percent_prot'] ) )
        list_label.append( "TAP Thres >= " + str( hash_filters['thres_percent_tap'] ) )
        list_label.append( "Endogen Freq Thres <= " + str( hash_filters['thres_endogenous_freq'] ) )

        return (' | ').join( list_label )




    def filter_neoeps( self, df_neoep_orig = None ):
        """
        Filters rows from dataframe of neoepitopes "self.df_neoep"
        -NOTE: note that "self.mhc_affinity" & "self.neoep_stat" may clash as they both filter for MHC-neoepitope affinity
        # """
        # frame_stat = hash_filters['frame_stat']
        # neoep_stat = hash_filters['neoep_stat']
        # check_NMD = hash_filters['check_NMD']
        # thres_percent_exp = hash_filters['thres_percent_exp']
        # thres_percent_prot = hash_filters['thres_percent_prot']
        # thres_percent_tap = hash_filters['thres_percent_tap']
        # thres_endogenous_freq = hash_filters['thres_endogenous_freq']

        #make deep copy of original dataframe
        if df_neoep_orig is None:
            df_neoep = self.df_neoep.copy( deep = True )        #make a copy of this as to not alter the original dataframe - deep parameter's default is True, just wrote this to make sure it is a deep copy anyways
        else:
            df_neoep = df_neoep_orig.copy( deep = True )        #make a copy of this as to not alter the original dataframe - deep parameter's default is True, just wrote this to make sure it is a deep copy anyways

        #filter based on MHC affinity
        if self.mhc_affinity == 1:
            df_neoep = df_neoep[ df_neoep[g_aff_neoep_alt].isin( ["HIGH", "MID", "LOW"] ) ]
        elif self.mhc_affinity == 2:
            df_neoep = df_neoep[ df_neoep[g_aff_neoep_alt].isin( ["HIGH", "MID"] ) ]
        elif self.mhc_affinity == 3:        #"LOW" only
            df_neoep = df_neoep[ df_neoep[g_aff_neoep_alt] == "LOW" ]
        elif self.mhc_affinity == 4:        #"MID" only
            df_neoep = df_neoep[ df_neoep[g_aff_neoep_alt] == "MID" ]
        elif self.mhc_affinity == 5:        #"HIGH" only
            df_neoep = df_neoep[ df_neoep[g_aff_neoep_alt] == "HIGH" ]

        # alt_column = 'alt_readframe_type' if with_rf else 'my_change_type'
        if self.frame_stat == 0:     #only look at genomic alterations (SNV, insertions, deletions)
            alt_column = 'my_change_type'
        elif self.frame_stat == 1:       #only look at in-frame or frameshift events
            alt_column = 'in_frame'
        elif self.frame_stat == 2:       #look at the combination of genomic alteration (SNV, insertions, deletions) + reading frame status (in-frame or frameshift)
            alt_column = 'alt_readframe_type'

        #select neoepitope based on affinity to MHC and relative performance (MHC affinity, proteasome cleavage, TAP transport) to the orignial, unaltered peptide
        if self.neoep_stat == 1:
            #quantify the number of unique alterations that produce high-affinity neoepitopes
            df_neoep = df_neoep[ df_neoep[g_aff_neoep_alt].isin( ["HIGH", "MID"] ) ]
        elif self.neoep_stat == 2:
            #METHOD 2: use (df_neoep[g_aff_neoep_alt].isin["HIGH", "MID"]) to find high to intermediate neoepitopes
            df_neoep = df_neoep[ 
            (df_neoep[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) & 
            (df_neoep["delta_mhc_score"].isin( ['higher', 'equal'] )) & 
            (df_neoep["delta_proteasome_score"].isin( ['higher', 'equal'] )) & 
            (df_neoep["delta_TAP_score"].isin( ['higher', 'equal'] )) ]


        #first retrieve the unique MHC alleles & genomic alteration types
        list_samples = df_neoep[g_sample_id].unique()
        list_alt_types = df_neoep[alt_column].unique()

        #determine if I need to remove events that are susceptible to NMD
        if self.check_NMD:
            df_neoep = df_neoep[ df_neoep['bool_NMD'].str.contains('FA', case = False) ]        #if bool_NMD = False, then this means the 

        #Gene Expression Percentile Threshold: apply threshold if threshold parameter > -1
        if self.thres_percent_exp > -1:
            df_neoep = df_neoep[ df_neoep['express_percentile'] >= self.thres_percent_exp ]

        #Proteasome Cleavage Percentile Threshold: apply threshold if threshold parameter > -1
        if self.thres_percent_prot > -1:
            df_neoep = df_neoep[ df_neoep['proteasome_percentile_2'] >= self.thres_percent_prot ]

        #TAP Transport Percentile Threshold: apply threshold if threshold parameter > -1
        if self.thres_percent_tap > -1:
            df_neoep = df_neoep[ df_neoep['tap_percentile_2'] >= self.thres_percent_tap ]

        #Threshold for Endogenous Frequency: apply threshold if threshold parameter > -1 (I will usually use 0 for this to find peptides that are unique)
        if self.thres_endogenous_freq > -1:
            df_neoep = df_neoep[ df_neoep['pep_matches_2'] <= self.thres_endogenous_freq ]

        return df_neoep


    ##DO I NEED THIS FUNCTION?? - I think this will be replaced by quant_neoeps_samples()
    # def quantify_num_neoepitopes( self ):
    #     """
    #     Quantifies the number of neoepitopes based on the filters "hash_filters"
    #     Protocol:
    #         -filter neoepitopes using self.filter_neoeps()
    #         -quantify the number of neoepitopes per sample
    #     """
    #     df_neoep_filt = self.filter_neoeps()

    #     #hash_all_alts will be used to quantify each type of alteration for each sample
    #     hash_all_alts = {k:[] for k in self.list_alt_types}

    #     #create a hash for each sample
    #     case_alt_count = {}      #reads genomic alteration frequency per sample, k = sample_barcode, value = hash where k2 = specific genomic alteration type & v2 = frequency of that alteration in sample_barcode
    #     for each_sample in list_samples:
    #         case_alt_count[each_sample] = dict( hash_all_alts )     #create hash for sample that will record alterations associated with sample

    #         for each_alt in self.list_alt_types:
    #             #quantify all types of alterations for each sample
    #             if neoep_stat == 0:         #quantify all alterations, whether they lead to high-affinity neoepitopes or not
    #                 list_alts = df_neoep[ (df_neoep[g_sample_id] == each_sample) & (df_neoep[alt_column].str.contains(each_alt, case = False) ) ]['genome_pos'].unique()
    #             elif neoep_stat == 1:
    #                 #quantify the number of unique alterations that produce high-affinity neoepitopes
    #                 list_alts = df_neoep[ (df_neoep[g_sample_id] == each_sample) & (df_neoep[alt_column].str.contains(each_alt, case = False) ) & ( df_neoep[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) ]['genome_pos'].unique()
    #             elif neoep_stat == 2:
    #                 #METHOD 2: use (df_neoep[g_aff_neoep_alt].isin["HIGH", "MID"]) to find high to intermediate neoepitopes
    #                 list_alts = df_neoep[ 
    #                 (df_neoep[g_sample_id] == each_sample) & 
    #                 (df_neoep[alt_column].str.contains( each_alt, case = False )) & 
    #                 (df_neoep[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) & 
    #                 (df_neoep["delta_mhc_score"].isin( ['higher', 'equal'] )) & 
    #                 (df_neoep["delta_proteasome_score"].isin( ['higher', 'equal'] )) & 
    #                 (df_neoep["delta_TAP_score"].isin( ['higher', 'equal'] )) ]['genome_pos'].unique()



    def quant_neoeps_samples( self, df_neoep_orig = None, neoep_stat = None, bool_uniq_neoep = True ):
        """
        Quantifies the number of unique neoepitopes produced per patient based on which types of neoepitopes are are looking for (based on neoep_stat -> All neoepitopes, High-affinity neoepitopes, or High-efficacy neoepitopes)
        Args:
            -df_neoep = dataframe that contains neoepitopes of interest (either filtered or not filtered)
            -neoep_stat = integer that will perform a specific type of analysis with neoepitopes. Default = 0
                -0 = will consider ALL genomic alterations that produce any neoepitope
                -1 = will only consider high-affinity MHC neoepitopes that result from genomic alterations
                -2 = will only consider high-efficacious neoepitopes that result from genomic alterations
            -bool_uniq_neoep = boolean where,
                -True = if True, will only count peptides with unique sequence
                -False = if False, will count all peptides, regardless if unique sequence or not
        Returns:
            -returns a hash where k = sample name, v = count of # of neoepitopes found
        """
        if df_neoep_orig is None:
            # df_neoep = self.df_neoep.copy( deep = True )
            df_neoep = self.filter_neoeps()
        if not neoep_stat:
            neoep_stat = self.neoep_stat

        #hash_all_alts will be used to quantify each type of alteration for each sample
        hash_all_alts = {k:[] for k in self.list_alt_types}

        #create a hash for each sample
        case_alt_count = {}      #reads genomic alteration frequency per sample, k = sample_barcode, value = hash where k2 = specific genomic alteration type & v2 = frequency of that alteration in sample_barcode
        for each_sample in self.list_samples:
            case_alt_count[each_sample] = dict( hash_all_alts )     #create hash for sample that will record alterations associated with sample

            for each_alt in self.list_alt_types:
                #quantify all types of alterations for each sample
                if neoep_stat == 0:         #quantify all alterations, whether they lead to high-affinity neoepitopes or not
                    df_neoep_type = df_neoep[ (df_neoep[g_sample_id] == each_sample) & (df_neoep[alt_column].str.contains(each_alt, case = False) ) ]
                elif neoep_stat == 1:
                    #quantify the number of unique alterations that produce high-affinity neoepitopes
                    case_alt_count[each_sample][each_alt] = df_neoep[ (df_neoep[g_sample_id] == each_sample) & (df_neoep[alt_column].str.contains(each_alt, case = False) ) & ( df_neoep[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) ]
                elif neoep_stat == 2:
                    #METHOD 2: use (df_neoep[g_aff_neoep_alt].isin["HIGH", "MID"]) to find high to intermediate neoepitopes
                    df_neoep_type = df_neoep[ 
                    (df_neoep[g_sample_id] == each_sample) & 
                    (df_neoep[alt_column].str.contains( each_alt, case = False )) & 
                    (df_neoep[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) & 
                    (df_neoep["delta_mhc_score"].isin( ['higher', 'equal'] )) & 
                    (df_neoep["delta_proteasome_score"].isin( ['higher', 'equal'] )) & 
                    (df_neoep["delta_TAP_score"].isin( ['higher', 'equal'] )) ]


                #calculate the number of neoepitopes 
                if bool_uniq_neoep:
                    case_alt_count[each_sample][each_alt] = df_neoep_type['peptide_2'].unique().count()
                else:
                    case_alt_count[each_sample][each_alt] = df_neoep_type['peptide_2'].count()

        #convert into a dataframe (& potential write it to a file)
        df_neoep_count = pd.DataFrame( case_alt_count.values(), index = case_alt_count.keys(), columns = case_alt_count.values()[0].keys() )

        return df_neoep_count


    def retrieve_mutation_count_all_samples( self ):
        """
        Similar to def retrieve_mutation_count(), but does this for all samples
        """
        case_alt_count = {}      #reads genomic alteration frequency per sample, k = sample_barcode, value = hash where k2 = specific genomic alteration type & v2 = frequency of that alteration in sample_barcode

        df_neoep = self.filter_neoeps()
        for each_sample in self.list_samples:
            case_alt_count[each_sample] = self.retrieve_mutation_count( each_sample, df_neoep, self.neoep_stat )

        #create dataframe and write dataframe to file
        df_alt_count = pd.DataFrame( case_alt_count.values(), index = case_alt_count.keys(), columns = case_alt_count.values()[0].keys() )

        return df_alt_count

    def create_chart_label( self ):
        #create label for title based on alteration, reading frame, or both
        if self.frame_stat == 0:
            label_alt_type = "Alterations Only"
        elif self.frame_stat == 1:
            label_alt_type = "Reading-Frame Only"
        elif self.frame_stat == 2:
            label_alt_type = "Alterations + Reading-Frame"

        #gene expression threshold
        label_thres = "Gene Exp Threshold " + str( self.thres_percent_exp ) if thres_percent_exp > -1 else "No Threshold Applied"
        #NMD filter label
        label_nmd = "NMD removed" if self.check_NMD else "May have NMD"
        label_prot = "ProtCleav Thres: " + str(self.thres_percent_prot) if thres_percent_prot > -1 else "NO ProtCleav"
        label_tap = "TAP Transport Thres: " + str(self.thres_percent_tap) if thres_percent_tap > -1 else "NO TAP"
        label_endog = "Endogenous Pep Freq Thres: " + str(self.thres_endogenous_freq) if thres_percent_tap > -1 else "NO Endogen Pep"

        #combine all the labels
        append_labels = label_alt_type + ", " + label_thres + ", " + label_nmd + ", " + label_prot + ", " + label_tap + ", " + label_endog

        #create title based on type of neoepitope (all, high-affinity, high-efficacious)
        if self.neoep_stat == 0:
            title = "Total Number Genomic " + label_alt_type + ": " + append_labels
        elif self.neoep_stat == 1:
            title = "Number Genomic " + label_alt_type + " Producing High-Affinity Neoepitopes: " + append_labels
        elif self.neoep_stat == 2:
            title = "Number Genomic " + label_alt_type + " Producing High-Efficacious Neoepitopes: " + append_labels

        return title



    def retrieve_mutation_count( self, sample_name, df_neoep_filt = None, neoep_stat = None ):
        """
        Quantifies the number of each type of alteration per sample
        """
        if df_neoep_filt is None:
            # df_neoep_filt = self.df_neoep.copy( deep = True )
            df_neoep_filt = self.filter_neoeps()
        if not neoep_stat:
            neoep_stat = self.neoep_stat


        #go through each sample & alteration and quantify their frequency
        alt_count = {}      #reads genomic alteration frequency per sample, k = sample_barcode, value = hash where k2 = specific genomic alteration type & v2 = frequency of that alteration in sample_barcode
        for each_alt in self.list_alt_types:
            #quantify all types of alterations for each sample
            if neoep_stat == 0:         #quantify all alterations, whether they lead to high-affinity neoepitopes or not
                list_alts = df_neoep_filt[ (df_neoep_filt[g_sample_id] == sample_name) & (df_neoep_filt[alt_column].str.contains(each_alt, case = False) ) ]['genome_pos'].unique()
            elif neoep_stat == 1:
                #quantify the number of unique alterations that produce high-affinity neoepitopes
                list_alts = df_neoep_filt[ (df_neoep_filt[g_sample_id] == sample_name) & (df_neoep_filt[alt_column].str.contains(each_alt, case = False) ) & ( df_neoep_filt[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) ]['genome_pos'].unique()
            elif neoep_stat == 2:       #look for high-efficacious neoepitopes (proteasome, TAP, MHC)
                #METHOD 2: use (df_neoep_filt[g_aff_neoep_alt].isin["HIGH", "MID"]) to find high to intermediate neoepitopes
                list_alts = df_neoep_filt[ 
                (df_neoep_filt[g_sample_id] == sample_name) & 
                (df_neoep_filt[alt_column].str.contains( each_alt, case = False )) & 
                (df_neoep_filt[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) & 
                (df_neoep_filt["delta_mhc_score"].isin( ['higher', 'equal'] )) & 
                (df_neoep_filt["delta_proteasome_score"].isin( ['higher', 'equal'] )) & 
                (df_neoep_filt["delta_TAP_score"].isin( ['higher', 'equal'] )) ]['genome_pos'].unique()

            #record the number of alterations for each sample
            alt_count[each_alt] = len( list_alts )

        return alt_count


    def retrieve_list_mutations( self, sample_name, df_neoep_filt = None, neoep_stat = None ):
        """
        determine all the mutations that occur in each sample
        Args:
            -sample_name = string that is sample name in the column that contains the sample name (e.g. "case_id")
            -df_neoep_filt = pandas dataframe that contains neoepitopes but perhaps has already been filtered
            -neoep_stat = integer that
        """
        if df_neoep_filt is None:
            # df_neoep_filt = self.df_neoep.copy( deep = True )
            df_neoep_filt = self.filter_neoeps()
        if not neoep_stat:
            neoep_stat = self.neoep_stat

        for each_alt in self.list_alt_types:
            #quantify all types of alterations for each sample
            if neoep_stat == 0:         #quantify all alterations, whether they lead to high-affinity neoepitopes or not
                list_alts = df_neoep_filt[ (df_neoep_filt[g_sample_id] == sample_name) & (df_neoep_filt[alt_column].str.contains(each_alt, case = False) ) ]['genome_pos'].unique()
            elif neoep_stat == 1:
                #quantify the number of unique alterations that produce high-affinity neoepitopes
                list_alts = df_neoep_filt[ (df_neoep_filt[g_sample_id] == sample_name) & (df_neoep_filt[alt_column].str.contains(each_alt, case = False) ) & ( df_neoep_filt[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) ]['genome_pos'].unique()
            elif neoep_stat == 2:       #look for high-efficacious neoepitopes (proteasome, TAP, MHC)
                #METHOD 2: use (df_neoep_filt[g_aff_neoep_alt].isin["HIGH", "MID"]) to find high to intermediate neoepitopes
                list_alts = df_neoep_filt[ 
                (df_neoep_filt[g_sample_id] == sample_name) & 
                (df_neoep_filt[alt_column].str.contains( each_alt, case = False )) & 
                (df_neoep_filt[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) & 
                (df_neoep_filt["delta_mhc_score"].isin( ['higher', 'equal'] )) & 
                (df_neoep_filt["delta_proteasome_score"].isin( ['higher', 'equal'] )) & 
                (df_neoep_filt["delta_TAP_score"].isin( ['higher', 'equal'] )) ]['genome_pos'].unique()

        return list_alts


    def count_neoepitope_per_alt( self, sample_name, df_neoep_filt = None, neoep_stat = None ):
        """
        Quantifies the number of neoepitopes per alteration
        """
        if not df_neoep_filt:
            # df_neoep_filt = self.df_neoep.copy( deep = True )
            df_neoep_filt = self.filter_neoeps()
        if not neoep_stat:
            neoep_stat = self.neoep_stat

        #retrieve the list of genomic mutations
        list_alts = self.retrieve_list_mutations( sample_name, df_neoep_filt = None, neoep_stat = None )


        #go through each sample and quantify the number of neoepitopes associated with genomic alteration position
        hash_sample_pep_freq = {}       #key = sample ID, value = array where each element is a unique peptide sequence
        if neoep_stat == 0:
            pep_freq = df_genome_pos[ df_genome_pos[g_sample_id] == each_sample ][g_alt_peptide].unique()
        elif neoep_stat == 1:
            pep_freq = df_genome_pos[ (df_genome_pos[g_sample_id] == each_sample) &
            (df_genome_pos[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) ][g_alt_peptide].unique()
        elif neoep_stat == 2:
            #METHOD 2: use (df_genome_pos[g_aff_neoep_alt].isin["HIGH", "MID"]) to find high to intermediate neoepitopes
            pep_freq = df_genome_pos[ 
            (df_genome_pos[g_sample_id] == each_sample) & 
            (df_genome_pos[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) & 
            (df_genome_pos["delta_mhc_score"].isin( ['higher', 'equal'] )) & 
            (df_genome_pos["delta_proteasome_score"].isin( ['higher', 'equal'] )) & 
            (df_genome_pos["delta_TAP_score"].isin( ['higher', 'equal'] )) ][g_alt_peptide].unique()

        #number of neoepitopes associated with this position at this sample
        return len( pep_freq )

    def avg_neoeps_sample_alt():
        """
        Calculates statistics for the number of neoepitopes produced per sample
        """
        #calculate mean & standard deviation for each sample
        case_alt_neoep_statistics = {}      #k = sample ID, v = hash where k2 = alteration type:mean/stdev & v2 = corresponding value
        for k2, v2 in case_alt_count.iteritems():        #k2 = sample ID, v2 = hash where k2_a = genomic alteration type based on 'frame_stat' & v2_a = array of integer that is the number of unique neoepitopes associated with all genomic alteration position for sample ID
            if not k2 in case_alt_neoep_statistics:
                case_alt_neoep_statistics[k2] = {}
            for k2_a, v2_a in v2.iteritems():

                ##TEST:: print "SUITE.Q_NUM_NEOEP: v2_a = ", v2_a

                list_neoep_len = [get_len for get_len in v2_a if get_len > 0]
                case_alt_neoep_statistics[k2][k2_a + ":avg"] = np.mean( list_neoep_len )
                case_alt_neoep_statistics[k2][k2_a + ":std"] = np.std( list_neoep_len )

    ##NEW FUNCTION - NEED TO TEST THIS FUNCTION
    def count_neoepitope_per_alt_v2( df_neoep, alt_genome_pos, hash_filters ):
        """
        Counts the number of neoepitopes associated with specific genomic position
        Args:
            -df_neoep = dataframe that is the file containing all original & altered neoepitopes
            -alt_genome_pos = string that is the genomic position
            -neoep_stat = integer that will perform a specific type of analysis with neoepitopes
                -0 = will consider ALL genomic alterations that produce any neoepitope
                -1 = will only consider high-affinity MHC neoepitopes that result from genomic alterations
                -2 = will only consider high-efficacious neoepitopes that result from genomic alterations
            -check_NMD = boolean that considers Nonsense-mediate Decay (NMD)
                -True = will only count events (genomic alterations) where NMD should not occur ()
                -False = will consider all events regardless if it is susceptible to NMD or not
        CAUTION: as a genomic alteration can be associated with multiple isoforms, need to make sure to count
        NOTE:
            -make sure the "bool_NMD" column is formatted as a string type (i.e. df[['bool_NMD']] = df[['bool_NMD']].astype(str) )
        """
        frame_stat = hash_filters['frame_stat']
        neoep_stat = hash_filters['neoep_stat']
        check_NMD = hash_filters['check_NMD']
        thres_percent_exp = hash_filters['thres_percent_exp']
        thres_percent_prot = hash_filters['thres_percent_prot']
        thres_percent_tap = hash_filters['thres_percent_tap']
        thres_endogenous_freq = hash_filters['thres_endogenous_freq']

        #first retrieve the unique MHC alleles & genomic alteration types - IMPORTANT TO DO THIS FIRST SO I CAN WHICH SAMPLES HAVE THIS MUTATION AND WHICH DO NOT
        list_samples = df_neoep[g_sample_id].unique()

        #retrieve all alterations with the same genomic positions across all samples
        df_genome_pos = df_neoep[ df_neoep['genome_pos'] == alt_genome_pos ]

        #if check_NMD is true, then keep any instance of transcripts where bool_NMD = False, meang these transcripts evades NMD (bool_NMD = True in my dataset means they are susceptible to NMD targeting)
        if check_NMD:
            df_genome_pos = df_genome_pos[ df_genome_pos['bool_NMD'].str.contains('FA', case = False) ]

        # alt_column = 'alt_readframe_type' if with_rf else 'my_change_type'
        if frame_stat == 0:     #only look at genomic alterations (SNV, insertions, deletions)
            alt_column = 'my_change_type'
        elif frame_stat == 1:       #only look at in-frame or frameshift events
            alt_column = 'in_frame'
        elif frame_stat == 2:       #look at the combination of genomic alteration (SNV, insertions, deletions) + reading frame status (in-frame or frameshift)
            alt_column = 'alt_readframe_type'

        #determine if I need to remove events that are susceptible to NMD
        if check_NMD:
            df_genome_pos = df_genome_pos[ df_genome_pos['bool_NMD'].str.contains('FA', case = False) ]

        #Gene Expression Percentile Threshold: apply threshold if threshold parameter > -1
        if thres_percent_exp > -1:
            df_genome_pos = df_genome_pos[ df_genome_pos['express_percentile'] >= thres_percent_exp ]

        #Proteasome Cleavage Percentile Threshold: apply threshold if threshold parameter > -1
        if thres_percent_prot > -1:
            df_genome_pos = df_genome_pos[ df_genome_pos['proteasome_percentile_2'] >= thres_percent_prot ]

        #TAP Transport Percentile Threshold: apply threshold if threshold parameter > -1
        if thres_percent_tap > -1:
            df_genome_pos = df_genome_pos[ df_genome_pos['tap_percentile_2'] >= thres_percent_tap ]

        #Threshold for Endogenous Frequency: apply threshold if threshold parameter > -1 (I will usually use 0 for this to find peptides that are unique)
        if thres_endogenous_freq > -1:
            df_genome_pos = df_genome_pos[ df_genome_pos['pep_matches_2'] <= thres_endogenous_freq ]

        #go through each sample and quantify the number of neoepitopes associated with genomic alteration position
        hash_sample_pep_freq = {}       #key = sample ID, value = array where each element is a unique peptide sequence
        for each_sample in list_samples:
            if neoep_stat == 0:
                pep_freq = df_genome_pos[ df_genome_pos[g_sample_id] == each_sample ][g_alt_peptide].unique()
            elif neoep_stat == 1:
                pep_freq = df_genome_pos[ (df_genome_pos[g_sample_id] == each_sample) &
                (df_genome_pos[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) ][g_alt_peptide].unique()
            elif neoep_stat == 2:
                #METHOD 2: use (df_genome_pos[g_aff_neoep_alt].isin["HIGH", "MID"]) to find high to intermediate neoepitopes
                pep_freq = df_genome_pos[ 
                (df_genome_pos[g_sample_id] == each_sample) & 
                (df_genome_pos[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) & 
                (df_genome_pos["delta_mhc_score"].isin( ['higher', 'equal'] )) & 
                (df_genome_pos["delta_proteasome_score"].isin( ['higher', 'equal'] )) & 
                (df_genome_pos["delta_TAP_score"].isin( ['higher', 'equal'] )) ][g_alt_peptide].unique()

            #number of neoepitopes associated with this position at this sample
            hash_sample_pep_freq[each_sample] = len( pep_freq )

        #return hash that contains number of neoepitopes per genomic alteration position for each sample
        return hash_sample_pep_freq


    def count_neoepitope_all_alts_v2( df_neoep, list_alts, hash_filters ):
        """
        Quantifies the # of neoepitopes generated from each genomic alteration
        Args:
            -df_neoep = dataframe that is the file containing all original & altered neoepitopes
            -list_alts = array of genomic alterations (string format), usually refers to positions of alterations
            -neoep_stat = integer that will perform a specific type of analysis with neoepitopes
                -0 = will consider ALL genomic alterations that produce any neoepitope
                -1 = will only consider high-affinity MHC neoepitopes that result from genomic alterations
                -2 = will only consider high-efficacious neoepitopes that result from genomic alterations
            -check_NMD = boolean that considers Nonsense-mediate Decay (NMD)
                -True = will only count events (genomic alterations) where NMD should not occur ()
                -False = will consider all events regardless if it is susceptible to NMD or not
        Returns:
            returns a hash of hashes 
        CAUTION: as a genomic alteration can be associated with multiple isoforms, need to make sure to count
        """
        hash_alt_neoep_count = {}       #k = genomic alteration position, v = another hash where k2 = sample ID & v2 = integer that is the number of neoepitopes associated with genomic alteration position for sample ID
        for each_alt_pos in list_alts:
             hash_alt_neoep_count[each_alt_pos] = count_neoepitope_per_alt( df_neoep, each_alt_pos, hash_filters )

        return hash_alt_neoep_count




    ##------THIS IS FOR REFERENCE
    def quantify_num_neoepitopes_v2_OLD( df_neoep_orig, file_summary, hash_filters ):
        """
        This quantifies the total number of neoepitopes originating from each alteration
        Args:
            -hash_filters = a hashtable that contains the following keys & values. NOTE: the filters for proteasome cleavage, TAP transport, & endogenous
                -frame_stat = integer with the following status: (NOTE: Default = 0)
                    -0 = will only consider both the alteration type (mutation, indel) (column = 'my_change_type')
                    -1 = will only consider the reading frame (column = 'in_frame')
                    -2 = will consider both the alteration type (mutation, indel) & frame-preservation (column = 'alt_readframe_type')
                -neoep_stat = integer that will perform a specific type of analysis with neoepitopes. Default = 0
                    -0 = will consider ALL genomic alterations that produce any neoepitope
                    -1 = will only consider high-affinity MHC neoepitopes that result from genomic alterations
                    -2 = will only consider high-efficacious neoepitopes that result from genomic alterations
                -check_NMD = boolean that considers Nonsense-mediate Decay (NMD) -> default is False
                    -True = will only count events (genomic alterations) where NMD should not occur ()
                    -False = will consider all events regardless if it is susceptible to NMD or not
                -thres_percent_exp = number (integer or float) that serves as gene expression threshold, but only if the expression level is greater than -1. Default = -1
                -thres_percent_prot = number (integer or float) that serves as the proteasome percentile threshold, but only if the expression level is greater than -1. Default = -1
                -thres_percent_tap = number (integer or float) that serves as the TAP percentile threshold, but only if the expression level is greater than -1. Default = -1
                -thres_endogenous_freq = the limit of the number times the peptide is found as an endogenous peptide -> I will usually use 0 for this to find peptides that are unique. Default = -1
        PROTOCOL:
            -retrieve the specific type of events (based on frame_stat & neoep_stat)
            -find the number of unique genomic events
            -for each genomic event, calculate how many neoepitopes are associated with each event
        """
        frame_stat = hash_filters['frame_stat']
        neoep_stat = hash_filters['neoep_stat']
        check_NMD = hash_filters['check_NMD']
        thres_percent_exp = hash_filters['thres_percent_exp']
        thres_percent_prot = hash_filters['thres_percent_prot']
        thres_percent_tap = hash_filters['thres_percent_tap']
        thres_endogenous_freq = hash_filters['thres_endogenous_freq']
        
        #make deep copy of original dataframe
        df_neoep = df_neoep_orig.copy( deep = True )        #deep parameter's default is True, just wrote this to make sure it is a deep copy anyways

        # alt_column = 'alt_readframe_type' if with_rf else 'my_change_type'
        if frame_stat == 0:     #only look at genomic alterations (SNV, insertions, deletions)
            alt_column = 'my_change_type'
        elif frame_stat == 1:       #only look at in-frame or frameshift events
            alt_column = 'in_frame'
        elif frame_stat == 2:       #look at the combination of genomic alteration (SNV, insertions, deletions) + reading frame status (in-frame or frameshift)
            alt_column = 'alt_readframe_type'


        #first retrieve the unique MHC alleles & genomic alteration types
        list_samples = df_neoep[g_sample_id].unique()
        list_alt_types = df_neoep[alt_column].unique()

        #determine if I need to remove events that are susceptible to NMD
        if check_NMD:
            df_neoep = df_neoep[ df_neoep['bool_NMD'].str.contains('FA', case = False) ]

        #Gene Expression Percentile Threshold: apply threshold if threshold parameter > -1
        if thres_percent_exp > -1:
            df_neoep = df_neoep[ df_neoep['express_percentile'] >= thres_percent_exp ]

        #Proteasome Cleavage Percentile Threshold: apply threshold if threshold parameter > -1
        if thres_percent_prot > -1:
            df_neoep = df_neoep[ df_neoep['proteasome_percentile_2'] >= thres_percent_prot ]

        #TAP Transport Percentile Threshold: apply threshold if threshold parameter > -1
        if thres_percent_tap > -1:
            df_neoep = df_neoep[ df_neoep['tap_percentile_2'] >= thres_percent_tap ]

        #Threshold for Endogenous Frequency: apply threshold if threshold parameter > -1 (I will usually use 0 for this to find peptides that are unique)
        if thres_endogenous_freq > -1:
            df_neoep = df_neoep[ df_neoep['pep_matches_2'] <= thres_endogenous_freq ]

        #hash_all_alts will be used to quantify each type of alteration for each sample
        hash_all_alts = {k:[] for k in list_alt_types}

        #create a hash for each sample
        case_alt_count = {}      #reads genomic alteration frequency per sample, k = sample_barcode, value = hash where k2 = specific genomic alteration type & v2 = frequency of that alteration in sample_barcode
        for each_sample in list_samples:
            case_alt_count[each_sample] = dict( hash_all_alts )     #create hash for sample that will record alterations associated with sample

        #go through each sample & alteration and quantify their frequency
        for each_sample in list_samples:
            case_alt_count[each_sample] = dict( hash_all_alts )     #create hash for sample that will record alterations associated with sample
            for each_alt in list_alt_types:
                #quantify all types of alterations for each sample
                if neoep_stat == 0:         #quantify all alterations, whether they lead to high-affinity neoepitopes or not
                    list_alts = df_neoep[ (df_neoep[g_sample_id] == each_sample) & (df_neoep[alt_column].str.contains(each_alt, case = False) ) ]['genome_pos'].unique()
                elif neoep_stat == 1:
                    #quantify the number of unique alterations that produce high-affinity neoepitopes
                    list_alts = df_neoep[ (df_neoep[g_sample_id] == each_sample) & (df_neoep[alt_column].str.contains(each_alt, case = False) ) & ( df_neoep[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) ]['genome_pos'].unique()
                elif neoep_stat == 2:
                    #METHOD 2: use (df_neoep[g_aff_neoep_alt].isin["HIGH", "MID"]) to find high to intermediate neoepitopes
                    list_alts = df_neoep[ 
                    (df_neoep[g_sample_id] == each_sample) & 
                    (df_neoep[alt_column].str.contains( each_alt, case = False )) & 
                    (df_neoep[g_aff_neoep_alt].isin( ["HIGH", "MID"] )) & 
                    (df_neoep["delta_mhc_score"].isin( ['higher', 'equal'] )) & 
                    (df_neoep["delta_proteasome_score"].isin( ['higher', 'equal'] )) & 
                    (df_neoep["delta_TAP_score"].isin( ['higher', 'equal'] )) ]['genome_pos'].unique()

                #retrieve the number of neoepitopes associated with each genomic alteration position 'alt_genome_pos' for each sample
                hash_alt_neoep_count = count_neoepitope_all_alts_v2( df_neoep, list_alts, hash_filters )

                #record the list of neoepitopes associated with each genomic alteration position 'alt_genome_pos' for each sample
                for k1,v1 in hash_alt_neoep_count.iteritems():      #k1 = string that is genomic alteration position, v1 = hash where k1_a = sample ID, v1_a = integer that is the number of unique neoepitopes associated with genomic alteration position for sample ID
                    for k1_a, v1_a in v1.iteritems():
                        case_alt_count[k1_a][each_alt] += [v1_a]

        #calculate mean & standard deviation for each sample
        case_alt_neoep_statistics = {}      #k = sample ID, v = hash where k2 = alteration type:mean/stdev & v2 = corresponding value
        for k2, v2 in case_alt_count.iteritems():        #k2 = sample ID, v2 = hash where k2_a = genomic alteration type based on 'frame_stat' & v2_a = array of integer that is the number of unique neoepitopes associated with all genomic alteration position for sample ID
            if not k2 in case_alt_neoep_statistics:
                case_alt_neoep_statistics[k2] = {}
            for k2_a, v2_a in v2.iteritems():

                ##TEST:: print "SUITE.Q_NUM_NEOEP: v2_a = ", v2_a

                list_neoep_len = [get_len for get_len in v2_a if get_len > 0]
                case_alt_neoep_statistics[k2][k2_a + ":avg"] = np.mean( list_neoep_len )
                case_alt_neoep_statistics[k2][k2_a + ":std"] = np.std( list_neoep_len )

        #create dataframe and write dataframe to file
        df_case_alt_neoep_statistics = pd.DataFrame( case_alt_neoep_statistics.values(), index = case_alt_neoep_statistics.keys(), columns = case_alt_neoep_statistics.values()[0].keys() )

        # ##FOR REFERENCE:
        # #create dataframe and write dataframe to file
        # df_alt_count = pd.DataFrame( case_alt_count.values(), index = case_alt_count.keys(), columns = case_alt_count.values()[0].keys() )

        #create label for title based on alteration, reading frame, or both
        if frame_stat == 0:
            label_alt_type = "Alterations Only"
        elif frame_stat == 1:
            label_alt_type = "Reading-Frame Only"
        elif frame_stat == 2:
            label_alt_type = "Alterations + Reading-Frame"

        #gene expression threshold
        label_thres = "Gene Exp Threshold " + str( thres_percent_exp ) if thres_percent_exp > -1 else "No Threshold Applied"
        #NMD filter label
        label_nmd = "NMD removed" if check_NMD else "May have NMD"
        label_prot = "ProtCleav Thres: " + str(thres_percent_prot) if thres_percent_prot > -1 else "NO ProtCleav"
        label_tap = "TAP Transport Thres: " + str(thres_percent_tap) if thres_percent_tap > -1 else "NO TAP"
        label_endog = "Endogenous Pep Freq Thres: " + str(thres_endogenous_freq) if thres_percent_tap > -1 else "NO Endogen Pep"

        #combine all the labels
        append_labels = label_alt_type + ", " + label_thres + ", " + label_nmd + ", " + label_prot + ", " + label_tap + ", " + label_endog

        #create title based on type of neoepitope (all, high-affinity, high-efficacious)
        if neoep_stat == 0:
            # title = "Total Number Genomic " + label_alt_type + ": " + label_thres + ", " + label_nmd 
            title = "Number of Neoepitopes per Alteration ALL: " + append_labels
        elif neoep_stat == 1:
            # title = "Number Genomic " + label_alt_type + " Producing High-Affinity Neoepitopes: " + label_thres + ", " + label_nmd
            title = "Number of Neoepitopes per Alteration Producing High-Affinity Neoepitopes: " + append_labels
        elif neoep_stat == 2:
            # title = "Number Genomic " + label_alt_type + " Producing High-Efficacious Neoepitopes: " + label_thres + ", " + label_nmd
            title = "Number of Neoepitopes per Alteration Producing High-Efficacy Neoepitopes: " + append_labels

        #write to the summary file
        write_df_to_file( df_case_alt_neoep_statistics, title, file_summary )