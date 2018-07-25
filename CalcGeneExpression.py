#/usr/bin/python

import scipy
from scipy import stats
import pandas as pd

import HTSeq


##NOTE: I HAVE NOT FULLY WRITTEN UP THIS CLASS YET
class CalcGeneExpression():

    def __init__( self ):
        pass

    """
    Calculate Gene Expression Percentile
    """
    def calculate_percentileofscore( list_gene_exp, curr_exp ):
        """
        Args:
            list_gene_exp = an array (list or numpy array) of gene expression values. As of now, this should not contain any duplicate values
            curr_exp = float value that is the gene expression value
        Function: determine the 
        """
        #Method 1
        return scipy.stats.percentileofscore( list_gene_exp, curr_exp, 'strict' )

    def checkIfFloat( x ):
        """
        Function: determines if the value 'x' can be converted from string to float (assumes x is a string, but doesn't have to be a string)
        """
        try:
            float( x )
            return True
        except ValueError:
            return False

    def checkIfFloatAboveZero( x ):
        """
        Function: determines if the value 'x' can be converted from string to float (assumes x is a string, but doesn't have to be a string)
        """
        try:
            float( x )
            if float(x) > 0:
                return True
            else:
                return False
        except ValueError:
            return False

    def retrieve_sample_expression( df_expression, sample_name, isoform_id, gene_sym ):
        """
        Retrieves the max gene expression for a specific 

        Args:
            -df_expression = dataframe of the expression profiles where the rows are the gene IDs & column = sample names. The rows are in the format "isoform_id;gene_sym;genomic_range" (e.g. NM_001184730;PIGT;chr20:44044706-44054885)
            -sample_name = string that is the sample of interest
            -isoform_id = string taht is the isoform ID. This is the first method to find 
        """
        key_isoform_id = isoform_id + ';'
        key_gene_sym = ';' + gene_sym + ';'


        ##TEST:: print "RSE 0: isoform_id = ", isoform_id, " & gene_sym = ", gene_sym

        df_gene_exp = df_expression.filter( like = key_isoform_id, axis = 0 )
        ##TEST: print "RSE 1: ", len( df_gene_exp ), " found for isoform = ", isoform_id
        if len( df_gene_exp ) == 0:
            df_gene_exp = df_expression.filter( like = key_gene_sym, axis = 0 )
            # print "RSE 2: ", len( df_gene_exp ), " found for gene_sym = ", gene_sym

        ##TEST::
        # print "length of df_expression = ", len( df_expression )
        # print "df_expression columns: "
        # print df_expression.columns.values

        # print "df_gene_exp: ", len( df_gene_exp )
        # print df_gene_exp


        return float(df_gene_exp[sample_name].max())


    def calculate_percentileofscore_allsamples( df_expression, hash_sample_allgenes_exp, list_samples, isoform_id, gene_sym ):
        """
        Calculates the gene expression percentile for all samples that contain aberrant splicing of gene
        
        Args:
            -df_expression = dataframe of the expression profiles for each sample. This has the RPKM expression for each gene for each sample
            -hash_sample_allgenes_exp = hash where k = sample name, v = list of expression values for the sample in the key. This will be used to calculate the percentile expression for a specific gene
            -list_samples = an array of sample names of interest
            -row_id = index for df_expression (this should be some sort of gene identification)

        Output:
            returns 3 hashes
                -hash_sample_percentile: key = sample names & value = corresponding percentile expression for specific gene 'isoform_id' (or if can't use isoform_id, then use 'gene_sym')
                -hash_sample_absolute_exp: key = sample names & value = corresponding absolute expression for specific gene 'isoform_id' (or if can't use isoform_id, then use 'gene_sym')
                -hash_sample_allgenes_exp: key = sample names & values = gene expression in RPKM (or whatever is reported in df_expression). This hash is updated if the sample does not contain a list of expressions
        """

        hash_sample_percentile = {}     #k = sample name, v = percentile of gene expressed
        hash_sample_absolute_exp = {}   #k = sample name, v = absolute gene expression (as reported by df_expression, which usually is in RPKM)
        for each_sample_name in list_samples:
            #if sample is already in hash_sample_allgenes_exp, that means it already contains a list of expressed genes
            if not each_sample_name in hash_sample_allgenes_exp:
                hash_sample_allgenes_exp[each_sample_name] = [ float(x) for x in df_expression[each_sample_name].unique() if checkIfFloatAboveZero(x) ]

            #record percentile expression for each sample
            # sample_curr_exp = df_expression[each_sample_name][row_id]
            sample_curr_exp = retrieve_sample_expression( df_expression, each_sample_name, isoform_id, gene_sym )
            hash_sample_absolute_exp[each_sample_name] = sample_curr_exp
            # hash_sample_percentile[each_sample_name] = scipy.stats.percentileofscore( hash_sample_allgenes_exp[each_sample_name], sample_curr_exp, 'strict' )
            hash_sample_percentile[each_sample_name] = calculate_percentileofscore( hash_sample_allgenes_exp[each_sample_name], sample_curr_exp )

        return [hash_sample_percentile, hash_sample_absolute_exp, hash_sample_allgenes_exp]

    def sj_read_support_variety_reads_multisample( list_samples, genomic_range ):
        """
        finds reads that support splice junctions by finding reads that uniquely map to splice junction position, and does this for multiple samples

        Args:
            -list_samples = array of sample names that contain the SJ position with the range 'genomic_range'
            -genomic_range = string in format chrom:start-end. If None, then uses the SJ position recorded in "self"
        Output:
            returns a hash where k = sample name, v = hash from def sj_read_support_variety_reads(). See sj_read_support_variety_reads() to see what each field means

        """
        hash_sample_sj_support = {}     #k = sample name, v = hash from def sj_read_support_variety_reads()
        for each_sample in list_samples:
            path_bam = DIR_RNASEQ + "/tophat_sample_" + each_sample + "/accepted_hits.bam"
            bam_reader = HTSeq.BAM_Reader( path_bam )
            hash_sample_sj_support[each_sample] = sj_read_support_variety_reads( bam_reader, genomic_range )

        return hash_sample_sj_support


    def sj_read_support_variety_reads( bam_reader, genomic_range ):
        """
        Args:
            -bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
            -genomic_range = string in format chrom:start-end. If None, then uses the SJ position recorded in "self"
        Function: finds reads that support splice junctions by finding reads that uniquely map to splice junction position
        Output:
            {"all_count": count, "unique_count": uniq_count, "unique_count_50": align_score_max_count, "count_all_variety_reads": len( list_variety_reads ), "count_unique_variety_reads": len( unique_variety_reads ) }
            returns a hash with the following values:
                -"all_count" = all reads that support splicing event with position 'genomic_range'
                -"unique_count" = uniquely-mapped reads supporting splicing event with position 'genomic_range'. This looks for NH == 1
                -"unique_count_50" = same as "unique_count", but looks at quality of read (aQual == 50), where 50 means uniquely mapped reads
                    -see "Protocol: 15.10.30 - Samtools" 
                        -50 (or 255): unique mapping (NH:i:1)
                        -3: maps to 2 locations in the target (NH:i:2, but I have also seen NH:i:3)
                        -2: maps to 3 locations
                        -1: maps to 4-9 locations (NH:i:4 or higher)
                        -0: maps to 10 or more locations
                        TO RETRIEVE UNIQUELY MAPPED READS: use command "samtools -q 5 file.bam" means any values above 4 are unique, where for tophat2 values 0 <= x <= 3 means multiple mapping at 50 means unique
                -"count_all_variety_reads" = total count of all reads with different end positions that support the splicing event. This is considering all reads, therefore there will be duplicates (meaning they will have the same start & end points - think of "thickBlocks" & "thinBlocks" for UCSC Genome Browser)
                -"count_unique_variety_reads" = count of # of reads with different end positions. This is the actually number of reads that support the splicing event & has different end points.
                    -The hypothesis is the reads with the same end points are just the same RNA fragments sequenced during RNA-sequencing, therefore the more different types of reads supporting a splicing event, the more convincine the support.
        """
        # if not genomic_range:
        #     genomic_range = self.chrom + ':' + str( self.start ) + '-' + str( self.end )

        hash_gr = Isoform.split_genome_pos( genomic_range )       #hash_gr = hash table of genomic range
        #these are the counts
        count = 0       #count all reads that map to "genomic_range"
        uniq_count = 0      #count all uniquely mapped reads to "genomic_range"
        align_score_max_count = 0       #this also counts all uniquely mapped reads to "genomic_range"
        list_variety_reads = []     #records the variety of reads associated with splicing event. (By "variety", I mean reads that with different "end positions")
        for i, a in enumerate( bam_reader.fetch( region = genomic_range ) ):
            # score = 0
            # if a.optional_field( "NH" ) == 1:
            #     score += 1

            """
            NOTE: a.cigar breaks down the meaning for the cigar, finding the end position for each cigar. 
            """
            for cigop in a.cigar:
                if cigop.type == 'N' and cigop.ref_iv.start == hash_gr['start'] and cigop.ref_iv.end == hash_gr['end']:

                    #record the end position & the CIGAR associated with each read that supports the splicing event
                    #I can split the information for each read by splitting by tab-delimiter '\t'
                    sam_line = a.get_sam_line()
                    list_sam = sam_line.split('\t')
                    str_read_info = list_sam[2] + "|" + list_sam[3] + "|" + list_sam[5]
                    list_variety_reads.append( str_read_info )
                    
                    #quantify the # of reads and unique reads that maps to splicing event
                    count += 1
                    if a.optional_field( "NH" ) == 1:
                        uniq_count += 1

                    #this also counts uniquely mapped reads
                    if a.aQual == 50:
                        align_score_max_count += 1

        #quantify # of unique "variety" of reads that support splicing
        unique_variety_reads = list( set(list_variety_reads) )

        print "show unique_variety_reads: "
        print unique_variety_reads


        print "list_variety_reads = ", len( list_variety_reads )
        print "# of unique variety reads = ", len( unique_variety_reads )

        return {"all_count": count, "unique_count": uniq_count, "unique_count_50": align_score_max_count, "count_all_variety_reads": len( list_variety_reads ), "count_unique_variety_reads": len( unique_variety_reads ) } 