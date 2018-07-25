#/usr/bin/python
import sys

from cruzdb import Genome
import HTSeq
import pysam

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
# from SVSv6 import IsoformSJ, Isoform, MultiIsoform, SpliceJunction
from SVSv6 import Isoform, MultiIsoform, SpliceJunction

#Constants: file paths
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv5"
DIR_DATA = DIR_CURR + "/TestData"
DIR_RESULTS = DIR_CURR + "/TestResults"
#get mapped reads
DIR_RNASEQ = DIR_PROJ + "/150802_TophatSamples"
DIR_TABIX = DIR_PROJ + '/160427_SJFalsePositives/Data/Tabix'
DIR_THRES = DIR_PROJ + "/160803_ReadCountToThreshold/Results"

def create_obj_sj( hash_sj_info ):
    """
    Args:
        hash_sj_info = a hash_sj_info from pandas Dataframe, where each hash_sj_info is indexed by the column labels
    Function: creates a SpliceJunction instance for the splice junction recorded in the file. Information about the splice junction is recorded in a hash_sj_info in the file contained in the variable "arr_rc"
    """
    sj_id = hash_sj_info['sj_id']
    hash_sj_pos = Isoform.split_genome_pos( hash_sj_info['sj_range'] )
    chrom = hash_sj_pos['chrom']
    start = hash_sj_pos['start']
    end = hash_sj_pos['end']
    strand = '-' if int( hash_sj_info['strand'] ) == -1 else '+'       #needs to be in string format
    read_count = hash_sj_info['read_count']
    gene_sym = hash_sj_info['gene_name']
    # isoform_id = hash_sj_info['isoform_id']
    isoform_id = None
    sample_prevalence = hash_sj_info['prevalence_all']
    control_prevalence = hash_sj_info['prevalence_control']
    bool_intronic = True
    # obj_sj = SpliceJunction( sj_id, chrom, start, end, strand, read_count, gene_sym, isoform_id = None, sample_prevalence = 0, control_prevalence = 0, bool_intronic = False )
    obj_sj = SpliceJunction( sj_id, chrom, start, end, None, None, strand, read_count, gene_sym, isoform_id, sample_prevalence, control_prevalence, bool_intronic )

    return obj_sj

def sj_read_support_TEST( bam_reader, genomic_range ):
    """
    Just playing around
    """
    hash_gr = Isoform.split_genome_pos( genomic_range )       #hash_gr = hash table of genomic range
    #these are the counts
    count = 0       #count all reads that map to "genomic_range"
    uniq_count = 0      #count all uniquely mapped reads to "genomic_range"
    align_score_max_count = 0       #this also counts all uniquely mapped reads to "genomic_range"
    for i, a in enumerate( bam_reader.fetch( region = genomic_range ) ):
        # score = 0
        # if a.optional_field( "NH" ) == 1:
        #     score += 1
        print i, "read ", i, " - ", a
        print i, "dir(a) = ", dir( a )
        print i, "a.aligned = ", a.aligned
        print i, "a.read = ", a.read
        print i, "a.read_as_aligned = ", a.read_as_aligned
        print i, "a._read = ", a._read
        print i, "a._read_as_sequenced = ", a._read_as_sequenced
        print i, "a.get_sam_line = ", a.get_sam_line()     #this retrieves the read information & presents it as a .sam file line
        print i, ": a.cigar = ", a.cigar        #this is an array that contains information for each cigar. For example, if the CIGAR is 99M, then there is only 1 element in array "a.cigar" [99M. If CIGAR is 2M4926N74M, then there are 3 elements in the array "a.cigar" [2M, 4926N, 74M]

        #I can split the information for each read by splitting by tab-delimiter '\t'
        sam_line = a.get_sam_line()
        list_sam = sam_line.split('\t')
        print i, "list_sam = ", list_sam

        print "----------------\n"

def sj_read_support( bam_reader, genomic_range ):
    """
    Args:
        bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
        genomic_range = string in format chrom:start-end. If None, then uses the SJ position recorded in "self"
    Function: finds reads that support splice junctions by finding reads that uniquely map to splice junction position
    """
    # if not genomic_range:
    #     genomic_range = self.chrom + ':' + str( self.start ) + '-' + str( self.end )

    hash_gr = Isoform.split_genome_pos( genomic_range )       #hash_gr = hash table of genomic range
    #these are the counts
    count = 0       #count all reads that map to "genomic_range"
    uniq_count = 0      #count all uniquely mapped reads to "genomic_range"
    align_score_max_count = 0       #this also counts all uniquely mapped reads to "genomic_range"
    for i, a in enumerate( bam_reader.fetch( region = genomic_range ) ):
        # score = 0
        # if a.optional_field( "NH" ) == 1:
        #     score += 1

        """
        NOTE: a.cigar breaks down the meaning for the cigar, finding the end position for each cigar. 
    
        """
        for cigop in a.cigar:
            if cigop.type == 'N' and cigop.ref_iv.start == hash_gr['start'] and cigop.ref_iv.end == hash_gr['end']:

                print i, ": a.get_sam_line = ", a.get_sam_line()        #this retrieves the read information & presents it as a sam line
                print i, ": cigop = ", cigop
                print i, ": a.cigar = ", a.cigar
                print i, ">>>>>>>>>>>>>>>>>>>>\n"
                
                # if a.optional_field( "NH" ) == 1:     #check if read is uniquely-mapped read
                ##TEST:: see the output of read
                print i, ": a = ", dir(a)       #see all possible properties of object
                print i, ": a.get_sam_line = ", a.get_sam_line()        #this retrieves the read information & presents it as a sam line
                print i, ": a.from_SAM_line = ", a.from_SAM_line
                print i, ": dir( a.from_SAM_line ) = ", dir( a.from_SAM_line )
                print i, ": a.aligned = ", a.aligned
                print i, ": a.flag = ", a.flag
                print i, ": a.get_sam_line = ", a.get_sam_line, " & dir = ", dir( a.get_sam_line )
                print i, ": a.read = ", a.read      #this is the nucleotide sequence
                print i, ": a.read_as_aligned = ", a.read_as_aligned      #this is the nucleotide sequence. "a.read_as_aligned" could be the reverse-complement to "a.read"
                print i, ": a.from_pysam_AlignedRead = ", a.from_pysam_AlignedRead, " & dir = ", dir( a.from_pysam_AlignedRead )
                print i, ": NH = ", a.optional_field( "NH" )
                print i, ": aligned = ", a.aligned, " & aQual = ", a.aQual, " & read quality = ", a.read.qual
                print i, ": dir( a.cigar ) = ", dir( a.cigar )
                print i, ": dir( cigop ) = ", dir( cigop )
                print i, ": cigop.ref_iv = ", cigop.ref_iv, " & chrom = ", cigop.ref_iv.chrom, " & start = ", cigop.ref_iv.start, " & end = ", cigop.ref_iv.end
                print i, ": cigop.size = ", cigop.size
                print i, ": cigop.type = ", cigop.type
                print i, ": cigop.check = ", cigop.check, " & dir = ", dir( cigop.check )
                print i, ": cigop.query_from = ", cigop.query_from        #CONJ: I think this refers the start of the range of nucleotides that map to the genome
                print i, ": cigop.query_to = ", cigop.query_to            #CONJ: I think this refers the end of the range of nucleotides that map to the genome

                #I can split the information for each read by splitting by tab-delimiter '\t'
                sam_line = a.get_sam_line()
                list_sam = sam_line.split('\t')
                print i, "list_sam = ", list_sam

                print '------------------\n\n'
                
                count += 1
                if a.optional_field( "NH" ) == 1:
                    uniq_count += 1

                #this also counts uniquely mapped reads
                if a.aQual == 50:
                    align_score_max_count += 1


    ##TEST::
    # print " | total reads = ", count,
    # print " | unique map = ", count_uniq,
    # print " | quality 50 count = ", align_score_max_count

    return {"all_count": count, "unique_count": uniq_count, "unique_count_50": align_score_max_count}

def sj_read_support_variety_reads( bam_reader, genomic_range ):
    """
    Args:
        bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
        genomic_range = string in format chrom:start-end. If None, then uses the SJ position recorded in "self"
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
                    TO RETRIEVE UNIQUELY MAPPED READS: use command "samtools -q 4 file.bam" means any values above 4 are unique, where for tophat2 values 0 <= x <= 3 means multiple mapping at 50 means unique
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
    
    list_variety_reads = []
    for i, a in enumerate( bam_reader.fetch( region = genomic_range ) ):
        # score = 0
        # if a.optional_field( "NH" ) == 1:
        #     score += 1

        """
        NOTE: a.cigar breaks down the meaning for the cigar, finding the end position for each cigar. 
        """
        for cigop in a.cigar:
            if cigop.type == 'N' and cigop.ref_iv.start == hash_gr['start'] and cigop.ref_iv.end == hash_gr['end']:

                #if splicing event matches end positions 
                #I can split the information for each read by splitting by tab-delimiter '\t'
                sam_line = a.get_sam_line()
                list_sam = sam_line.split('\t')
                str_read_info = list_sam[2] + "|" + list_sam[3] + "|" + list_sam[5]
                list_variety_reads.append( str_read_info )
                
                
                count += 1
                if a.optional_field( "NH" ) == 1:
                    uniq_count += 1

                #this also counts uniquely mapped reads
                if a.aQual == 50:
                    align_score_max_count += 1

    unique_variety_reads = list( set(list_variety_reads) )

    print "show unique_variety_reads: "
    print unique_variety_reads


    print "list_variety_reads = ", len( list_variety_reads )
    print "# of unique variety reads = ", len( unique_variety_reads )



    ##TEST::
    # print " | total reads = ", count,
    # print " | unique map = ", count_uniq,
    # print " | quality 50 count = ", align_score_max_count

    # return {"all_count": count, "unique_count": uniq_count, "unique_count_50": align_score_max_count}: 


def sj_canonical_elsewhere( obj_sj ):
    """
    Determines if a splicing event is canonical across other genes

    """
    ##PROTOCOL: retrieve all genes & isoforms associated with position -> check if it is canonical for each isoform, if so, then record

    all_canon_isoforms = []     #records all isoforms
    # all_isoforms = Isoform.obj_cruzdb.bin_query( 'refGene', obj_sj.chrom, obj_sj.start, obj_sj.end ).all()
    all_isoforms = Isoform.obj_cruzdb.bin_query( 'knownGene', obj_sj.chrom, obj_sj.start, obj_sj.end ).all()

    print "all_isoforms = ", all_isoforms, " & type = ", type( all_isoforms )

    for i1, each_isoform in enumerate( all_isoforms ):

        # print i1, ": visit isoform = ", each_isoform.name2, " -> ", each_isoform.name
        print i1, ": visit isoform = ", each_isoform.name

        #retrieve all genes associated with position
        obj_mi = MultiIsoform( obj_sj.chrom, obj_sj.start, obj_sj.end, each_isoform.name2, each_isoform.name )
        other_annots = True         #this will check across other genomic consortiums (specifically, Ensembl)
        canon_sj = obj_mi.is_canon_sj_all( obj_sj.start, obj_sj.end, other_annots )     #return all isoforms that have canonical splicing behavior
        all_canon_isoforms += [k for k in canon_sj if canon_sj[k]]

    print "all_canon_isoforms = ", all_canon_isoforms

    rel_pos_start = "exonLeft"
    rel_pos_end = "exonRight"
    list_isoform_feat = []       #list that will record the isoform & ligated features
    for i2, each_canon_iso in enumerate( all_canon_isoforms ):
        get_gene_sym = obj_mi.hash_isoforms[k].gene_sym
        iso_feats = obj_mi.sj_get_ligated_exons_isoform( each_canon_iso, obj_sj.start, obj_sj.end, rel_pos_start, rel_pos_end )
        list_isoform_feat.append( get_gene_sym + ',' + each_canon_iso + '=' + iso_feats )

    #combine all the features 
    str_isoform_feat = " || ".join( list_isoform_feat )

    print "str_isoform_feat: "
    print str_isoform_feat

    # return str_isoform_feat


print "------------ TDD: 170720_SJ_Other_Canons.py ------------"


g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )

#example 2
# sample_name = 'yufulo'
# gene_sym = 'MITF'
# sj_pos = 'chr3:70014008-70014109'

#example 3
# sample_name = 'yuww165'
# gene_sym = 'WASH7P'
# sj_pos = 'chr1:17055-17605'

#example 4
# sample_name = 'yuhimo'
# # sample_name = 'yukadi'
# gene_sym = 'LOC100132287'
# sj_pos = 'chr1:27623647-27624428'

#example 5
# sample_name = 'gapi'
# gene_sym = 'PSMC6'
# sj_pos = 'chr14:53185756-53190682'

#example 6
#generic info for each SJ instance
sj_id = 'SJ_1'
read_count = 10
sample_prevalence = 1
control_prevalence = 0
#specific information for each SJ instance
sample_name = 'gopik'
sj_range = 'chr16:89986616-89986997'
strand = 1
gene_sym = 'MC1R'
isoform_id = 'NM_002386' 

#retrieve the mapped reads to see which reads uniquely mapped
path_bam = DIR_RNASEQ + "/tophat_sample_" + sample_name + "/accepted_hits.bam"
bam_reader = HTSeq.BAM_Reader( path_bam )
pysam_file = pysam.AlignmentFile( path_bam, 'rb' )


#show the types of reads that map to splicing event
# hash_htseq = sj_read_support( bam_reader, sj_range )
# print "see read counts: ", hash_htseq

# sj_read_support_TEST( bam_reader, sj_range )

#calculate the # of different types of reads mapping to a splicing event
# sj_read_support_variety_reads( bam_reader, sj_range )


#create SpliceJunction instance
hash_sj_info = {'sj_id': sj_id,
'sj_range': sj_range,
'strand': strand,
'read_count': read_count,
'gene_name': gene_sym,
'isoform_id': isoform_id,
'prevalence_all': sample_prevalence,
'prevalence_control': control_prevalence
}
obj_sj = create_obj_sj( hash_sj_info )

sj_canonical_elsewhere( obj_sj )

tubb3 = Isoform.obj_cruzdb.refGene.filter_by( name2 = 'TUBB3' ).all()
for i, each_tubb in enumerate( tubb3 ):
    print "each tubb ", i, " --> ", each_tubb


print "------------ TDD Completed: 170720_SJ_Other_Canons.py ------------"