#/usr/bin/python

import numpy as np

#import python classes
from Isoform import Isoform
from MultiIsoform import MultiIsoform
from KinaseFusion import KinaseFusion


#Constants - fusions.out columns once split by '\t'
COL_CHROM = 0
COL_START = 1
COL_END = 2
COL_ORIENTATION = 3
COL_READSPAN = 4
COL_MATEPAIR = 5        #mate pair reads that map to around fusion junction
COL_MATEPAIR_BREAK = 6  #mate pair reads that map to around fusion junction and one of the mate pairs spans the fusion junction
COL_CONTRADICT = 7      #number of reads that contradict the fusion (can contradict by read mapping elsewhere in genome)
COL_BP_LEFT = 8         #basepair left = number of bases covered on left side by spanning reads
COL_BP_RIGHT = 9        #basepair right = number of bases covered on right side by spanning reads

#Constants - minimum threshold for MultiIsoformFusions
min_readcount = 3
thres_readspan = 15
thres_matepair = 15
thres_matepair_break = 15
##TEST::
# min_readcount = 0
# thres_readspan = 0
# thres_matepair = 0
# thres_matepair_break = 0

"""
Requirements
-build off of SVSv4
-need cruzdb object to identify each gene
-finding duplicates
    -need to measure distance between fusions, set some sort of threshold (perhaps based on Tophat threshold)
    -I think I need Tabix to retrieve fusions that could be duplicates - BUT HOW DO I USE TABIX FOR FUSIONS!! --> I need a way to retrieve fusions that contain the same WHAT?? isoform IDs OR within the same range?


KinaseFusionV4
"""

"""
Things to add:
-if one the genes in IsoformFusion is a kinase, then need to create a KinaseFusion instance!
"""

class IsoformFusion():
    """
    Class: records only 2 isoforms (one from each part of the fusion) associated with the fusion
    """
    def __init__( self, db_type, isoform1, pos1, isoform2, pos2, orientation ):
        """
        Args:
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
        """
        #retrieve chromosome & position for isoform - this will be used to find the isoform version closest to the fusion break point
        chrom_iso1 = Isoform.obj_cruzdb.refGene.filter_by( name = isoform1 ).first()
        hash_pos1 = None if not chrom_iso1 else {'chrom': chrom_iso1.chrom, 'pos_oi': pos1}
        chrom_iso2 = Isoform.obj_cruzdb.refGene.filter_by( name = isoform2 ).first()
        hash_pos2 = None if not chrom_iso2 else {'chrom': chrom_iso2.chrom, 'pos_oi': pos2}

        self.pos1 = pos1
        self.isoform1 = Isoform( db_type, isoform1, hash_pos1 )
        self.kinase1 = None     #if not None, then will record an array of hashes ( see KinaseFusion.locate_isoform_feature() )
        
        self.pos2 = pos2
        self.isoform2 = Isoform( db_type, isoform2, hash_pos2 )
        self.kinase2 = None     #if not None, then will record an array of hashes ( see KinaseFusion.locate_isoform_feature() )
        self.orientation = orientation      #is either 'ff', 'rr', 'rf', 'fr'

        #determine kinase status
        self.kinase_stat = self.determine_kinase()

        if self.kinase_stat == 1:       #if 1st gene is a kinase
            self.kinase1 = self.retrieve_kinase_info( self.pos1, self.isoform1 )
        elif self.kinase_stat == 2:     #if 2nd gene is a kinase
            self.kinase2 = self.retrieve_kinase_info( self.pos2, self.isoform2 )
        elif self.kinase_stat == 3:     #if both genes are kinases
            self.kinase1 = self.retrieve_kinase_info( self.pos1, self.isoform1 )
            self.kinase2 = self.retrieve_kinase_info( self.pos2, self.isoform2 )

    def determine_compatibility( self ):
        """ 
        Args:
            orientation = 2-character string with one of the forms: ff, rr, fr, rf
        Function: checks if the orientation & the gene strands are compatible
        """
        strand1 = self.isoform1.strand
        strand2 = self.isoform2.strand
        if not strand1 or not strand2:
            return 'Unknown'

        #if orientation is the same (ff, rr) & gene strands are the same (+ & +, - & -), then they are compatible
        if ( strand1 == strand2 ) and ( self.orientation[0].lower() == self.orientation[1].lower() ):
            return 'Yes'
        #if orientation is different (fr, rf) & gene strands are different (+ & -, - & +), then they are compatible
        elif ( strand1 != strand2 ) and ( self.orientation[0].lower() != self.orientation[1].lower() ):
            return 'Yes'
        else:
            return 'No'

    def determine_kinase( self ):
        """
        Method: Determines if gene fusion contains a kinase gene
        """
        kinase_stat = 0     #records if kinase present in fusion, where 0 = no kinase, 1 = first gene kinase, 2 = second gene kinase, 3 = both genes are kinases
        if KinaseFusion.kinase_sym( self.isoform1.gene_sym ):
            kinase_stat += 1
        if KinaseFusion.kinase_sym( self.isoform2.gene_sym ):
            kinase_stat += 2

        return kinase_stat

    def retrieve_kinase_info( self, position, obj_isoform ):
        """
        Method: retrieves kinase information (isoform, features containing position), and returns an array of hashes
        """
        #retrieve an array of hashes, where each hash is information about a kinase
        chrom = obj_isoform.chrom
        isoform_id = obj_isoform.isoform_id
        list_kinase = KinaseFusion.kinase_info( chrom, position, isoform_id, 1 )

        #return only the first entry, which should be a h
        return list_kinase[0] if list_kinase else {}

    def isoform_info( self, bool_start = True ):
        """
        Args:
            bool_start = boolean that, if True, means look up feature information for first position, else if False, then look up information for 2nd position
        Function: returns a hash of information about the isoform (assuming the isoform is not associated with a kinase, else use KinaseFusion.kinase_info to retrieve information)
        """
        obj_isoform = self.isoform1 if bool_start else self.isoform2
        position = self.pos1 if bool_start else self.pos2

        obj_exon = self.isoform_info_canon_exon( bool_start )
        if obj_exon:        #means position is within exon
            feat_name = 'exon' + str( obj_exon.exonNum )
            rel_pos = obj_isoform.in_exon( position )
        else:       #means position is not within exon (within the intron)
            obj_intron = obj_isoform.get_intron( position )
            feat_name = 'intron' + str( obj_intron[0].exonNum ) if obj_intron else 'Unknown'
            rel_pos = None
            # rel_pos = obj_isoform.in_intron( position )

        #hash is the same type as in KinaseFusion.kinase_info & KinaseFusion.locate_isoform_feature
        return {'gene_sym': obj_isoform.gene_sym, 'isoform': obj_isoform.isoform_id, 'feature_name': feat_name, 'relative_pos': rel_pos, 'kinase_domain': None, 'strand': obj_isoform.strand}

    def isoform_info_canon_exon( self, bool_start ):
        """
        Args:
            bool_start = boolean that, if True, means look up feature information for first position, else if False, then look up information for 2nd position
            position = integer that is the genomic position
        Function: retrieves the canonical exon for a fusion depending on the isoform strand & relative
        Rules:
            if gene on +:
                -bool_start = True -> look at right side of exon
                -bool_start = False -> look at left side of exon
            if gene on -:
                -bool_start = True -> look at left side of exon
                -bool_start = False -> look at right side of exon
        """
        obj_isoform = self.isoform1 if bool_start else self.isoform2
        position = self.pos1 if bool_start else self.pos2

        obj_exon = obj_isoform.get_exon( position )

        #retrieve the exon based on strand sign & bool_start (if the gene is the start or end of the fusion)
        if obj_isoform.strand > 0:      #plus strand gene
            obj_exon = obj_isoform.get_element( position, True, 0, False, 'exonRight' ) if bool_start else  obj_isoform.get_element( position, True, 0, False, 'exonLeft' )
        else:       #minus strand gene
            obj_exon = obj_isoform.get_element( position, True, 0, False, 'exonLeft' ) if bool_start else  obj_isoform.get_element( position, True, 0, False, 'exonRight' )

        if not obj_exon:
            obj_exon = obj_isoform.get_element( position )

        return obj_exon


    def return_feature( self, bool_start = True ):
        """
        Args:
            bool_start = boolean that, if True, means look up feature information for first position, else if False, then look up information for 2nd position
        Function: returns the feature information (exon, intron, relative position within feature) for start or end position of fusion depending on 'bool_start'
        """
        kinase_info = self.kinase1 if bool_start else self.kinase2
        kin_num = [1,3] if bool_start else [2,3]

        #check if gene is a kinase
        feature_info = {}
        if self.kinase_stat in kin_num and kinase_info:
            feature_info = kinase_info
            feature_info['kinase'] = True
        else:
            #retrieve the position information for the isoform
            feature_info = self.isoform_info( bool_start )
            feature_info['kinase'] = False

        return feature_info

    """
    Functions: determine which exons are included/excluded from gene fusion
    """

    def exons_included_fusion( self, first_gene = True, inclusion = True, pos_oi = None ):
        """
        WARNING: The orientation doesn't matter (e.g. 'r' or 'f') for each gene?? I feel like I should have to incorporate orientation to see if the gene is flipped or not (flipped by 'r')
        Args:
            -first_gene = boolean where
                -True = retrieves information about the first gene in the fusion
                -False = retrieves information about the second gene in the fusion
            -inclusion = boolean where,
                -True = returns an array of exons included in the gene fusion
                -False = returns an array of exons excluded in the gene fusion
            -pos_oi = the "end" position for the gene fusion, which will be used to extract exons of interest. If it is None, then will use either self.pos1 or self.pos2, depending on boolean 'first_gene'
                -! if pos_oi is defined, make sure it is within the range of the appropriate gene ()
        Function: determines which exons are included in the fusion based on the fusion break & the orientation of the gene
        """
        #if pos_oi is None, then need to retrieve the correct breakpoint based on boolean 'first_gene'
        if not pos_oi:
            pos_oi = self.pos1 if first_gene else self.pos2

        #retrieve the genomic range based on gene fusion orientation ('f' or 'r') and boolean 'first_gene'
        if first_gene:
            fusion_range = (pos_oi, self.isoform1.boundary[1]) if self.orientation[0].lower() == 'r' else (self.isoform1.boundary[0], pos_oi)
        else:
            fusion_range = (self.isoform2.boundary[0], pos_oi) if self.orientation[1].lower() == 'r' else (pos_oi, self.isoform2.boundary[1])
        #retrieve the genomic range
        # fusion_range = (self.isoform1.boundary[0], pos_oi) if first_gene else (pos_oi, self.isoform2.boundary[1])
        
        #determine if should return exons included in gene fusion or excluded
        if inclusion:
            if first_gene:      #look at isoform1
                all_exons_included = [x for x in self.isoform1.hashExonList.values() if fusion_range[0] <= x.exonPos.location.end <= fusion_range[1]]
            else:               #ELSE look at isoform2
                all_exons_included = [x for x in self.isoform2.hashExonList.values() if fusion_range[0] <= x.exonPos.location.start <= fusion_range[1]]
            exons_oi = all_exons_included
        else:       #look for exons excluded
            if first_gene:      #look at isoform1
                all_exons_excluded = [x for x in self.isoform1.hashExonList.values() if not (fusion_range[0] <= x.exonPos.location.end <= fusion_range[1]) ]
            else:               #ELSE look at isoform2
                all_exons_excluded = [x for x in self.isoform2.hashExonList.values() if not (fusion_range[0] <= x.exonPos.location.start <= fusion_range[1]) ]
            exons_oi = all_exons_excluded
        #sort by starting position
        exons_oi.sort( key = lambda x: x.exonPos.location.start )

        return exons_oi

    def eio_find_in_out( self, bool_first_gene ):        #eio_find_in_out = Exons In-Out (included & excluded from fusion) find exons In & Out of fusion
        """
        Calculates exon expression (e.g. RPKM) for exons included and excluded from gene fusion
        Args:
            bool_first_gene = boolean where
                -True = will look at first gene in gene fusion (i.e. obj_isofuse.isoform1)
                -False = will look at first gene in gene fusion (i.e. obj_isofuse.isoform2)
        Returns:
            array of arrays, where array 1 = expression of exons included gene fusion & array 2 = expression of exons excluded from gene fusion
        """
        # exons_included_fusion( first_gene = True, inclusion = True, pos_oi = None )
        xi = self.exons_included_fusion( bool_first_gene, True, None )        #xi = exons included
        xx = self.exons_included_fusion( bool_first_gene, False, None )       #xx = exons excluded

        return [xi, xx]

    def calc_fusion_rpkm_exons_only( self, bam_reader, library_size, unique_reads = False, inclusion = True ):
        """
        Args:
            -bam_reader = HTSeq.BAM_Reader instance
            -library_size = integer that is the total number of mapped exons in the library
            -unique_reads = boolean
                -True = only consider uniquely mapped reads that map to genomic_range (i.e. a.optional_field( "NH" ) == 1)
                -False = consider all reads that map to genomic_range (unique + multimapped)
            -inclusion = boolean where,
                -True = returns an array of exons included in the gene fusion
                -False = returns an array of exons excluded in the gene fusion
        Function: calculates the RPKM of an isoform, but only considers the total length of the exons when calculate the RPKM (instead of the entire length of the gene)
        """
        #sum all the read counts for all exons in an isoform -> sum total length of all exons for a gene -> calculate RPKM by 

        #retrieve all exons of interest (exons included if inclusion = True, else exons excluded if inclusion = False)
        gene_1_exons = self.exons_included_fusion( True, inclusion )
        gene_2_exons = self.exons_included_fusion( True, inclusion )
        all_exons = gene_1_exons + gene_2_exons
        #calculate number of reads per exon
        sum_count = 0       #sum the number of reads
        exon_len = []
        for x in all_exons:
            genomic_range = x.str_genomic_pos( True )
            sum_count += Isoform.quant_genes_rpkm( bam_reader, genomic_range, unique_reads )
            exon_len.append( x.exonPos.location.end - x.exonPos.location.start )

        #sum length of all exons
        exon_len_sum = sum( exon_len )

        if exon_len_sum > 0 and library_size > 0:
            return ( 10**9 * float( sum_count ) ) / ( library_size * exon_len_sum )
        else:
            return -1       #reason I used -1 as to be able to distinguish between no expression (0) & error because of division by 0

    """
    Functions: calculate exon gene expression 
    """
    @staticmethod
    def calc_exon_express( list_exons, bam_reader, library_size, bool_rpkm = True, unique_reads = False ):
        """
        Args:
            -bam_reader = HTSeq.BAM_Reader instance, i.e. HTSeq.BAM_Reader( bam_path )
            -list_exons = array of Exon instances
            -library_size = integer that is the total number of reads in library. This is calculated by Isoform.total_mapped_reads( bam_path )
            -unique_reads = boolean
                -True = only consider uniquely mapped reads that map to genomic_range (i.e. a.optional_field( "NH" ) == 1)
                -False = consider all reads that map to genomic_range (unique + multimapped)
        Function: calculates the RPKM for each exon
        PROTOCOL: go through each exon -> calculate the # of reads mapping to exon
        """
        hash_erd = {}       #hash_erd = hash exon read density, where key = exon range, value = read density (i.e. exon RPKM or FPKM)
        for each_exon in list_exons:
            exon_pos = each_exon.str_genomic_pos( True )
            hash_erd[ exon_pos ] = each_exon.calc_read_density( bam_reader, library_size, bool_rpkm , unique_reads )
        
        return hash_erd






class MultiIsoformFusion():
    """
    Class: records all isoforms associated with each fusion 
    """
    obj_cruzdb = None       #the cruzdb object that will query the UCSC genome browser
    #kinase annotations within SVS folder
    # kinase_index = pd.read_csv( '/home/mokha/Documents/Krauthammer_Lab/PythonClasses/SVSv3/150908_AllKinase_GEPD.txt', sep = '\t', header = 0 )
    #kinase annotations elsewhere in Atlas

    #thresholds
    thres_gene_dist = 10**5
    
    def __init__( self, hash_fusion ):
        """
        Args:
            hash_fusion = dictionary with the following keys
                -chrom_start = chromosome number for starting chromosome (format: 'chr9')
                -chrom_end = chromosome number for end chromosome (format: 'chr12')
                -pos_start = integer that is one side of the fusion break
                -pos_end = integer that is the other end fusion break
                -orientation = one of the following: 'ff', 'rr', 'fr', 'rf', where r = reverse & f = forward
                -read_span = integer that is read count supporting the fusion break
                -read_matepair = integer that is the read count that maps to both sides of the fusion break (but neither read maps across the fusion break)
                -read_matepair_break
        """

        #split row into individual elements
        self.chrom_start = hash_fusion['chrom_start']
        self.chrom_end = hash_fusion['chrom_end']
        self.start = int( hash_fusion['pos_start'] )
        self.end = int( hash_fusion['pos_end'] )
        self.orientation = hash_fusion['orientation']

        self.mi_start = MultiIsoform( self.chrom_start, self.start, self.start )
        self.mi_end = MultiIsoform( self.chrom_end, self.end, self.end )

        #record read support
        self.read_span = hash_fusion['read_span']
        self.read_matepair = hash_fusion['read_matepair']
        self.read_matepair_break = hash_fusion['read_matepair_break']

        #MAY NEED TO DELETE
        # self.isoforms1 = MultiIsoform( arr_row[COL_CHROM].split('-')[0], int( arr_row[COL_START] ), int( arr_row[COL_START] ) )
        # self.isoforms2 = MultiIsoform( arr_row[COL_CHROM].split('-')[1], int( arr_row[COL_END] ), int( arr_row[COL_END] ) )

        #MAY NEED TO DELETE
        # self.hash_isoforms1 = create_isoform_obj( self.chrom_start, self.start, self.start )     #key = isoform id, value = Isoform class
        # self.hash_isoforms2 = create_isoform_obj( self.chrom_end, self.end, self.end )     #key = isoform id, value = IsoformFusion class

        self.isoform_fusions = self.create_isoform_fusion()        #key = isoformIDs in fusion (format: isoformID_1:isoformID_2), value = IsoformFusion instance

        self.gene_fusions = self.create_gene_fusion()

        self.viable = self.check_fusion()      #property that determines if a fusion is indeed a viable fusion (True = it is a viable fusion otherwise False)


    def check_fusion( self ):
        """
        Function: determines if a fusion can truly be considered a fusion
        """
        #check 1: make sure both genes are not the same
        gene_sym1 = MultiIsoform.get_gene_syms( self.chrom_start, self.start, self.start )
        gene_sym2 = MultiIsoform.get_gene_syms( self.chrom_end, self.end, self.end )
        if set( gene_sym1 ).intersection( gene_sym2 ):
            return False

        #check 2: see if there is a enough space between both genes (for intrachromosomal fusions only)
        if self.chrom_start == self.chrom_end:
            if np.absolute( self.start - self.end ) < MultiIsoformFusion.thres_gene_dist:
                return False

        #check 3: see if there is enough reads to support fusion
        if not self.threshold_pass_v2():
            return False

        #return True if all criteria have been passed
        return True

    def create_gene_fusion( self ):
        """
        Function: create IsoformFusion instance for each pair of gene
        """
        #retrieve all gene for both positions - returns a list of strings that are the isoform IDs
        ig_1 = MultiIsoform.get_isoform_gene( self.chrom_start, self.start, self.start )    #ig = isoform + gene symbol
        ig_2 = MultiIsoform.get_isoform_gene( self.chrom_end, self.end, self.end )      #ig = isoform + gene symbol

        #if nothing present, then return empty hash
        if not ig_1 or not ig_2:
            return {}

        #go through each isoform ID, recording each pair as an IsoformFusion instance
        hash_gene_fusion = {}   #key = isoform IDs for both genes, value = IsoformFusion instance
        for k1, v1 in ig_1.iteritems():     #k = isoform ID, v = gene symbol
            for k2, v2 in ig_2.iteritems():     #k = isoform ID, v = gene symbol
                key = v1 + ':' + v2         #key is a combination of both gene symbols (not isoform ID)
                hash_gene_fusion[key] = IsoformFusion( k1, self.start, k2, self.end, self.orientation )

        return hash_gene_fusion

    def create_isoform_fusion( self ):
        """
        Function: create IsoformFusion instance for each pair of isoforms
        """
        #retrieve all isoforms for both positions - returns a list of strings that are the isoform IDs
        isoforms_1 = MultiIsoform.get_isoforms( self.chrom_start, self.start, self.start )
        isoforms_2 = MultiIsoform.get_isoforms( self.chrom_end, self.end, self.end )

        #if nothing present, then return empty hash
        if not isoforms_1 or not isoforms_2:
            return {}

        #go through each isoform ID, recording each pair as an IsoformFusion instance
        hash_isoforms_fusion = {}   #key = isoform IDs for both genes, value = IsoformFusion instance
        for iso1 in isoforms_1:
            for iso2 in isoforms_2:
                key = iso1 + ':' + iso2         #to save all isoform versions
                hash_isoforms_fusion[key] = IsoformFusion( iso1, self.start, iso2, self.end, self.orientation )

        return hash_isoforms_fusion

    @staticmethod
    def get_isoforms( chrom, start, end ):
        all_isoforms = Isoform.obj_cruzdb.bin_query( 'refGene', chrom, start, end ).all()
        #Method 1: If I only want to consider mRNA protein-coding isoforms
        # return [x.name for x in all_isoforms if "NM_" in x.name]

        #Method 2: if I want all isoforms
        """
        IMPORTANT: the reason I decided to retrieve all isoforms & not just protein-coding isoforms is because there could a splice junction that could correctly splice a non-mRNA transcript, therefore reducing erroneously stating an SJ is aberrant when indeed it maps to a non-mRNA transcript. So need to consider reads that map to non-mRNA transcripts as well.
        """
        return [x.name for x in all_isoforms]

    def threshold_pass_v2( self ):
        """
        Function: determine if the row recording the fusion is sufficient to be recorded
        """
        #if any of the read counts are below a specific minimum threshold, then skip
        if ( self.read_span < min_readcount ) or ( self.read_matepair < min_readcount ) or ( self.read_matepair_break < min_readcount ):
            return False

        #apply threshold V1 - if does not pass threshold then skip this fusion
        if ( self.read_span < thres_readspan ) and ( self.read_matepair < thres_matepair ) and ( self.read_matepair_break < thres_matepair_break ):
            return False

        #if all thresholds are pass
        return True


    @staticmethod
    def threshold_pass( row ):
        """
        Function: determine if the row recording the fusion is sufficient to be recorded
        """
        arr_row = row.split( '\t' )
        #if any of the read counts are below a specific minimum threshold, then skip
        if ( int( arr_row[COL_READSPAN] ) < min_readcount ) or ( int( arr_row[COL_MATEPAIR] ) < min_readcount ) or ( int( arr_row[COL_MATEPAIR_BREAK] ) < min_readcount ):
            return False

        #apply threshold V1 - if does not pass threshold then skip this fusion
        if ( int( arr_row[COL_READSPAN] ) < thres_readspan ) and ( int( arr_row[COL_MATEPAIR] ) < thres_matepair ) and ( int( arr_row[COL_MATEPAIR_BREAK] ) < thres_matepair_break ):
            return False

        #if all thresholds are pass
        return True


    @staticmethod
    def chrom_int( str_chrom ):
        """
        Function: returns the integer form of the chromosom
        """
        if str_chrom.replace('chr', '').isdigit():      #format: chr#, e.g. chr9, chr12, chr7
            return int( str_chrom.replace('chr', '') )
        elif str_chrom.split('_')[0].replace('chr', '').isdigit():      #format: chr10_GL383546v1_alt or something similar
            return int( str_chrom.split('_')[0].replace('chr', '') )
        else:
            return None


    @staticmethod
    def sort_gene_order( row ):
        """
        Function: sorts order of fusions based on their position
        """
        arr_row = row.split( '\t' )

        #retrieve chromsome numbers
        chroms = arr_row[COL_CHROM].split('-')
        chrom_start = MultiIsoformFusion.chrom_int( chroms[0] )
        chrom_end = MultiIsoformFusion.chrom_int( chroms[1] )

        reverse_order = False       #if True, then reverse order of positions in fusion, else if False keep order
        if not chrom_start or not chrom_end:        #if either is None, then do not reverse order
            reverse_order = False
        elif chrom_start > chrom_end:
            reverse_order = True
        elif chrom_start == chrom_end:
            if int( arr_row[COL_START] ) > int( arr_row[COL_START] ):
                reverse_order = True

        if reverse_order:
            chrom_start = arr_row[COL_CHROM].split('-')[1]
            chrom_end = arr_row[COL_CHROM].split('-')[0]
            pos_start = int( arr_row[COL_END] )
            pos_end = int( arr_row[COL_START] )
            #reverse orientation of orientation
            reverse_orientation = {'ff': 'rr', 'rr': 'ff', 'rf': 'rf', 'fr': 'fr'}
            orientation = reverse_orientation[ arr_row[COL_ORIENTATION] ]
        else:
            chrom_start = arr_row[COL_CHROM].split('-')[0]
            chrom_end = arr_row[COL_CHROM].split('-')[1]
            pos_start = int( arr_row[COL_START] )
            pos_end = int( arr_row[COL_END] )
            orientation = arr_row[COL_ORIENTATION]

        return {'chrom_start': chrom_start, 'chrom_end': chrom_end, 'pos_start': pos_start, 'pos_end': pos_end, 'orientation': orientation}

    @staticmethod
    def row_to_hash( row, sort_order = True ):
        """
        Args:
            row = tab-delimited string where, once split by '/t', will have rows corresponding to fusions.out column constants defined above
            sort_order = boolean
                -True = sort from least to greatest
                -False = keep order as is
        Function: retrieve a row of fusion data from tophat fusion's fusions.out
        """
        arr_row = row.split( '\t' )
        hash_fusion = {}
        if sort_order:
            hash_fusion = MultiIsoformFusion.sort_gene_order( row )
        else:
            hash_fusion['chrom_start'] = arr_row[COL_CHROM].split('-')[0]
            hash_fusion['chrom_end'] = arr_row[COL_CHROM].split('-')[1]
            hash_fusion['pos_start'] = int( arr_row[COL_START] )
            hash_fusion['pos_end'] = int( arr_row[COL_END] )
            hash_fusion['orientation'] = arr_row[COL_ORIENTATION]

        hash_fusion['read_span'] = arr_row[COL_READSPAN]
        hash_fusion['read_matepair'] = arr_row[COL_MATEPAIR]
        hash_fusion['read_matepair_break'] = arr_row[COL_MATEPAIR_BREAK]

        return hash_fusion


