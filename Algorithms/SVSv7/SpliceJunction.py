#/usr/bin/python
import numpy as np

import HTSeq
import pysam

from Exon import Exon
from Isoform import Isoform
from MultiIsoform import MultiIsoform


class SpliceJunction( MultiIsoform ):
    inst_mi = {}       #instances of MultiIsoform class, where key = isoform_id, value = MultiIsoform instance 

    def __init__( self, db_type, sj_id, chrom, start, end, leftmost, rightmost, strand, read_count, gene_sym, isoform_id = None, bool_record_canon_isoforms = False, sample_prevalence = 0, control_prevalence = 0, bool_intronic = False ):
        """
        Args:
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
            -sj_id = string that is the splice junction ID (e.g. JUNC00003 or anything that is a string)
            -chrom = string that is the chromosome of the splice junction (format: chr#, e.g. chr2, chr5, chr12)
            -start & end = integers that are the start & end genomic positions of the splicing event. Regardless of the strand sign, start < end
            -leftmost & rightmost = integers that are the leftmost & rightmost positions that map to this splice junction. Regardless of the strand sign, leftmost < rightmost
                -Basically some reads have gapped alignments (not continuously mapped), yet these reads have a lower & upper bound to where the map in the genome. The idea behind using leftmost & rightmost is to determine how this splicing event will be incorporated into a reconstructed transcript (transcript reconstruction is in IsoformSJ)
            -strand = string in the format '+' or '-'
            -gene_sym = string that is the name of the gene symbol. Need this so I can find the constitutive exons and constitutive introns (aka constitutive splice junctions).
            -isoform_id = optional parameter, if this is provided then will only record one isoform. If gene_sym is also defined, then isoform_id needs to be an isoform of gene_sym else it will not be created
                -NOTE: assigning isoform_id creates only Isoform instance for each SJ -> this can save a lot of time and memory when generating each SpliceJunction instance (one Isoform instance for MultiIsoform as opposed to 70 isoforms. I think gene TTN has this issue)
            -bool_record_canon_isoforms = boolean where
                -True = will record all isoforms where this SJ is canonical. Use this when I want to find constitutive exons or see how prevalent an exon is across isoforms
                -False = will only record isoform in 'isoform_id', assuming it is given.
            -bool_intronic = boolean
                -True = will record aberrants for SJ if the SJ lands within the intron
                -False = will NOT record aberrants for SJ if the SJ lands within the intron
        NOTE: Don't want to use isoform as the SJ could be canonical for one isoform but not another.
        NOTE: with super( className ), the classes inherited need to have new-style class definitions
        New Style: class SomeClass( object ): --> this requires the (object) element
        Old Style: class OldClass():
        """
        #if gene_sym = None (no gene symbol is assigned), then retrieve the gene symbol where the SJ is the closest
        if not gene_sym:
            list_gene_syms = SpliceJunction.get_gene_syms( chrom, start, end, db_type )
            gene_sym = list_gene_syms[0]

        ##TEST:: print "SpliceJunction.isoform_id = ", isoform_id, " & chrom = ", chrom, " | start = ", start, " | end = ", end, " & gene_sym = ", gene_sym

        super( SpliceJunction, self ).__init__( db_type, chrom, start, end, gene_sym, isoform_id )     #allows for multiple inheritance - inherits everything in parenthesis
        # MultiIsoform.__init__( chrom, start, end )      #inheritance - only inherits MultiIsoform class
        self.sj_id = sj_id
        self.chrom = chrom
        self.start = int( start )
        self.end = int( end )
        #NOTE: new in SVSv6: This is considering "leftmost" & "rightmost"
        self.leftmost = int( leftmost ) if leftmost else int( start )
        self.rightmost = int( rightmost ) if rightmost else int( start )
        convert_ss = { '+' : 1, '-' : -1 }
        self.sj_strand = convert_ss.get( strand, 0 )
        # self.sj_strand = strand
        self.read_count = int( read_count )

        ##TEST:: print "SpliceJunction - MI instance: ", self.gene_sym_oi, " & self.hash_isoforms = ", self.hash_isoforms

        #need the strand sign for the gene - assuming all isoforms are the same 
        # self.strand = self.hash_isoforms.values()[0].strand       #I don't think I need this

        #determine if SJ is canonical
        if bool_record_canon_isoforms:
            self.assigned_isoform = self.assign_canon_isoform()      #should be an array of isoforms to assign
            if self.assigned_isoform:       #if canonical isoform is assigned, then SJ is true
                self.isoform_aberrants = {}
                self.canon = True
            else:       #if no canonical isoforms assigned, look for isoforms that experience aberrant splicing
                self.isoform_aberrants = self.assign_aberrant_isoform( bool_intronic )      #hash where keys = isoform id, value = aberrations occurring to isoform (exon skips, frameshifts, etc.)
                self.assigned_isoform = self.isoform_aberrants.keys() if self.isoform_aberrants else []
                self.canon = False
        else:
            self.isoform_aberrants = self.assign_aberrant_isoform( bool_intronic )      #hash where keys = isoform id, value = aberrations occurring to isoform (exon skips, frameshifts, etc.)
            self.assigned_isoform = [isoform_id] if isoform_id != None else []
            #determine if the SJ is canonical or not in "isoform_id"
            self.canon = self.hash_isoforms[isoform_id].is_canon_sj( self.start, self.end, False )      #boolean that, if True, means SJ is canonical, else if False, means not canonical splice junction
            

        ##NOTE: maybe should be in 
        self.sample_prevalence = sample_prevalence
        self.control_prevalence = control_prevalence

    def __eq__( self, obj_sj ):
        """ Function: compares two Splice Junction objects to see if they are equivalent """
        if not obj_sj:
            return False

        if self.chrom == obj_sj.chrom and self.start == obj_sj.start and self.end == obj_sj.end:
            return True
        else:
            return False
    

    def __str__( self ):
        """ Function: return string representation of class """
        pos = self.chrom + ':' + str( self.start ) + '-' + str( self.end )
        return self.sj_id + ", " + pos + ", read = " + str( self.read_count ) + ", strand sign = " + str( self.sj_strand ) + ", canonical = " + str( self.canon ) + ", assigned isoforms = " + " | ".join( self.assigned_isoform ) + ", gene symbol = " + self.gene_sym_oi

    def str_genomic_pos( self, str_only = True ):
        """
        Args:
            str_only = boolean that:
                -True = will return in string format (e.g chrom:start-end)
                -False = returns in array format, where [chrom, start, end], where chrom = string, start & end = int
        Function: returns position of gene in string form
        """
        if str_only:
            return str( self.chrom ) + ":" + str( self.start ) + "-" + str( self.end )
        else:
            return {'chrom': self.chrom, 'start': self.start, 'end': self.end}

    #NOTE: new in SVSv6: This is considering "leftmost" & "rightmost"
    def determine_sj_issues_all( self ):
        """
        Determines if there are any concerns with the splicing event. Issues include:
            -is self.leftmost & self.start in the same feature (i.e. in the same exon or intron)
            -is self.rightmost & self.end in the same feature (i.e. in the same exon or intron)

        Args:
            -pos_1 & pos_2 = integers that are genomic positions. Assumes both positions are from the same chromosome
            -isoform_id = isoform_id that will be used to retrieve 
        Output:
            returns boolean, where it returns True if the features that contain both positions are the same, else return False
        """
        hash_feat_start = self.positions_same_feature_isoform( self.start, self.leftmost )
        hash_feat_end = self.positions_same_feature_isoform( self.end, self.rightmost )

        return [same_feat_start, same_feat_end]

    #NOTE: new in SVSv6: This is considering "leftmost" & "rightmost"
    def determine_sj_issues_isoform( self, isoform_id ):
        """
        Determines if there are any concerns with the splicing event. Issues include:
            -is self.leftmost & self.start in the same feature (i.e. in the same exon or intron)
            -is self.rightmost & self.end in the same feature (i.e. in the same exon or intron)

        Args:
            -pos_1 & pos_2 = integers that are genomic positions. Assumes both positions are from the same chromosome
            -isoform_id = isoform_id that will be used to retrieve 
        Output:
            returns boolean, where it returns True if the features that containt both positions are the same, else return False
        """
        bool_feat_start = self.positions_same_feature_isoform( self.start, self.leftmost, isoform_id )
        bool_feat_end = self.positions_same_feature_isoform( self.end, self.rightmost, isoform_id )

        return [bool_feat_start, bool_feat_end]

    #NOTE: new in SVSv6: This is considering "leftmost" & "rightmost"
    def positions_same_feature_all( self, pos_1, pos_2 ):
        """
        same as def positions_same_feature_isoform, but does this for all isoforms
        """
        isoform_same_feat = {}      #k = isoform_id, v = boolean, where True means pos_1 & pos_2 are in the same feature & False means they are not in the same feature
        for k in self.hash_isoforms:      #k = isoform id, self.hash_isoforms[k] = Isoform instance
            isoform_same_feat[k] = self.positions_same_feature_isoform( pos_1, pos_2, k )

        return isoform_same_feat

    #NOTE: new in SVSv6: This is considering "leftmost" & "rightmost"
    def positions_same_feature_isoform( self, pos_1, pos_2, isoform_id ):
        """
        determines if both positions pos_1 & pos_2 are on the same feature in a specific isoform

        Args:
            -pos_1 & pos_2 = integers that are genomic positions. Assumes both positions are from the same chromosome
            -isoform_id = isoform_id that will be used to retrieve 
        Output:
            returns boolean, where it returns True if the features that containt both positions are the same, else return False
        """
        feat_1 = self.hash_isoforms[isoform_id].get_element_v2( pos_1, 1, -1, True, None )     #parameters for get_element_v2( self, position, feat_search = 1, direction = 0, bool_strand = True, rel_pos = None )
        if not feat_1:      #if not found in exon, look for it in intron
            feat_1 = self.hash_isoforms[isoform_id].get_element_v2( pos_1, 3, -1, True, None )

        feat_2 = self.hash_isoforms[isoform_id].get_element_v2( pos_2, 1, -1, True, None )     #parameters for get_element_v2( self, position, feat_search = 1, direction = 0, bool_strand = True, rel_pos = None )
        if not feat_2:      #if not found in exon, look for it in intron
            feat_2 = self.hash_isoforms[isoform_id].get_element_v2( pos_2, 3, -1, True, None )

        #if both variables are 'None'
        if feat_1 == None and feat_2 == None:
            return None

        #check if both features are equivalent or not
        if feat_1 == None or feat_2 == None:
            return False
        elif feat_1 == feat_2:
            return True
        else:
            return False




    @staticmethod
    def get_gene_syms( chrom, start, end, db_type ):
        """
        Args:
            -chrom = string in the format chr# (e.g. chr9, chr12)
            -start & end = integers that are positions of interest (ends of the SJ)
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
        Function: finds the gene symbol associated with the SJ position
        Returns: returns a list of all gene symbols associated with position
        """
        #retrieve the gene symbol that is the closest in terms of position
        gene_sym = None
        hash_all_isoforms = Isoform.get_isoforms_by_pos_db_all( chrom, start, end, db_type )
        if hash_all_isoforms:
            gene_sym = hash_all_isoforms.name2
            gene_sym = np.unique( [x.name for x in isoform_variants] ) if db_type == 3 else np.unique( [x.name2 for x in isoform_variants] )
        # else:
        #     #if no gene symbol is found, then find the first gene symbol occurrence that is protein-coding, else just find the first occurrence regardless if it is protein-coding or not.
        #     all_isoforms = Isoform.obj_cruzdb.bin_query( 'refGene', chrom, start, end ).all()
        #     isoforms_pc = [x for x in all_isoforms if "NM_" in x]       #isoforms_pc = isoforms protein coding
        #     try:
        #         gene_sym = isoforms_pc[0].name2 if isoforms_pc else all_isoforms[0].name2
        #     except IndexError:  #Error: no gene symbol is associated with position 'chrom:start-end'
        #         # print "SpliceJunction.get_gene_sym: no gene symbol found in position ", chrom, ":", start, "-", end
        #         gene_sym = None

        return gene_sym

    def str_aberrations( self ):
        """ Function: returns string of all aberrations associated with  """
        #For each isoform, extract aberration_score, exon skip, frameshift, reading frame, effects, sample prevalence, control prevalence
        if not self.isoform_aberrants:
            return None

        record_aberrants = self.sj_id + '::'
        for i, (k,v) in enumerate( self.isoform_aberrants.iteritems() ):      #k = isoform ID, v = hash where k2 = categories of aberrations (e.g. exon skip, frameshift, reading frame), & v2 = value of each category
            if i > 0:
                record_aberrants += ">>"

            record_aberrants += "isoform_id = " + k + ';'
            record_aberrants += "exon_skip = " + ','.join( [str(x) for x in v['exon_skip']] ) + ';'
            record_aberrants += "frame_preserved = " + str( v['frame_preserved'] ) + ';'
            record_aberrants += "calc_rf = " + str( v['calc_rf'][0] ) + '-' + str( v['calc_rf'][1] ) + ';'
            record_aberrants += "effects = " + ','.join( v['effects'] )
        
        return record_aberrants

    """
    Functions: assign isoform to splice junction
    """
    def check_canon_sj_all_db( self, sj_pos, list_db = [1,2,3,4] ):
        """
        Uses def check_canon_sj_each_db() and goes through all database annotations 
        Args:
            -sj_pos = string that is the format "chrNum:start-end" (e.g. chr9:5-89). Note that "chr" needs to be present. This should be the start & end position for a splicing event to determine if this SJ position occurs across different 
            -list_db = an array that contains all numbers referring to annotation databases "db_type"
                -db_type = integer that chooses the genomic database of interest
                    -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                    -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                    -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                    -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        """
        test_list_isoforms = []
        for db_type in list_db:
            test_list_isoforms += self.check_canon_sj_each_db( sj_pos, db_type )

        return list( set( test_list_isoforms ) )

    def check_canon_sj_each_db( self, sj_pos, db_type ):
        """
        this retrieves all annotations for a given position range 'sj_pos' for a specific annotation database
        Args:
            -obj_sj = instance of SpliceJunction class
            -sj_pos = string that is the format "chrNum:start-end" (e.g. chr9:5-89). Note that "chr" needs to be present. This should be the start & end position for a splicing event to determine if this SJ position occurs across different 
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        Returns: returns a list of gene symbols that contain the sj position "sj_pos"
        """
        hash_pos = Isoform.split_genome_pos( sj_pos )
        #retrieve all annotatations from specific database
        isoform_info = Isoform.get_isoforms_by_pos_db_all( hash_pos['chrom'], hash_pos['start'], hash_pos['end'], db_type )

        #if no isoform found, then no need to search if "sj_pos" is in isoform since there are no isoforms
        if not isoform_info:
            return []

        #retrieve all gene symbols & go through each one
        if db_type == 3:
            all_gene_syms = list( set( [ str(x.name) for x in isoform_info ] ) )
        else:
            all_gene_syms = list( set( [ str(x.name2) for x in isoform_info ] ) )

        
        #go through each gene symbol to see if 'sj_pos' (genomic range of SJ) 
        test_list_isoforms = []
        for each_sym in all_gene_syms:
            check_canon_sym = self.is_canon_sj_other_annots( hash_pos['start'], hash_pos['end'], db_type, each_sym )

            if check_canon_sym:
                # return True
                test_list_isoforms.append( each_sym )

        # return False
        return test_list_isoforms

    def assign_canon_isoform( self ):
        """ Function: returns a list of isoforms that have canonical splicing behavior (i.e. splice at donor-acceptor sites) """
        # canon_sj = self.is_canon_sj_all( self.start, self.end, False )     #returns hash where k = isoform id, v = boolean where True = canonical & False = non-canonical

        ##TEST:: print "\n----SJ START: class_SJ.assign_canon_isoform = ", self.str_genomic_pos()

        #Version 1 - will compare SJ to only known gene models based on self.db_type
        canon_sj = super( SpliceJunction, self ).is_canon_sj_all( self.start, self.end, False )     #returns hash where k = isoform id, v = boolean where True = canonical & False = non-canonical

        # #Version 2 - will compare SJ to only known gene models to a genomic annotation database different from sefl.db_type, but I need to specify which database I'm using (e.g. db_type = 4 will use GENCODE)
        # other_db_type = 4       #4 --> this is for GENCODE
        # canon_sj = super( SpliceJunction, self ).is_canon_sj_all( self.start, self.end, True, other_db_type )     #returns hash where k = isoform id, v = boolean where True = canonical & False = non-canonical

        ##TEST::  print "----SJ END: class_SJ.assign_canon_isoform = ", self.str_genomic_pos(), " & canon_isoforms = ", [k for k in canon_sj if canon_sj[k]], " & # of isoforms = ", len(canon_sj), " & # of isoforms that are true = ", len([k for k in canon_sj if canon_sj[k]]), "\n"

        return [k for k in canon_sj if canon_sj[k]]     #return all isoforms that have canonical splicing behavior


    def assign_aberrant_isoform( self, intronic = False ):
        """
        Args:
            intronic = boolean
                -True = return an isoform ID even if it does land in the intronic space
                -False = do not return isoform ID if only lands in intronic space
        Function: if no canonical isoforms detected, then find the aberrations for each isoform & return list of isoforms with the minimum amount of defects
        """
        #need to adjust the start & end from 0-based to 1-based. Note that start position is less than end position regardless of strand sign
        pos_start = self.adjust_sj_to_exon( True )
        pos_end = self.adjust_sj_to_exon( False )
        isoform_aberrations = self.is_aberrant_sj_all( pos_start, pos_end, intronic )     #key = isoform ID, value = hash with elements 'exon_start', 'exon_end', 'exon skip', 'preserved_frame'

        ##TEST:: print "Class SJ.AssAbIso 1: isoform_aberrations = ", isoform_aberrations

        if not isoform_aberrations:
            return {}

        #severity of defect (greatest to least): intronic, frameshift, exonic (not on donor-acceptor end), exon skip
        # isoform_scores = { k:v['score'] for k,v in isoform_aberrations.iteritems() if v }
        isoform_scores = [ v['score'] for v in isoform_aberrations.values() if v ]
        # #if return a blank hash, then return None
        # if not isoform_scores:
        #     return None
        min_score = min( isoform_scores )
        select_isoforms = { str(k):v for k,v in isoform_aberrations.iteritems() if v['score'] == min_score }

        ##TEST:: print "Class SJ.AssAbIso 1: select_isoforms with MIN SCORE = ", select_isoforms

        ##TEST::
        # print "assign isoform 1: isoform_aberrations = ", isoform_aberrations
        # print self.sj_id, " - assign isoform 2: isoform_scores = ", isoform_scores
        # print self.sj_id, " - assign isoform 3: min isoform = ", select_isoforms

        return select_isoforms
        

    def isoform_check_before_coding( self, isoform_id ):
        """
        Determines if the start of the SJ (5' end of the SJ) is before the coding region of of the isoform
        Args:
            isoform_id = string that is the isoform ID, will be used to retrieve info for specific isoform ID in self.hash_isoforms
        Returns:
            Boolean where True means the start of the SJ comes before the start of the coding region, & False means it does not come before the start of the coding region.
        """
        #check if isoform ID is present
        if not isoform_id in self.hash_isoforms:
            return False
        #if gene is a non-coding gene (e.g. NR_), then this will not have "pos_start_codon" (because doesn't code proteins)
        if not self.hash_isoforms[isoform_id].pos_start_codon:
            return False

        if self.hash_isoforms[isoform_id].strand < 0:
            return True if self.hash_isoforms[isoform_id].pos_start_codon < self.end else False
        else:
            return True if self.hash_isoforms[isoform_id].pos_start_codon > self.start else False

    def all_check_before_coding( self, isoform_id ):
        """
        Determines if the start of the SJ (5' end of the SJ) is before the coding region of of the isoform
        Args:
            isoform_id = string that is the isoform ID, will be used to retrieve info for specific isoform ID in self.hash_isoforms
        Returns:
            returns a hash where the key = isoform ID, value = boolean (True = SJ start comes before coding region, False = SJ starts is after the coding region)
        """
        hash_all_isoforms = {}      #key = isoform ID, value = boolean (True = SJ start comes before coding region, False = SJ starts is after the coding region)
        for k,v in self.hash_isoforms.iteritems():      #k = isoform ID, v = Isoform instance
            hash_all_isoforms[k] = self.isoform_check_before_coding( k )

        return hash_all_isoforms

    """
    Function: find next exon
    """
    def spliced_elems( self, isoform_id, canon = False ):       #formerly known as spliced_exons()
        """ 
        Args:
            isoform_id = string that is the isoform ID
            canon = boolean:
                -if True, only look for exons where splicing occurs at the exon ends (aka donor-acceptor sites) - therefore "exonLeft" for starting position, & "exonRight" for ending position
                -if False, then look for all exons that are spliced together ('exonLeft', 'exonRight', 'withinExon')
        Function: retrieves the exons spliced together by splice junction, returns a tuple where [0] = exon that contains start splice site & [1] = exon that contains end splice site
        """

        ##TEST:: print "SJ - splice_exons: self.hash_isoforms = ", self.hash_isoforms
        
        if not isoform_id in self.hash_isoforms:
            return ( None, None )

        #if canon is True, then look for splicing events that are canonical
        canon_start = 'exonRight' if canon else None
        canon_end = 'exonLeft' if canon else None

        #retrieve exon - new & improved b/c allows to choose direction (lower or higher exon), choose CDS or exon, & take into account the strand sign
        # x_start = self.hash_isoforms[isoform_id].get_element( self.start, True, -1, True, canon_start )     #parameters for get_element( self, position, bool_exon = True, direction = 0, bool_strand = True, rel_pos = None )
        # x_end = self.hash_isoforms[isoform_id].get_element( self.end, True, 1, True, canon_end )

        #retrieve the exons (or introns) that contain the start & end position 
        x_start = self.hash_isoforms[isoform_id].get_element_v2( self.start, 1, -1, True, canon_start )     #parameters for get_element_v2( self, position, feat_search = 1, direction = 0, bool_strand = True, rel_pos = None )
        if not x_start:     #if not found, then look in intron
            x_start = self.hash_isoforms[isoform_id].get_element_v2( self.start, 3, -1, True, canon_start )
       
        x_end = self.hash_isoforms[isoform_id].get_element_v2( self.end, 1, 1, True, canon_end )
        if not x_end:     #if not found, then look in intron
            x_end = self.hash_isoforms[isoform_id].get_element_v2( self.end, 3, -1, True, canon_end )
 
        return ( x_start, x_end )

    def elem_nearest_exon( self, elem, sj_pos, which_exon, isoform_id ):
        """
        Args:
            elem = Exon instance that is either type 'exon' or 'intron'. If type 'intron', then will look for nearest exon
            sj_pos = integer that is one of the ends of the splice junction (either self.start or self.end)
            which_exon = integer that is either 0 or 1, where:
                -0 = means this is the ligated exon with the lower numerical genomic position (for + strand: lower exon number; for - strand: higher exon number)
                -1 = means this is the ligated exon with the higher numerical genomic position (for + strand: higher exon number; for - strand: lower exon number)
            isoform_id = string that is the isoform ID. USE = used for self.hash_isoforms
        Function: finds the nearest exon in isoform ID based on Exon instance 'elem'. This function is primarily used by def spliced_elems_position_range()
        NOTE: when grouping parameters for the Exon instance & SJ position (i.e. elem, sj_pos):
            -if starting position of SJ, then sj_pos = self.start & elem = elems[0]
            -if ending position of SJ, then sj_pos = self.end & elem = elems[1]
        """
        if not elem:
            str_sj_pos = "None"
            exon_canon = False         #True if previous exon keeps it canonical form (original start & end position), else False (meaning SJ spliced within exon)
            new_elem = None
        else:        #make sure it is not None 
            #first check if element is an exon or intron
            if elem.exonPos.type.lower() == 'intron':       #if intron, then get the exon
                if which_exon == 0:     #for the lower genome-positioned exon, retrieve the nearest exon that has a lower genomic position than the intron 
                    calc_exon_num = elem.exonNum if self.hash_isoforms[isoform_id].strand > 0 else elem.exonNum + 1       #get the exon that is before the SJ start position, differs for different gene strand signs
                else:       #for the higher genome-positioned exon, retrieve the nearest exon that has a higher genomic position than the intron 
                    calc_exon_num = elem.exonNum + 1 if self.hash_isoforms[isoform_id].strand > 0 else elem.exonNum
                get_exon = self.hash_isoforms[isoform_id].get_exon_num( calc_exon_num )
            else:
                get_exon = elem

            if not get_exon:
                str_sj_pos = "None"
                exon_canon = False         #True if previous exon keeps it canonical form (original start & end position), else False (meaning SJ spliced within exon)
                new_elem = None
            else:
                exon_str_pos = get_exon.str_genomic_pos( False )   #returns array where [0] = chrom, [1] = start pos, & [2] = end pos
                str_sj_pos = exon_str_pos['chrom'] + ':' + str( exon_str_pos['start'] ) + '-' + str( sj_pos )
                exon_canon = True if sj_pos == exon_str_pos['end'] else False     #True if previous exon keeps it canonical form (original start & end position), else False (meaning SJ spliced within exon)
                new_elem = get_exon

        return {'str_sj_pos': str_sj_pos, 'exon_canon': exon_canon, 'new_elem': new_elem}

    
    def spliced_elems_position_range( self, isoform_id, canon = False, allow_none = True ):     #formerly known as spliced_exons_position_range
        """ 
        Args:
            isoform_id = string that is the isoform ID. USE = used for self.hash_isoforms
            canon = boolean:
                -if True, only look for exons where splicing occurs at the exon ends (aka donor-acceptor sites) - therefore "exonLeft" for starting position, & "exonRight" for ending position
                -if False, then look for all exons that are spliced together ('exonLeft', 'exonRight', 'withinExon')
            allow_none = boolean:
                -if True, then will still output string even if either of the exons from self.spliced_elems() returns None
                -if False, then will only output string only if both exons are not None, else this function will return None
        Function: similar to def spliced_elems(), but will return the position from the start & end exon connected by the SJ, but will also will incorporate positions of SJ. Returns string in format 'chrA:startA-SJ_A*chrB:SJ_B-endB' 
        """
        elems = self.spliced_elems( isoform_id, canon )       #returns tuple where [0] = preceding exon & [1] = proceeding exon connected by SJ
        if not allow_none:
            if not elems[0] or not elems[1]:
                return {}

        hash_exon_start = self.elem_nearest_exon( elems[0], self.start, 0, isoform_id )
        hash_exon_end = self.elem_nearest_exon( elems[1], self.end, 1, isoform_id )

        # #get the lower exon (b/c lower genomic position) ligated by the splice junction
        # new_elems = {}      #this will record the exon elements ligated by 
        # str_sj_pos = ''
        # if elems[0]:        #make sure it is not None 
        #     #first check if element is an exon or intron
        #     if elems[0].exonPos.type.lower() == 'intron':       #if intron, then get the exon
        #         calc_exon_num = elems[0].exonNum if self.hash_isoforms[isoform_id].strand > 0 else elems[0].exonNum + 1       #get the exon that is before the SJ start position, differs for different gene strand signs
        #         get_exon = self.hash_isoforms[isoform_id].get_exon_num( calc_exon_num )
        #     else:
        #         get_exon = elems[0]

        #     if not get_exon:
        #         str_sj_pos = "None"
        #         prev_exon_canon = False         #True if previous exon keeps it canonical form (original start & end position), else False (meaning SJ spliced within exon)
        #         new_elems['exon_start'] = None
        #     else:
        #         pos_0 = get_exon.str_genomic_pos( False )   #returns array where [0] = chrom, [1] = start pos, & [2] = end pos
        #         str_sj_pos = pos_0['chrom'] + ':' + str( pos_0['start'] ) + '-' + str( self.start )
        #         prev_exon_canon = True if self.start == pos_0['end'] else False     #True if previous exon keeps it canonical form (original start & end position), else False (meaning SJ spliced within exon)
        #         new_elems['exon_start'] = get_exon
        # else:
        #     str_sj_pos = "None"
        #     prev_exon_canon = False         #True if previous exon keeps it canonical form (original start & end position), else False (meaning SJ spliced within exon)
        #     new_elems['exon_start'] = None
        
        # str_sj_pos+= '*'

        # #get the higher exon (b/c higher genomic position) ligated by the splice junction
        # if elems[1]:        #make sure it is not None
        #     #first check if element is an exon or intron
        #     if elems[1].exonPos.type.lower() == 'intron':       #if intron, then get the exon
        #         calc_exon_num = elems[1].exonNum + 1 if self.hash_isoforms[isoform_id].strand > 0 else elems[1].exonNum       #get the exon that is after the SJ start position, differs for different gene strand signs
        #         get_exon = self.hash_isoforms[isoform_id].get_exon_num( calc_exon_num )
        #         # get_exon = self.hash_isoforms[isoform_id].get_exon_num( int( elems[1].exonNum ) )
        #     else:
        #         get_exon = elems[1]

        #     pos_1 = get_exon.str_genomic_pos( False )  #returns array where [0] = chrom, [1] = start pos, & [2] = end pos
        #     str_sj_pos+= pos_1['chrom'] + ':' + str( self.end ) + '-' + str( pos_1['end'] )
        #     next_exon_canon = True if self.end == pos_1['start'] else False       #True if subsequent exon keeps it canonical form (original start & end position), else False (meaning SJ spliced within exon)
        #     new_elems['exon_end'] = get_exon
        # else:
        #     str_sj_pos+= "None"
        #     next_exon_canon = False         #True if subsequent exon keeps it canonical form (original start & end position), else False (meaning SJ spliced within exon)
        #     new_elems['exon_end'] = None

        # """
        # *NOTE: prev_elem < next_elem in terms of numerical position
        # -str_sj_pos = string that is in the format 'chrA:startA-SJ_A*chrB:SJ_B-endB' from lower to higher position
        # -prev_elem = the previous element (exon or intron) spliced by the SJ, contains SJ's starting position
        # -next_elem = the next element (exon or intron) spliced by the SJ, contains SJ's ending position
        # -prev_exon_canon = boolean that, if True, means the starting element is canonical, else if False means it is non-canonical (i.e. the splicing event spliced somewhere other than the canonical splice site)
        # -next_exon_canon = same as 'prev_exon_canon', but for the element (exon or intron) contains SJ's ending position
        # """
        # return {'str_sj_pos': str_sj_pos, 'prev_elem': new_elems['exon_start'], 'next_elem': new_elems['exon_end'], 'prev_exon_canon': prev_exon_canon, 'next_exon_canon': next_exon_canon}

        str_sj_pos = hash_exon_start['str_sj_pos'] + '*' + hash_exon_end['str_sj_pos']
        return {'str_sj_pos': str_sj_pos, 'prev_elem': hash_exon_start['new_elem'], 'next_elem': hash_exon_end['new_elem'], 'prev_exon_canon': hash_exon_start['exon_canon'], 'next_exon_canon': hash_exon_end['exon_canon']}


    def make_sj_modified_exon( self, isoform_id, bool_start ):
        """
        Args:
            isoform_id = string that is the isoform ID
            bool_start = boolean that
                -True = uses the start position of the SJ
                -False = uses the end position of the SJ
        Function: create a new Exon based on the position of the splicing, 
        """
        isoform_ss = self.hash_isoforms[isoform_id].strand      #isoform_ss = isoform strand sign

        #retrieve the exons
        elems = self.spliced_elems( isoform_id, False )       #returns tuple where [0] = preceding exon & [1] = proceeding exon connected by SJ
        elem_oi = elems[0] if bool_start else elems[1]
        sj_curr_pos = self.start if bool_start else self.end 

        if elem_oi:        #make sure it is not None 
            #first check if element is an exon or intron
            if elem_oi.exonPos.type.lower() == 'intron':       #if intron, then get the exon
                calc_exon_num = elem_oi.exonNum if self.hash_isoforms[isoform_id].strand > 0 else elem_oi.exonNum + 1       #get the exon that is before the SJ start position, differs for different gene strand signs
                get_exon = self.hash_isoforms[isoform_id].get_exon_num( calc_exon_num )
            else:
                calc_exon_num = elem_oi.exonNum
                get_exon = elem_oi

            #retrieve the exon position
            exon_gp = get_exon.str_genomic_pos( False )   #exon_gp = exon genomic position, returns array where [0] = chrom, [1] = start pos, & [2] = end pos
            new_exon_gp = {'chrom': exon_gp['chrom'], 'start': exon_gp['start'], 'end': sj_curr_pos} if bool_start else {'chrom': exon_gp['chrom'], 'start': sj_curr_pos, 'end': exon_gp['end']}

            #Check if SJ position is at exon end - True if previous exon keeps it canonical form (original start & end position), else False (meaning SJ spliced within exon)
            if bool_start:      #canonical if SJ starting position matches with the end of the previous exon (5' splice site for "+" strand gene & 3' splice site for "-" strand gene)
                stat_exon_canon = True if sj_curr_pos == exon_gp['end'] else False
                new_exon_gp = {'chrom': exon_gp['chrom'], 'start': exon_gp['start'], 'end': sj_curr_pos}
            else:       #canonical if the SJ ending position matches the start position of the next exon (3' end for "+" strand gene & 5' end for "-" strand gene)
                stat_exon_canon = True if sj_curr_pos == exon_gp['start'] else False
                new_exon_gp = {'chrom': exon_gp['chrom'], 'start': sj_curr_pos, 'end': exon_gp['end']}
            
            #retrieve the 5' splice site for the exon -> this only depends on the strand sign, not on bool_start
            exon_pos_five_prime = new_exon_gp['end'] if isoform_ss > 0 else new_exon_gp['start']
            #get reading frame - retrieve the 5' position of the exon as this position will be assigned the reading frame
            get_rf = self.hash_isoforms[isoform_id].calc_reading_frame( exon_pos_five_prime, isoform_ss )

            ##TEST::print "SJ MODIFIED EXONS: new_exon_gp = ", new_exon_gp, " & SJ = ", str( self )

            #create an Exon instance
            hash_exon_info = self.hash_isoforms[isoform_id].get_feat_info( new_exon_gp['start'], new_exon_gp['end'], get_rf, calc_exon_num, stat_exon_canon, 'exon' )
            #create Exon instance 
            obj_exon = Exon( hash_exon_info )
            return obj_exon
        else:
            return None


    def spliced_elems_position_range_v2( self, isoform_id, canon = False, allow_none = True ):     #formerly known as spliced_exons_position_range
        """ 
        Args:
            isoform_id = string that is the isoform ID
            canon = boolean:
                -if True, only look for exons where splicing occurs at the exon ends (aka donor-acceptor sites) - therefore "exonLeft" for starting position, & "exonRight" for ending position
                -if False, then look for all exons that are spliced together ('exonLeft', 'exonRight', 'withinExon')
            allow_none = boolean:
                -if True, then will still output string even if either of the exons from self.spliced_elems() returns None
                -if False, then will only output string only if both exons are not None, else this function will return None
        Function: similar to def spliced_elems(), but will return the position from the start & end exon connected by the SJ, but will also will incorporate positions of SJ. Returns a hash of both elements
        """
        isoform_ss = self.hash_isoforms[isoform_id].strand      #isoform_ss = isoform strand sign
        elems = self.spliced_elems( isoform_id, canon )       #returns tuple where [0] = preceding exon & [1] = proceeding exon connected by SJ
        if not allow_none:
            if not elems[0] or not elems[1]:
                return {}

        #retrieve each exon - each one will return an Exon instance
        exon_start = self.make_sj_modified_exon( isoform_id, True )
        exon_end = self.make_sj_modified_exon( isoform_id, False )

        return {'exon_start': exon_start, 'exon_end': exon_end}



    #MIGHT DELETE: May not need this function because of def spliced_elems_position_range()
    def missing_exon_partial( self, isoform_id, bool_prev_exon = True ):
        """
        Function: retrieves the part of the exon that is spliced out from aberrant splicing. Returns a hash that contains all information about the exon spliced into
        """
        if bool_prev_exon:
            sj_pos = self.start
            exon_direction = -1     #will be used for Isoform.get_element(), direction parameter
        else:
            sj_pos = self.end
            exon_direction = 1     #will be used for Isoform.get_element(), direction parameter

        #find the exon where the splicing is within the exon
        #NOTE: did bool_strand = False because I'm looking for exons based on absolute position (lower exon = lower genomic position, higher exon = higher genomic position), not by it's relation to the gne strand
        partial_exon = self.hash_isoforms[isoform_id].get_element( sj_pos, 1, exon_direction, False, 'withinElem' )       #( position, bool_exon = True, direction = 0, bool_strand = True, rel_pos = None )

        if not partial_exon:
            return None
        else:       #retrieve the exonic position removed by the aberrant splice junction
            if bool_prev_exon:      #if previous exon, then exon end to start of SJ
                missing_exon_info = partial_exon.get_exon_info()
                missing_exon_info['start'] = sj_pos
            else:       #if subsequent exon, then exon start to end of SJ
                missing_exon_info = partial_exon.get_exon_info()
                missing_exon_info['end'] = sj_pos

        #return the missing exon segment as an Exon instance
        missing_exon = Exon( missing_exon_info )
        return missing_exon

    def get_sj_exon_skips( self ):
        """
        Function: use in conjuction with def missing_exons() to find SJs that experience exon skips, full or partial
        """
        #record all isoforms & associated aberrations in the hash
        hash_es = {k:v for k,v in self.isoform_aberrants.iteritems() if v['exon_skip']}     #k = isoform ID, v = defect categories (e.g. exon_skips, frame_preserved) - see function Isoform.is_aberrant_sj()
        if hash_es:
            isoform_max_es = max( self.isoform_aberrants, key = lambda k: len(hash_es) )
            return isoform_max_es


        hash_es_partial = {k:v for k,v in self.isoform_aberrants.iteritems() if 'exonic' in v['effects']}
        if hash_es_partial:
            return hash_es_partial.keys()[0]

        #if no full exon skips or partial exon skips are found, then return None
        return None


    def missing_exons( self, isoform_id = None ):
        """
        Args:
            isoform_id = string that is the isoform ID. If not specified, then will look at isoform that experiences that most full exon skips 
        Function: returns a list of Exon instances of exons skipped or lost due to the splicing event. This includes exons that are fully skipped 
        """
        #if isoform_id is None, then find isoform with the maximum number of exon skips
        if not isoform_id:
            isoform_id = self.get_sj_exon_skips()
            #if no isoforms found, then this means no full or partial exon skips are found
            if not isoform_id:
                return []
            
        #retrieve all exons that are fully skipped -> retrieve all regions that are partially skipped
        es_full = []      #records all positions of skipped exons
        for x in self.isoform_aberrants[isoform_id]['exon_skip']:
            es_full.append( self.hash_isoforms[isoform_id].get_exon_num(x) )

        #get partial exon skipping (basically aberrant SJ splices into exon)
        es_partial = []     #record exonic regions missing due to exonic splicing
        if 'exonic' in self.isoform_aberrants[isoform_id]['effects']:
            ligated_exon_prev = self.missing_exon_partial( isoform_id, True )
            if ligated_exon_prev:
                es_partial.append( ligated_exon_prev )

            ligated_exon_next = self.missing_exon_partial( isoform_id, False )
            if ligated_exon_next:
                es_partial.append( ligated_exon_next )

        #return all exonic regions missing, full skips & partial skips, where each is an Exon instance
        missing_es_all = es_full + es_partial
        return missing_es_all

    def get_retained_intronic_region( self, isoform_id, bool_start = True ):
        """
        Args:
            isoform_id = string that is the 
            bool_start = boolean that
                -True = checks if the start position of the SJ is intronic
                -False = checks if the end position of the SJ is intronic
        Function: Used in conjuction with "def intron_retained()", retrieve the intronic region that is retained by a splicing event. If an intronic region is found, then an Exon instance with the retained intron is returned, else None
        """
        #check if should use start or end position
        sj_pos = self.start if bool_start else self.end

        #NOTE: I used 'withinElem' because 
        sj_intronic = self.hash_isoforms[isoform_id].get_element_v2( sj_pos, 3, 0, False, 'withinElem' )
        #see if the start position of the splice junction is within the intron, and retrieve the intronic region retained
        if sj_intronic:
            hash_intron = sj_intronic.get_exon_info()
            #determine the adjusted position of the intron based on the SJ 
            if bool_start:  #if bool_start is True, then the range of the intron is (intron start) - (sj_start)
                hash_intron['end'] = sj_pos
            else:           #if bool_start is False, then the range of the intron is (sj_end) - (intron end)
                hash_intron['start'] = sj_pos
            obj_retained_intron = Exon( hash_intron )    #create Exon instance that is the retained intronic region
        else:
            obj_retained_intron = None         #this means the SJ start position is not intronic

        return obj_retained_intron


    def intron_retained( self, isoform_id = None ):
        """
        Function: returns retained introns due to the aberrant splicing event, returns it as an Exon instance
        """
        #if isoform_id is None, then retrieve the first isoform that contains an intronic effect
        isoform_intron = [k for k,v in self.isoform_aberrants.iteritems() if 'intronic' in v['effects']]
        if not isoform_intron:
            return None
        else:
            isoform_id = isoform_intron[0]

        #need to find element in hashIntronList that contains 
        if 'intronic' in self.isoform_aberrants[isoform_id]['effects']:
            #NOTE: I used 'withinElem' 
            sj_start_intron = self.hash_isoforms[isoform_id].get_element_v2( self.start, 3, 0, False, 'withinElem' )
            #see if the start position of the splice junction is within the intron, and retrieve the intronic region retained
            obj_intron_1 = self.get_retained_intronic_region( isoform_id, True )
            obj_intron_2 = self.get_retained_intronic_region( isoform_id, False )
        else:
            obj_intron_1 = None
            obj_intron_2 = None

        #return Exon instances of both introns retained by start & end position (if no intron retained, then that value will be None)
        return { "retained_intron_1": obj_intron_1, "retained_intron_2": obj_intron_2 }


    def spliced_elems_all( self, canon = False ):
        """
        Args:
            isoform_id = string that is the isoform ID
            canon = boolean:
                -if True, only look for exons where splicing occurs at the exon ends (aka donor-acceptor sites) - therefore "exonLeft" for starting position, & "exonRight" for ending position
                -if False, then look for all exons that are spliced together ('exonLeft', 'exonRight', 'withinExon')
        Function: retrieves the exons spliced together by splice junction for all isoforms, returns a tuple where [0] = exon that contains start splice site & [1] = exon that contains end splice site """
        #go through all isoforms assigned to SJ -> retrieve the exons associated
        isoform_exons = {}      #k = isoform id, v = tuple where [0] = exon containing start position & [1] = exon containing end position

        for isoform_id in self.assigned_isoform:
            isoform_exons[isoform_id] = self.spliced_elems( isoform_id, canon )

        return isoform_exons

    def has_exon( self, exon, isoform_id, direction = 0, canon = False ):
        """ 
        Args:
            exon = instance of Exon that will be used to see if it is present in splice junction
            isoform_id = string for isoform to look into
            direction = integer that will check if exon is in a specific position, i.e. start or end
                -0 = is if exon is present, position doesn't matter
                -1 = checks if exon is the start exon
                -2 = checks if exon is the end exon
            canon = boolean:
                -if True, only look for exons where splicing occurs at the exon ends (aka donor-acceptor sites) - therefore "exonLeft" for starting position, & "exonRight" for ending position
                -if False, then look for all exons that are spliced together ('exonLeft', 'exonRight', 'withinExon')
        Function: check if SJ contains exon for specific isoform 'isoform_id'
        """
        #if exon is not an Exon object but instead is None, return False
        if exon == None:
            return False

        sj_exons = self.spliced_elems( isoform_id, canon )
        #if isoform_id not in self.hash_isoforms, then return None
        if not sj_exons:
            return None
    
        #check if exon is present depending on direction
        if direction == 0:
            return True if exon in sj_exons else False 
        else:
            return True if exon == sj_exons[ direction - 1 ] else False

    def has_exon_modified( self, exon, isoform_id, direction = 0, canon = False ):
        """
        Args:
            exon = instance of Exon that will be used to see if it is present in splice junction
            isoform_id = string for isoform to look into
            direction = integer that will check if exon is in a specific position, i.e. start or end
                -0 = is if exon is present, position doesn't matter
                -1 = checks if exon is the start exon
                -2 = checks if exon is the end exon
            canon = boolean:
                -if True, only look for exons where splicing occurs at the exon ends (aka donor-acceptor sites) - therefore "exonLeft" for starting position, & "exonRight" for ending position
                -if False, then look for all exons that are spliced together ('exonLeft', 'exonRight', 'withinExon')
        Function: similar to SJ.has_exon(), however searches for exons that could have been modified, especially if one end of the SJ splices into an intron (effectively "stretching the exon" as the exon will include part of the intronic region)
        """
        #if exon is not an Exon object but instead is None, return False
        if exon == None:
            return False

        hash_exon_ligated = self.spliced_elems_position_range( isoform_id, canon )
        sj_exons = ( hash_exon_ligated['prev_elem'], hash_exon_ligated['next_elem'] )
        #if isoform_id not in self.hash_isoforms, then return None
        if not sj_exons:
            return None
    
        #check if exon is present depending on direction
        if direction == 0:
            return True if exon in sj_exons else False 
        else:
            return True if exon == sj_exons[ direction - 1 ] else False


    def has_exon_all( self, exon, direction ):
        """ 
        Args:
            exon = instance of Exon that will be used to see if it is present in splice junction
            isoform_id = string for isoform to look into
            direction = integer that will check if exon is in a specific position, i.e. start or end
                -0 = is if exon is present, position doesn't matter
                -1 = checks if exon is the start exon
                -2 = checks if exon is the end exon
        Function: check if SJ contains exon across all isoforms
        """
        isoform_bool = {}        #k = isoform id, v = boolean that, if True, means the exon is present for that isoform, otherwise False
        for isoform_id in self.assigned_isoform:
            isoform_bool[isoform_id] = self.has_exon( exon, direction )

        return isoform_bool

    """
    Functions: calculate reading frame
    """

    def adjust_sj_to_exon( self, bool_start ):
        """
        Args:
            bool_start = boolean that, if True, will look at the start position, else will look at the end position
        Function: this will adjust the splice junction end to the exon position, basically switching from 0-based to 1-based (where, for + strand, use -1 & for - strand use +1)
        
        NOTE: end > start just based on numerical position
        NOTE: reason for -1 & +1 - shifting from splice junction to exon is almost like shifting from 0-based to 1-based, therefore need to move back 1 base for + gene & move forward 1 base for - gene
        """
        return self.start if bool_start else self.end + 1

    ##DELETE? I DON'T THINK THIS WORKS PROPERLY AND SHOULD PERHAPS DISCARD IT
    # def adjust_sj_to_exon_backup( self, bool_start ):
    #     """
    #     Args:
    #         bool_start = boolean that, if True, will look at the start position, else will look at the end position
    #     Function: this will adjust the splice junction end to the exon position, almost like switching from 0-based to 1-based (where, for + strand, use -1 & for - strand use +1)
        
    #     NOTE: end > start just based on numerical position
    #     NOTE: reason for -1 & +1 - shifting from splice junction to exon is almost like shifting from 0-based to 1-based, therefore need to move back 1 base for + gene & move forward 1 base for - gene
    #     """
    #     if bool_start:
    #         position = self.start - 1 if self.strand > 0 else self.start
    #     else:
    #         position = self.end if self.strand > 0 else self.end + 1

    #     return position

    def sj_reading_frame( self, bool_start = True ):
        """
        Args:
            bool_start = boolean that, if True, will look at the start position, else will look at the end position
        Function: calculates the reading frame for one end (depends on bool_start) of the SJ and returns a hash where k = isoform ID & v = integer that is a frame number (0, 1, or 2) or None, depends on if position is found in exon
        """
        isoform_frame = {}      #k = isoform id (string), v = integer that is a frame number (0, 1, or 2) or None, depends on if position is found in exon
        
        #readjust positions from 0-based to 1-based & calculate reading frame 
        position = self.adjust_sj_to_exon( bool_start ) if bool_start else self.adjust_sj_to_exon( bool_end )
        isoform_frame = self.calc_reading_frame_all( position )     #isoform_frame = hash where key = isoform ID & v = integer that is a frame number (0, 1, or 2) or None, depends on if position is found in exon

        return isoform_frame

    def sj_frame_preserved( self, preserved_only = False ):
        """
        Args: 
            preserved_only = boolean:
                -True = only return hash of isoforms where reading frame is preserved
                -False = return hash for all isoforms, whether reading frame is preserved or not
        Function: determines if the SJ preserves the reading frame - returns a hash where key = isoform ID & value = boolean where True = frame preserved & False = frame not preserved
        """
        #reason for '-1': shifting from splice junction to exon is almost like shifting from 0-based to 1-based, therefore need to move back 1 base to be at  
        pos_start = self.adjust_sj_to_exon( True )
        pos_end = self.adjust_sj_to_exon( False )
        isoform_frame_preserved = self.frame_preserved_all( pos_start, pos_end )        #k = isoform ID, v = boolean where True = frame preserved & False = frame not preserved
        if preserved_only:
            isoform_frame_preserved = {k:v for k,v in isoform_preserved.iteritems() if v['rf_preserved']}

        return isoform_frame_preserved

    def compare_overlap_sj( self, obj_sj ):
        """
        Function: see if 'self' (current SJ) overlaps the entire or partial SpliceJunction 'obj_sj'. Will output a number that has the following meanings:
            - -1 = self & obj_sj do not reside 
        NOTE: I think the funciotn SpliceJunction.canonical_overlapping_sj() is similar to this function as it retrieves all canonical SJs that are overlapped by self
        """
        if self.chrom != obj_sj.chrom:
            return -1

        overlap_score = 0
        #check if starting point of other SJ is contained by self SJ
        if self.start < obj_sj.start < self.end:
            overlap_score += 1
        #check if ending point of other SJ is contained by self SJ
        if self.start < obj_sj.end < self.end:
            overlap_score += 2
        
        return 

    def compare_overlap_sj_v2( self, obj_sj, consider_most = False ):
        """
        see if 'self' (current SJ) overlaps the entire or partial SpliceJunction 'obj_sj'. Will output a number that has the following meanings:
            - -1 = self & obj_sj do not reside 

        Args:
            -obj_sj = a different SpliceJunction instance (not self). Variable 'obj_sj' will be used to determine if it is within the range of 'self'
            -consider_most = boolean where
                -True = Will consider self.leftmost & self.rightmost as endpoints to see if 'obj_sj' is contained 
                -False = will considered self.start & self.end as the endpoints to see if 'obj_sj' is contained
        Output:
            outputs an integer with the following meaning:
                - -1 = self & obj_sj do not reside on the same chromosome
                -0 = neither end of obj_sj is contained in self
                -1 = the start position of obj_sj is the range of self
                -2 = the end position of obj_sj is the range of self
                -3 = both start & end position of obj_sj is in the range of self
        NOTE: I think the funciotn SpliceJunction.canonical_overlapping_sj() is similar to this function as it retrieves all canonical SJs that are overlapped by self
        """
        if self.chrom != obj_sj.chrom:
            return -1

        #determine
        if consider_most:
            start = self.leftmost
            end = self.rightmost
        else:
            start = self.start
            end = self.end

        overlap_score = 0
        #check if starting point of other SJ is contained by self SJ
        if start <= obj_sj.start < end:
            overlap_score += 1
        #check if ending point of other SJ is contained by self SJ
        if start < obj_sj.end <= end:
            overlap_score += 2
        
        return overlap_score
        

    def canonical_overlapping_sj( self, bool_local = True ):
        """
        Args:
            bool_local = boolean where:
                -True = only look at elements that overlap genome_pos
                -False = find longest element that overlaps genome_pos, then find all elements
        Function: retrieves all canonical SJ overlapping the splice junction in "self", returns a hash of all canonical SJs overlapping the SJ of interest, where key = genomic range of element & value = Exon object
        """
        genome_pos = self.chrom + ':' + str( self.start ) + '-' + str( self.end )
        return self.find_overlapping_elements( genome_pos, self.gene_sym_oi, bool_local, True )     #bool_intron is True because I'm looking for overlapping canonical splice junctions, key = genomic range of element & value = Exon object (which are actually intronic elements that are overlapped as the introns are the boundaries of the splice site)

    def get_canon_splice_sites( self, position, which_ss = None ):
        """
        Args:
            position = integer that is the position of interest -> will use this to find the exon that contains this position
            which_ss = string that can have the following values:
                -None = keeps both splice sites (lower_ss & higher_ss, where ss = splice site).
                -'donor' = keeps the splice site on the donor side (this depends on strand sign)
                -'acceptor' = keeps the splice site on the acceptor side (this depends on strand sign)
        Function: retrieves all the canonical splice sites (either 5' or 3' splice site) for each isoform and returns the splice site position for each isoform
        """
        return self.get_exon_splice_sites_all_isoforms( position, which_ss )

    def sj_read_support( self, bam_reader, genomic_range = None ):
        """
        Args:
            bam_reader = HTSeq.BAM_Reader instance, used to quantify the number of reads that map to SJ genomic range 'genomic_range' (command: bam_reader = HTSeq.BAM_Reader(path_to_bam_file) )
            genomic_range = string in format chrom:start-end. If None, then uses the SJ position recorded in "self"
        Function: finds reads that support splice junctions by finding reads that uniquely map to splice junction position
        """
        if not genomic_range:
            genomic_range = self.chrom + ':' + str( self.start ) + '-' + str( self.end )

        hash_gr = Isoform.split_genome_pos( genomic_range )       #hash_gr = hash table of genomic range
        #these are the counts
        count = 0       #count all reads that map to "genomic_range"
        uniq_count = 0      #count all uniquely mapped reads to "genomic_range"
        align_score_max_count = 0       #this also counts all uniquely mapped reads to "genomic_range"
        for i, a in enumerate( bam_reader.fetch( region = genomic_range ) ):
            # score = 0
            # if a.optional_field( "NH" ) == 1:
            #     score += 1
            for cigop in a.cigar:
                if cigop.type == 'N' and cigop.ref_iv.start == hash_gr['start'] and cigop.ref_iv.end == hash_gr['end']:
                    ##TEST:: see the output of read
                    # print i, ": NH = ", a.optional_field( "NH" )
                    # print i, ": aligned = ", a.aligned, " & aQual = ", a.aQual, " & read quality = ", a.read.qual
                    # print i, ": cigop.ref_iv = ", cigop.ref_iv, " & chrom = ", cigop.ref_iv.chrom, " & start = ", cigop.ref_iv.start, " & end = ", cigop.ref_iv.end
                    # print i, ": cigop.size = ", cigop.size
                    # print i, ": cigop.type = ", cigop.type
                    # print '\n'
                    
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


    def sj_read_support_pysam( self, pysam_file, genomic_range = None, uniq_only = False ):
        """
        Args:
            pysam_file = pysam.AlignmentFile that opens up the mapped reads bam file (e.g. accepted_hits.bam)
            genomic_range = string that is the position of interest (format = chrom:start-end)
            uniq_only = boolean
                -True = will only quantify the uniquely-mapped gapped reads that map to 'genomic_range'
                -False = will quantify the uniquely & non-uniquely-mapped gapped reads that map to 'genomic_range'
        Function: this function retrieves gapped reads that supports range 'genomic_range'. NOTE that this function is much faster (maybe 6x faster) than SpliceJunction's def sj_read_support()
        """
        if not genomic_range:
            genomic_range = self.chrom + ':' + str( self.start ) + '-' + str( self.end )

        all_count = 0
        unique_count = 0
        # hash_gr = Isoform.split_genome_pos( genomic_range )
        hash_gr = Isoform.split_genome_pos( genomic_range )
        for read in pysam_file.fetch( hash_gr['chrom'], hash_gr['start'], hash_gr['end'] ):
            #this if if I only want to quantify uniquely-mapped reads -> slightly faster
            if uniq_only:
                if read.mapq != 50:
                    continue

                if not any(x for x in read.blocks if hash_gr['start'] in x):
                    continue
                if not any(x for x in read.blocks if hash_gr['end'] in x):
                    continue

                unique_count += 1
            else:       #this if if I only want to quantify non-unique & uniquely-mapped reads -> slightly slower
                if not any(x for x in read.blocks if hash_gr['start'] in x):
                    continue
                if not any(x for x in read.blocks if hash_gr['end'] in x):
                    continue
                all_count += 1

                if read.mapq != 50:
                    continue
                unique_count += 1

        # return [{'all_count': all_count, 'unique_count': unique_count}, hash_query_test, hash_query_pos]
        return {'all_count': all_count, 'unique_count': unique_count}


    @staticmethod
    def calc_sjpos(line):
        """
        Args:
            line = a string that is a row in the file, where each column is separated by '\t', else could be an array (as this is the case when retrieving with Tabix)
        Function:
            this function will calculate the end positions, start & end 
        Output:
            will return an array where [0] = true left position and [1] = true right position
        """
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

        #STEP: split line into an array, and calculate the true left & right positions
        arrRowCol = line.split('\t') if type( line ) is str else line
        #split the columns with 2 numbers
        arrBlockLen = arrRowCol[COLSJ_BLOCKLEN].split(",")
        arrBlockStart = arrRowCol[COLSJ_BLOCKSTART].split(",")
        #Find true left position  =  COLSJ_LEFTMOST + COLSJ_BLOCKLEN[0], where COLSJ_BLOCKLEN is 2 numbers separated by ",", where [0]  =  length of first block (see UCSC Genome Browser), [1]  =  length of second block (see UCSC Genome Browser)
        trueLeftPos = int( arrRowCol[COLSJ_LEFTMOST] ) + int( arrBlockLen[0] )
        #Find true right position  =  COLSJ_RIGHTMOST - COLSJ_BLOCKLEN[1], where COLSJ_BLOCKLEN is 2 numbers separated by ",", where [0]  =  first number, [1]  =  second number
        trueRightPos = int( arrRowCol[COLSJ_RIGHTMOST] ) - int( arrBlockLen[1] )

        return [ trueLeftPos, trueRightPos ]
