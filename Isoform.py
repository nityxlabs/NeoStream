#/usr/bin/python
import numpy as np
import tabix
import pysam

from Exon import Exon

exon_base = 1           #if 0, then exons are 0-based, if 1, then exons are 1-based

class Isoform( object ):
    obj_cruzdb = None
    #make exon range string ( format: chrom:start-end )
    make_key_range = staticmethod( lambda chrom, exon_start, exon_end : str(chrom) + ':' + str(exon_start) + '-' + str(exon_end) )

    def __init__( self, db_type, isoform_id, hash_pos = None ):
        """
        Args:
            -isoform_id = string that is the isoform ID of interest. This isoform_id can
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
            -hash_pos = hash that will be used to determine the closest position of interest. It will have 2 keys:
                -"chrom" = string that is the chromosome of interest (format: chr# (like chr9, chr11))
                -"pos_oi" = integer that suggest which isoform to select by selecting the isoform version that is closest to this position 
        Class: this class will retrieve all information for an associated isoform
        """

        # isoform_info = Isoform.get_isoform( isoform_id )
        info_isoform = None
        if hash_pos:
            # isoform_variants = Isoform.obj_cruzdb.refGene.filter_by( name = isoform_id ).all()
            isoform_variants = self.get_isoforms_db_all( isoform_id, db_type )
            info_isoform = Isoform.closest_isoform_variant_isoform( isoform_id, db_type, hash_pos['chrom'], hash_pos['pos_oi'] )

        #if nothing assigned to info_isoform (either by hash_pos = None OR nothing found for 'info_isoform')
        if not info_isoform:
            # info_isoform = Isoform.obj_cruzdb.refGene.filter_by( name = isoform_id ).first()
            info_isoform = self.get_isoforms_db_first( isoform_id, db_type )

        self.hash_pos = hash_pos            #has 2 elements, "chrom" & "pos"
        self.db_type = db_type
        self.isoform_id = str( info_isoform.name )
        self.gene_sym = str( info_isoform.name2 )
        self.is_coding = str( info_isoform.is_coding )
        # self.chrom = str( info_isoform.chrom )
        self.chrom = str( info_isoform.chrom ) if not '_' in info_isoform.chrom else str( info_isoform.chrom.split('_')[0] )

        convert_ss = { '+' : 1, '-' : -1 }
        self.strand = convert_ss.get( info_isoform.strand, 0 )
        self.boundary = ( info_isoform.txStart, info_isoform.txEnd )       #tuple that records the gene boundary
        #get boundary for coding region -> this may not exist if self.is_coding is False
        if info_isoform.cds:
            cds_start = min( [x[0] for x in info_isoform.cds] )
            cds_end = max( [x[1] for x in info_isoform.cds] )
            self.coding_boundary = ( cds_start, cds_end )
        else:
            self.coding_boundary = ( None, None )
        self.tss = info_isoform.tss()       #a tuple with 2 elements, the position of the transcription start site. Usually both positions are the same number. NOTE: this is the NOT the start of coding region, just the start of trancription of the mRNA seqeunce --> self.coding_boundary gives the start & end position of the coding region!
        self.pos_start_codon = None         #this will stay None for non-coding genes (e.g. NR_, anything with blank self.cds)
        if info_isoform.cds:
            self.pos_start_codon = info_isoform.cds[-1][1] if self.strand < 0 else info_isoform.cds[0][0]

        self.utr_3 = info_isoform.utr3      #a tuple with 2 elements, the start & end position of the 3' UTR 
        self.utr_5 = info_isoform.utr5      #a tuple with 2 elements, the start & end position of the 5' UTR 

        # self.hashExonList = self.organize_exons( info_isoform )     #key = string that is exon range (chrom:start-end), value = Exon object

        self.hashExonList = self.organize_features( info_isoform, True )     #key = string that is exon range (chrom:start-end), value = Exon object
        self.hashIntronList = self.organize_features( info_isoform, False )     #key = string that is exon range (chrom:start-end), value = Exon object

        # self.last_exon_num = len( self.hashExonList.keys() ) - 1      #for 0-based exons, -1 because exon starts at 0, not 1
        self.last_exon_num = len( self.hashExonList.keys() )            #for 1-based exons (i.e. first exon is exon 1, NOT exon 0) 

        self.donor_acceptor_sites = self.get_donor_acceptor_sites()

    @staticmethod
    def get_isoforms_db_first( isoform_id, db_type ):
        """
        Returns the first isoforms associated with specific isoform_id
        Args:
            -isoform_id = string that is the isoform ID of interest. This isoform_id can
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        """
        if db_type == 2:
            isoform_variant = Isoform.obj_cruzdb.ensGene.filter_by( name = isoform_id ).first()
        elif db_type == 3:
            isoform_variant = Isoform.obj_cruzdb.knownGene.filter_by( name = isoform_id ).first()
        elif db_type == 4:
            if not "ENST" in isoform_id and not "." in isoform_id:
                print "Error in Isoform.get_isoforms_db_first() - Gene isoform not in GENCODE format (i.e. no 'ENST' and no '.', e.g. ENST00000319349.5)"
            isoform_variant = Isoform.obj_cruzdb.wgEncodeGencodeBasicV19.filter_by( name = isoform_id ).first()
        else:
            isoform_variant = Isoform.obj_cruzdb.refGene.filter_by( name = isoform_id ).first()

        return isoform_variant

    @staticmethod
    def get_isoforms_db_all( isoform_id, db_type ):
        """
        Returns all isoforms associated with specific isoform_id
        Args:
            -isoform_id = string that is the isoform ID of interest. This isoform_id can
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        """
        #NOTE: for db_type = 4, the isoform ID needs to be in Ensembl form, and even then it needs a "." as I think this signifies some sort of version for the gene.
        if db_type == 2:
            isoform_variants = Isoform.obj_cruzdb.ensGene.filter_by( name = isoform_id ).all()
        elif db_type == 3:
            isoform_variants = Isoform.obj_cruzdb.knownGene.filter_by( name = isoform_id ).all()
        elif db_type == 4:
            if not "ENST" in isoform_id and not "." in isoform_id:
                print "Error in Isoform.get_isoforms_db_first() - Gene isoform not in GENCODE format (i.e. no 'ENST' and no '.', e.g. ENST00000319349.5)"
            isoform_variants = Isoform.obj_cruzdb.wgEncodeGencodeBasicV19.filter_by( name = isoform_id ).all()
        else:
            isoform_variants = Isoform.obj_cruzdb.refGene.filter_by( name = isoform_id ).all()

        #retrieve the all the isoforms associated with gene symbol -> I found that I can retrieve more isoforms when using the gene symbol as a filter as opposed to using .bin_query
        return Isoform.get_all_isoforms( isoform_variants, db_type )


    @staticmethod
    def get_isoforms_by_pos_db_all( chrom, start, end, db_type ):
        """
        Function: based on the position, will find and return the closest isoform ID. Returns hash that contains gene symbol, isoform ID, and distance

        Args:
            -isoform_id = string that is the isoform ID
            -chrom = string in the format 'chr#' (chr9, chr12)
            -start & end = integer that is position to where other positions will be looked for
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        """
        if db_type == 2:
            isoform_variants = Isoform.obj_cruzdb.bin_query( 'ensGene', chrom, start, end ).all()
        elif db_type == 3:
            isoform_variants = Isoform.obj_cruzdb.bin_query( 'knownGene', chrom, start, end ).all()
        elif db_type == 4:
            isoform_variants = Isoform.obj_cruzdb.bin_query( 'wgEncodeGencodeBasicV19', chrom, start, end ).all()
        else: 
            isoform_variants = Isoform.obj_cruzdb.bin_query( 'refGene', chrom, start, end ).all()

        #retrieve the all the isoforms associated with gene symbol -> I found that I can retrieve more isoforms when using the gene symbol as a filter as opposed to using .bin_query
        return Isoform.get_all_isoforms( isoform_variants, db_type )

    @staticmethod
    def get_all_isoforms( isoform_variants, db_type ):
        """
        get all isoforms based on "snapshot" of isoforms in 'isoform_variants'
        Args:
            -isoform_variants = information from cruzdb, usually from get_isoforms_db_all() or get_isoforms_by_pos_db_all()
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        """
        if not isoform_variants:
            return None

        #retrieve the all the isoforms associated with gene symbol -> I found that I can retrieve more isoforms when using the gene symbol as a filter as opposed to using .bin_query
        list_gene_sym = np.unique( [x.name for x in isoform_variants] ) if db_type == 3 else np.unique( [x.name2 for x in isoform_variants] )
        gene_info = []
        for each_sym in list_gene_sym:
            add_info = Isoform.get_gene_sym_db_all( each_sym, db_type )
            if add_info:
                gene_info += add_info

        ##TEST:: print "\t\tIsoform.GAI - gene_info = ", gene_info

        return gene_info

    @staticmethod
    def get_gene_sym_db_all( gene_sym, db_type ):
        """
        Returns all isoforms associated with specific gene symbol. NOTE: for UCSC, there is no "gene symbol", so need to use the same ID as would be used for isoform
        Args:
            -gene_sym = string that is the gene symbol (e.g. BRAF, RAF1)
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene) --> e.g. ETV4, BRAF, AKT1
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene) --> needs to being with "ENSG"
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene) --> e.g. 'uc010tyk.2', 'uc002idv.4'
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....) --> can use the same gene symbol as RefSeq
        """
        if db_type == 2:
            if not "ENSG" in gene_sym:
                print "Error in Isoform.get_gene_sym_db_all() - Gene symbol not in Ensembl format (i.e. no 'ENSG')"
            isoform_variants = Isoform.obj_cruzdb.ensGene.filter_by( name2 = gene_sym ).all()
        elif db_type == 3:
            if not "uc" in gene_sym:
                print "Error in Isoform.get_gene_sym_db_all() - Gene symbol not in UCSC format (i.e. no 'uc')"
            isoform_variants = Isoform.obj_cruzdb.knownGene.filter_by( name = gene_sym ).all()
        elif db_type == 4:
            isoform_variants = Isoform.obj_cruzdb.wgEncodeGencodeBasicV19.filter_by( name2 = gene_sym ).all()
        else:
            isoform_variants = Isoform.obj_cruzdb.refGene.filter_by( name2 = gene_sym ).all()

        return isoform_variants

    @staticmethod
    def closest_isoform_variant_pos( chrom, start, end, db_type ):
        """
        Args:
            -isoform_id = string that is the isoform ID
            -chrom = string in the format 'chr#' (chr9, chr12)
            -start & end = integer that is position to where other positions will be looked for
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        Function: based on the position, will find and return the closest isoform ID. Returns hash that contains gene symbol, isoform ID, and distance
        """
        isoform_variants = Isoform.get_isoforms_by_pos_db_all( chrom, start, end, db_type )

        get_closest_variant = Isoform.closest_isoform_variant_process( isoform_variants, chrom, start, db_type )
        return get_closest_variant


    @classmethod
    def is_obj_possible( cls_obj, isoform_id, pos_range, db_type ):
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
        hash_genomic_pos = Isoform.split_genome_pos( pos_range )
        # isoform_variants = Isoform.obj_cruzdb.refGene.filter_by( name = isoform_id ).all()    ##DELETE THIS LINE - this is inflexible as it only access RefSeq instead of being open to different databases
        isoform_variants = Isoform.get_isoforms_db_all( isoform_id, db_type )
        #check starting position
        get_isoform_variant = Isoform.closest_isoform_variant_check( isoform_variants, hash_genomic_pos['chrom'], hash_genomic_pos['start'], db_type )
        if get_isoform_variant:
            return {'chrom': hash_genomic_pos['chrom'], 'pos_oi': hash_genomic_pos['start']}
        
        #Else, if starting position doesn't work, check ending position
        get_isoform_variant = Isoform.closest_isoform_variant_check( isoform_variants, hash_genomic_pos['chrom'], hash_genomic_pos['end'], db_type )

        if get_isoform_variant:
            return {'chrom': hash_genomic_pos['chrom'], 'pos_oi': hash_genomic_pos['end']}

        #if no position is found to be in the confines of the isoform (start to end position of isoform), then return an empty hash, which means it is not possible to make an IsoformSJ instance as the SJ is out of the isoform's range
        return {}

    @classmethod
    def closest_isoform_variant_isoform( cls_obj, isoform_id, db_type, chrom, pos_oi ):
        """
        Args:
            isoform_id = string that is the isoform ID
            chrom = string in the format 'chr#' (chr9, chr12)
            pos_oi = integer that is position to where other positions will be looked for
        Function: finds the isoform that is the closest to returns the isoform that is closest
        """
        # #depending on self.db_type, retrieve annotations from either RefSeq or Ensembl
        # if db_type == 2:        #retrieve annotations from Ensembl
        #     # isoform_variants = cls_obj.get_ensembl_gene_obj( isoform_id )
        #     #need to retrieve gene symbol based on isoform ID
        #     isoform_variants = cls_obj.get_ensembl_isoforms( isoform_id )       
        #     #if no isoform variant is found, then return None
        #     if not isoform_variants:
        #         return None
        # else:       #default --> RefSeq
        #     isoform_variants = cls_obj.obj_cruzdb.refGene.filter_by( name = isoform_id ).all()

        isoform_variants = cls_obj.get_isoforms_db_all( isoform_id, db_type )

        return cls_obj.closest_isoform_variant_process( isoform_variants, chrom, pos_oi, db_type )

    @classmethod
    def closest_isoform_variant_process( cls_obj, isoform_variants, chrom, pos_oi, db_type ):
        """
        Args:
            -isoform_variants = array of all isoforms found using cruzdb.Genome.bin_query or cruzdb.Genome.refGene.filter_by
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        Function: this will find the isoform variant with the closest distance to integer position "pos_oi"
        """
        get_isoform_variant = cls_obj.closest_isoform_variant_check( isoform_variants, chrom, pos_oi, db_type )

        if get_isoform_variant:
            return get_isoform_variant
        else:       #else if no elements are found within the position of interest, then return the first instance
            try:
                return isoform_variants[0]
            except IndexError:      #Error meaning no isoform associated with position 'chrom:pos_oi'
                # print "Error with Isoform ", chrom, ":", pos_oi, " - no isoform found"
                return None

    @classmethod
    def closest_isoform_variant_check( cls_obj, isoform_variants, chrom, pos_oi, db_type ):
        """
        Args:
            isoform_variants = array of all isoforms found using cruzdb.Genome.bin_query or cruzdb.Genome.refGene.filter_by
            -pos_oi = position of interest to see which isoform variant in 'isoform_variants' is closest to this position
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        Function: this will find the isoform variant with the closest distance to integer position "pos_oi". Use this function in conjunction with Isoform.closest_isoform_variant_process().
        Output: Will output a specific refGene object from cruzdb (or ensGene, depends on gene_db) that is closest to 'chrom:pos_oi'. However, if position 'chrom:pos_oi' is not located in any of the isoforms, then return 
            -NOTE: This function only returns one of the isoforms that is considered the "closest", even if there are multiple isoforms that are equally as close
        """
        ##CAN PROBABLY DELETE THIS
        # #depending on "gene_db", retrieve annotations from either RefSeq or Ensembl
        # if gene_db == 2:        #retrieve annotations from Ensembl
        #     # isoform_variants = cls_obj.get_ensembl_gene_obj( isoform_id )
        #     #need to retrieve gene symbol based on isoform ID
        #     isoform_variants = cls_obj.get_ensembl_isoforms( isoform_id )       
        #     #if no isoform variant is found, then return None
        #     if not isoform_variants:
        #         return None
        # else:       #default --> RefSeq
        #     isoform_variants = cls_obj.obj_cruzdb.refGene.filter_by( name = isoform_id ).all()


        hash_closest_dist = {}      #k = array index (from loop's enumerate), v = minimum distance
        hash_gene_sym = {}          #k = array index (from loop's enumerate), v = gene symbol
        for i, each_isoform in enumerate( isoform_variants ):
            #retrieve chromosome - if not the same chromosome then skip
            iso_chrom = each_isoform.chrom if not '_' in each_isoform.chrom else each_isoform.chrom.split('_')[0]
            if chrom != iso_chrom:
                continue
            #if position is not in isoform, then skip
            if not pos_oi in range( each_isoform.start, each_isoform.end ):     #can also use each_isoform.txStart & each_isoform.txEnd - same thing
                continue

            #flatten all the exons, and find the smallest distance between "pos_oi" and all the exon ends
            flatten_exons = list( sum(each_isoform.exons, ()) )
            min_diff = min( [abs(x - pos_oi) for x in flatten_exons] )
            hash_closest_dist[i] = min_diff

            #record gene symbol - DID I DO THIS RIGHT FOR ENSEMBL???
            hash_gene_sym[i] = each_isoform.name if db_type == 3 else each_isoform.name2

        #find the isoform with smallest distance
        if hash_closest_dist:
            isoform_index = min( hash_closest_dist, key = hash_closest_dist.get )
            return isoform_variants[isoform_index]
        else:
            return None

    
    def __str__( self ):
        pos = self.chrom + ':' + str( self.boundary[0] ) + '-' + str( self.boundary[1] )
        strand_sign = '+' if self.strand > 0 else '-'
        return self.isoform_id + " (" + self.gene_sym + "): " + pos + ", (" + strand_sign + ")"

    def display_exon_list( self ):
        """
        Function: displays the exons in the isoform, returns string that gives exon number & the genomic range (e.g. chrom:start-end)
        """
        exon_info = {}      #hash where key = exon number & value = Exon Object
        for k,v in self.hashExonList.iteritems():       #key = string that is exon range (chrom:start-end), value = Exon object
            exon_info[ v.exonNum ] = v

        exon_order = exon_info.keys()
        exon_order.sort( reverse = False )

        str_exons = ""
        for i in exon_order:
            str_exons += "exon " + str( i ) + ": " + exon_info[i].str_genomic_pos() + "\n"

        # str_exons = ""
        # for i in range( 1, len( self.hashExonList.keys() ) + 1 ):
        #     str_exons += "exon " + str( i ) + ": ", 

        # for k,v in self.hashExonList.iteritems():       #key = string that is exon range (chrom:start-end), value = Exon object
        #     str_exons += "exon " + str( v.exonNum ) + ": " + k + "\n"

        return str_exons

    def json_isoform_info( self ):
        """
        Function: returns a hash that contains information about the isoform - gene name, isoform id, genomic position, exons & introns, etc. 
        """
        info = {}
        info['isoform_id'] = self.isoform_id
        info['gene_sym'] = self.gene_sym
        info['chrom'] = self.chrom
        info['pos_start'] = self.boundary[0]
        info['pos_end'] = self.boundary[1]
        info['strand_sign'] = self.strand
        info['last_exon_num'] = self.last_exon_num

        info['features'] = {}

        #get all the exons in the isoform
        for i in range(exon_base, self.last_exon_num):
            obj_feat_num = self.get_exon_num( i )
            obj_feat = {
            'type': 'exon',
            'chrom': obj_feat_num.chrom,
            'start': obj_feat_num.exonPos.location.start,
            'end': obj_feat_num.exonPos.location.end,
            'strand': obj_feat_num.exonPos.strand,
            'featureNum': obj_feat_num.exonNum
            }

            #add object to hash of features
            feat_name = 'exon' + str( i )
            info['features'][feat_name] = obj_feat

        #get all the introns in the isoform
        for i in range(exon_base, self.last_exon_num - 1):
            obj_feat_num = self.get_intron_num( i )
            obj_feat = {
            'type': 'intron',
            'chrom': obj_feat_num.chrom,
            'start': obj_feat_num.exonPos.location.start,
            'end': obj_feat_num.exonPos.location.end,
            'strand': obj_feat_num.exonPos.strand,
            'featureNum': obj_feat_num.exonNum
            }

            #add object to hash of features
            feat_name = 'intron' + str( i )
            info['features'][feat_name] = obj_feat

        return info


    """
    Function: organizing genomic information
    """
    def organize_features( self, isoform, bool_exon = True ):
        """
        Args:
            isoform = cruzdb object that contains information for a specific isoform
        Function: this will organize the feats, cds, & reading frame for a specific gene
        """
        #hash_feats = will contain feats sorted by start position (key = index, value = feat)
        hash_feats = {}       #key = string that is feat range (chrom:start-end), value = feat object

        #get reading frames
        

        ##TEST:: print "organize_feats: isoform = ", isoform.feats
        if bool_exon:
            list_features = isoform.exons
            feat_frames = map( int, [x for x in isoform.exonFrames.split(',') if x] )
        else:
            list_features = isoform.introns

        for i, feat in enumerate( list_features ):      #i = feat number
            #check to see if the key_range exists 
            key_range = Isoform.make_key_range( self.chrom, feat[0], feat[1] )

            #calculate the feat number
            # feat_num = i if self.strand == 1 else ( len(list_features) - i - 1 )        #for 0-based exons, use this
            # feat_num = i + 1 if self.strand == 1 else ( len(list_features) - i )        #for 1-based exons, use this
            feat_num = i + exon_base if self.strand == 1 else ( len(list_features) - i - ( exon_base - 1 ) )        #this handles both 0-based & 1-based

            #NOTE: UCSC has 0-based genome, meaning the first position of the feat is actually the last position in the previous intron, that is why I add '+1'
            if bool_exon:
                feat_info = self.get_feat_info( feat[0], feat[1], feat_frames[i], feat_num, True, 'exon' )
            else:
                feat_info = self.get_feat_info( feat[0], feat[1], None, feat_num, True, 'intron' )
            hash_feats[ key_range ] = Exon( feat_info )

        #go through all cds, assign CDS to each feat
        if bool_exon:
            for each_cds in isoform.cds:     #each_cd = tuple where [0] = start position & [1] = end position
                hash_feats = self.organize_exons_cds( each_cds[0], each_cds[1], hash_feats )

        return hash_feats


    
    ##MAY DELETE THIS BECAUSE NOW I HAVE - organize_features
    # def organize_exons( self, isoform ):
    #     """
    #     Args:
    #         isoform = cruzdb object that contains information for a specific isoform
    #     Function: this will organize the exons, cds, & reading frame for a specific gene
    #     """
    #     #hashExonList = will contain exons sorted by start position (key = index, value = Exon)
    #     hashExonList = {}       #key = string that is exon range (chrom:start-end), value = Exon object

    #     #get reading frames
    #     exon_frames = map( int, [x for x in isoform.exonFrames.split(',') if x] )

    #     ##TEST:: print "organize_exons: isoform = ", isoform.exons

    #     for i, exon in enumerate( isoform.exons ):      #i = exon number
    #         #check to see if the key_range exists 
    #         key_range = Isoform.make_key_range( self.chrom, exon[0], exon[1] )

    #         #calculate the exon number
    #         exon_num = i if self.strand == 1 else ( len(isoform.exons) - i - 1 )

    #         #NOTE: UCSC has 0-based genome, meaning the first position of the exon is actually the last position in the previous intron, that is why I add '+1'
    #         exon_info = self.get_feat_info( exon[0], exon[1], exon_frames[i], exon_num, True )
    #         hashExonList[ key_range ] = Exon( exon_info )

    #     #go through all cds, assign CDS to each exon
    #     for each_cds in isoform.cds:     #each_cd = tuple where [0] = start position & [1] = end position
    #         hashExonList = self.organize_exons_cds( each_cds[0], each_cds[1], hashExonList )

    #     return hashExonList


    def organize_exons_cds( self, start_cds, end_cds, hashExonList ):
        """
        Args:
            start_cds & end_cds = the starting (lower numerical position) & ending CDS positions (higher numerical position) that is fused by the splice junction.
            tuple_cds = tuple where [0] = start position of cds & [1] = 
            hashExonList = hash that records the exon for the gene, where key = string that is exon range (chrom:start-end), value = Exon object
        Function: maps the cds exons to the original exons in hashExonList
        """
        for k,v in hashExonList.iteritems():        #key = exon range (chrom:start-end), value = Exon object
            v.add_cds( start_cds, end_cds )

        return hashExonList

    """
    Function: get exons based on exon number
    """
    def get_exon_num( self, exon_num ):
        """ 
        Args:
            exon_num = integer that is the exon number of interest
        Function: returns Exon object with exon number
        """
        for k,v in self.hashExonList.iteritems():       #k = string (chrom:start-end), v = Exon object
            if v.exonNum == exon_num:
                return v

        #if no exon found with that exon number, then return None
        return None

    def get_intron_num( self, intron_num ):
        """ 
        Args:
            intron_num = integer that is the 
        Function: returns intron object with intron number
        """
        for k,v in self.hashIntronList.iteritems():       #k = string (chrom:start-end), v = Exon object
            if v.exonNum == intron_num:
                return v

        #if no intron found with that intron number, then return None
        return None

    #MAY DELETE THIS - DO NOT NEED THIS BECAUSE I HAVE def get_exon_splice_sites()
    # def get_canon_splice_site( self, position, bool_donor ):
    #     """
    #     Args:
    #         position = integer that is the position of interest -> will use this to find the exon that contains this position
    #         bool_donor = boolean
    #             -True = will retrieve the 5' splice site (aka splice donor)
    #             -False = will retrieve the 3' splice site (aka splice acceptor)
    #     Function: this will retrieve the canonical splice sites associated with an exon (or intron), based on the parameter 'position' given.
    #     """
    #     #find which exons contain this position
    #     get_elem = self.get_element_v2( position, 1 )
    #     if get_elem:
    #         if bool_donor:
    #             get_pos = get_elem.exonPos.location.end if self.strand > 0 else get_elem.exonPos.location.start
    #         else:
    #             get_pos = get_elem.exonPos.location.start if self.strand > 0 else get_elem.exonPos.location.end
    #     else:       #else if position is not found in an exon, then look into introns
    #         get_elem = self.get_element_v2( position, 3 )
    #         #if nothing found, then return None
    #         if not get_elem:
    #             return None

    #         #retrieve the splice site position from the intron
    #         if bool_donor:      #retrieve the 5' splice site (splice donor)
    #             get_pos = get_elem.exonPos.location.start if self.strand > 0 else get_elem.exonPos.location.end
    #         else:               #retrieve the 3' splice site (splice acceptor)
    #             get_pos = get_elem.exonPos.location.end if self.strand > 0 else get_elem.exonPos.location.start

    #     #return the splice site position
    #     return get_pos


    def get_exon_splice_sites( self, position ):
        """
        Args:
            position = integer that is a genomic position
        Function: retrieves the splice sites adjacent to the exon of interest, returns a hash that contains 2 elements, lower_ss & higher_ss, where each value is an Intron object (instance of class Exon) that refers the splice site, else if no splice site is found, then value is None
        """
        #find which exons contain this position
        get_elem = self.get_element_v2( position, 1 )
        if get_elem:        #this means the position lands within an exon
            #if position lands in first exon - there will not be splice site before the exon
            if get_elem.exonNum == exon_base:
                higher_ss = self.get_intron_num( get_elem.exonNum )     #higher_ss = higher splice site
                hash_ss = {'lower_ss': None, 'higher_ss': higher_ss }   #hash_ss = hash splice site
            #if position lands in last exon - there will not be splice site after the exon
            elif get_elem.exonNum == self.last_exon_num:
                lower_ss = self.get_intron_num( get_elem.exonNum - 1 )  #lower_ss = lower splice site
                hash_ss = {'lower_ss': lower_ss, 'higher_ss': None }    #hash_ss = hash splice site
            #else position lands in one of the interior exons
            else:
                lower_ss = self.get_intron_num( get_elem.exonNum - 1 )  #lower_ss = lower splice site
                higher_ss = self.get_intron_num( get_elem.exonNum )     #higher_ss = higher splice site
                hash_ss = {'lower_ss': lower_ss, 'higher_ss': higher_ss }   #hash_ss = hash splice site
        #else the position may have landed in the intron
        else:           #else this means look into introns to find where the position is located
            get_elem = self.get_element_v2( position, 3 )
            if not get_elem:
                hash_ss = {'lower_ss': None, 'higher_ss': None }        #hash_ss = hash splice site
            else:
                hash_ss = {'lower_ss': get_elem, 'higher_ss': get_elem }    #hash_ss = hash splice site. I recorded an element for both so it will be selected for when sifting for acceptor or donor splice sites (as can be seen in MultiIsoform's def get_exon_splice_sites_all_isoforms() )


        ##TEST::
        # print "Isoform.gess - isoform = ", self.isoform_id, " & position = ", position, " & strand = ", self.strand
        # for k,v in hash_ss.iteritems():
        #     print "k = ", k, " & v = ", v

        #return the hash splice sites
        return hash_ss


    """
    Function Start: finding CDS that contain position
    """
    def get_element_v2( self, position, feat_search = 1, direction = 0, bool_strand = True, rel_pos = None ):
        """
        Args:
            position = integer that is the genomic position of interest
            feat_search = integer where:
                -2 = retrieve cds features
                -3 = retrieve intronic features
                -else if any other number (0,1) = retrieve exon features
            direction = integer that is the direction to find the next exon
                - 0 = direction does not matter
                - 1 = find next exon to the right (highest position - higher exon number)
                - -1 = find next exon to the left (lowest position - lower exon number)
            bool_strand = boolean (Note: direction = 0 will not pick a specific exon based on direction )
                -True = will consider strand sign to find exon
                    -if direction = 1, then will find the higher exon number (higher position for + & lower position for -)
                    -if direction = -1, then will find lower exon number (lower position for + & higher position for -)
                -False = will consider finding the higher-positioned exon regardless of exon number & strand sign
                    -if direction = 1, then will just find exon with higher position
                    -if direction = -1, then will just find exon will lower position
            rel_pos = string that is relative position of the element. Can be the following values: (look at Exon.in_exon() for appropriate labels)
                -None = position can be anywhere within exon
                -"withinElem" = found within the exon
                -"exonLeft" = left side of exon 
                -"exonRight" = right side of exon
        Function: retrieve exons from all isoforms, returns the first Exon object that contains integer 'position' 
        Assumption: assumes that exons are not overlapping (therefore only 1 exon, if any, will contain position)
        Output: Returns an instance of the class "Exon" (could be an exon or intron object)
        """
        #retrieve the element of interest
        if feat_search == 2:
            elems = self.get_cds( position, rel_pos )
        elif feat_search == 3:
            elems = self.get_intron( position, rel_pos )
        else:
            elems = self.get_exon( position, rel_pos )


        #if considering gene strand, then configure direction to take into account the strand sign
        if bool_strand:
            direction *= self.strand

        #determine which exon to return based on if exons are found & the direction
        if not elems:
            elem = None
        elif direction > 0:      #find the elem with the highest numerical position
            elem = max( elems, key = lambda x: x.exonPos.location.start )
        elif direction < 0:     #find the elem with the lowest numerical position
            elem = min( elems, key = lambda x: x.exonPos.location.end )
        else:
            elem = elems[0]

        ##TEST:: print "Isoform.get_element: gene = ", self.gene_sym," & position = ", position, " & elems = ", elems

        return elem

    ##NEED TO DELETE THIS --> replaced by "def get_element_v2()"
    def get_element( self, position, bool_exon = True, direction = 0, bool_strand = True, rel_pos = None ):
        """
        Args:
            position = integer that is the genomic position of interest
            bool_exon = boolean
                -True = retrieve exon features
                -False = retrieve cds features
            direction = integer that is the direction to find the next exon
                - 0 = direction does not matter
                - 1 = find next exon to the right (highest position - higher exon number)
                - -1 = find next exon to the left (lowest position - lower exon number)
            bool_strand = boolean (Note: direction = 0 will not pick a specific exon based on direction )
                -True = will consider strand sign to find exon
                    -if direction = 1, then will find the higher exon number (higher position for + & lower position for -)
                    -if direction = -1, then will find lower exon number (lower position for + & higher position for -)
                -False = will consider finding the higher-positioned exon regardless of exon number & strand sign
                    -if direction = 1, then will just find exon with higher position
                    -if direction = -1, then will just find exon will lower position
            rel_pos = string that is relative position of the element. Can be the following values: (look at Exon.in_exon() for appropriate labels)
                -None = position can be anywhere within exon
                -"withinElem" = found within the exon
                -"exonLeft" = left side of exon 
                -"exonRight" = right side of exon
        Function: retrieve exons from all isoforms, returns the first Exon object that contains integer 'position' 
        Assumption: assumes that exons are not overlapping (therefore only 1 exon, if any, will contain position)
        """
        elems = self.get_exon( position, rel_pos ) if bool_exon else self.get_cds( position, rel_pos )

        #if considering gene strand, then configure direction to take into account the strand sign
        if bool_strand:
            direction *= self.strand

        #determine which exon to return based on if exons are found & the direction
        if not elems:
            elem = None
        elif direction > 0:      #find the elem with the highest numerical position
            elem = max( elems, key = lambda x: x.exonPos.location.start )
        elif direction < 0:     #find the elem with the lowest numerical position
            elem = min( elems, key = lambda x: x.exonPos.location.end )
        else:
            elem = elems[0]

        ##TEST:: print "Isoform.get_element: gene = ", self.gene_sym," & position = ", position, " & elems = ", elems

        return elem

    def get_exon( self, position, rel_pos = None ):
        """
        Args:
            position = integer that is the genomic position of interest
            rel_pos = string that is relative position of the element. Can be the following values: (look at Exon.in_exon() for appropriate labels)
                -None = position can be anywhere within exon
                -"withinElem" = found within the exon
                -"exonLeft" = left side of exon 
                -"exonRight" = right side of exon
        Function: retrieve exons from all isoforms, returns the first Exon object that contains integer 'position' 
        Assumption: assumes that exons are not overlapping (therefore only 1 exon, if any, will contain position)
        """
        if not rel_pos:
            exons = [self.hashExonList[k] for k in self.hashExonList if self.hashExonList[k].in_exon( position )]
        else:
            exons = [self.hashExonList[k] for k in self.hashExonList if self.hashExonList[k].in_exon( position ) == rel_pos]

        #Assumption: that the exons are not overlapping, therefore only 1 exon will be retrieved
        # if len( exon ) > 1:
        #     print "!!Isoform.get_exon Warning: More than 1 exon found. Returning only first occurrence!!"

        return exons

    ##MAY NEED TO DELETE
    # def get_exon_backup( self, position, rel_pos = None ):
    #     """
    #     Args:
    #         position = integer that is the genomic position of interest
    #         rel_pos = string that is relative position of the element. Can be the following values:
    #             -None = position can be anywhere within exon
    #             -"exonLeft" = left side of exon 
    #             -"exonRight" = right side of exon
    #     Function: retrieve exons from all isoforms, returns the first Exon object that contains integer 'position' 
    #     Assumption: assumes that exons are not overlapping (therefore only 1 exon, if any, will contain position)
    #     """
    #     if not rel_pos:
    #         exon = [self.hashExonList[k] for k in self.hashExonList if self.hashExonList[k].in_exon( position )]
    #     else:
    #         exon = [self.hashExonList[k] for k in self.hashExonList if self.hashExonList[k].in_exon( position ) == rel_pos]

    #     #Assumption: that the exons are not overlapping, therefore only 1 exon will be retrieved
    #     # if len( exon ) > 1:
    #     #     print "!!Isoform.get_exon Warning: More than 1 exon found. Returning only first occurrence!!"

    #     return exon[0] if exon else None


    def get_cds( self, position, rel_pos = None ):
        """
        Args:
            position = integer that is the genomic position of interest
            rel_pos = string that is relative position of the element. Can be the following values:
                -None = position can be anywhere within exon
                -"exonLeft" = left side of exon 
                -"exonRight" = right side of exon
        Function: retrieve exons from all isoforms, returns an array of Exon objects that contain integer 'position'
        Assumption: assumes that exons are not overlapping (therefore only 1 exon, if any, will contain position)
        """

        #NOTE: cdss is just the plural of CDS (Coding Sequence)
        if not rel_pos:
            cdss = [self.hashExonList[k] for k in self.hashExonList if self.hashExonList[k].in_cds( position )]
        else:
            cdss = [self.hashExonList[k] for k in self.hashExonList if self.hashExonList[k].in_cds( position ) == rel_pos]
        
        # #Assumption: that the exons are not overlapping, therefore only 1 exon will be retrieved
        # if len( cds ) > 1:
        #     print "!!Isoform.get_cds Warning: More than 1 cds found. Returning only first occurrence!!"

        return cdss


    #MAYBE DELETE THIS AS THIS DOESN'T USE "self.hashIntronList"
    # def get_intron( self, position ):
    #     """ Function: determines if position is within intron """
    #     #as the donor-acceptor sites span the intron, use it to see if position is within introns
    #     da_sites = self.get_donor_acceptor_sites()
    #     intron_num = None
    #     for i, intron in enumerate( da_sites ):
    #         if intron[0] < position <= intron[1]:
    #             intron_num = i
    #             break

    #     return intron_num

    def get_intron( self, position, rel_pos = None ):
        """
        Args:
            position = integer that is the genomic position of interest
            rel_pos = string that is relative position of the element. Can be the following values: (look at Exon.in_intron() for appropriate labels)
                -None = position can be anywhere within intron
                -"withinElem" = found within the intron
                -"exonLeft" = left side of intron 
                -"exonRight" = right side of intron
        Function: retrieve introns from all isoforms, returns the first Exon object that contains integer 'position' 
        Assumption: assumes that introns are not overlapping (therefore only 1 intron, if any, will contain position)
        """
        #NOTE: cdss is just the plural of CDS (Coding Sequence)
        if not rel_pos:
            obj_intron = [self.hashIntronList[k] for k in self.hashIntronList if self.hashIntronList[k].in_exon( position )]
        else:
            obj_intron = [self.hashIntronList[k] for k in self.hashIntronList if self.hashIntronList[k].in_exon( position ) == rel_pos]


        ##TEST:: print "Isoform.get_intron: gene = ", self.gene_sym," & position = ", position, " & obj_intron = ", obj_intron
        
        # #Assumption: that the exons are not overlapping, therefore only 1 exon will be retrieved
        # if len( cds ) > 1:
        #     print "!!Isoform.get_cds Warning: More than 1 cds found. Returning only first occurrence!!"

        #get the first element
        # return obj_intron[0] if obj_intron else None
        #get all possible elements
        return obj_intron

    def in_exon( self, position ):
        """ Function: returns the relative position within the exon ('exonLeft', 'exonRight', 'withinElem') """
        rel_pos = None      #records the position relative in the exon ('exonLeft', 'exonRight', 'withinElem')
        for k,v in self.hashExonList.iteritems():       #k = string (chrom:start-end), v = Exon object
            #if position is found, then return relative position
            rel_pos = v.in_exon( position )
            if rel_pos:
                break

        return rel_pos

    def in_cds( self, position ):
        """ Function: returns the relative position within the exon ('exonLeft', 'exonRight', 'withinElem') """
        rel_pos = None      #records the position relative in the exon ('exonLeft', 'exonRight', 'withinElem')
        for k,v in self.hashExonList.iteritems():       #k = string (chrom:start-end), v = Exon object
            #if position is found, then return relative position
            rel_pos = v.in_cds( position )
            if rel_pos:
                break

        return rel_pos

    def in_intron( self, position ):
        """ Function: similar to def in_exon(), but does this for introns. returns the relative position within the exon ('exonLeft', 'exonRight', 'withinExon') """
        rel_pos = None      #records the position relative in the exon ('exonLeft', 'exonRight', 'withinExon')
        for k,v in self.hashIntronList.iteritems():       #k = string (chrom:start-end), v = Exon object
            #if position is found, then return relative position
            rel_pos = v.in_exon( position )
            if rel_pos:
                break

        return rel_pos


    def find_containing_feature_type( self, position ):
        """
        Function: determines if position is in CDS, exon, or intron
        """
        str_containing_elem = None
        get_elem = self.in_exon( position )
        if get_elem:
            str_containing_elem = "exon: " + get_elem
        else:       #else check intron
            get_elem = self.in_intron( position )
            if get_elem:
                str_containing_elem = "intron: " + get_elem

        return str_containing_elem


    def find_closest_element( self, pos, direction, feat_search = 1 ):
        """
        Args:
            pos = integer that is the genomic position
            direction = integer that is the direction to find the next exon
                - 0 = direction does not matter
                - 1 = find next exon to the right (highest position - higher exon number)
                - -1 = find next exon to the left (lowest position - lower exon number)
            feat_search = integer where:
                -2 = retrieve cds features
                -3 = retrieve intronic features
                -else if any other number (0,1) = retrieve exon features
        Function: finds the closest element (exon or intron) based on the position given, the direction (either towards lower)
        """
        #retrieve the element of interest
        if feat_search == 2:
            list_elems = [x for x in self.hashExonList.values() if x.cdsPos != None]
        elif feat_search == 3:
            list_elems = self.hashIntronList.values()
        else:
            list_elems = self.hashExonList.values()

        #find the nearest element
        if direction == -1:
            list_elems = [x for x in list_elems if x.exonPos.location.end <= pos]
        elif direction == 1:
            list_elems = [x for x in list_elems if x.exonPos.location.start >= pos]

        hash_closest_dist = {}      #k = array index (from loop's enumerate), v = minimum distance
        for i, each_elem in enumerate( list_elems ):
            hash_closest_dist[i] = min( [ abs(pos - each_elem.exonPos.location.start), abs(pos - each_elem.exonPos.location.end) ] )

        #find the isoform with smallest distance
        if hash_closest_dist:
            elem_index = min( hash_closest_dist, key = hash_closest_dist.get )
            return list_elems[elem_index]
        else:       #else if no elements are found within the position of interest, then return the first instance
            try:
                return list_elems[0]
            except IndexError:      #Error meaning no elements (exons, introns, CDS) in list_elems
                # print "Error with Isoform ", chrom, ":", pos_oi, " - no isoform found"
                return None


    """
    Function Start: calculate reading frame
    """
    def calc_reading_frame( self, position, direction ):
        """
        Args:
            position = integer that is the genomic position of interest
            direction = integer that is the direction to find the next exon (NOTE: I feel like this parameter is confusing and redundant)
                - 0 = direction does not matter
                - 1 = find next exon to the right (highest position - higher exon number)
                - -1 = find next exon to the left (lowest position - lower exon number)
        Function: return reading frame for isoform. Returns integer 0, 1, or 2. If position not in any exons, then return None, else if no reading frame for exon then return -1
        """
        #parameters for get_element( self, position, bool_exon = True, direction = 0, bool_strand = True, rel_pos = None ) --> don't consider strand sign & retrieves CDS, not exon
        bool_strand = False
        rel_pos = None
        # bool_exon = False
        # cds = self.get_element( position, bool_exon, direction, bool_strand, rel_pos )    #retrieves Exon object
        feat_type = 2       #means look into CDS in def get_element_v2()
        cds = self.get_element_v2( position, feat_type, direction, bool_strand, rel_pos )    #retrieves Exon object

        if cds:
            return cds.calc_reading_frame( position )
        else:
            return None

    ##DELETE?? - 
    # def calc_reading_frame( self, position, isoform ):
    #     """ Function: return reading frame for isoform. Returns integer 0, 1, or 2 """
    #     #get exons that contain position
    #     cds_exon = self.get_cds( position, isoform )       #k = exon number, v = Exon Model that contains position & is from specific isoform

    #     #if no exons contain position, then return
    #     if not cds_exon:
    #         return None

    #     #if -1, this means there is no reading frame
    #     if cds_exon.read_frame == -1:
    #         return -1

    #     zero_base = position - cds_exon.cdsPos.location.start
    #     if self.strand == -1:      #minus genes, 5' end on higher numerical position of exon (right side of exon)
    #         dist = cds_exon.cdsPos.location.end - position
    #         if zero_base == 0:       #similar to dist = exons[cds_num].exonPos.location.end - ( position + 1 )
    #             dist -= 1       #CORRECTION (ONLY A TEMPORARY FIX): Since UCSC (and as a result cruzdb) is 0-based, the first nucleotide position for each exon is actually the last position nucleotide of the previous intron, so need to add +1 in order to be at the correct position

    #         ##TEST:: print "CRF minus gene --> dist = ", dist, " & position = ", position, " & end_pos = ", exons[cds_num].exonPos.location.end, " & cds_num = ", cds_num
    #     else:       #plus genes: 5' end on lower numerical position of exon (left side of exon)
    #         dist = position - cds_exon.cdsPos.location.start
    #         if zero_base > 0:       #similar to dist = position - ( exons[cds_num].exonPos.location.start + 1 )
    #             dist -= 1       #CORRECTION (ONLY A TEMPORARY FIX): Since UCSC (and as a result cruzdb) is 0-based, the first nucleotide position for each exon is actually the last position nucleotide of the previous intron, so need to add +1 in order to be at the correct position

    #     return (dist + cds_exon.read_frame) % 3

    """
    Function Start: determine if reading frame is preserved
    """
    ##MAY NEED TO DELETE THIS FUNCTION
    # def sj_position_reading_frame( self, position, bool_start = True ):
    #     """
    #     Args:
    #         position = integer that is the genomic position of interest - will be used to find the reading frame at that position
    #         bool_start = is this the starting or ending position of a splice junction 
    #     Helpful when: this function is helping when considering the position of a splice junction
    #     """

    #     if bool_start:      #if this is the starting position
    #         if self.strand > 0:     #if plus strand
    #             start_direction = -1        #retrieve the lowest-positioned cds (left-most cds)
    #         else:                   #if minus strand
    #             start_direction = 1     #retrieve the highest-positioned cds (right-most cds)
    #     else:           #else if this is the ending position
    #         if self.strand > 0:     #if plus strand
    #             end_direction = 1           #retrieve the highest-positioned cds (right-most cds)
    #         else:                   #if minus strand
    #             end_direction = -1      #retrieve the lowest-positioned cds (left-most cds)

    #     return self.calc_reading_frame( position, direction )


    def frame_preserved( self, pos_a, pos_b ):
        """ 
        Args:
            pos_a & pos_b = integer the is the genomic position of interest, where pos_a < pos_b
            isoform_id: string that isoform id (RefSeq ID)
        Function: determines if the reading frame is preserved. Returns True if reading frame is preserved, else returns false
        """
        #need to take strand sign into account (if plus then reading frame goes left to right, else if minus then goes right to left)
        if self.strand > 0:     #if plus strand
            start = pos_a
            end = pos_b
            start_direction = -1        #retrieve the lowest-positioned cds (left-most cds)
            end_direction = 1           #retrieve the highest-positioned cds (right-most cds)
        else:                   #if minus strand
            start = pos_b
            end = pos_a
            start_direction = 1     #retrieve the highest-positioned cds (right-most cds)
            end_direction = -1      #retrieve the lowest-positioned cds (left-most cds)

        #calculate the reading frame for position A & B
        frame_a = self.calc_reading_frame( start, start_direction )
        frame_b = self.calc_reading_frame( end, end_direction )

        #check if frame found for both positions
        if frame_a == None or frame_b == None:
            return { 'rf_start': frame_a, 'rf_end': frame_b, 'rf_preserved': None }

        #calculate next frame
        next_frame = ( frame_a + 1 ) % 3
        if next_frame == frame_b:
            stat = True
        else:
            stat = False

        return { 'rf_start': frame_a, 'rf_end': frame_b, 'rf_preserved': stat }

    """
    Function: check if exons skipped in isoform
    """

    def get_donor_acceptor_sites( self, bool_strand = False ):
        """ 
        Args:
            bool_strand = boolean
                -True = will consider the strand sign when return the donor-acceptor sites (if plus strand, then will return from least to greatest genomic position, else if minus strand then will return from greatest to least genomic position)
                -False = will return donor-acceptor pairs from least to greatest genomic position regardless of gene strand
        Function: returns a tuple of exon ends 
        Notes: for each tuple: for + genes, [0] = 5' splice site (splice donor) & [1] = 3' splice site (splice acceptor). For - gene, [0] = 3' splice site (splice acceptor) & [1] = 5' splice site (splice donor).
        """
        sj_sites = []
        #METHOD 1: use the start & end position of the introns to record splice sites
        introns = sorted( self.hashIntronList.values(), key = lambda k: k.exonPos.location.start, reverse = False )
        for each_ss in introns:     #each_ss = each splice site
            sj_sites.append( (each_ss.exonPos.location.start, each_ss.exonPos.location.end) )
        #METHOD 2: use exons to determine splice site positions
        # exons = sorted( self.hashExonList.values(), key = lambda k: k.exonPos.location.start, reverse = False )
        # for i in range( 0, len( exons ) - 1 ):
        #     sj_sites.append( ( exons[i].exonPos.location.end, exons[i + 1].exonPos.location.start )  )


        ##TEST::
        # print "!!!!!!!!!!In class ISOFORM, see the Introns making up SJs: - ", self.isoform_id
        # for each_i in introns:
        #     print "\tIntron = ", each_i, " & ", each_i.arrAllIsoforms

        #see if strand sign matters
        if bool_strand:
            if self.strand < 0:
                sj_sites = sj_sites[::-1]

        return sj_sites

    def get_donor_sites_only( self, bool_strand = False ):
        """
        Args:
            bool_strand = boolean
                -True = will consider the strand sign when return the donor-acceptor sites (if plus strand, then will return from least to greatest genomic position, else if minus strand then will return from greatest to least genomic position)
                -False = will return donor-acceptor pairs from least to greatest genomic position regardless of gene strand
        Function: only retrieve the 5' splice site (splice donors) and returns an array of splice donor positions based on the strand sign. If 'bool_strand' is True, then will return positions from least to greatest if gene is plus strand, else return False if gene is minus strand.
        """
        da_sites = self.get_donor_acceptor_sites( bool_strand )
        donor_sites = [ x[0] for x in da_sites ] if self.strand > 0 else [ x[1] for x in da_sites ]
        return donor_sites

    def get_acceptor_sites_only( self, bool_strand = False ):
        """
        Args:
            bool_strand = boolean
                -True = will consider the strand sign when return the donor-acceptor sites (if plus strand, then will return from least to greatest genomic position, else if minus strand then will return from greatest to least genomic position)
                -False = will return donor-acceptor pairs from least to greatest genomic position regardless of gene strand
        Function: only retrieve the 3' splice site (splice acceptors) and returns an array of splice acceptor positions based on the strand sign. If 'bool_strand' is True, then will return positions from least to greatest if gene is plus strand, else return False if gene is minus strand.
        """
        da_sites = self.get_donor_acceptor_sites( bool_strand )
        acceptor_sites = [ x[1] for x in da_sites ] if self.strand > 0 else [ x[0] for x in da_sites ]
        return acceptor_sites


    def ensembl_get_exons( self ):
        """
        Function: retrieves all exons associated with Ensembl and returns an array of tuples where each tuple is the start & end position of the exon ([0] = start position of exon, [1] = end position of exon)
        """
        isoform_ensembl = Isoform.closest_isoform_variant_isoform( self.isoform_id, self.db_type, self.hash_pos['chrom'], self.hash_pos['pos_oi'], 2 )

        if not isoform_ensembl:
            return None

        return isoform_ensembl.exons 

    def ensembl_get_donor_acceptor_sites( self ):
        """
        Function: retrieve donor & acceptor sites based on end points of exon
        """
        sj_sites = []       #array of tuples, where for each tuple contains ends of exon
        exons = self.ensembl_get_exons()
        if not exons:
            return None

        for i in range( 0, len( exons ) - 1 ):
            sj_sites.append( ( exons[i][1], exons[i + 1][0] )  )

        return sj_sites


    ##BM1
    def retrieve_all_isoform_exons( self, gene_sym = None, db_type = None ):
        """
        Retrieves all exons for all isoforms based on 
        Args:
            -gene_sym = string that is the gene symbol - NOTE: look at def get_gene_sym_db_all() to see what gene_sym should be based on "db_type"
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        """
        if not gene_sym:
            gene_sym = self.gene_sym
        if not db_type:
            db_type = self.db_type

        all_isoforms = Isoform.get_gene_sym_db_all( gene_sym, db_type )

        if not all_isoforms:
            return {}

        hash_exons = {}     #hash of all exons for each isoform k = isoform ID, v = list of exons
        for each_iso in all_isoforms:
            hash_exons[ each_iso.name ] = each_iso.exons

        return hash_exons

    def annots_get_donor_acceptor_sites( self, gene_sym = None, db_type = None ):
        """
        Function: will retrieve donor & acceptor sites (end positions of introns) based on end points of exon and specific annotations retrieved from database based on "db_type"
        """
        if not gene_sym:
            gene_sym = self.gene_sym
        if not db_type:
            db_type = self.db_type

        hash_exons = self.retrieve_all_isoform_exons( gene_sym, db_type )
        if not hash_exons:
            return None

        hash_introns = {}
        for k,v in hash_exons.iteritems():      #k = isoform ID, v = array of tuples that refer to exons positions, where [0] = is the start position & [1] = is the end position, note that end position > start position regardless of strand sign
            hash_introns[k] = []
            for i in range( 0, len( v ) - 1 ):
                hash_introns[k].append( ( v[i][1], v[i + 1][0] ) )

        return hash_introns

    def is_canon_sj_other_annots( self, start, end, db_type, gene_sym = None ):
        """
        determines if a splice junction is canonical in a specific database (given by 'db_type') with the given positions 'start' and 'end'
        Args:
            -start & end = integers that are start & end position, where start < end regardless of strand sign
            -other_annots = boolean:
                -True = will see if SJ position is canonical in other databases. CAUTION: if going to make true, then know that I will be mixing different database annotations up (may not be a good idea)
                -False = will only use RefSeq to see if SJ is canonical (documented position)
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        """
        sj_pos = (start, end)

        #if gene_sym is not set (i.e. gene_sym = None), then just use self.gene_sym
        if not gene_sym:
            gene_sym = self.gene_sym

        hash_da_annots = self.annots_get_donor_acceptor_sites( gene_sym, db_type )        #hash_da_annots = Hash Donor-Acceptor Sites for annotations of interest

        # #Alternative - return the isoform IDs that contain the splice junction position 'sj_pos'
        # list_isoforms = []
        # for k,v in hash_da_annots.iteritems():
        #     if sj_pos in v:
        #         list_isoforms.append( k )
        # return list_isoforms

        for k,v in hash_da_annots.iteritems():
            if sj_pos in v:
                return True

        return False


    def sj_get_ligated_exons( self, start, end, rel_pos_start = None, rel_pos_end = None ):
        """
        returns a string of the gene features (exon or intron) ligated by a splicing event

        Args:
            -start & end: genomic position of splicing event, where start < end. This assumes the splicing event is on the same chromosome as the Isoform.
            -rel_pos_start & rel_pos_end = string that is relative position of the element. Can be the following values: (look at Exon.in_exon() for appropriate labels)
                -None = position can be anywhere within exon
                -"withinElem" = found within the exon
                -"exonLeft" = left side of exon 
                -"exonRight" = right side of exon
        Outputs:
            a string of the ligated exons for the splicing event (format exonA:exonB)
        """
        bool_strand = False             #as of now, I don't consider the strand sign

        direction_minus = -1
        start_num = self.get_element_v2( start, 1, direction_minus, bool_strand, rel_pos_start )    #retrieve the left-most exon (exon with lowest numerical position)
        direction_plus = 1
        end_num = self.get_element_v2( end, 1, direction_plus, bool_strand, rel_pos_end )     #retrieve the right-most exon (exon with highest numerical position)

        str_connected_exons = start_num.exonPos.type + '_' + str( start_num.exonNum ) + ":" + end_num.exonPos.type + '_' + str( end_num.exonNum )
        return str_connected_exons


    def is_exon_skip( self, start, end ):
        """
        Args:
            start & end = genomic positions of start & end position of SJ, where start < end
        Function: returns an array of integers, where each integer represents the exons skipped. Blank array means no exons were skipped 
        NOTE: the exon numbers are zero-based, so exon 0 = 1st exon, exon 1 = 2nd exon
        """
        bool_strand = False
        rel_pos = None

        direction_minus = -1
        start_num = self.get_element_v2( start, 1, direction_minus, bool_strand, rel_pos )    #retrieve the left-most exon (exon with lowest numerical position)
        direction_plus = 1
        end_num = self.get_element_v2( end, 1, direction_plus, bool_strand, rel_pos )     #retrieve the right-most exon (exon with highest numerical position)

        #MAY DELETE - I am using get_element_v2() now
        # #retrieve Exon objects that contain position
        # start_num = self.get_element( start, True, -1, False, None )    #retrieve the left-most exon (exon with lowest numerical position)
        # end_num = self.get_element( end, True, 1, False, None )     #retrieve the right-most exon (exon with highest numerical position)

        if start_num == None or end_num == None:
            return None

        #determine the lower & higher exon number
        skip_range = range( start_num.exonNum + 1, end_num.exonNum ) if start_num.exonNum < end_num.exonNum else range( end_num.exonNum + 1, start_num.exonNum )

        return [i for i in skip_range]


    """
    Function: determine the effect of a splice junction in the isoform
    """
    # def is_canon_sj( self, start, end ):
    #     """ Function: checks if positions of splice sites land in canonical donor-acceptor sites - returns True if position exists array of donor-acceptor sites, else returns False """
    #     #check start & end position of splice junction to see if it is canonical
    #     sj_pos = (start, end)
    #     if sj_pos in self.donor_acceptor_sites:
    #         return True
    #     else:
    #         return False

    def is_canon_sj( self, start, end, other_annots = False, db_type = 4 ):
        """ 
        Args:
            -start & end = integers that are start & end position, where start < end regardless of strand sign
            -other_annots = boolean:
                -True = will see if SJ position is canonical in other databases. CAUTION: if going to make true, then know that I will be mixing different database annotations up (may not be a good idea)
                -False = will only use RefSeq to see if SJ is canonical (documented position)
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....) - I use this as the default as of now - CONJ: the reason is because I THINK it is the most comprehensive genomic annotated database as of now
        Function: checks if positions of splice sites land in canonical donor-acceptor sites - returns True if position exists array of donor-acceptor sites, else returns False """
        #check start & end position of splice junction to see if it is canonical
        sj_pos = (start, end)
        #get all donor-acceptor sites from both annotations
        da_annot = self.donor_acceptor_sites       #da_annot = Donor-Acceptor Sites for RefSeq

        ##TEST:: print "\t\tClass Isoform - ", self.isoform_id, " & sj_pos = ", sj_pos, " & da_sites = ", da_annot, " & sj_pos in da_sites? - ", sj_pos in da_annot

        if sj_pos in da_annot:
            sj_canon =  True
        else:
            sj_canon = False

        #if want to see if SJ is canonical in other databases
        if other_annots and not sj_canon:
            sj_canon = self.is_canon_sj_other_annots( start, end, db_type )

        return sj_canon

    #I THINK I CAN DELETE THIS...
    # def is_canon_sj_other_annots( self, start, end ):
    #     """
    #     Function: use this function to see if SJ is canonical in other annotation databases
    #     """
    #     da_ensembl = self.ensembl_get_donor_acceptor_sites()        #da_ensembl = Donor-Acceptor Sites for Ensembl
    #     if da_ensembl and sj_pos in da_ensembl:
    #         return True
    #     else:
    #         return False



    def is_aberrant_sj( self, start, end, intronic = False ):
        """
        Args:
            start & end = integers that are the start & end genomic position 
            intronic = boolean
                -True = return an isoform ID even if it does land in the intronic space
                -False = do not return isoform ID if only lands in intronic space
        Function: retrieve information about start & end splicing to see if how the splice junction is aberrant, including whether splicing is in a exon (donor site, acceptor site, in exon) or intronic, any exons skipped, and any frameshifts """
        #see where position of start & end positions is located
        exon_start = self.in_exon( start )
        exon_start = 'exon_canon' if exon_start == 'exonRight' else exon_start      #'exon_canon' means splicing occurring at canonical position
        exon_end = self.in_exon( end )
        exon_end = 'exon_canon' if exon_end == 'exonLeft' else exon_end         #'exon_canon' means splicing occurring at canonical position

        #if intronic == False, then check if SJ is intronic - if so then discard SJ
        if not intronic:
            #if no exon (or intron) is returned, return None
            if exon_start == None or exon_end == None:
                return None
        
        exon_skips = self.is_exon_skip( start, end )        #check for exons skipped
        check_frame = self.frame_preserved( start, end )    #check to see if reading frame is preserved

        #exon_skips = array of integers referring to exons skipped
        #frame_preserved = True (SJ preserves frame), False (SJ does not preserve frame), or None (no frame, probably b/c SJ is intronic)
        #intronic (when I do implement it) = True (means SJ is intronic), False = (SJ is not intronic)
        calc_rf = ( check_frame['rf_start'], check_frame['rf_end'] )
        aberrations = { 'exon_start': exon_start, 'exon_end': exon_end, 'exon_skip': exon_skips, 'frame_preserved': check_frame['rf_preserved'], 'calc_rf': calc_rf }
        aberrations.update( self.score_aberrations( aberrations ) )

        return aberrations

    def score_aberrations( self, aberrations ):
        """ 
        Args:
            aberrations = hash of all aberrations, keys are 'exon_start', 'exon_end', 'exon_skips', 'frame_preserved'
        Function: scores the aberrations associated with the abnormalities associated with splice junction
        Note on scoring: 0 = no aberrations, 1 = exon skip, 100 = exonic (not on exon ends), 200 = frameshift, 500 = intronic
        """
        score = 0
        effects = []        #record the effects
        #quantify number of exon skips - make sure not aberrations['exon_skip'] is not 'NoneType'
        score += len( aberrations['exon_skip'] ) if aberrations['exon_skip'] is list else 0
        if score > 0:
            effects.append( 'exon_skip: ' + ','.join( aberrations['exon_skip'] ) )

        #check if splice junction is intronic (one of the positions lands in the intron) or exonic (within exon body, not at ends of exon)
        if aberrations['exon_start'] == None or aberrations['exon_end'] == None:        #splicing in intron
            score += 500
            effects.append( 'intronic' )
        elif aberrations['exon_start'] != 'exon_canon' or aberrations['exon_end'] != 'exon_canon':      #splicing in exon (not at donor-acceptor sites)
            score += 100
            effects.append( 'exonic' )
        #check if reading frame is preserved
        if not aberrations['frame_preserved']:
            score += 200
            effects.append( 'frameshift' )

        return { 'score': score, 'effects': effects }
    
    """
    Functions: Retrieve information for gene & isoform, including start & end position 
    """

    ##DELETE? - MAY DELETE THIS LATER AS I THINK I CAN RETRIEVE THE ISOFORM JUST BY USING CRUZDB
    @classmethod
    def get_isoform( cls_obj, isoform_id, db_type ):
        """ retrieves all isoforms that have the same gene id 'isoform_id' """
        #retrieve information from cruzdb
        #NOTE: I am using "refGene.filter_by" because for some reason more isoform IDs retrieved by it than by .bin_query
        # cruzdbInfo = ExonList.objCruzDB.bin_query( "refGene", chrom, pos_start, pos_end ).all()       ##DELETE THIS LINE - can't restrict myself to just "RefSeq"
        # info = Isoform.obj_cruzdb.refGene.filter_by( name = isoform_id ).all()        ##DELETE THIS LINE - can't restrict myself to just "RefSeq"
        info = Isoform.get_isoforms_db_all( isoform_id, db_type )

        #save all isoforms with geneSym into array
        isoform_info = None
        for x in info:
            #Method 1: If I only want to consider mRNA protein-coding isoforms
            # if x.name2 == isoform_id and "NM_" in x.name:        #if want to consider only mRNA isoforms
            #     isoform_info = x
            
            #Method 2: if I want all isoforms
            """
            IMPORTANT: the reason I decided this is because there could a splice junction could correctly splice a non-mRNA transcript, therefore reducing erroneously stating an SJ is aberrant when indeed it maps to a non-mRNA transcript. So need to consider reads that map to non-mRNA transcripts as well.
            """
            if x.name == isoform_id:      #if want to consider all isoforms (mRNA, ncRNA, etc.)
                isoform_info = x
                break

        return isoform_info

    @classmethod
    def set_cruzdb( cls_obj, obj_cruzdb ):
        """ Function: assigns cruzdb object to class variable """
        cls_obj.obj_cruzdb = obj_cruzdb

    def get_feat_info( self, pos_start, pos_end, exon_frame, exon_num, canonical, feat_type ):
        """
        Args:
            pos_start & pos_end = integer that is the lower & higher nucleotide position for the exon, respectively
            exon_frame = integer that is 0, 1, or 2. Refers to the reading frame of the 5' end of exon (lower position for plus genes & higher position for minus genes). It is -1 if the exon is non-coding & None if reading frame is unknown
            exon_num = exon number in gene. For plus genes, exon number increases from lower to higher position. For minus genes, exon number decreases from lower to higher position. If None, means unknown what the exon number is (could be non-canonical exon number)
            canonical = boolean that, if True, means this exons is canonical (known to exist) in the gene, else if False, then is not a canonical exon in the gene
            feat_type = string that conveys the type of feature the element is, either 'exon' or 'intron'
        Function: creates & returns hash to create Exon
        """
        #prepare hash that will record elements into exon list
        exon_info = {
        "pos_start": int( pos_start ),
        "pos_end": int( pos_end ),
        "strand_sign": self.strand,
        "chrom": self.chrom,
        "isoform_id": self.isoform_id,
        "exon_frame": exon_frame,
        "exon_num": exon_num,
        "canonical": canonical,
        "feat_type": feat_type }

        return exon_info


    """
    Calculating exon & gene expression
    """
    def calc_isoform_rpkm_exons_only( self, bam_reader, library_size, unique_reads = True ):
        """
        Args:
            bam_reader = HTSeq.BAM_Reader instance
            library_size = integer that is the total number of mapped exons in the library
            unique_reads = boolean
                -True = only consider uniquely mapped reads that map to genomic_range (i.e. a.optional_field( "NH" ) == 1)
                -False = consider all reads that map to genomic_range (unique + multimapped)
        Function: calculates the RPKM of an isoform, but only considers the total length of the exons when calculate the RPKM (instead of the entire length of the gene)
        """
        #sum all the read counts for all exons in an isoform -> sum total length of all exons for a gene -> calculate RPKM by 

        #calculate number of reads per exon
        sum_count = 0
        exon_len = []
        for k,v in self.hashExonList.iteritems():      #k = string that is the genomic position, v = Exon instance
            genomic_range = k
            sum_count += Isoform.quant_genes_rpkm( bam_reader, genomic_range, unique_reads )
            exon_len.append( v.exonPos.location.end - v.exonPos.location.start )

        #sum length of all exons
        exon_len_sum = sum( exon_len )

        return ( 10**9 * float( sum_count ) ) / ( library_size * exon_len_sum )

    def calc_isoform_rpkm_exons_only_pysam( self, pysam_file, library_size, unique_reads = True ):
        """
        Args:
            pysam_file = pysam.AlignmentFile that opens up the mapped reads bam file (e.g. accepted_hits.bam)
            library_size = integer that is the total number of mapped exons in the library
            unique_reads = boolean
                -True = only consider uniquely mapped reads that map to genomic_range (i.e. a.optional_field( "NH" ) == 1)
                -False = consider all reads that map to genomic_range (unique + multimapped)
        Function: calculates the RPKM of an isoform, but only considers the total length of the exons when calculate the RPKM (instead of the entire length of the gene)
        """
        #sum all the read counts for all exons in an isoform -> sum total length of all exons for a gene -> calculate RPKM by 

        #calculate number of reads per exon
        sum_count = 0
        exon_len = []
        for k,v in self.hashExonList.iteritems():      #k = string that is the genomic position, v = Exon instance
            genomic_range = k
            sum_count += Isoform.quant_genes_rpkm( pysam_file, genomic_range, unique_reads )
            exon_len.append( v.exonPos.location.end - v.exonPos.location.start )

        #sum length of all exons
        exon_len_sum = sum( exon_len )

        return ( 10**9 * float( sum_count ) ) / ( library_size * exon_len_sum )

    @staticmethod
    def total_mapped_reads( path_bam ):
        """
        Args:
            path_bam = string that is the path to the .bam file
        Function: calculates library size - quantifies total number of mapped reads in library & returns array of number of reads that map to 
        NOTE: I think this is considering total reads, not total fragments, so perhaps this is calculating RPKM instead of FPKM. How do I determine the number of fragments
        """
        # #Method 1: the original method - but for some reason this is not working?!?!
        # mapped_reads = [ int( x.split('\t')[2] ) for x in pysam.idxstats( path_bam ) if len( x.split('\t') ) >= 2 and 'chr' in x.split('\t')[0] ]

        #Method 2: I want to get read of this method though...
        mapped_reads = []
        for i, x in enumerate( pysam.idxstats( path_bam ).split('\n') ):
            if len( x.split('\t') ) >= 2 and 'chr' in x.split('\t')[0]:
                mapped_reads.append( int( x.split('\t')[2] ) )

        return sum( mapped_reads )

    @staticmethod
    def quant_genes_rpkm( bam_reader, genomic_range, unique_reads = True ):
        """
        Args:
            bam_reader = HTSeq.BAM_Reader instance
            genomic_range = string with the format "chr_num:start-end"
            unique_reads = boolean
                -True = only consider uniquely mapped reads that map to genomic_range (i.e. a.optional_field( "NH" ) == 1)
                -False = consider all reads that map to genomic_range (unique + multimapped)
        Function: quantify the number of reads that mapes to genomic range 'genomic_range'
        """
        #determine if only unique reads should be considered (unique_reads = True) or all reads (unique_reads = False)
        if unique_reads:
            unique_count = 0
            for i, a in enumerate( bam_reader.fetch( region = genomic_range ) ):
                if a.optional_field( "NH" ) == 1:
                    unique_count += 1
            
            return unique_count
        else:
            return sum( 1 for x in bam_reader.fetch( region = genomic_range ) )


    @staticmethod
    def quant_genes_rpkm_pysam( pysam_file, genomic_range, unique_reads = True ):
        """
        Args:
            pysam_file = pysam.AlignmentFile that opens up the mapped reads bam file (e.g. accepted_hits.bam)
            genomic_range = string with the format "chr_num:start-end"
            unique_reads = boolean
                -True = only consider uniquely mapped reads that map to genomic_range (i.e. a.optional_field( "NH" ) == 1)
                -False = consider all reads that map to genomic_range (unique + multimapped)
        Function: quantify the number of reads that mapes to genomic range 'genomic_range'
        """
        #determine if only unique reads should be considered (unique_reads = True) or all reads (unique_reads = False)
        if unique_reads:
            # unique_count = 0
            # for i, a in enumerate( pysam_file.fetch( region = genomic_range ) ):
            #     if a.mapq == 50:
            #         unique_count += 1
            
            # return unique_count

            return sum( 1 for x in pysam_file.fetch( region = genomic_range ) if x.mapq == 50 )
        else:
            return sum( 1 for x in pysam_file.fetch( region = genomic_range ) )

    @staticmethod
    def quant_genes_rpkm_v2( bam_reader, genomic_range, count_NH = 1 ):
        """
        Args:
            bam_reader = HTSeq.BAM_Reader instance
            genomic_range = string with the format "chrom:start-end"
            unique_reads = boolean
                -True = only consider uniquely mapped reads that map to genomic_range (i.e. a.optional_field( "NH" ) == 1)
                -False = consider all reads that map to genomic_range (unique + multimapped)
            count_NH = integer that is the number of positions where a read maps to (less than or equal to)
                -1 = Uniquely-mapped read. This read maps to 1 position.
                -2 = Read maps to 2 or less different positions
                -3 = Read maps to 3 or less positions, etc.
        Function: quantify the number of reads that mapes to genomic range 'genomic_range'
        """
        #determine if only unique reads should be considered (unique_reads = True) or all reads (unique_reads = False)
        read_count = 0
        for i, a in enumerate( bam_reader.fetch( region = genomic_range ) ):
            if a.optional_field( "NH" ) == count_NH:
                read_count += 1
        
        return read_count


    @staticmethod
    def quant_genes_fpkm( bam_reader, genomic_range, unique_reads = True ):
        """
        Args:
            bam_reader = HTSeq.BAM_Reader instance
            genomic_range = string with the format "chr_num:start-end"
            unique_reads = boolean
                -True = only consider uniquely mapped reads that map to genomic_range (i.e. a.optional_field( "NH" ) == 1)
                -False = consider all reads that map to genomic_range (unique + multimapped)
        Function: same as def quant_genes_rpkm(), but uses pair-end reads (FPKM) instead of single-end reads (RPKM)
        """
        #retrieve the number of reads that mape to genomic region
        if unique_reads:
            unique_count = 0
            pair_reader = HTSeq.pair_SAM_alignments( bam_reader.fetch( region = genomic_range ) )
            for i, (a, b) in enumerate( pair_reader ): 
                if a and b:     #make sure a & b are not NONE
                    if a.optional_field( "NH" ) == 1 and b.optional_field( "NH" ) == 1:
                        unique_count += 1
            
            return unique_count
        else:
            pair_reader = HTSeq.pair_SAM_alignments( bam_reader.fetch( region = genomic_range ) )
            return sum( 1 for a,b in pair_reader if a and b )       #make sure both reads in pair are not NONE

    @staticmethod
    def calc_read_density( count, genomic_range, library_size ):
        """
        Args:
            count = number of reads mapped to region, calculated by quant_exons()
            genomic_range = string in the format "chrom:start-end"
            library_size = integer that is the total number of reads in library
        Function: calculates RPKM or FPKM, depending on how parameter 'count' was generated (def quant_genes_rpkm() or def quant_genes_fpkm())
        """
        #calculate RPKM = 10^9 * num_of_reads / (region_length * library_size )
        #split genomic range 
        num_pos = genomic_range.split(':')[1].split('-')

        exon_len = int( num_pos[1] ) - int( num_pos[0] )
        return ( 10**9 * float( count ) ) / ( library_size * exon_len )


    """
    Retrieve information from Ensembl database
    """
    @staticmethod
    def get_ucsc_gene_name( isoform_id ):
        """
        Args:
            -isoform_id = string that is the isoform ID, CONJ: I think this should be the RefSeq form of the isoform ID  
            -gene_sym = string that is the gene symbol (e.g. BRAF, RAF1, MUC5B) - this is property "name2" in the refGene table [cruzdb.Genome(db = "hg19").refGene.filter_by(name2 = gene_sym)]
        Function:
            this will retrieve the UCSC gene name. This can be used to:
            -retrieve the mRNA sequence
            -retrieve the associate Ensembl ID
        """
        # #get the gene symbol - will need this to retrieve the UCSC gene info
        # gene_refseq = Isoform.obj_cruzdb.refGene.filter_by( name = isoform_id ).first()
        # # gene_sym = get_gene.name2

        # # gene1 = Isoform.obj_cruzdb.refGene.filter_by(name2 = geneName)     #if did not use ".first()", will need to do gene1.all()[0] to retrieve the properties

        # gene_ucsc = Isoform.obj_cruzdb.knownToRefSeq.filter_by( value = gene_refseq.name ).first()     #this will give me another ID for the gene that I can use to retrieve the mRNA sequence

        gene_ucsc = Isoform.obj_cruzdb.knownToRefSeq.filter_by( value = isoform_id ).first()     #this will give me another ID for the gene that I can use to retrieve the mRNA sequence

        return None if gene_ucsc is None else gene_ucsc.name

    @staticmethod
    def get_ensembl_gene_name( isoform_id ):
        """
        Args:
            -isoform_id = string that is the isoform ID, CONJ: I think this should be the RefSeq form of the isoform ID  
            -Isoform.obj_cruzdb = cruzDB object (e.g. Isoform.obj_cruzdb = cruzdb.Genome(db = "hg19"))
            -gene_sym = string that is the gene symbol (e.g. BRAF, RAF1, MUC5B) - this is property "name2" in the refGene table [cruzdb.Genome(db = "hg19").refGene.filter_by(name2 = gene_sym)]
        Function:
            this function will return the Ensembl gene ID based on the gene symbol (e.g. BRAF, MUC5B)
        """
        #STEP: get the UCSC gene name & then retrieve the gene mRNA sequence
        ucsc_gene_name = Isoform.get_ucsc_gene_name( isoform_id )
        if not ucsc_gene_name:
            return None

        gene_ensembl = Isoform.obj_cruzdb.knownToEnsembl.filter_by( name = str(ucsc_gene_name) ).first()       #this will return 2 properties, "name" & "seq"

        #return the Ensembl Transcript ID (e.g. ENST, such as ENST00000458694)
        return None if gene_ensembl is None else gene_ensembl.value

 
    @staticmethod
    def get_ensembl_isoforms( isoform_id ):
        """
        Args:
            -isoform_id = string that is the isoform ID, CONJ: I think this should be the RefSeq form of the isoform ID           
            -gene_sym = string that is the gene symbol (e.g. BRAF, RAF1, MUC5B) - this is property "name2" in the refGene table [cruzdb.Genome(db="hg19").refGene.filter_by(name2=gene_sym)]
        Function:
            this will retrieve the multiple Ensembl transcript IDs associated with 
        """
        #retrieve the ensembl transcript ID, find the ensembl gene ID, and then find all transcript IDs associated with the ensembl gene ID
        ensembl_transcript_id = Isoform.get_ensembl_gene_name( isoform_id )     #in the format "ENST" (e.g. ENST00000458694)
        if not ensembl_transcript_id:
            return None

        ensembl_gene_id = Isoform.obj_cruzdb.ensGene.filter_by( name = ensembl_transcript_id ).first()      #retrieve ensembl gene ID, in the format ENSG (e.g. ENSG00000124593)
        ensembl_isoforms = Isoform.obj_cruzdb.ensGene.filter_by( name2 = ensembl_gene_id.name2 ).all()

        return ensembl_isoforms

    """
    Other functions: static
    """
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

