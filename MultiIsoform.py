#/usr/bin/python
from Isoform import Isoform

class MultiIsoform( Isoform ):
    """ record all isoforms associated with position """
    def __init__( self, db_type, chrom, start, end, gene_sym = None, isoform_id = None ):
        """
        Args:
            -chrom = string that is in the format "chrNum" (e.g. chr2, chr14)
            -start, end = integers that are the start & end position. Also, the position "start" will be used to determine with variant of the isoform will be used (this is because cruzdb will retrieve an isoform that will have multiple positions recorded for 1 isoform)
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
            -gene_sym = string, if defined, will only record isoforms that are assigned to any of the gene symbols recorded in this array
            -isoform_id = optional parameter, if this is provided then will only record one isoform. If gene_sym is also defined, then isoform_id needs to be an isoform of gene_sym else it will not be created
        """
        if gene_sym:
            all_isoform_info = Isoform.get_gene_sym_db_all( gene_sym, db_type )
        elif isoform_id:
            all_isoform_info = Isoform.get_isoforms_db_all( isoform_id, db_type )
        else:
            all_isoform_info = MultiIsoform.get_isoform_gene( chrom, start, end, db_type )

        if all_isoform_info:
            all_isoforms = [x.name for x in all_isoform_info]
            all_gene_sym = [x.name2 for x in all_isoform_info] if db_type != 3 else [x.name for x in all_isoform_info]
            gene_sym = max( set(all_gene_sym), key = all_gene_sym.count )     #retrieve the most common gene symbol
        else:
            all_isoforms = []
            gene_sym = None

        if isoform_id in all_isoforms:
            all_isoforms = [isoform_id]
        # all_isoforms_filter = [x for x in all_isoforms if x == isoform_id]
        # if all_isoforms_filter:
        #     all_isoforms = all_isoforms_filter
        
        hash_pos = {'chrom': chrom, 'pos_oi': start}        #this hash will be used with the Isoform class to find isoform variant closest to this position

        self.db_type = db_type
        self.gene_sym_oi = gene_sym
        self.hash_isoforms = {}     #key = isoform id, value = Isoform class
        for each_isoform in all_isoforms:
            #if gene_sym is defined, then only record isoforms that are assigned to this gene symbol
            if self.gene_sym_oi:
                gene_info = Isoform.get_isoforms_db_first( each_isoform, self.db_type )
                # gene_info = Isoform.obj_cruzdb.refGene.filter_by( name = each_isoform ).first()

                ##TEST:: print "\tMULTIISOFORM.__init__() - gene_info.name2 = ", gene_info.name2, " & ", self.gene_sym_oi, " & each_isoform = ", each_isoform


                if str( gene_info.name2 ) == self.gene_sym_oi:
                    # self.hash_isoforms[each_isoform] = Isoform( self.db_type, each_isoform, hash_pos )
                    self.hash_isoforms[each_isoform] = Isoform( self.db_type, each_isoform, None )

                    ##TEST:: print "\tMULTIISOFORM.__init__ with GENE_SYM: each_isoform = ", each_isoform, " & ", self.hash_isoforms[each_isoform], " & gene_sym = ", self.gene_sym_oi

            #else, just record all isoforms
            else:
                # self.hash_isoforms[each_isoform] = Isoform( self.db_type, each_isoform, hash_pos )
                self.hash_isoforms[each_isoform] = Isoform( self.db_type, each_isoform, None )

                ##TEST:: print "\tMULTIISOFORM.__init__: isoform = ", each_isoform, " & ", self.hash_isoforms[each_isoform]

        #for each isoform, determine if each element (exon & intron) are constitutive elements or not
        #record all constitutive exons
        c_exons = self.find_constitutive_element( self.gene_sym_oi, False )     #hash where k = string element pos (chrom:start-end), v = Exon objects
        #record all constitutive introns (also means constitutive splice junctions)
        c_introns = self.find_constitutive_element( self.gene_sym_oi, True )     #hash where k = string element pos (chrom:start-end), v = Exon object  
        for each_isoform in self.hash_isoforms:
            for k,v in self.hash_isoforms[each_isoform].hashExonList.iteritems():     #key = string that is exon range (chrom:start-end), value = Exon object
                stat_constitutive = True if k in c_exons.keys() else False
                v.set_constitutive( stat_constitutive )

            for k,v in self.hash_isoforms[each_isoform].hashIntronList.iteritems():     #key = string that is intron range (chrom:start-end), value = Exon object
                stat_constitutive = True if k in c_introns.keys() else False
                v.set_constitutive( stat_constitutive )

        self.mi_strand = self.get_gene_strand()     #mi_strand = MultiIsoform Strand, this is the predominant strand sign for the gene 


        ##TEST::
        # for k,v in self.hash_isoforms.iteritems():
        #     print "\tMI - show isoform ", k, " & v = ", v


        #get the lowest & highest boundary for all genes
        # self.boundary_all = self.get_gene_range_isoforms()

    
    def get_gene_strand( self ):
        """
        Function: retrieves the predominant strand sign
        """
        all_ss = [v.strand for k,v in self.hash_isoforms.iteritems()]
        all_ss = [x for x in all_ss if x == -1 or x == 1]
        return 5 if not all_ss else max( set(all_ss), key = all_ss.count )      #I return "5" because this is stating there is no gene strand associated with this gene, don't want to assign "None" because it won't pass numerical comparison operators (e.g. self.mi_strand == 1)


    @staticmethod
    def get_isoform_gene( chrom, start, end, db_type ):
        """
        Function: retrieve all isoforms associated with position chrom:start-end
        Args:
            -chrom = string in the format 'chr#' (chr9, chr12)
            -start & end = integer that is position to where other positions will be looked for
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        """
        # all_isoforms = Isoform.get_isoforms_by_pos_db_all( chrom, start, end, db_type )
        # # all_isoforms = Isoform.obj_cruzdb.bin_query( 'refGene', chrom, start, end ).all()
        # #Method 1: If I only want to consider mRNA protein-coding isoforms
        # # return [x.name for x in all_isoforms if "NM_" in x.name]

        # #Method 2: if I want all isoforms
        # """
        # IMPORTANT: the reason I decided to retrieve all isoforms & not just protein-coding isoforms is because there could a splice junction that could correctly splice a non-mRNA transcript, therefore reducing erroneously stating an SJ is aberrant when indeed it maps to a non-mRNA transcript. So need to consider reads that map to non-mRNA transcripts as well.
        # """
        # return {x.name.upper():x.name2.upper() for x in all_isoforms} if db_type != 3 else {x.name.upper():x.name.upper() for x in all_isoforms}     #return hash where key = isoform ID & value = gene symbol

        all_isoforms = Isoform.get_isoforms_by_pos_db_all( chrom, start, end, db_type )
        return all_isoforms

    @staticmethod
    def get_isoforms( chrom, start, end, db_type ):
        """
        Function: retrieve all isoforms associated with position chrom:start-end
        Args:
            -chrom = string in the format 'chr#' (chr9, chr12)
            -start & end = integer that is position to where other positions will be looked for
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        """
        
        all_isoforms = Isoform.get_isoforms_by_pos_db_all( chrom, start, end, db_type )
        # all_isoforms = Isoform.obj_cruzdb.bin_query( 'refGene', chrom, start, end ).all()
        #Method 1: If I only want to consider mRNA protein-coding isoforms
        # return [x.name for x in all_isoforms if "NM_" in x.name]

        #Method 2: if I want all isoforms
        """
        IMPORTANT: the reason I decided to retrieve all isoforms & not just protein-coding isoforms is because there could a splice junction that could correctly splice a non-mRNA transcript, therefore reducing erroneously stating an SJ is aberrant when indeed it maps to a non-mRNA transcript. So need to consider reads that map to non-mRNA transcripts as well.
        """
        return [str(x.name) for x in all_isoforms]

    @classmethod
    def get_gene_syms( cls_obj, chrom, start, end ):
        """ Function: retrieve all isoforms associated with position chrom:start-end """
        
        all_isoforms = Isoform.get_isoforms_by_pos_db_all( chrom, start, end, db_type )
        # all_isoforms = Isoform.obj_cruzdb.bin_query( 'refGene', chrom, start, end ).all()
        #Method 1: If I only want to consider mRNA protein-coding isoforms
        # return [x.name for x in all_isoforms if "NM_" in x.name]

        #Method 2: if I want all isoforms
        """
        IMPORTANT: the reason I decided to retrieve all isoforms & not just protein-coding isoforms is because there could a splice junction that could correctly splice a non-mRNA transcript, therefore reducing erroneously stating an SJ is aberrant when indeed it maps to a non-mRNA transcript. So need to consider reads that map to non-mRNA transcripts as well.
        """
        gene_syms = [x.name2 for x in all_isoforms]

        return list( set(gene_syms) )

    ##MAY DELETE & REPLACE WITH get_gene_range_v2()
    @staticmethod
    def get_gene_range( gene_sym, db_type ):
        """
        Function: retrieves the min & max range of a gene in 'gene_sym'
        Args:
            -gene_sym = string that is the gene symbol (e.g. BRAF, RAF1)
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        """
        all_isoforms = Isoform.get_gene_sym_db_all( gene_sym, db_type )
        # all_isoforms = Isoform.obj_cruzdb.refGene.filter_by( name2 = gene_sym ).all()
        if db_type == 3:
            start_min = min( [x.txStart for x in all_isoforms if x.name == gene_sym] )
            end_max = max( [x.txEnd for x in all_isoforms if x.name == gene_sym] )
        else:
            start_min = min( [x.txStart for x in all_isoforms if x.name2 == gene_sym] )
            end_max = max( [x.txEnd for x in all_isoforms if x.name2 == gene_sym] )

        return ( start_min, end_max )

    @staticmethod
    def get_gene_range_v2( gene_sym, db_type ):
        """
        Function: same as get_gene_range(), but retrieves the chrom position as well as the min & max range of a gene in 'gene_sym'
        Args:
            -gene_sym = string that is the gene symbol (e.g. BRAF, RAF1)
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
                -4 = uses GENCODE (Isoform.obj_cruzdb.wgEncodeGencodeBasicV19 --> CONJ: this is only for hg19 build I think....)
        """
        all_isoforms = Isoform.get_gene_sym_db_all( gene_sym, db_type )
        #if no isoform information found, then return nothing
        if not all_isoforms:
            return {}

        #find chromosome
        chrom = all_isoforms[0].chrom
        for x in all_isoforms:
            if x.name == gene_sym:
                chrom = x.chrom
                break
        chrom = chrom if not '_' in chrom else chrom.split('_')[0]      #sometimes chromosome shows up with strange name (e.g. chr9_asfwgr or something, so I only want the chromosome position before the "_")

        # all_isoforms = Isoform.obj_cruzdb.refGene.filter_by( name2 = gene_sym ).all()
        if db_type == 3:
            start_min = min( [x.txStart for x in all_isoforms if x.name == gene_sym] )
            end_max = max( [x.txEnd for x in all_isoforms if x.name == gene_sym] )
        else:
            start_min = min( [x.txStart for x in all_isoforms if x.name2 == gene_sym] )
            end_max = max( [x.txEnd for x in all_isoforms if x.name2 == gene_sym] )

        return {'chrom': chrom, 'start_min': start_min, 'end_max': end_max}

    def get_gene_range_isoforms( self ):
        """
        Function: retrieves the genomic range for all isoforms, finding the lowest & highest positions
        """
        list_start = [v.boundary[0] for k,v in self.hash_isoforms.iteritems()]
        list_end = [v.boundary[1] for k,v in self.hash_isoforms.iteritems()]
        #if no position recorded for either list, then return an empty tuple
        if not list_start or not list_end:
            return ()

        boundary_start = min( list_start )
        boundary_end = max( list_end )

        return ( boundary_start, boundary_end )

    """
    Function: find the prevalence of each exon for a gene (i.e. number of isoforms that contain an isoform)
    """
    def all_exon_prevalence( self, gene_sym, bool_intron = True ):
        """
        Args:
            gene_sym = string that is the gene symbol, only record elements that belong to the same gene
            bool_intron = boolean where:
                -True = will look for constitutive introns in the gene (i.e. splice junctions that are present across all isoforms)
                -False = will look for constitutive exons in the gene (exons present across all isoforms of the gene)
        Function: determines the prevalence of each elem (exon or intron), returns 2 hashes - both hashes keys are string of their position (format: chrom:start-end), elem_prevalence: value is the count of the number of isoforms that contains the elemnt
        """
        #calculate the total number of isoforms
        # total_isoform_count = len( [k for k,v in self.hash_isoforms.iteritems() if v.gene_sym == gene_sym] )

        elem_prevalence = {}      #k = genomic position (chrom:start-end), v = number of isoforms that have contain the isoform
        elem_obj = {}   #k = genomic position (chrom:start-end), v = Exon object associated with genomic position
        for k,v in self.hash_isoforms.iteritems():      #k = string that is isoform_id, v = Isoform instance
            hash_elem = v.hashIntronList if bool_intron else v.hashExonList

            #if not the same gene symbol, then skip
            if v.gene_sym != gene_sym:
                continue

            for k2,v2 in hash_elem.iteritems():    #k2 = string element pos (chrom:start-end), v2 = Exon object
                if not k2 in elem_prevalence:
                    elem_prevalence[k2] = 1
                    elem_obj[k2] = v2
                else:
                    elem_prevalence[k2]+= 1

        return [elem_prevalence, elem_obj]


    """
    Function: find position within all isoforms
    """
    def find_constitutive_element( self, gene_sym, bool_intron = True ):
        """
        Args:
            gene_sym = string that is the gene symbol, only record elements that belong to the same gene
            bool_intron = boolean where:
                -True = will look for constitutive introns in the gene (i.e. splice junctions that are present across all isoforms)
                -False = will look for constitutive exons in the gene (exons present across all isoforms of the gene)
        Function: find the elements that are constitutive elements (i.e. present in all isoforms) for a given gene 'gene_sym' - returns a list of Exon objects that are constitutive elements in the gene
        """
        #go through all elements, and record count & the Exon instance -> only record elements that are present 
        #retrieve all isoforms
        # gene = Isoform.obj_cruzdb.refGene.filter_by( name2 = gene_sym ).all()
        # total_isoform_count = len( [x for x in gene if x.name2 == gene_sym ] )
        total_isoform_count = len( [k for k,v in self.hash_isoforms.iteritems() if v.gene_sym == gene_sym] )

        # elem_isoform_presence = {}      #k = genomic position (chrom:start-end), v = number of isoforms that have contain the isoform
        # elem_objs = {}   #k = genomic position (chrom:start-end), v = Exon object associated with genomic position
        # for k,v in self.hash_isoforms.iteritems():      #k = string that is isoform_id, v = Isoform instance
        #     hash_elem = v.hashIntronList if bool_intron else v.hashExonList

        #     if v.gene_sym != gene_sym:
        #         continue

        #     for k2,v2 in hash_elem.iteritems():    #k2 = string element pos (chrom:start-end), v2 = Exon object
        #         if not k2 in elem_isoform_presence:
        #             elem_isoform_presence[k2] = 1
        #             elem_objs[k2] = v2
        #         else:
        #             elem_isoform_presence[k2]+= 1

        #get information about exon prevalence
        [elem_prevalence, elem_obj] = self.all_exon_prevalence( gene_sym, bool_intron )     #exon_prevalence = hash where key = string, genomic position (format = chrom:start-end), v = integer, number of isoforms with exon, elem_obj = hash where key = string, genomic position (format = chrom:start-end), v = Exon instance

        #determine which elements are constitutive by comparing the isoform count to the 
        hash_constitutive_elems = {}            #k = string that is genomic position of element, v = Exon instance
        for k2,v2 in elem_prevalence.iteritems():       #k = string that is genomic position of element, v = integer that is the number of elements containing this isoform
            if v2 == total_isoform_count:
                hash_constitutive_elems[k2] = elem_obj[k2]

        return hash_constitutive_elems


    def find_alternative_element( self, gene_sym, bool_intron = True ):
        """
        Args:
            gene_sym = string that is the gene symbol, only record elements that belong to the same gene
            bool_intron = boolean where:
                -True = will look for alternative introns in the gene (i.e. splice junctions that are present across all isoforms)
                -False = will look for alternative exons in the gene (exons present across all isoforms of the gene)
        Function: find the elements that are alternative elements (i.e. present in all isoforms) for a given gene 'gene_sym'
        """
        #go through all elements, and record count & the Exon instance -> only record elements that are present 
        #retrieve all isoforms
        # gene = Isoform.obj_cruzdb.refGene.filter_by( name2 = gene_sym ).all()
        # total_isoform_count = len( [x for x in gene if x.name2 == gene_sym ] )
        total_isoform_count = len( [k for k,v in self.hash_isoforms.iteritems() if v.gene_sym == gene_sym] )

        # elem_isoform_presence = {}      #k = genomic position (chrom:start-end), v = number of isoforms that have contain the isoform
        # elem_objs = {}   #k = genomic position (chrom:start-end), v = Exon object associated with genomic position
        # for k,v in self.hash_isoforms.iteritems():      #k = string that is isoform_id, v = Isoform instance
        #     hash_elem = v.hashIntronList if bool_intron else v.hashExonList

        #     if v.gene_sym != gene_sym:
        #         continue

        #     for k2,v2 in hash_elem.iteritems():    #k2 = string element pos (chrom:start-end), v2 = Exon object
        #         if not k2 in elem_isoform_presence:
        #             elem_isoform_presence[k2] = 1
        #             elem_objs[k2] = v2
        #         else:
        #             elem_isoform_presence[k2]+= 1

        #get information about exon prevalence
        [elem_prevalence, elem_obj] = self.all_exon_prevalence( gene_sym, bool_intron )     #exon_prevalence = hash where key = string, genomic position (format = chrom:start-end), v = integer, number of isoforms with exon, elem_obj = hash where key = string, genomic position (format = chrom:start-end), v = Exon instance

        #determine which elements are alternative by comparing the isoform count to the 
        hash_alternative_elems = {}         #k = string that is genomic position of element, v = Exon instance
        for k2,v2 in elem_prevalence.iteritems():       #k = string that is genomic position of element, v = integer that is the number of elements containing this isoform
            if v2 < total_isoform_count:
                hash_alternative_elems[k2] = elem_obj[k2]

        return hash_alternative_elems

    def find_overlapping_elements( self, genome_pos, gene_sym, bool_local = True, bool_intron = False ):
        """
        Args:
            genome_pos = string that is the genomic position (in the format chrom:start-end)
            gene_sym = string that is the gene symbol, only record elements that belong to the same gene
            bool_local = boolean where:
                -True = only look at elements that overlap genome_pos
                -False = find longest element that overlaps genome_pos, then find all elements
            bool_intron = boolean where:
                -True = will look for constitutive introns in the gene (i.e. splice junctions that are present across all isoforms)
                -False = will look for constitutive exons in the gene (exons present across all isoforms of the gene)
        Function: find overlapping elements (e.g. exons, introns) and returns hash of all overlapping elements, where key = genomic range & value = Exon object.
        NOTE: This only retrieves canonical elements (e.g. canonical exons, canonical SJ from canonical introns)
        NOTE: REASON THIS IS NOT in Isoform is because I want to use this for alternative splicing instances. Basically I want to see overlap across all isoforms for a given gene.
        """
        #get the element of interest -> find elements that overlap this element (this element overlaps a part of another element)
        hash_genome_pos = Isoform.split_genome_pos( genome_pos )
        pos_range = set( range(hash_genome_pos['start'], hash_genome_pos['end']) )

        hash_overlapped_elem = {}       #records all elements (introns if bool_intron is True, else exons) that overlap, where key = genomic range of element & value = Exon object
        for k,v in self.hash_isoforms.iteritems():      #k = string that is isoform_id, v = Isoform instance
            hash_elem = v.hashIntronList if bool_intron else v.hashExonList
            #if not same gene symbol, then skip
            if v.gene_sym != gene_sym:
                continue

            #check if element is within
            for k2,v2 in hash_elem.iteritems():    #k2 = string element pos (chrom:start-end), v2 = Exon object
                elem_range = range( v2.exonPos.location.start, v2.exonPos.location.end )

                ##TEST:: print "MultiIsoform.overlap_elem: genome_pos = ", genome_pos, " & elem_range = ", v2, " & intersection length = ", len( pos_range.intersection( elem_range ) )

                if pos_range.intersection( elem_range ):
                    #record overlapping element
                    key = v2.str_genomic_pos( True )
                    hash_overlapped_elem[key] = v2      #this ensures that no duplicate Exon instances (intron or exon) are in the hash

        #find longest element that has overlap, and then return the elements that overlap with the longest element
        if not bool_local:
            max_key = max( hash_overlapped_elem, key = lambda k: hash_overlapped_elem[k] )
            # max( hash_overlapped_elem, key = hash_overlapped_elem.get )       #this returns max key & value
            return self.find_overlapping_elements( max_key, gene_sym, True, bool_intron  )

            
        #return all elements that overlap genome_pos
        return hash_overlapped_elem

    def get_exon_splice_sites_all_isoforms( self, position, which_ss = None ):
        """
        Args:
            position = integer that is the position of interest -> will use this to find the exon that contains this position
            which_ss = string that can have the following values:
                -None = keeps both splice sites (lower_ss & higher_ss, where ss = splice site).
                -'donor' = keeps the splice site on the donor side (this depends on strand sign)
                -'acceptor' = keeps the splice site on the acceptor side (this depends on strand sign)
        Function: retrieves all the canonical splice sites (either 5' or 3' splice site) for each isoform and returns the splice site position for each isoform
        """
        hash_isoform_splice_sites = {}      #k = isoform ID, v = integer that is the posiiton
        for k,v in self.hash_isoforms.iteritems():      #k = isoform ID, v = Isoform instance
            splice_sites = v.get_exon_splice_sites( position )
            if 'donor' in which_ss.lower():
                hash_isoform_splice_sites[k] = splice_sites['higher_ss'] if v.strand > 0 else splice_sites['lower_ss']
            elif 'acceptor' in which_ss.lower():
                hash_isoform_splice_sites[k] = splice_sites['lower_ss'] if v.strand > 0 else splice_sites['higher_ss']
            else:
                hash_isoform_splice_sites[k] = splice_sites     #pos_ss = position splice site

        return hash_isoform_splice_sites


    
    ##MAYBE I CAN DELETE THIS - find_constitutive_element() takes care of both constitutive exons & introns
    def find_constitutive_exons( self, gene_sym ):
        """
        Function: find the exons that are constitutive exons (i.e. present in all isoforms) for a given gene 'gene_sym'
        """
        #go through all exons, and record count & the Exon instance -> only record exons that are present 
        #retrieve all isoforms
        # gene = Isoform.obj_cruzdb.refGene.filter_by( name2 = gene_sym ).all()
        # total_isoform_count = len( [x for x in gene if x.name2 == gene_sym ] )
        total_isoform_count = len( [k for k,v in self.hash_isoforms.iteritems() if v.gene_sym == gene_sym] )

        exon_isoform_presence = {}      #k = genomic position (chrom:start-end), v = number of isoforms that have contain the isoform
        exon_objs = {}   #k = genomic position (chrom:start-end), v = Exon object associated with genomic position
        for k,v in self.hash_isoforms.iteritems():      #k = string that is isoform_id, v = Isoform instance
            if v.gene_sym != gene_sym:
                continue

            for k2,v2 in v.hashExonList.iteritems():    #k2 = string exon pos (chrom:start-end), v2 = Exon object
                if not k2 in exon_isoform_presence:
                    exon_isoform_presence[k2] = 1
                    exon_objs[k2] = v2
                else:
                    exon_isoform_presence[k2]+= 1

        #determine which exons are constitutive by comparing the isoform count to the 
        list_constitutive_exons = []
        for k2,v2 in exon_isoform_presence.iteritems():       #k = string that is genomic position of exon, v = integer that is the number of exons containing this isoform
            if v2 == total_isoform_count:
                list_constitutive_exons.append( exon_objs[k2] )

        return list_constitutive_exons



    def find_longest_isoform( self ):
        """
        Function: find the longest isoform based on the number of exons present. Returns an array of Isoform objects that have the most number of exons
        Assumption: the isoform with the most exons is the longest, though that may not be the case
        """
        #get all isoforms associated with gene_sym -> count the numnber of exons for each isoform -> return the 
        max_exons = 0
        list_isoforms_max_exons = []

        for k,v in self.hash_isoforms.iteritems():      #k = string that is the isoform ID, v = Isoform object
            if max_exons < len( v.hashExonList ):
                list_isoforms_max_exons = [v]
            elif max_exons == len( v.hashExonList ):
                list_isoforms_max_exons.append( v )

        return list_isoforms_max_exons


    def in_exon_isoform( self, position, isoform_id ):
        """ Function: returns Exon object of element that contains position in specific isoform, else returns 'None' if no exons found """
        return self.hash_isoforms[isoform].in_exon( position )

    def in_exon_all( self, position ):
        """ Function: returns Exon object of element that contains position in specific isoform, else returns 'None' if no exons found """
        isoform_exon = {}       #k = isoform id (string), v = Exon object or None, depends on if position is found in exon
        for k,v in self.hash_isoforms.iteritems():      #k = isoform id (string), v = Isoform object
            isoform_exon[k] = self.in_exon_isoform( position, k )

        return isoform_exon

    def in_cds_isoform( self, position, isoform_id ):
        """ Function: returns Exon object of element that contains position in specific isoform, else returns 'None' if no exons found """
        return self.hash_isoforms[isoform_id].in_cds( position )

    def in_cds_all( self, position ):
        """ Function: returns Exon object of element that contains position in specific isoform, else returns 'None' if no exons found """
        isoform_cds = {}       #k = isoform id (string), v = Exon object or None, depends on if position is found in cds
        for k,v in self.hash_isoforms.iteritems():      #k = isoform id (string), v = Isoform object
            isoform_cds[k] = self.in_cds_isoform( position, k )
        return isoform_cds

    def in_intron_isoform( self, position, isoform_id ):
        """ Function: similar to def in_exon(), but does this for introns. returns the relative position within the exon ('exonLeft', 'exonRight', 'withinExon') """
        return self.hash_isoforms[isoform_id].in_intron( position )

    def in_intron_all( self, position ):
        """ Function: returns Exon object of element that contains position in specific isoform, else returns 'None' if no exons found """
        isoform_intron = {}       #k = isoform id (string), v = Exon object or None, depends on if position is found in cds
        for k,v in self.hash_isoforms.iteritems():      #k = isoform id (string), v = Isoform object
            isoform_intron[k] = self.in_intron_isoform( position, k )

        return isoform_intron

    """
    Function Start: calculate reading frame
    """
    def calc_reading_frame( self, position, isoform_id ):
        """ Function: return the reading frame for a specific position """
        return self.hash_isoforms[isoform_id].calc_reading_frame( position )        

    def calc_reading_frame_all( self, position ):
        """ Function: returns Exon object of element that contains position in specific isoform, else returns 'None' if no exons found """
        isoform_rf = {}       #k = isoform id (string), v = integer that is a frame number (0, 1, or 2) or None, depends on if position is found in exon
        for k,v in self.hash_isoforms.iteritems():      #k = isoform id (string), v = Isoform object
            isoform_rf[k] = self.calc_reading_frame( position, k )

        return isoform_rf

    """
    Function Start: determine if reading frame is preserved
    """
    def frame_preserved( self, pos_a, pos_b, isoform_id ):
        """ 
        Args:
            pos_a & pos_b = integer the is the genomic position of interest, where pos_a < pos_b
            isoform_id: string that isoform id (RefSeq ID)
        Function: determines if the reading frame is preserved. Returns True if reading frame is preserved, else returns false
        """
        return self.hash_isoforms[isoform_id].frame_preserved( pos_a, pos_b )

    def frame_preserved_all( self, pos_a, pos_b ):
        """ 
        Args:
            pos_a & pos_b = integer the is the genomic position of interest, where pos_a < pos_b
            isoform_id: string that isoform id (RefSeq ID) 
        Function: same as def frame_preserved() except retrieves the reading frame for all isoforms. Returns hash where key = isoform ID & value = boolean (True if reading frame preserved, else False if reading frame not preserved) 
        """
        isoform_frame = {}        #key = isoform ID, value = reading frame position (0, 1, or 2) based on the position
        for isoform in self.hash_isoforms:
            isoform_frame[isoform] = self.frame_preserved( pos_a, pos_b, isoform )

        return isoform_frame


    """
    Function: determine if canonical splicing, and if not then what type of splicing has occurred (e.g. exon skip, frameshift, intronic)
    """
    def is_canon_sj( self, pos_a, pos_b, isoform_id, other_annots = False, db_type = None ):
        """
        Args:
            start & end = integers that are start & end position, where start < end regardless of strand sign
            other_annots = boolean:
                -True = will see if SJ position is canonical in other databases. CAUTION: if going to make true, then know that I will be mixing different database annotations up (may not be a good idea)
                -False = will only use RefSeq to see if SJ is canonical (documented position)
        Function: returns True if position exists in array of donor-acceptor sites for isoform_id, else returns False """
        ##TEST:: print "\tMultiIsoform.Is_canon_sj: isoform = ", isoform_id, " & self.hash_isoforms[isoform_id] = ", self.hash_isoforms[isoform_id]
        if not db_type:
            db_type = self.db_type
        return self.hash_isoforms[isoform_id].is_canon_sj( pos_a, pos_b, other_annots, db_type )

    def is_canon_sj_all( self, pos_a, pos_b, other_annots = False, db_type = None ):
        """ 
        Args:
            start & end = integers that are start & end position, where start < end regardless of strand sign
            other_annots = boolean:
                -True = will see if SJ position is canonical in other databases. CAUTION: if going to make true, then know that I will be mixing different database annotations up (may not be a good idea)
                -False = will only use RefSeq to see if SJ is canonical (documented position)
        Function: returns hash where key = isoform id, value = boolean where True means canonical SJ & False means aberrant SJ
        """
        if not db_type:
            db_type = self.db_type

        isoform_canon_sj = {}        #key = isoform ID, value = boolean, where True means canonical SJ & False means aberrant SJ
        ##TEST:: print "\tMultiIsoform.ICSALL: self.hash_isoforms.keys() -> ", self.hash_isoforms.keys(), " & len = ", len( self.hash_isoforms.keys() )
        for isoform_id in self.hash_isoforms.keys():
            isoform_canon_sj[isoform_id] = self.is_canon_sj( pos_a, pos_b, isoform_id, other_annots, db_type )

            ##TEST:: print "\tMultiIsoform.ICSALL 2: isoform = ", isoform_id, " & isoform_canon_sj[isoform_id] = ", isoform_canon_sj[isoform_id], " for pos = ", pos_a, "-", pos_b, " & self.hash_isoforms[isoform_id] = ", self.hash_isoforms[isoform_id]

        return isoform_canon_sj

    def sj_get_ligated_exons_isoform( self, isoform_id, start, end, rel_pos_start = None, rel_pos_end = None ):
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
        return self.hash_isoforms[isoform_id].sj_get_ligated_exons( start, end, rel_pos_start = None, rel_pos_end = None )

    def sj_get_ligated_exons_all( self, start, end, rel_pos_start = None, rel_pos_end = None ):
        """
        Same as function sj_get_ligated_exons_isoform(), but will do it across all isoforms
        """
        isoform_ligated_exons = {}        #key = isoform ID, value = boolean, where True means canonical SJ & False means aberrant SJ
        for isoform in self.hash_isoforms:
            isoform_ligated_exons[isoform] = self.hash_isoforms[isoform_id].sj_get_ligated_exons( start, end, rel_pos_start = None, rel_pos_end = None )

        return isoform_ligated_exons


    def is_exon_skip( self, pos_a, pos_b, isoform_id ):
        """ Function: returns an array of all skipped exons """
        return self.hash_isoforms[isoform_id].is_exon_skip( pos_a, pos_b )

    def is_exon_skip_all( self, pos_a, pos_b ):
        """ Function: returns hash where key = isoform id, value = boolean where True means canonical SJ & False means aberrant SJ """
        isoform_exon_skip = {}        #key = isoform ID, value = boolean, where True means canonical SJ & False means aberrant SJ
        for isoform in self.hash_isoforms:
            isoform_exon_skip[isoform] = self.is_exon_skip( pos_a, pos_b, isoform )

        return isoform_exon_skip

    def is_aberrant_sj( self, pos_a, pos_b, isoform_id, intronic = False ):
        """ 
        Args:
            pos_a & pos_b = integers that are the start & end genomic position 
            isoform_id = string that is the isoform ID
            intronic = boolean
                -True = return an isoform ID even if it does land in the intronic space
                -False = do not return isoform ID if only lands in intronic space
        Function: retrieve information about start & end splicing to see if how the splice junction is aberrant, including whether splicing is in a exon (donor site, acceptor site, in exon) or intronic, any exons skipped, and any frameshifts
        """
        return self.hash_isoforms[isoform_id].is_aberrant_sj( pos_a, pos_b, intronic )


    def is_aberrant_sj_all( self, pos_a, pos_b, intronic = False ):
        """
        Args:
            pos_a & pos_b = integers that are the start & end genomic position 
            isoform_id = string that is the isoform ID
            intronic = boolean
                -True = return an isoform ID even if it does land in the intronic space
                -False = do not return isoform ID if only lands in intronic space 
        Function: returns hash where key = isoform id, value = boolean where True means canonical SJ & False means aberrant SJ """
        isoform_aberrant_sj = {}        #key = isoform ID, value = hash with elements 'exon_start', 'exon_end', 'exon skip', 'preserved_frame'
        for isoform in self.hash_isoforms:
            isoform_aberrant = self.is_aberrant_sj( pos_a, pos_b, isoform, intronic )
            if isoform_aberrant:        #if a hash is returned and not 'None', then record value
                isoform_aberrant_sj[str(isoform)] = self.is_aberrant_sj( pos_a, pos_b, isoform, intronic )

        return isoform_aberrant_sj

    """
    Function: find overlapped splice junction based given position(s)
    """
    def find_overlapped_canon_sj_isoform( self, pos_a, pos_b, isoform_id ):
        """
        Args:
            pos_a & pos_b = integers that are the start & end genomic position, respectively (usually these positions refer to the start & end of a splice junction)
            isoform_id = string that is the isoform ID
        Function: finds the nearest canonical splice junction (including overlapped) based on specific isoform. Returns the intron (Exon instance) that is nearest to one of the positions "pos_a" or "pos_b"
        NOTE: this assumes the strand sign as assigned to "self"
        """
        pos_oi = pos_b if self.mi_strand == -1 else pos_a
        #retrieve all introns for isoform
        list_introns = self.hash_isoforms[isoform_id].hashIntronList.values()
        if self.mi_strand == -1:
            #get all introns whose 3' end (location.end) is greater than the 5' end of the SJ
            select_introns = [x for x in list_introns if x.exonPos.location.start < pos_oi]
            #if no introns found, then return None, else return Intron with lowest numerical genomic position (as this will be the closest to the SJ)
            if not select_introns:
                return None
            return max( select_introns, key = lambda y: y.exonPos.location.end )
        else:
            #get all introns whose 3' end (location.end) is greater than the 5' end of the SJ. The reason it is "> pos_oi" and not ">= pos_oi" is because if
            select_introns = [x for x in list_introns if x.exonPos.location.end > pos_oi]
            #if no introns found, then return None, else return Intron with lowest numerical genomic position (as this will be the closest to the SJ)
            if not select_introns:
                return None
            return min( select_introns, key = lambda y: y.exonPos.location.start )


    def find_overlapped_canon_sj_all( self, pos_a, pos_b ):
        """
        Function: same as def find_overlapped_canon_sj_isoform(), but will return a hash where the key = isoform_id, value = Exon object which is the closest intron 
        """
        hash_overlapped_sj = {}     #k = isoform_id, v = Exon object which is the closest intron 
        for k,v in self.hash_isoforms.iteritems():      #k = isoform_id, v = Isoform instance of isoform
            hash_overlapped_sj[k] = self.find_overlapped_canon_sj_isoform( pos_a, pos_b, k )

        return hash_overlapped_sj


    def calc_gene_rpkm( self, bam_reader, library_size ):
        """
        Args:
            bam_reader = HTSeq.BAM_Reader instance
            library_size = integer that is the total number of mapped exons in the library
        Function: calculates the expression of the gene (takes the longest isoform, i.e. the isoform with the most exons)
        """
        #get the longest isoform
        list_isoform_id = self.find_longest_isoform()
        return self.hash_isoforms[ list_isoform_id[0].isoform_id ].calc_isoform_rpkm_exons_only( bam_reader, library_size ) if list_isoform_id else None



    """
    Other functions: static
    """
    # @staticmethod
    # def split_genome_pos( str_gene_range ):
    #     """
    #     Args:
    #         str_gene_range = string in the format (chrom:start-end)
    #     Function: splits string 'str_gene_range' (format: chrom:start-end) into individual elements
    #     """
    #     num_pos = str_gene_range.split(':')[1].split('-')
    #     exon_len = int( num_pos[1] ) - int( num_pos[0] )
    #     return {'chrom': str_gene_range.split(':')[0], 'start': int( num_pos[0] ), 'end': int( num_pos[1] ), 'exon_len': exon_len }

