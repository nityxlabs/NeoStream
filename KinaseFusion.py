#/usr/bin/python
import pandas as pd

from Isoform import Isoform

#columns for kinase index
c_gene_sym = 'geneName_short'
c_chrom = 'chrNum'
c_exon = 'exonNum'
c_start = 'exon_posStart'
c_end = 'exon_posEnd'
c_strand = 'strandSign'
c_isoform_refseq = 'geneID (RefSeq)'
c_isoform_ccds = 'geneID (CCDS_ID)'
c_isoform_uniprot = 'geneID (Uniprot)'
c_kinase_domain = 'exon_statKinaseDomain'


class KinaseFusion():
    """
    Class: records all gene fusions that contain at least 1 kinase
    Notes on KinaseFusion:
    -id_type: column for identify isoform id
        -1 = RefSeq isoform column
        -2 = CCDS isoform column
        -3 = Uniprot isoform column
    """
    kinase_index = None
    def __init__( self, hash_fusion, obj_cruzdb ):
        """
        Args:
            hash_fusion = hash that contains the following keys to populate the object
                orientation = string that is the orientation of the fusion
                isoform_1 & isoform_2 = string that is the isoform_id part of the fusion
                chrom_1 & chrom_2 = string in the format 'chrNum' (e.g. chr9, chr12)
                pos_1 & pos_2 = integer that is the genomic position of the gene fusion
            obj_cruzdb = cruzdb.Genome instance of the genome build of interest (e.g. hg19, hg38)
        """
        self.orientation = hash_fusion['orientation']
        self.orientation_gene_range = self.get_orientation_gene_range()     #will be use to determine which part of the gene to retrieve
        self.kinase_num = hash_fusion['kinase_num']     #1 = first gene is kinase, 2 = second gene is kinase, 3 = both genes are kinases

        self.isoform_1 = hash_fusion['isoform_1']
        self.gene_sym_1 = KinaseFusion.get_gene_sym( self.isoform_1 )
        self.isoform_2 = hash_fusion['isoform_2']
        self.gene_sym_2 = KinaseFusion.get_gene_sym( self.isoform_2 )

        self.chrom_1 = hash_fusion['chrom_1']
        self.chrom_2 = hash_fusion['chrom_2']
        self.pos_1 = int( hash_fusion['pos_1'] )
        self.pos_2 = int( hash_fusion['pos_2'] )

        #retrieve dataframe that contains the isoform ID
        self.df_1 = KinaseFusion.isoform_dataframe( self.isoform_1 )
        self.df_2 = KinaseFusion.isoform_dataframe( self.isoform_2 )

        #set cruzdb Genome database instance
        Isoform.set_cruzdb( obj_cruzdb )


    @classmethod
    def isoform_dataframe( cls_obj, isoform_id ):
        """
        Args:
            isoform_id = string that is the isoform ID (either RefSeq, CCDS, or Uniprot)
        Function: return a subset of the DataFrame that contains the isoform ID of interest
        """
        df = cls_obj.kinase_index

        col_isoform = cls_obj.isoform_id_get_column( isoform_id )
        return df[ df[col_isoform] == isoform_id ]

    def get_orientation_gene_range( self ):
        """
        Function: based on fusion orientation, this function returns True or False for whether the genomic region before the fusion break is retrieve for both genes. True = retrieve exons between the lowest genomic position (start exon for + gene & end exon for - gene) & the fusion break point, & False = retrieve exons between the fusion break point & the high genomic position (end exon for + gene & start exon for - gene)
        """
        if self.orientation == 'ff':
            return ( True, False )      #this means for first gene = retrieve genomic position before fusion break, second gene = retrieve genomic position after fusion break
        elif self.orientation == 'rr':
            return ( False, True )
        elif self.orientation == 'fr':
            return ( True, True )
        elif self.orientation == 'rf':
            return ( False, False )

    @classmethod
    def set_kinasefile( cls_obj, path_file ):
        """
        Args:
            cls_obj: should be 'KinaseFusion'
            path_file: absolute path to kinase file
        Function: sets the file path for the panda dataframe to open a file that contains all annotations associated with kinases
        NOTE: usually the path to the kinase annotation file is /home/mokha/Documents/Krauthammer_Lab/160510_GeneFusions/Data/160910_KinaseAnnots_hg38_Final.txt
        """
        cls_obj.kinase_index = pd.read_csv( path_file, sep = '\t', header = 0 )

    @staticmethod
    def detect_isoform_id_type( isoform_id ):
        """
        Arg:
            isoform_id = string that is the isoform ID, should be one of the follwoign of the following
        Function: returns the type of isoform ID (CCDS, RefSeq, Uniprot)
        NOTE: id_type = integer that designates the column to retrieve (see top for notes on id_type)
            -1 = RefSeq isoform column
            -2 = CCDS isoform column
            -3 = Uniprot isoform column
        """
        if 'NM_' in isoform_id or 'NR_' in isoform_id:      #RefSeq
            return 1
        elif 'CCDS' in isoform_id:      #CCDS ID
            return 2
        else:       #usually Uniprot ID starts with 'P' or 'Q'
            return 3


    @staticmethod
    def get_column_isoform( id_type = 1 ):
        """
        Arg:
            gene_sym = string that is the gene symbol (e.g. BRAF, RAF1)
            id_type = integer that designates the column to retrieve (see top for notes on id_type)
        Function: retrieves all isoform IDs (CCDS IDs) associated with gene symbol 'gene_sym'
        """
        if id_type == 2:
            return c_isoform_ccds
        elif id_type == 3:
            return c_isoform_uniprot
        else:
            return c_isoform_refseq

    @staticmethod
    def isoform_id_get_column( isoform_id ):
        """
        Function: returns the column based on the isoform ID (basically combining def detect_isoform_id_type() & get_column_isoform() )
        """
        id_type = KinaseFusion.detect_isoform_id_type( isoform_id )
        return KinaseFusion.get_column_isoform( id_type )


    """
    Get kinase information
    """
    @classmethod
    def check_kinase( cls_obj, chrom, position ):
        """
        Args:
            chrom: string in the format 'chr#' (e.g chr2, chr9)
        Function: determines if genomic position contains kinase gene(s) """
        if not Isoform.obj_cruzdb:
            print "Error: need to set cruzdb object."
        
        genes = Isoform.obj_cruzdb.bin_query( 'refGene', chrom, position, position ).all()
        
        kinase_genes = []
        for gene in genes:
            if cls_obj.kinase_sym( gene.name2 ):
                kinase_genes.append( str( gene.name2 ) )

        #return a unique list of kinase names
        return list( set( kinase_genes ) )

    @classmethod
    def kinase_sym( cls_obj, gene_sym ):
        """ Function: returns true if gene fusion contains kinase, else returns false """
        #find all the gene symbols
        if gene_sym in cls_obj.kinase_index[c_gene_sym].unique():
            return True 
        else:
            return False


    def retrieve_kinase_exons( self ):
        """
        Function: returns all rows in dataframe exons coding for the kinase domain, returns tuple where [0] = # of kinase-coding exons for first gene & [1] = # of kinase-coding exons for second gene
        """
        df_kinase_exons_1 = self.df_1[ self.df_1[c_kinase_domain].astype(int) == 0 ]
        # count_2 = len( self.df_2[c_kinase_domain].astype(int) == 0 )
        df_kinase_exons_2 = self.df_2[ self.df_2[c_kinase_domain].astype(int) == 0 ]

        return ( df_kinase_exons_1, df_kinase_exons_2 )

    def count_kinase_exons( self ):
        """
        Function: quantifies the number of exons that code for the kinase domain for both genes
        """
        ( df_kinase_exons_1, df_kinase_exons_2 ) = self.retrieve_kinase_exons()
        return ( len( df_kinase_exons_1 ), len( df_kinase_exons_2 ) )


    @classmethod
    def get_isoforms( cls_obj, gene_sym, id_type = 1 ):
        """
        Arg:
            gene_sym = string that is the gene symbol (e.g. BRAF, RAF1)
            id_type = integer that designates the column to retrieve (see top for notes on id_type)
        Function: retrieves all isoform IDs (CCDS IDs) associated with gene symbol 'gene_sym'
        """ 
        col_isoform = cls_obj.get_column_isoform( id_type )
        return cls_obj.kinase_index[ cls_obj.kinase_index[c_gene_sym] == gene_sym ][col_isoform].unique()

    @classmethod
    def get_gene_sym( cls_obj, isoform_id ):
        """
        Arg:
            gene_sym = string that is the gene symbol (e.g. BRAF, RAF1)
        Function: retrieves all isoform IDs (CCDS IDs) associated with gene symbol 'gene_sym'
        """
        id_type = cls_obj.detect_isoform_id_type( isoform_id )     #get the ID type (CCDS, RefSeq, or Uniprot)
        col_isoform = cls_obj.get_column_isoform( id_type )
        # gene_sym = cls_obj.kinase_index[ cls_obj.kinase_index[col_isoform] == isoform_id ][c_gene_sym].iloc[0]
        gene_sym = cls_obj.kinase_index[ cls_obj.kinase_index[col_isoform] == isoform_id ][c_gene_sym].unique()

        return gene_sym[0] if len( gene_sym ) > 0 else None


    """
    Find location within gene
    """
    @classmethod
    def in_elem( cls_obj, elem, position ):
        """
        Args:
            elem = a hash that has information about the elem (exon or intron)
        Function: checks if position is within elem
        """
        if position == elem['start']:
            return "onLeft"
        elif position == elem['end']:
            return "onRight"
        elif elem['start'] < position < elem['end']:
            return "withinElem"
        else:
            return None


    def isoform_exon_range( self, gene_num ):
        """
        Args:
            gene_num = if 1, then will look into 
            bool_before = boolean
                -True = retrieve all positions before number 'pos'
                -False = retrieve all positions after number 'pos'
        Function: retrieves the exons within the range of start & end
        """
        if gene_num == 2:
            df = self.df_2
            pos = self.pos_2
            orientation = list( self.orientation )[1]
            bool_before = self.orientation_gene_range[1]
        else:
            df = self.df_1
            pos = self.pos_1
            orientation = list( self.orientation )[0]
            bool_before = self.orientation_gene_range[0]

        #retrieve exons associated with start & end
        df_sub = df[ df[c_start].astype(int) < pos ] if bool_before else df[ df[c_end].astype(int) > pos ]

        return df_sub

    def exon_range_kinase_coding( self, gene_num ):
        #retrieve exons within range of position (depends on orientation, start OR end position, & position of fusion break)
        df = self.exon_range( gene_num )
        return df[ df[c_kinase_domain].astype(int) == 0 ]


    @classmethod
    def get_position_in_kinase( cls_obj, gene_sym, position ):
        """
        Function: finds the position of the 
        """
        df = cls_obj.kinase_index
        #retrieve all entries with the same gene symbol
        df_gene = df[ df['c_gene_sym'] == gene_sym ]

        #if df_gene is empty then return None
        if df_gene.empty():
            return None

        #go through each position and see if position is contained within start & end
        df_gene_exons = df_gene[ df_gene[c_start].astype('int') <= position & df_gene[c_end].astype('int') >= position ]
        
        return df_gene_exons if not df_gene_exons.empty() else None


    """
    Function: retrieve information about kinases
    """
    @classmethod
    def isoform_info( cls_obj, isoform, position, db_type = 1, id_type = 1 ):
        """
        Arg:
            gene_sym = string that is the isoform ID
            position = integer that is the genomic position
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
            id_type = integer that designates the column to retrieve. This should be in line with the string 'isoform'. (see top for notes on id_type)
        Function: retrieve information about the isoform based on the position
        """

        #find the kinase feature information relative to the position
        hash_feature = cls_obj.locate_isoform_feature( isoform, position, db_type, id_type )
        
        #make sure feature_name (exon/intron) is not None
        hash_feature['gene_sym'] = cls_obj.get_gene_sym( isoform ) if hash_feature['feature_name'] else None

        return hash_feature

    @classmethod
    def kinase_info( cls_obj, chrom, position, isoform_id = None, db_type = 1, id_type = 1 ):
        """ 
        Args:
            chrom = string that is the chromosome (format: chr#, e.g. chr2, chr17)
            position = integer that is the genomic position (usually the fusion point)
            isoform_id = string that can specify the isoform ID, but if None then will look at all isoforms
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
            id_type = integer that designates the column to retrieve. (see top for notes on id_type)
        Function: returns an array of kinases associated with position 'chrom:position', where each element will contain information about kinase gene symbol, isoform ID, feature (exon/intron) where position falls into, 
        """
        kinases = cls_obj.check_kinase( chrom, position )
        kinase_features = []       #array of hashes - records information about isoform, including which exon/intron & relative position of kinase domain

        #if no kinases found, then exit function
        if not kinases:
            return kinase_features

        for kinase in kinases:      #kinase = each kinase gene symbol, kinases = array of kinase gene symbols
            #go through each isoform  - find the features (exon/intron) where position falls within each isoform isoform
            isoforms = cls_obj.get_isoforms( kinase, id_type )
            for isoform in isoforms:
                
                #if specific isoform ID is specified, then process this isoform
                if isoform_id and isoform == isoform_id:
                    hash_feature = cls_obj.isoform_info( isoform, position, db_type )
                    if hash_feature:        #if not blank hash, then 
                        kinase_features.append( hash_feature )
                elif not isoform_id:
                    hash_feature = cls_obj.isoform_info( isoform, position, db_type )
                    if hash_feature:        #if not blank hash, then 
                        kinase_features.append( hash_feature )

        return kinase_features

    """
    Functions: retrieve list of features (e.g. exon, intron) 
    """
    @classmethod
    def kinase_isoform_exons( cls_obj, isoform ):
        """ 
        Args:
            isoform = string that is the isoform id (specifically the CCDS ID)
            id_type = integer that designates the column to retrieve (see top for notes)
        Function: returns an array of start and end positions & relative position of kinase domain for each exon
        Note: for kinase domain, -1 = before kinase domain, 0 = within kinase domain, 1 = after kinase domain
        """
        col_isoform = cls_obj.isoform_id_get_column( isoform )
        df_isoforms = cls_obj.kinase_index[ cls_obj.kinase_index[col_isoform] == isoform ]

        exons = {}      #key = exon number, value = hash that contains chromosome, start, end, & relative kinase domain position
        for i, row in df_isoforms.iterrows():
            #get strand, assign with number
            convert_ss = { '+' : 1, '-' : -1 }
            strand = convert_ss.get( row[c_strand], 0 )

            exons[ int( row[c_exon] ) ] = { 
            'chrom': row[c_chrom],
            'start': int( row[c_start] ), 
            'end': int( row[c_end] ), 
            'strand': strand,
            'kinase_domain': int( row[c_kinase_domain] ) if row[c_kinase_domain].lstrip('+-').isdigit() else row[c_kinase_domain],
            }

        return exons

    @classmethod
    def isoform_introns( cls_obj, isoform, id_type = 1 ):
        """ 
        Args:
            isoform = string that is the isoform id (specifically the CCDS ID)
            id_type = integer that designates the column to retrieve (see top for notes on id_type)
        Function: returns an array of start and end positions & relative position of kinase domain for each intron
        Note: for kinase domain, -1 = before kinase domain, 0 = within kinase domain, 1 = after kinase domain
        """
        col_isoform = cls_obj.get_column_isoform( id_type )
        df_isoforms = cls_obj.kinase_index[ cls_obj.kinase_index[col_isoform] == isoform ]

        introns = {}      #key = exon number, value = hash that contains chromosome, start, end, & relative kinase domain position
        num_rows = range( 0, len( df_isoforms.index ) )
        for i in num_rows:

            if i == max( num_rows ):
                break

            #sort the start & end position from least to greatest
            # intron_pos = [ int(df_isoforms.iloc[i][c_end]), int(df_isoforms.iloc[i + 1][c_start]) ].sort( reverse = False )
            intron_pos = sorted( [ int(df_isoforms.iloc[i][c_end]), int(df_isoforms.iloc[i + 1][c_start]) ] )
            # intron_pos = intron_pos.sort()
            introns[ int( df_isoforms.iloc[i][c_exon] ) ] = { 
            'chrom': df_isoforms.iloc[i][c_chrom],
            'start': intron_pos[0], 
            'end': intron_pos[1],
            'strand': df_isoforms.iloc[i][c_strand], 
            'kinase_domain': int( df_isoforms.iloc[i][c_kinase_domain] ) if df_isoforms.iloc[i][c_kinase_domain].lstrip('+-').isdigit() else df_isoforms.iloc[i][c_kinase_domain],       #relative position of kinase domain, where -1 = before kinase domain, 0 = on kinase domain, 1 = after kinase domain
            }

        return introns


    @classmethod
    def locate_isoform_feature( cls_obj, isoform, position, db_type = 1, id_type = 1 ):
        """
        Args:
            isoform = string that is the isoform id (specifically the CCDS ID)
            position = integer that is the genomic location of interest (usually the position where the fusion occurs)
            -db_type = integer that chooses the genomic database of interest
                -1 = uses RefSeq (Isoform.obj_cruzdb.refGene)
                -2 = uses Ensembl (Isoform.obj_cruzdb.ensGene)
                -3 = uses UCSC (CONJ: I think it is Isoform.obj_cruzdb.knownGene)
            id_type = integer that designates the column to retrieve (see top for notes on id_type)
        Function: finds the gene feature (e.g. exon, intron) that contains the position of interest
        """
        #see if position is within exon
        feature_name = None
        feature_position = None     #string that is t
        feature_kd = None           #feature kinase domain - records
        relative_pos = None
        strand = None 

        hash_cds = cls_obj.kinase_isoform_exons( isoform )
        for k,v in hash_cds.iteritems():      #k = exon number, v = hash (chrom, start, end, kinase domain)
            #check where position falls within exon
            rel_pos = cls_obj.in_elem( v, position )
            if rel_pos:     #if not 'None' (meaning position is within exon)
                # feature_name = 'cds' + str(k)       #this is the real notation for CDS
                feature_name = 'exon_' + str(k)         #TEMPORARY using this for now
                relative_pos = rel_pos
                feature_position = v['chrom'] + ':' + str( v['start'] ) + '-' + str( v['end'] )
                strand = v['strand']
                feature_kd = v['kinase_domain']     #relative position of kinase domain, where -1 = before kinase domain, 0 = on kinase domain, 1 = after kinase domain
                break

        #if exon feature_name found (not None), then return feature_name
        if feature_name:
            return {'isoform': isoform, 'feature_name': feature_name, 'feature_position': feature_position, 'relative_pos': relative_pos, 'kinase_domain': feature_kd, 'strand': strand}
        else:       #else if not found in CDS, then look into exons
            #create an isoform object and see which position 
            obj_isoform = Isoform( db_type, isoform )
            get_exon = obj_isoform.get_element( position )

            if get_exon:        #if found in exon, then record information about 
                feature_num = get_exon.exonNum
                feature_name = "exon" + str( feature_num )
                # feature_num = get_exon.exonNum + 1          #as exons in Isoform start with 0, need to add +1 to make it 1-based exon

                feature_position =  get_exon.chrom + ":" + str( get_exon.exonPos.location.start ) + "-" + str( get_exon.exonPos.location.end )
                strand = get_exon.exonPos.strand
                #find relative position
                hash_elem_pos = {'start': get_exon.exonPos.location.start, 'end': get_exon.exonPos.location.end }
                relative_pos = cls_obj.in_elem( hash_elem_pos, position )


                col_isoform = KinaseFusion.get_column_isoform( isoform )
                find_match = cls_obj.kinase_index[ (cls_obj.kinase_index[col_isoform] == isoform) & (cls_obj.kinase_index[c_exon] == str(feature_num) ) ]
                feature_kd = str( find_match[c_kinase_domain].values[0] ) if not find_match.empty else 'no_kinase_domain'
                
                return {'isoform': isoform, 'feature_name': feature_name, 'feature_position': feature_position, 'relative_pos': relative_pos, 'kinase_domain': feature_kd, 'strand': strand}
            else:       #if no exons found, then look into introns
                all_introns = obj_isoform.get_intron( position )

                #if nothing found for intron, then just return blank hash, meaning no kinase information found
                if not all_introns:
                    return {}
                get_intron = all_introns[0]     #extract the first intron

                feature_num = get_intron.exonNum
                # feature_num = get_intron.exonNum + 1          #as exons in Isoform start with 0, need to add +1 to make it 1-based exon

                feature_name = "intron" + str( feature_num )
                feature_position =  get_intron.chrom + ":" + str( get_intron.exonPos.location.start ) + "-" + str( get_intron.exonPos.location.end )
                strand = get_intron.exonPos.strand
                #find relative position
                # hash_elem_pos = {'start': get_intron.exonPos.location.start, 'end': get_intron.exonPos.location.end }
                # relative_pos = cls_obj.in_elem( hash_elem_pos, position )
                relative_pos = None
                
                #retrieve the feature number (exon12, intron5) & the relative kinase domain
                col_isoform = KinaseFusion.get_column_isoform( isoform )

                ##TEST::
                print "KinaseFusion.locate_isoform_feature: get_intron.exonNum = ", get_intron.exonNum, " & featureNum = ", feature_num

                find_match = cls_obj.kinase_index[ (cls_obj.kinase_index[col_isoform] == isoform) & (cls_obj.kinase_index[c_exon] == str(feature_num) ) ]
                feature_kd = str( find_match[c_kinase_domain].values[0] ) if not find_match.empty else 'no_kinase_domain'

                return {'isoform': isoform, 'feature_name': feature_name, 'feature_position': feature_position, 'relative_pos': relative_pos, 'kinase_domain': feature_kd, 'strand': strand}

        #if nothing found, then return blank hash
        return {}


    ##MAY DELETE THIS FUNCTION SINCE I UPDATED IT (see locate_isoform_feature() )
    # @classmethod
    # def locate_isoform_feature_BACKUP( cls_obj, isoform, position, id_type = 1 ):
    #     """
    #     Args:
    #         isoform = string that is the isoform id (specifically the CCDS ID)
    #         position = integer that is the genomic location of interest (usually the position where the fusion occurs)
    #         id_type = integer that designates the column to retrieve (see top for notes on id_type)
    #     Function: finds the gene feature (e.g. exon, intron) that contains the position of interest
    #     """
    #     #see if position is within exon
    #     feature_name = None
    #     feature_position = None     #string that is t
    #     feature_kd = None           #feature kinase domain - records
    #     relative_pos = None
    #     strand = None 

    #     hash_exons = cls_obj.kinase_isoform_exons( isoform )
    #     for k,v in hash_exons.iteritems():      #k = exon number, v = hash (chrom, start, end, kinase domain)
    #         #check where position falls within exon
    #         rel_pos = cls_obj.in_elem( v, position )
    #         if rel_pos:     #if not 'None' (meaning position is within exon)
    #             feature_name = 'exon' + str(k)
    #             relative_pos = rel_pos
    #             feature_position = v['chrom'] + ':' + str( v['start'] ) + '-' + str( v['end'] )
    #             strand = v['strand']
    #             feature_kd = v['kinase_domain']     #relative position of kinase domain, where -1 = before kinase domain, 0 = on kinase domain, 1 = after kinase domain
    #             break

    #     #if exon feature_name found (not None), then return feature_name
    #     if feature_name:
    #         return {'isoform': isoform, 'feature_name': feature_name, 'feature_position': feature_position, 'relative_pos': relative_pos, 'kinase_domain': feature_kd, 'strand': strand}
    #     else:       #else if no feature_name found, then look into introns
    #         hash_introns = cls_obj.isoform_introns( isoform, id_type )

    #         for k,v in hash_introns.iteritems():

    #             rel_pos = cls_obj.in_elem( v, position )
    #             if rel_pos:     #if not 'None' (meaning position is within exon)
    #                 feature_name = 'intron' + str(k)
    #                 relative_pos = rel_pos
    #                 feature_position = v['chrom'] + ':' + str( v['start'] ) + '-' + str( v['end'] )
    #                 strand = v['strand']
    #                 feature_kd = v['kinase_domain']     #relative position of kinase domain, where -1 = before kinase domain, 0 = on kinase domain, 1 = after kinase domain
    #                 break

    #         return {'isoform': isoform, 'feature_name': feature_name, 'feature_position': feature_position, 'relative_pos': relative_pos, 'kinase_domain': feature_kd, 'strand': strand}

    #     #if nothing found, then return blank hash
    #     return {}