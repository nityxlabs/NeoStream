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
    """
    kinase_index = None
    def __init__():
        pass

    @classmethod
    def set_kinasefile( cls_obj, path_file ):
        """
        Args:
            cls_obj: should be 'KinaseFusion'
            path_file: absolute path to kinase file
        Function: sets the file path for the panda dataframe to open """
        cls_obj.kinase_index = pd.read_csv( path_file, sep = '\t', header = 0 )

    @staticmethod
    def get_column_isoform( id_type = 1 ):
        """
        Arg:
            gene_sym = string that is the gene symbol (e.g. BRAF, RAF1)
            id_type = integer that designates the column to retrieve
                -1 = get CCDS
                -2 = get RefSeq
                -3 = get Uniprot
        Function: retrieves all isoform IDs (CCDS IDs) associated with gene symbol 'gene_sym'
        """
        if id_type == 2:
            return c_isoform_ccds
        elif id_type == 3:
            return c_isoform_uniprot
        else:
            return c_isoform_refseq

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

    @classmethod
    def get_isoforms( cls_obj, gene_sym, id_type = 1 ):
        """
        Arg:
            gene_sym = string that is the gene symbol (e.g. BRAF, RAF1)
            id_type = integer that designates the column to retrieve
                -1 = get CCDS
                -2 = get RefSeq
                -3 = get Uniprot
        Function: retrieves all isoform IDs (CCDS IDs) associated with gene symbol 'gene_sym'
        """ 
        col_isoform = KinaseFusion.get_column_isoform( id_type )
        return cls_obj.kinase_index[ cls_obj.kinase_index[c_gene_sym] == gene_sym ][col_isoform].unique()

    @classmethod
    def get_gene_sym( cls_obj, isoform_id, id_type ):
        """
        Arg:
            gene_sym = string that is the gene symbol (e.g. BRAF, RAF1)
            id_type = integer that designates the column to retrieve
                -1 = get CCDS
                -2 = get RefSeq
                -3 = get Uniprot
        Function: retrieves all isoform IDs (CCDS IDs) associated with gene symbol 'gene_sym'
        """
        col_isoform = KinaseFusion.get_column_isoform( id_type )
        # gene_sym = cls_obj.kinase_index[ cls_obj.kinase_index[col_isoform] == isoform_id ][c_gene_sym].iloc[0]
        gene_sym = cls_obj.kinase_index[ cls_obj.kinase_index[col_isoform] == isoform_id ][c_gene_sym].unique()

        return gene_sym[0] if gene_sym else None

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


    """
    Function: retrieve information about kinases
    """
    @classmethod
    def isoform_info( cls_obj, isoform, position, id_type = 1 ):
        """
        Arg:
            gene_sym = string that is the isoform ID
            position = integer that is the genomic position
            id_type = integer that designates the column to retrieve. This should be in line with the string 'isoform'
                -1 = get CCDS
                -2 = get RefSeq
                -3 = get Uniprot
        Function: retrieve information about the isoform based on the position
        """

        #find the kinase feature information relative to the position
        hash_feature = cls_obj.locate_isoform_feature( isoform, position, id_type )
        
        #make sure feature_name (exon/intron) is not None
        hash_feature['gene_sym'] = KinaseFusion.get_gene_sym( isoform, id_type ) if hash_feature['feature_name'] else None

        return hash_feature

    @classmethod
    def kinase_info( cls_obj, chrom, position, isoform_id = None, id_type = 1 ):
        """ 
        Args:
            chrom = string that is the chromosome (format: chr#, e.g. chr2, chr17)
            position = integer that is the genomic position (usually the fusion point)
            isoform_id = string that can specify the isoform ID, but if None then will look at all isoforms
            id_type = integer that designates the column to retrieve
                -1 = get CCDS
                -2 = get RefSeq
                -3 = get Uniprot
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
                    hash_feature = KinaseFusion.isoform_info( isoform, position )
                    if hash_feature:        #if not blank hash, then 
                        kinase_features.append( hash_feature )
                elif not isoform_id:
                    hash_feature = KinaseFusion.isoform_info( isoform, position )
                    if hash_feature:        #if not blank hash, then 
                        kinase_features.append( hash_feature )

        return kinase_features

    """
    Functions: retrieve list of features (e.g. exon, intron) 
    """
    @classmethod
    def kinase_isoform_exons( cls_obj, isoform, id_type = 1 ):
        """ 
        Args:
            isoform = string that is the isoform id (specifically the CCDS ID)
            id_type = integer that designates the column to retrieve
                -1 = get CCDS
                -2 = get RefSeq
                -3 = get Uniprot
        Function: returns an array of start and end positions & relative position of kinase domain for each exon
        Note: for kinase domain, -1 = before kinase domain, 0 = within kinase domain, 1 = after kinase domain
        """
        col_isoform = KinaseFusion.get_column_isoform( id_type )
        df_isoforms = cls_obj.kinase_index[ cls_obj.kinase_index[col_isoform] == isoform ]

        exons = {}      #key = exon number, value = hash that contains chromosome, start, end, & relative kinase domain position
        for i, row in df_isoforms.iterrows():
            exons[ int( row[c_exon] ) ] = { 
            'chrom': row[c_chrom],
            'start': int( row[c_start] ), 
            'end': int( row[c_end] ), 
            'strand': row[c_strand],
            'kinase_domain': int( row[c_kinase_domain] ) if row[c_kinase_domain].isdigit() else row[c_kinase_domain],
            }

        return exons

    @classmethod
    def isoform_introns( cls_obj, isoform, id_type = 1 ):
        """ 
        Args:
            isoform = string that is the isoform id (specifically the CCDS ID)
            id_type = integer that designates the column to retrieve
                -1 = get CCDS
                -2 = get RefSeq
                -3 = get Uniprot
        Function: returns an array of start and end positions & relative position of kinase domain for each intron
        Note: for kinase domain, -1 = before kinase domain, 0 = within kinase domain, 1 = after kinase domain
        """
        col_isoform = KinaseFusion.get_column_isoform( id_type )
        df_isoforms = cls_obj.kinase_index[ cls_obj.kinase_index[col_isoform] == isoform ]

        introns = {}      #key = exon number, value = hash that contains chromosome, start, end, & relative kinase domain position
        num_rows = range( 0, len( df_isoforms.index ) )
        for i in num_rows:
            if i == max( num_rows ):
                break

            introns[ int( df_isoforms.iloc[i][c_exon] ) ] = { 
            'chrom': df_isoforms.iloc[i][c_chrom],
            'start': int( df_isoforms.iloc[i][c_end] ), 
            'end': int( df_isoforms.iloc[i + 1][c_start] ),
            'strand': df_isoforms.iloc[i][c_strand], 
            'kinase_domain': int( df_isoforms.iloc[i][c_kinase_domain] ) if df_isoforms.iloc[i][c_kinase_domain].isdigit() else df_isoforms.iloc[i][c_kinase_domain],       #relative position of kinase domain, where -1 = before kinase domain, 0 = on kinase domain, 1 = after kinase domain
            }

        return introns

    @classmethod
    def locate_isoform_feature( cls_obj, isoform, position, id_type = 1 ):
        """
        Args:
            isoform = string that is the isoform id (specifically the CCDS ID)
            position = integer that is the genomic location of interest (usually the position where the fusion occurs)
            id_type = integer that designates the column to retrieve
                -1 = get CCDS
                -2 = get RefSeq
                -3 = get Uniprot
        Function: finds the gene feature (e.g. exon, intron) that contains the position of interest
        """
        #see if position is within exon
        feature_name = None
        feature_position = None     #string that is t
        feature_kd = None           #feature kinase domain - records
        relative_pos = None
        strand = None 

        hash_exons = cls_obj.kinase_isoform_exons( isoform, id_type )
        for k,v in hash_exons.iteritems():      #k = exon number, v = hash (chrom, start, end, kinase domain)
            #check where position falls within exon
            rel_pos = cls_obj.in_elem( v, position )
            if rel_pos:     #if not 'None' (meaning position is within exon)
                feature_name = 'exon' + str(k)
                relative_pos = rel_pos
                feature_position = v['chrom'] + ':' + str( v['start'] ) + '-' + str( v['end'] )
                strand = v['strand']
                feature_kd = v['kinase_domain']     #relative position of kinase domain, where -1 = before kinase domain, 0 = on kinase domain, 1 = after kinase domain
                break

        #if exon feature_name found (not None), then return feature_name
        if feature_name:
            return {'isoform': isoform, 'feature_name': feature_name, 'feature_position': feature_position, 'relative_pos': relative_pos, 'kinase_domain': feature_kd, 'strand': strand}
        else:       #else if no feature_name found, then look into introns
            hash_introns = cls_obj.isoform_introns( isoform, id_type )
            for k,v in hash_introns.iteritems():

                rel_pos = cls_obj.in_elem( v, position )
                if rel_pos:     #if not 'None' (meaning position is within exon)
                    feature_name = 'intron' + str(k)
                    relative_pos = rel_pos
                    feature_position = v['chrom'] + ':' + str( v['start'] ) + '-' + str( v['end'] )
                    strand = v['strand']
                    feature_kd = v['kinase_domain']     #relative position of kinase domain, where -1 = before kinase domain, 0 = on kinase domain, 1 = after kinase domain
                    break

            return {'isoform': isoform, 'feature_name': feature_name, 'feature_position': feature_position, 'relative_pos': relative_pos, 'kinase_domain': feature_kd, 'strand': strand}

        #if nothing found, then return blank hash
        return {}