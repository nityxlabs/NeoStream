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


"""
NOTE: THIS CLASS IS NOT DONE YET. NEED THE FOLLOWING:
-determine if aberrant splicing occurs before, within, or after the exons that encode the kinase domain
"""

class KinaseSplice( Isoform ):
    """
    Class: records all gene fusions that contain at least 1 kinase
    Notes on KinaseFusion:
    -id_type: column for identify isoform id
        -1 = RefSeq isoform column
        -2 = CCDS isoform column
        -3 = Uniprot isoform column
    """
    kinase_index = None
    def __init__( self, isoform_id, obj_sj ):
        """
        Args:
            hash_fusion = hash that contains the following keys to populate the object
                orientation = string that is the orientation of the fusion
                isoform_1 & isoform_2 = string that is the isoform_id part of the fusion
                chrom_1 & chrom_2 = string in the format 'chrNum' (e.g. chr9, chr12)
                pos_1 & pos_2 = integer that is the genomic position of the gene fusion
            obj_cruzdb = cruzdb.Genome instance of the genome build of interest (e.g. hg19, hg38)
        """
        #make sure a cruzdb instance is set for the Isoform class
        hash_pos = {'chrom': obj_sj.chrom, 'pos_oi': obj_sj.start}
        super( KinaseSplice, self ).__init__( isoform_id, hash_pos )

        #retrieve dataframe that contains the isoform ID
        self.df_kinase = KinaseSplice.isoform_dataframe( self.isoform_id )

    @classmethod
    def set_kinasefile( cls_obj, path_file ):
        """
        Args:
            cls_obj: should be 'KinaseSplice'
            path_file: absolute path to kinase file
        Function: sets the file path for the panda dataframe to open a file that contains all annotations associated with kinases
        NOTE: usually the path to the kinase annotation file is /home/mokha/Documents/Krauthammer_Lab/160510_GeneFusions/Data/160910_KinaseAnnots_hg38_Final.txt
        """
        cls_obj.kinase_index = pd.read_csv( path_file, sep = '\t', header = 0 )

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


    """
    Functions: find gene features associated with position & kinase domain-coding exons
    """
    def find_relative_kd( self, feature_num ):       #relative_kd_pos = relative kinase domain position
        """
        Args:
            feature_num = integer that is the feature number
        Function: returns a string that signifies where the integer 'feature_num' is with respect to the exons the encode the kinase domain, this is retrieved from the file set by def set_kinasefile(). Returns a number (-1, 0, 1) that signifies if position is before (-1), within (0), or after (1) the exons that encode the kinase domain, else returns 'no_kinase_domain' if no kinase domain found
        NOTE: the file set by def set_kinasefile() has the kinase domain position labeled with respect to the gene strand sign
        """

        col_isoform = KinaseSplice.get_column_isoform( isoform )
        find_match = KinaseSplice.kinase_index[ (KinaseSplice.kinase_index[col_isoform] == self.isoform) & (KinaseSplice.kinase_index[c_exon] == str(feature_num) ) ]
        #feature_kd = feature kinase domain, this will contain the value that determines if position is before (-1), within (0), or after (1) the exons that encode the kinase domain.
        feature_kd = str( find_match[c_kinase_domain].values[0] ) if not find_match.empty else 'no_kinase_domain'

        return feature_kd

    def locate_isoform_feature( self, position ):
        """
        Args:
            position = integer that is the genomic position
        Function: retrieve the gene feature (exon or intron) that contains "position", returns an Exon object that contains integer 'position' 
        """
        #see if position is found in exon
        containing_elem = self.get_element_v2( position, 1, 0 )
        if containing_elem:
            return containing_elem

        #else if position is not found in exon, look in intron
        containing_elem = self.get_element_v2( position, 3, 0 )
        if containing_elem:
            return containing_elem

        #if no gene feature is associated with position, then return "None"
        return None

    def pos_relative_kd( self, position ):
        """
        Args:
            position = integer that is the genomic position
        Function: uses a combination of both def locate_isoform_feature() & def find_relative_kd() to retrieve where the integer 'position'
        """
        get_elem = self.locate_isoform_feature( position )
        if not get_elem:
            return None

        return self.find_relative_kd( get_elem.exonNum )

