#/usr/bin/python

#import libraries
import tabix
from cruzdb import Genome

from Bio.SeqFeature import SeqFeature, FeatureLocation
# from Bio.SeqFeature import SeqFeature, FeatureLocation, BeforePosition, AfterPosition

class Transcript:
    def __init__( self ):
        pass

    def reconstruct_transcript( self ):
        """ Function: reconstructs transcript based on WHAT? - exons or splice junctions?? Or exon objects connect by splice junctions?? """
        pass

    def quantify_penalty( self, max_penalty, isoform_id ):
        """ Function: calculates number of splice junctions that violate the canonical form of the transcript """
        pass

class CollectionSJ( SpliceJunction ):
    """ retrieves all splice junctions associated with a gene """

class Genes( Gene ):
    """ records all genes """
    def __init__():
        pass





class ClassName(object):
    """docstring for ClassName"""
    def __init__(self, arg):
        super(ClassName, self).__init__()
        self.arg = arg
        

