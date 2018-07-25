#/usr/bin/python
#__init.py__ file for package SVSv4 (Splice Variant System version 4)

from Exon import Exon
from Isoform import Isoform
from MultiIsoform import MultiIsoform
from SpliceJunction import SpliceJunction
from IsoformSJ import IsoformSJ
# from TranscribeTranslate import TranscribeTranscript, TranslateTranscript
# from TranscribeTranslate_V2 import TranscribeTranscript, TranslateTranscript
# from TranscribeTranslate_V3 import TranscribeTranscript, TranslateTranscript
# from TranscribeTranslate_V4 import TranscribeTranscript, TranslateTranscript
from TranscribeTranslate_V5 import TranscribeTranscript, TranslateTranscript
from SimpleTT import SimpleTranscribe, SimpleTranslate

from IsoformFusion import IsoformFusion, MultiIsoformFusion, KinaseFusion

#KinaseSplice is similar to KinaseFusion, but is used to find aberrant splicing affecting kinases
from KinaseSplice import KinaseSplice

#Determine prevalence of splicing events - currently used in IsoformSJ as well
from SJPrevalence import SJPrevalence


from GeneSNV import GeneSNV
from GenomicVariant import GenomicVariant

#Determine gene expression & gene expression percentile
from CalcGeneExpression import CalcGeneExpression

#Import Library for Ensembl's REST API - access VEP (Variant Effect Predictor)
from EnsemblVEP import EnsemblVEP
# from VEPIsoform import VEPIsoform
# from VEPIsoformSJ import VEPIsoformSJ
# from VEPTranscribeTranslate_V5 import VEPTranscribeTranscript, VEPTranslateTranscript

#Epitope-MHC prediction algorithms
from NeoepitopeMHC import NeoepitopeMHC
from MHC_IEDB import MHC_IEDB
from MHC_IEDB_V2 import MHC_IEDB_V2

from SimpleNeoepitope import SimpleNeoepitopeAll, SimpleNeoepitopeIsoform
from SimpleNeoepitopeV2 import SimpleNeoepitopeAllV2, SimpleNeoepitopeIsoformV2       ##NEED TO TEST BEFORE USING
from SimpleNeoepitopeIsoformSJ import SimpleNeoepitopeIsoformSJ         #this the the SJ of SimpleNeoepitope, so where SimpleNeoepitopeV2 = mutations for neoepitopes (SNVs, indels), SimpleNeoepitopeIsoformSJ = aberrant SJs for neoepitopes

#Checking TAP & proteasomal scoring percentile
from CompareProteasomeTAP import CompareProteasomeTAP

#checking frequency of peptide occurrence in endogenous genes
from ProtInfoResource import ProtInfoResource

#this is only for files in path in Atlas: /home/mokha/Documents/Krauthammer_Lab/160427_SJFalsePositives/Data/Tabix_Unique_Count & algorithm 161130_ComparativeSJ/Algorithms/161204_AberrantSJs_TopXGenes_V4.py
from IsoformSJ_UniqFile import IsoformSJ_UniqFile

#this will retrieve neoepitope statistics from file (used for Anti-PD1, Veliparib dataset, TCGA, etc.)
from SuiteNeoepStatistics import NeoepStatistics

"""
-SVS Improvement 1: do not assign a MultiIsoform instance for each SJ instance --> too much memory
    -Idea 1: create a MultiIsoform instance for each gene symbol -> but only make 1 copy (perhaps this should be a class variable in SJ class)
    -Idea 2: create an instance for MultiIsoform that records all the hash_isoforms recorded for all MultiIsoform instances -> then I can refer to which isoform is referring to which array element in SpliceJunction.inst_mi
    -Needs:
        -need to record the isoform where the SJ is canonical (just need string, don't need object)
        -need a way to refer to MultiIsoform instance that contains isoform ID: can retrieve gene_sym, other isoforms, etc.
        -MultiIsoform will still record all isoforms associated with gene symbol
"""



"""
Exon
-perhaps need to save annotations associated with exon - see which types of exons are affected
-Q: should I have a way to save or consider mutations, indels, etc.?
"""