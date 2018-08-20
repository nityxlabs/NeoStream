#/usr/bin/python
import sys

from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv5 import MultiIsoformFusion, IsoformFusion, KinaseFusion, Isoform
from mokhaPy import mokhaPy

DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv5"
DIR_DATA = DIR_CURR + "/TestData"
DIR_RESULTS = DIR_CURR + "/TestResults"
# DIR_RESULTS = DIR_CURR + "/Results/160729_Analyze_KF"
# DIR_RESULTS = DIR_CURR + "/Results/160731_Analyze_KF"
# DIR_RESULTS = DIR_CURR + "/Results/160909_Analyze_KF"

DIR_FUSION = DIR_PROJ + "/160510_GeneFusions"

print "------------ TDD: 161002_IsoformFusion_1.py ------------"

#set kinase gene annotation file
KinaseFusion.set_kinasefile( DIR_FUSION + "/Data/160910_KinaseAnnots_hg38_Final.txt" )
obj_cruzdb = Genome( 'sqlite:////tmp/hg38_v2.db' )
#set cruzdb Genome database instance
Isoform.set_cruzdb( obj_cruzdb )

#CASE: This returns "None" for the kinase domain for the kinase gene (TLK2 - NM_001284363)
#Fusion - ASIC2:TLK2
hash_multi_isoform = { "orientation": 'fr',
"chrom_start": 'chr17',
"chrom_end": 'chr17',
"pos_start": 34038904,
"pos_end": 62565136,
"read_span": 5, 
"read_matepair": 5,
"read_matepair_break": 5 }

obj_mif = MultiIsoformFusion( hash_multi_isoform )      #MIF = MultiIsoform Fusion instance

for i, (k,v) in enumerate( obj_mif.isoform_fusions.iteritems() ):         #k = isoformIDs in fusion (format: isoformID_1:isoformID_2), v = IsoformFusion instance
    feat_info_1 = v.return_feature( True )
    feat_info_2 = v.return_feature( False )

    print "num ", i, " & kinase stat = ", v.kinase_stat
    print "v.kinase2 = ", v.kinase2, " & v.isoform2 = ",
    print "fusion gene 1 = ", feat_info_1
    print "fusion gene 2 = ", feat_info_2
    print '\n\n'


print "------------ TDD Completed: 161002_IsoformFusion_1.py ------------"

