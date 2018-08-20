#/usr/bin/python
import sys

from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv5 import KinaseFusion, Isoform
from mokhaPy import mokhaPy

DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv5"
DIR_DATA = DIR_CURR + "/TestData"
DIR_RESULTS = DIR_CURR + "/TestResults"
# DIR_RESULTS = DIR_CURR + "/Results/160729_Analyze_KF"
# DIR_RESULTS = DIR_CURR + "/Results/160731_Analyze_KF"
# DIR_RESULTS = DIR_CURR + "/Results/160909_Analyze_KF"

DIR_FUSION = DIR_PROJ + "/160510_GeneFusions"

#columns from kinase fusion (this is my random format, not an establish format)
COL_SAMPLE = 0
COL_LINE_ID = 2
COL_CHROM = 3
COL_POS_START = 4
COL_POS_END = 5
COL_FUSION_ID = 9
COL_ORIENTATION = 10
COL_COMPATIBLE = 11
COL_KINASE_POS = 12     #0 = no kinase in fusion, 1 = first gene is kinase, 2 = 2nd gene is kinase, 3 = both are kinases
COL_FUSION_NAME = 13
COL_GENE1 = 14
COL_GENE2 = 20
COL_ISO1 = 15
COL_ISO2 = 21
COL_STRAND1 = 16
COL_STRAND2 = 22
COL_KINASE_DOMAIN_1 = 19
COL_KINASE_DOMAIN_2 = 25

print "------------ TDD: 160929_KinaseFusion_2.py ------------"

#CASE 1: fusion in intron
#CASE 2: fusion where in intron in one isoform but outside in different isoform - look at gene NEK2


#set kinase gene annotation file
KinaseFusion.set_kinasefile( DIR_FUSION + "/Data/160910_KinaseAnnots_hg38_Final.txt" )


#Fusion - NEK4:AKTIP
# hash_fusion = { "orientation": 'fr',
# "isoform_1": 'NM_001193533',
# "isoform_2": 'NM_001012398',
# "kinase_num": '1',
# "chrom_1": 'chr3',
# "chrom_2": 'chr16',
# "pos_1": '52754588',
# "pos_2": '53492518' }

#Fusion - FRK:TPI1
# hash_fusion = { "orientation": 'ff',
# "isoform_1": 'NM_002031',
# "isoform_2": 'NM_001159287',
# "kinase_num": '1',
# "chrom_1": 'chr6',
# "chrom_2": 'chr12',
# "pos_1": '116039565',
# "pos_2": '6870454' }

#CASE: Position is in one isoform but outside in the other isoforms
#Fusion - FDPS:NEK2
# hash_fusion = { "orientation": 'fr',
# "isoform_1": 'NM_002004',
# "isoform_2": 'NM_001204182',
# "kinase_num": '2',
# "chrom_1": 'chr1',
# "chrom_2": 'chr1',
# "pos_1": '155318275',
# "pos_2": '211660781' }

#CASE: Position is outside of the CDS but in the exon (not in CDS of CDK16, but in non-coding exon of CDK16)
#Fusion - RHOA:CDK16
# hash_fusion = { "orientation": 'ff',
# "isoform_1": 'NM_001664',
# "isoform_2": 'NM_001170460',
# "kinase_num": '2',
# "chrom_1": 'chr1',
# "chrom_2": 'chr10',
# "pos_1": '92301475',
# "pos_2": '86856237' }

#CASE: Position is outside of the CDS but in the exon (not in CDS of CDK16, but in non-coding exon of CDK16)
#Fusion - RAB3GAP2:AURKA
# hash_fusion = { "orientation": 'rr',
# "isoform_1": 'NM_012414',
# "isoform_2": 'NM_001323303',
# "kinase_num": '2',
# "chrom_1": 'chr1',
# "chrom_2": 'chr20',
# "pos_1": '220267207',
# "pos_2": '56373550' }

#CASE: This returns "None" for the kinase domain for the kinase gene (TLK2 - NM_001284363)
#Fusion - ASIC2:TLK2
hash_fusion = { "orientation": 'fr',
"isoform_1": 'NM_001094',
"isoform_2": 'NM_001284363',
"kinase_num": '2',
"chrom_1": 'chr17',
"chrom_2": 'chr17',
"pos_1": '34038904',
"pos_2": '62565136' }

print "hash_fusion = ", hash_fusion


obj_cruzdb = Genome( 'sqlite:////tmp/hg38_v2.db' )
obj_kf = KinaseFusion( hash_fusion, obj_cruzdb )

##TEST:: see pandas Dataframe
# print "show df_1:\n", obj_kf.df_1
# print "show df_2:\n", obj_kf.df_2
print "\n-------------------\n"

##TEST:: see exon range based on orientation of fusion
# print "exons for gene 1 = ", obj_kf.isoform_exon_range( 1 ), "\n"
# print "exons for gene 2 = ", obj_kf.isoform_exon_range( 2 ), "\n"

print "number of exons coding for kinase = ", obj_kf.count_kinase_exons()

print "BEFORE kinase search: isoform = ", obj_kf.isoform_1, " & pos = ", obj_kf.pos_1
hash_feature_1 = KinaseFusion.locate_isoform_feature( obj_kf.isoform_1, obj_kf.pos_1, 1 )
print "hash_feature for isoform 1 = ", hash_feature_1

hash_feature_2 = KinaseFusion.locate_isoform_feature( obj_kf.isoform_2, obj_kf.pos_2, 1 )
print "hash_feature for isoform 2 = ", hash_feature_2


##TDD 2: Test access
# #columns for kinase index
# c_gene_sym = 'geneName_short'
# c_chrom = 'chrNum'
# c_exon = 'exonNum'
# c_start = 'exon_posStart'
# c_end = 'exon_posEnd'
# c_strand = 'strandSign'
# c_isoform_refseq = 'geneID (RefSeq)'
# c_isoform_ccds = 'geneID (CCDS_ID)'
# c_isoform_uniprot = 'geneID (Uniprot)'
# c_kinase_domain = 'exon_statKinaseDomain'


# df_ki = KinaseFusion.kinase_index
# df_1 = df_ki[ df_ki[c_isoform_refseq] == 'NM_002031' ]
# # df_1 = df_ki[ df_ki[c_isoform_refseq] == 'NM_0017' ]

# if not df_1.empty:
#     # print "df_1 = ", df_1
#     for x in df_1[c_exon]:
#         print "x = ", x, " & type = ", type( x )
# else:
#     print "df_1 has NOTHING!"

# df_2 = df_ki[ (df_ki[c_isoform_refseq] == 'NM_002031') & (df_ki[c_exon] == '3') ]
# # df_2 = df_ki[ (df_ki[c_isoform_refseq] == 'NM_002031') & (df_ki[c_exon].astype('int') == 3) ][c_kinase_domain]
# if not df_2.empty:
#     print "df_2 = ", df_2
# else:
#     print "df_2 has NOTHING!"

# df_3 = df_ki[ (df_ki[c_isoform_refseq] == 'NM_002031') & (df_ki[c_exon] == '34') ][c_kinase_domain]
# print "df_3:\n", df_3
# print "is df_3 empty?", df_3.empty

print "------------ TDD Completed: 160929_KinaseFusion_2.py ------------"

