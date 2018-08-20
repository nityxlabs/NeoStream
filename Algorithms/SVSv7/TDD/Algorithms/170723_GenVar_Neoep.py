#/usr/bin/python
import sys
import time

import scipy
from scipy import stats
import pandas as pd

from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv6 import Isoform, GenomicVariant
from mokhaPy import mokhaPy

#Constants - directories
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv6"
DIR_DATA = DIR_CURR + "/TestData"
DIR_RESULTS = DIR_CURR + "/TestResults"

DIR_DATA_GENEVAR = DIR_PROJ + "/170304_NeoantigenExtract/Data"
DIR_DATA_VELIP = DIR_DATA_GENEVAR + "/Velip"
DIR_RESULTS = DIR_CURR + "/Results"
#directory for genome
DIR_GENOME = DIR_PROJ + '/ArchiveData/hg19.fa'      #directory for samtool-indexed genome


print "------------ TDD: 170723_GenVar_Neoep.py ------------"

g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )

df_express = pd.read_csv( DIR_DATA_VELIP + "/170602_Velip_GeneCounts.txt", sep = '\t')

df_muts = pd.read_csv( DIR_DATA_VELIP + "/170602_Velip_WGS.txt", skiprows = 1, sep = '\t' )      #whole genome sequencing

for i, (i_row, row) in enumerate( df_muts.iterrows() ):
    count_skip = 3153
    if i < count_skip:
        continue
    if i > (count_skip + 1):
        break

    #Hugo_Symbol, Chromosome, Start_Position, End_Position, Strand, Reference_Allele, Tumor_Seq_Allele1   Tumor_Seq_Allele2, MAYBE: Variant_Type (this describes the type of mutation)

    print "row info = ", row
    print "gene_sym = ", row['Hugo_Symbol'], " | mutation = ", row['HGVSc'], " | codon = ", row['Codons'], " | AA change = ", row['HGVSp_Short']

    #Tumor_Seq_Allele1
    # obj_gv = GenomicVariant( row['Chromosome'], row['Start_Position'], row['End_Position'], row['Strand'], row['Reference_Allele'], row['Tumor_Seq_Allele1'], DIR_GENOME, row['Hugo_Symbol'], None )
    #Tumor_Seq_Allele2
    # obj_gv = GenomicVariant( row['Chromosome'], row['Start_Position'], row['End_Position'], row['Strand'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type'], DIR_GENOME, row['Hugo_Symbol'], None )
    obj_gv = GenomicVariant( row['Chromosome'], row['Start_Position'], row['End_Position'], row['Strand'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type'], DIR_GENOME, None, None )

    ##TEST::
    print "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "STANDARD: ROW ", i, " - row = ", row['Hugo_Symbol'], " & strand = ", row['Strand']
    print "obj_gv = ", obj_gv.obj_mi.hash_isoforms.keys(), " & AA change = ", row['HGVSp_Short']

    obj_gv.determine_aa_change()

##TEST::
# get_gene = Isoform.obj_cruzdb.bin_query( 'refGene', 'chr10', 30316219, 30316219 ).all()
# get_gene = Isoform.obj_cruzdb.bin_query( 'ensGene', 'chr10', 30316219, 30316219 ).all()

# print "get_gene = ", get_gene




print "------------ TDD Completed: 170723_GenVar_Neoep.py ------------"