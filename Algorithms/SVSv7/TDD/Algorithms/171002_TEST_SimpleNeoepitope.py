#/usr/bin/python
import sys

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
# from SVSv7 import Isoform, IsoformSJ, TranscribeTranscript, TranslateTranscript, EnsemblVEP
from SVSv7 import Isoform, SimpleNeoepitopeAll, SimpleNeoepitopeIsoform
from mokhaPy import mokhaPy

#Constants - directories
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv7"
DIR_DATA = DIR_CURR + "/TestData"
DIR_RESULTS = DIR_CURR + "/TestResults"

DIR_GENOME = DIR_PROJ + '/ArchiveData/hg19.fa'      #directory for samtool-indexed genome

print "------------ TDD: 171002_TEST_SimpleNeoepitope.py ------------"

g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )

##IMPORTANT TEST - THE "X" AMINO ACID
# #simulate point mutation - the genomic range for this "X" amino acid is "X:153146127-153146128"
# isoform_id = "ENST00000452593"
# # genomic_range = "X:153146127-153146128"     #this is the range of the "X" amino acid
# genomic_range = "X:153146127-153146127"
# orig = None
# alt = "T"

##MUTATIONS - minus gene
# #simulate point mutation - RESULT: codon =  aAt/at & amino acids =  N/X
# isoform_id = "ENST00000376887"
# # genomic_range = "13:95815411-95815411"
# genomic_range = "13:95953564-95953564"
# orig = None
# alt = "G"

##DELETION - minus gene
# #simulate deletion
# isoform_id = "ENST00000376887"
# # genomic_range = "13:95815411-95815411"
# genomic_range = "13:95953555-95953564"
# alt = "-"

##INSERTION
#simulate insertion
# isoform_id = "ENST00000376887"
# # genomic_range = "13:95815411-95815411"
# genomic_range = "13:95953564-95953564"
# orig = '-'
# alt = "GGGGG"

##INSERTION 2 - plus isoform of gene
# isoform_id = "ENST00000452593"
# genomic_range = "X:153151280-153151281"
# orig = '-'
# alt = "CC"


##DELETIONS
# #simulate single deletion - RESULT: codon =  aAt/at & amino acids =  N/X
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765077-17765077"
# alt = "-"

# #simulate single deletion - RESULT: codon =  Aat/at & amino acids =  N/X
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765076-17765076"
# alt = "-"

# #simulate 2-base deletions - RESULT: codon =  AAt/t  & amino acids =  N/X
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765077-17765078"
# orig = None
# alt = "-"

# #simulate 3-base deletions - RESULT: codon =  taCAAt/tat  & amino acids =  YN/Y
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765077-17765079"
# orig = None
# alt = "-"

# #simulate single deletion - 'codon_change': 'cTg/cg' & 'aa_change': 'L/X' (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046908"
# orig = None
# alt = "-"     #this is just a SNV, no insertion

# #simulate 2 base deletion - 'codon_change': 'cTG/c' & 'aa_change': 'L/X' (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046909"
# orig = None
# alt = "-"     #this is just a SNV, no insertion

# #simulate 3 base deletion - 'codon_change': 'cTGGag/cag' & 'aa_change': 'LE/Q' (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046910"
# alt = "-"     #this is just a SNV, no insertion

# #simulate 4 base deletion - 'codon_change': 'cTGGAg/cg' & 'aa_change': 'LE/X' (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046911"
# alt = "-"     #this is just a SNV, no insertion

# #simulate 5 base deletion - 'codon_change': 'cTGGAG/c' & 'aa_change': 'LE/X' (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046912"
# alt = "-"     #this is just a SNV, no insertion

# #simulate 6 base deletion - 'codon_change': 'cTGGAGAcg/ccg' & 'aa_change': 'LET/P' (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046913"
# alt = "-"     #this is just a SNV, no insertion

##MUTATION - plus strand genes
#simulate single point mutation - RESULT: codon =  aAt/aCt  & amino acids =  N/T (minus strand gene)
isoform_id = "ENST00000502297"
genomic_range = "6:17765077-17765077"
orig = None
alt = "G"     #this is just a SNV, no insertion

# #simulate single point mutation - codons = cTg/cGg & amino acid consequence =  L/R (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046908"
# alt = "G"     #this is just a SNV, no insertion

# #simulate 2 base mutation - 'codon_change': 'cTG/cGG' & amino acid consequence =  L/R (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046909"
# alt = "GG"     #this is just a SNV, no insertion

# #simulate 3 base mutation - 'codon_change': 'cTGGag/cGGGag' & 'aa_change': 'LE/RE' (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046910"
# alt = "GGG"     #this is just a SNV, no insertion

# #simulate 4 base mutation - 'codon_change': 'cTGGAg/cAAAAg' & 'aa_change': 'LE/QK' (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046911"
# alt = "AAAA"     #this is just a SNV, no insertion

# #simulate 6 base mutation - 'codon_change': 'taCGAGAAg/taAAAAAAg' & 'aa_change': 'YEK/*KK' (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046900-67046905"
# alt = "AAAAAA"     #this is just a SNV, no insertion

##INSERTIONS - plus strand gene
# #simulate single insertion - RESULT: codon =  aAt/aCCt  & amino acids =  N/TX
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765077-17765077"
# alt = "GG"      #one base for replace nucleotide at position, and +1 for single insertion

# #simulate double insertion - RESULT: codon =  aAt/aCCCt  & amino acids =  N/TX
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765077-17765077"
# alt = "GGG"     #one base for replace nucleotide at position, and +2 for double insertion

# #simulate triple insertion - RESULT: codon =  aAt/aCCCCt  & amino acids =  N/TP 
# isoform_id = "ENST00000502297"
# genomic_range = "6:17765077-17765077"
# alt = "GGGG"     #one base for replace nucleotide at position, and +3 for double insertion

# #simulate 2 insertion - 'codon_change': 'cTg/cGGg' & 'aa_change': 'L/RX'(plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046908"
# alt = "GG"

# #simulate 3 insertion - 'codon_change': 'cTg/cGGGg' & 'aa_change': 'L/RX'(plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046908"
# alt = "GGG"

# #simulate 4 insertion - 'codon_change': 'cTg/cGGGGg' & 'aa_change': 'L/RG'(plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046908"
# alt = "GGGG"

# #simulate 5 insertion - 'codon_change': 'cTg/cGGGGGg' & 'aa_change': 'L/RGX' (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046908"
# alt = "GGGGG"

# #simulate 7 insertion - 'codon_change': 'cTg/cGGGGGGGg' & 'aa_change': 'L/RGG' (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046908"
# alt = "GGGGGGG"

# #simulate 1 insertion - 'codon_change': 'ctg/cGtg' & 'aa_change': 'L/RX'(plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046907"
# alt = "G"

# #simulate 2 insertion - 'codon_change': 'ctg/cGGtg' & 'aa_change': 'L/RX'(plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046907"
# alt = "GG"

# #simulate 3 insertion - 'codon_change': 'ctg/cGGGtg' & 'aa_change': 'L/RV'(plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046907"
# alt = "GGG"

# #simulate 4 insertion - 'codon_change': 'ctg/cGGGGtg' & 'aa_change': 'L/RGX'(plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046907"
# alt = "GGGG"

##DELETIONS - plus strand gene
# #single deletion - 'codon_change': 'cTg/cg' & aa_change = L/X(plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046908"
# alt = "-"     

# #1-base deletion - 'codon_change': 'cTG/c' & aa_change = L/X(plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046909"
# alt = "-"     #this is just a SNV, no insertion

# #2-base deletion - 'codon_change': 'cTG/c' & aa_change = L/X(plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046909"
# alt = "-"     #this is just a SNV, no insertion

# #3-base deletion - 'codon_change': 'cTGGag/cag' & 'aa_change': 'LE/Q' (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046910"
# alt = "-"     #this is just a SNV, no insertion

# #4-base deletion - 'codon_change': 'cTGGAg/cg' & 'aa_change': 'LE/X (plus strand gene)
# isoform_id = "NM_001619.3"
# genomic_range = "11:67046908-67046911"
# alt = "-"     #this is just a SNV, no insertion

strand = 1
if "NM_" in isoform_id:
    opt_param = "variant_class=1&hgvs=1&refseq=1"
else:
    opt_param = "variant_class=1&hgvs=1"
build_hg38 = False
# specific_isoform = None
specific_isoform = isoform_id
obj_sna = SimpleNeoepitopeAll( genomic_range, orig, alt, DIR_GENOME, strand, opt_param, build_hg38, specific_isoform )


print "number of obj_sna.list_isoforms = ", len( obj_sna.list_isoforms )

for i, each_iso in enumerate(obj_sna.list_isoforms):        #each_iso = SimpleNeoepitopeIsoform instance
    print "mRNA: ----------", i, "----------"
    print "string = ", each_iso, " & len = ", len( each_iso.mRNA )

    mrna_orig = each_iso.mRNA
    mrna_alt = each_iso.mRNA_alt

    # [subseq_orig, subseq_alt] = each_iso.retrieve_mRNA_neoep_subseq( mrna_orig, mrna_alt, 1 )
    # # print "LEN subseq_org = ", len(subseq_orig)
    # # print "LEN subseq_alt = ", len(subseq_alt)
    # print "str subseq_org = ", ''.join(subseq_orig)
    # print "str subseq_alt = ", ''.join(subseq_alt)

    [subseq_orig, subseq_alt] = each_iso.retrieve_mRNA_neoep_subseq_v2( mrna_orig, mrna_alt, 5 )
    print "LEN subseq_org = ", len(subseq_orig)
    print "LEN subseq_alt = ", len(subseq_alt)
    str_subseq_orig =  ''.join(subseq_orig)
    str_subseq_alt =  ''.join(subseq_alt)
    # print "str subseq_org = ", str_subseq_orig
    # print "str subseq_alt = ", str_subseq_alt
    #look at the reverse complement - this is helpful for validating results for minus strand genes on UCSC Genome Browser
    print "str subseq_org = ", str_subseq_orig, " & Reverse Complement = ", str( Seq(str_subseq_orig.upper(), IUPAC.unambiguous_dna).reverse_complement() )
    print "str subseq_alt = ", str_subseq_alt, " & Reverse Complement = ", str( Seq(str_subseq_alt.upper(), IUPAC.unambiguous_dna).reverse_complement() )

    [aa_orig, aa_alt, aa_orig_2, aa_alt_2]  = each_iso.retrieve_comparative_neoeps( mrna_orig, mrna_alt, 5 )
    print "aa_org = ", aa_orig
    print "aa_alt = ", aa_alt
    print "aa_org_2 = ", aa_orig_2
    print "aa_alt_2 = ", aa_alt_2

    [seq_orig, seq_alt] = each_iso.extract_changed_codon_self()
    print "seq_orig = ", seq_orig
    print "seq_alt = ", seq_alt

    [codon_orig, codon_alt] = each_iso.compare_codon_orig_alt( True )
    print "codon_orig = ", codon_orig
    print "codon_alt = ", codon_alt

    # [aa_orig, aa_alt] = SimpleNeoepitopeIsoform.retrieve_comparative_neoeps( mrna_orig, mrna_alt, 9 )
    # print i2, " - aa_orig = ", aa_orig
    # print i2, " - aa_alt = ", aa_alt

    # for i2 in range( 0, len( each_iso.list_mRNA ) ):
    #     mrna_orig = each_iso.list_mRNA[i2]
    #     mrna_alt = each_iso.list_mRNA_alt[i2]

    #     #retrieve nucleotide
    #     # [subseq_orig, subseq_alt] = each_iso.retrieve_mRNA_neoep_subseq( mrna_orig, mrna_alt )
    #     # print "subseq_orig = ", subseq_orig
    #     # print "subseq_alt = ", subseq_alt

    #     #retrieve amino acid sequence        
    #     [aa_orig, aa_alt] = SimpleNeoepitopeIsoform.retrieve_comparative_neoeps( mrna_orig, mrna_alt, 9 )
    #     print i2, " - aa_orig = ", aa_orig
    #     print i2, " - aa_alt = ", aa_alt

    print "mRNA: ********", i, "********\n"


print "------------ TDD Completed: 171002_TEST_SimpleNeoepitope.py ------------"