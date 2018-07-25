#/usr/bin/python
import sys

from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv6 import Isoform, IsoformSJ, GenomicVariant, TranscribeTranscript, TranslateTranscript
from mokhaPy import mokhaPy

#Constants - directories
DIR_PROJ = "/home/mokha/Documents/Krauthammer_Lab"
DIR_CURR = DIR_PROJ + "/PythonClasses/SVSv6"
DIR_DATA = DIR_CURR + "/TestData"
DIR_RESULTS = DIR_CURR + "/TestResults"

DIR_GENOME = DIR_PROJ + '/ArchiveData/hg19.fa'      #directory for samtool-indexed genome

def find_containing_isoform( check_pos, genomic_range, all_isoforms ):
    """
    Finds the isoform that contains the position 'check_pos'
    """
    list_containing_indices = []        #array of the index of isoforms that contain isoforms that position 'check_pos'
    for i, each_iso in enumerate( all_isoforms ):
        isoform_id = each_iso.name
        obj_tt = create_tt_instance( isoform_id, genomic_range )

        if check_pos in obj_tt.arr_genome_pos:
            i_pos = obj_tt.arr_genome_pos.index( check_pos )
            print i, ": i_pos = ", i_pos, " in isoform ", isoform_id, " & strand = ", each_iso.strand, " & obj_tt.strand = ", obj_tt.iso_sj.strand, " contains position ", check_pos
            list_containing_indices.append( i )
        else:
            print "Position ", check_pos, " not found in isoform ", isoform_id
            # if not check_pos in obj_tt.arr_genome_pos:
                # print "----- the position ", check_pos, " is not in isoform ", isoform_id, " -----"


        # try:
        #     #Testing new function
        #     i_pos = obj_tt.arr_genome_pos.index( check_pos )
        #     print i, ": i_pos = ", i_pos, " in isoform ", isoform_id
        # except ValueError:
        #     print "Position not found in isoform ", isoform_id
        #     # if not check_pos in obj_tt.arr_genome_pos:
        #     #     print "----- the position ", check_pos, " is not in isoform ", isoform_id, " -----"

    return list_containing_indices



def create_tt_instance( isoform_id, genomic_range ):
    """
    Create instance of class TranslateTranscript
    """

    #retrieve position info
    split_genome_pos = Isoform.split_genome_pos( genomic_range )
    hash_pos = { 'chrom': split_genome_pos['chrom'], 'pos_oi': split_genome_pos['start'] }      #used to find the closest position for an isoform

    #create instance of IsoformSJ
    bool_simulant_sj = False
    group_sj = 0        #this means splicing events will NOT be grouped into 5' or 3' competitive splicing
    iso_sj = IsoformSJ( isoform_id, [], -10, hash_pos, bool_simulant_sj, group_sj )

    #create instance of TranslateTranscript
    canon_transcript = iso_sj.create_canon_transcript()
    #reconstruct canonical transcript
    obj_tt = TranslateTranscript( canon_transcript, iso_sj, DIR_GENOME, {} )

    return obj_tt


print "------------ TDD: 170728_List_Genomic_Alterations.py ------------"

g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )


#example 1 - nonsynonymous mutation - correct answer: AA (Original/mutated) = P/S; Codon (original/mutated) = Cct/Tct
# gene_sym = 'EPB41L3'
# genomic_range = 'chr18:5416374-5416374'
# strand = 1
# base_orig = 'G'
# base_alt = 'A'
# change_type = 0         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - nonsynonymous mutation - correct answer: AA (Original/mutated) =  A/D; Codon (original/mutated) = gCt/gAt
# gene_sym = 'CDON'
# genomic_range = 'chr11:125830916-125830916'
# strand = 1
# base_orig = 'G'
# base_alt = 'T'
# change_type = 0         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - nonsynonymous mutation - correct answer: AA (Original/mutated) =  E/Q; Codon (original/mutated) = Gag/Cag
# gene_sym = 'SLC6A13'
# genomic_range = 'chr12:369101-369101'
# strand = 1
# base_orig = 'C'
# base_alt = 'G'
# change_type = 0         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - nonsynonymous mutation - correct answer: AA (Original/mutated) =  N/K; Codon (original/mutated) = aaT/aaG
# gene_sym = 'RIMKLB'
# genomic_range = 'chr12:8902471-8902471'
# strand = 1
# base_orig = 'T'
# base_alt = 'G'
# change_type = 0         #the type of alteration - mutation = 0, insertion = 1, deletion = 2


#example - DNP (2 nonsynonymous mutation) - correct answer: AA (Original/mutated) =  S/F; Codon (original/mutated) = tCC/tTT
# gene_sym = 'CHRNA4'
# genomic_range = 'chr20:61981956-61981957'
# strand = 1
# base_orig = 'GG'
# base_alt = 'AA'
# change_type = 0         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - DNP (2 nonsynonymous mutation) - correct answer: AA (Original/mutated) =  QL; Codon (original/mutated) = caGTtg/caACtg
# gene_sym = 'CUL3'
# genomic_range = 'chr2:225339081-225339082'
# strand = 1
# base_orig = 'AC'
# base_alt = 'GT'
# change_type = 0         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - DNP (2 nonsynonymous mutation) - correct answer: AA (Original/mutated) =  YD/YY; Codon (original/mutated) = taCGac/taTTac
# gene_sym = 'TRPV5'
# genomic_range = 'chr7:142625257-142625258'
# strand = 1
# base_orig = 'CG'
# base_alt = 'AA'
# change_type = 0         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - DNP (2 nonsynonymous mutation) - correct answer: AA (Original/mutated) =  R/L; Codon (original/mutated) = cGC/cTT
# gene_sym = 'WASHC2C'            #reported with gene_sym 'FAM21C'
# genomic_range = 'chr10:46222947-46222948'
# strand = 1
# base_orig = 'GC'
# base_alt = 'TT'
# change_type = 0         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

# example 2 - insertion - correct answer: AA (Original/mutated) = G/AGGFGAGFGTGGFGG, codon (original/mutated) = ggt/gCCGGCGGCTTCGGAGCTGGTTTCGGCACTGGTGGCTTTGGTGgt
# gene_sym = 'KRT4'
# genomic_range = 'chr12:53207583-53207584'
# strand = 1
# base_orig = '-'
# base_alt = 'CACCAAAGCCACCAGTGCCGAAACCAGCTCCGAAGCCGCCGG'
# change_type = 1         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example 2 -  insertion - correct answer: AA (Original/mutated) = -/A, codon (original/mutated) = -/GCC
# gene_sym = 'RBM23'
# genomic_range = 'chr14:23371255-23371256'
# strand = 1
# base_orig = '-'
# base_alt = 'GGC'
# change_type = 1         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - insertion - correct answer: AA (Original/mutated) = A/DA, codon (original/mutated) = gcg/gATGcg
# gene_sym = 'GOLGA6L2'
# genomic_range = 'chr15:23685004-23685005'
# strand = 1
# base_orig = '-'
# base_alt = 'CAT'
# change_type = 1         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - insertion - correct answer: AA (Original/mutated) = V/VL, codon (original/mutated) = gtg/gTGCtg
# gene_sym = 'CNDP1'
# genomic_range = 'chr18:72223591-72223592'
# strand = 1
# base_orig = '-'
# base_alt = 'TGC'
# change_type = 1         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - insertion - correct answer: AA (Original/mutated) = V/DV, codon (original/mutated) = gtc/gATGtc
# gene_sym = 'HRC'
# genomic_range = 'chr19:49657710-49657711'
# strand = 1
# base_orig = '-'
# base_alt = 'CAT'
# change_type = 1         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - insertion - correct answer: AA (Original/mutated) = -/X, codon (original/mutated) = -/C
# gene_sym = 'LCA10'
# genomic_range = 'chrX:153149708-153149709'
# strand = 1
# base_orig = '-'
# base_alt = 'C'
# change_type = 1         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - insertion - correct answer: AA (Original/mutated) = N/KX, codon (original/mutated) = aac/aaAc
# gene_sym = 'RNF145'
# genomic_range = 'chr5:158630629-158630630'
# strand = 1
# base_orig = '-'
# base_alt = 'T'
# change_type = 1         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - insertion - correct answer: AA (Original/mutated) = -/DX, codon (original/mutated) = -/GATG
# gene_sym = 'LFNG'
# genomic_range = 'chr7:2552881-2552882'
# strand = 1
# base_orig = '-'
# base_alt = 'GATG'
# change_type = 1         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - insertion - correct answer: AA (Original/mutated) = C/WGX, codon (original/mutated) = tgc/tGGGGgc
# gene_sym = 'AGAP3'
# genomic_range = 'chr7:150783922-150783923'
# strand = 1
# base_orig = '-'
# base_alt = 'GGGG'
# change_type = 1         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - insertion - correct answer: AA (Original/mutated) = K/KLEAALX, codon (original/mutated) = aag/aAGCTTGAGGCAGCCCTag
# gene_sym = 'LMNA'           #THIS IS A PLUS STRAND GENE!!
# genomic_range = 'chr1:156100562-156100563'
# strand = 1
# base_orig = '-'
# base_alt = 'AGCTTGAGGCAGCCCT'
# change_type = 1         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - insertion - correct answer: AA (Original/mutated) = S/SRERERERX, codon (original/mutated) = agt/agTCGTGAGAGAGAGAGAGAAAGACt
# gene_sym = 'CCAR1'        #THIS IS A PLUS STRAND GENE!!
# genomic_range = 'chr10:70509022-70509023'
# strand = 1
# base_orig = '-'
# base_alt = 'TCGTGAGAGAGAGAGAGAAAGAC'
# change_type = 1         #the type of alteration - mutation = 0, insertion = 1, deletion = 2


#example - insertion - correct answer: AA (Original/mutated) = G/GGLTGDGRSRX, codon (original/mutated) = ggg/ggAGGGCTTACGGGTGATGGGAGAAGTCGGg
# gene_sym = 'ACLY'
# genomic_range = 'chr17:40055037-40055038'
# strand = 1
# base_orig = '-'
# base_alt = 'CCGACTTCTCCCATCACCCGTAAGCCCT'
# change_type = 1         #the type of alteration - mutation = 0, insertion = 1, deletion = 2


#example - insertion - correct answer: AA (Original/mutated) = -/QQQQQQX, codon (original/mutated) = -/CAGCAGCAGCAGCAGCAGCA
# gene_sym = 'ATXN3'
# genomic_range = 'chr14:92537354-92537355'
# strand = 1
# base_orig = '-'
# base_alt = 'TGCTGCTGCTGCTGCTGCTG'
# change_type = 1         #the type of alteration - mutation = 0, insertion = 1, deletion = 2


#example - frameshift deletion - correct answer: AA (Original/mutated) = GPG/G, codon (original/mutated) = ggTCCCGGc/ggc
# gene_sym = 'UBXN11'
# genomic_range = 'chr1:26608878-26608883'
# strand = 1
# base_orig = 'CCGGGA'
# base_alt = '-'
# change_type = 2         #the type of alteration - mutation = 0, insertion = 1, deletion = 2


#example - frameshift deletion - correct answer: AA (Original/mutated) = EPGCTKVP/-, codon (original/mutated) = GAGCCAGGCTGTACCAAGGTCCCT/-
# gene_sym = 'SPRR3'
# genomic_range = 'chr1:152975659-152975682'
# strand = 1
# base_orig = 'GAGCCAGGCTGTACCAAGGTCCCT'
# base_alt = '-'
# change_type = 2         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - frameshift deletion - correct answer: AA (Original/mutated) = KE/K, codon (original/mutated) = aAAGaa/aaa
# gene_sym = 'AXDND1'
# genomic_range = 'chr1:179504026-179504028'
# strand = 1
# base_orig = 'AAG'
# base_alt = '-'
# change_type = 2         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - frameshift deletion - correct answer: AA (Original/mutated) = PSVQ/Q, codon (original/mutated) = cCCAGCGTGCag/cag
# gene_sym = 'DUSP16'
# genomic_range = 'chr12:12630666-12630674'
# strand = 1
# base_orig = 'GCACGCTGG'
# base_alt = '-'
# change_type = 2         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - frameshift deletion - correct answer: AA (Original/mutated) = KKK/K, codon (original/mutated) = aaGAAAAAa/aaa
# gene_sym = 'USP36'
# genomic_range = 'chr17:76798549-76798554'
# strand = 1
# base_orig = 'TTTTTC'
# base_alt = '-'
# change_type = 2         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - frameshift deletion - correct answer: AA (Original/mutated) = EAAPRGR/-, codon (original/mutated) = GAGGCGGCGCCCCGGGGGAGA/-
# gene_sym = 'SRRD'
# genomic_range = 'chr22:26879947-26879967'
# strand = 1
# base_orig = 'GAGGCGGCGCCCCGGGGGAGA'
# base_alt = '-'
# change_type = 2         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - frameshift deletion - correct answer: AA (Original/mutated) = PGPS/X, codon (original/mutated) = CCCGGCCCCAgt/gt
# gene_sym = 'UBXN11'
# genomic_range = 'chr1:26608849-26608858'
# strand = 1
# base_orig = 'TGGGGCCGGG'
# base_alt = '-'
# change_type = 2         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - frameshift deletion - correct answer: AA (Original/mutated) = GPG/X, codon (original/mutated) = GGTCCCGGt/t
# gene_sym = 'UBXN11'
# genomic_range = 'chr1:26608860-26608867'
# strand = 1
# base_orig = 'CCGGGACC'
# base_alt = '-'
# change_type = 2         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - frameshift deletion - correct answer: AA (Original/mutated) = E/X, codon (original/mutated) = gaA/ga
# gene_sym = 'PDE4DIP'
# genomic_range = 'chr1:144923729-144923729'
# strand = 1
# base_orig = 'T'
# base_alt = '-'
# change_type = 2         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - frameshift deletion - correct answer: AA (Original/mutated) = GS/GX, codon (original/mutated) = ggGTca/ggca
# gene_sym = 'SARM1'
# genomic_range = 'chr17:26708303-26708304'
# strand = 1
# base_orig = 'GT'
# base_alt = '-'
# change_type = 2         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - frameshift deletion - correct answer: AA (Original/mutated) = QC/X, codon (original/mutated) = caGTGT/ca
# gene_sym = 'DTNA'       #THIS IS A PLUS STRAND GENE
# genomic_range = 'chr18:32398198-32398201'
# strand = 1
# base_orig = 'GTGT'
# base_alt = '-'
# change_type = 2         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - frameshift deletion - correct answer: AA (Original/mutated) = C/X, codon (original/mutated) = Tgt/gt --> PROBLEM: THIS POSITION IS NOT FOUND IN THE GENE
# gene_sym = 'SHANK3'
# genomic_range = 'chr22:51135969-51135969'
# strand = 1
# base_orig = 'T'
# base_alt = '-'
# change_type = 2         #the type of alteration - mutation = 0, insertion = 1, deletion = 2

#example - frameshift deletion - correct answer: AA (Original/mutated) = K/X, codon (original/mutated) = Aaa/aa
gene_sym = 'OR6C76'
genomic_range = 'chr12:55820959-55820959'
strand = 1
base_orig = 'A'
base_alt = '-'
change_type = 2         #the type of alteration - mutation = 0, insertion = 1, deletion = 2


all_isoforms = Isoform.obj_cruzdb.refGene.filter_by( name2 = gene_sym ).all()

#check which isoform contains this position  - ONLY USE THIS WHEN I DON'T KNOW WHICH ISOFORM TO SELECT
split_genome_pos = Isoform.split_genome_pos( genomic_range )
isoform_indices = find_containing_isoform( split_genome_pos['start'], genomic_range, all_isoforms )


##TEST::
test_get_gene = Isoform.obj_cruzdb.bin_query( 'refGene', split_genome_pos['chrom'], split_genome_pos['start'], split_genome_pos['end'] )

print "See which genes are here"
for i_test, each_iso in enumerate( test_get_gene ):
    print i_test, ": ", each_iso, " ||||| gene_sym = ", each_iso.name2, " & iso = ", each_iso.name


for each_iso_i in isoform_indices:
    #retrieve isoform
    isoform_id = all_isoforms[each_iso_i].name
    #create instance TranslateTranscript
    obj_tt = create_tt_instance( isoform_id, genomic_range )

    print each_iso_i, " - ", all_isoforms[each_iso_i].name, ": obj_tt.iso_sj.strand = ", obj_tt.iso_sj.strand, " & strand = ", strand

    hash_alt_conseq = obj_tt.alteration_consequence( base_alt, genomic_range, strand, change_type )

    print "############################\n\n"



# print "\n\n>>>>>>>>>>>>> Fix Error with positions"


# hash_pos = Isoform.split_genome_pos( genomic_range )

# hash_pos_nuc = obj_tt.retrieve_genomic_subseq( genomic_range, True, False )
# print "hash_pos_nuc = ", hash_pos_nuc
# pos_start = obj_tt.find_codon_beginning_nearest( hash_pos['start'], True )
# pos_end = obj_tt.find_codon_rf_nearest_TEST( hash_pos['end'], 2, False )

# print "pos_start = ", pos_start
# print "pos_end = ", pos_end



print "------------ TDD Completed: 170728_List_Genomic_Alterations.py ------------"