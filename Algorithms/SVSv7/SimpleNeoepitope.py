#/usr/bin/python
import subprocess
import requests
import re
import copy

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from cruzdb import Genome

#need Isoform class for obj_cruzdb (object cruzDB)
from Isoform import Isoform
#import my own class functions
from EnsemblVEP import EnsemblVEP
from MHC_IEDB import MHC_IEDB

"""
Things I need to to
-import EnsemblVEP class
-use the HGVSc to make mutation
-Q: Can I just retrieve the amino acid change caused by the mutation?

PROTOCOL:
    -retrieve the nucleotide alteration & genomic range
    -retrieve information from VEP, specifically HGVSp & HGVSc
    -retrieve information for isoform using cruzdb
    -using samtool's faidx, retrieve the nucleotide sequence
    -translate the nucleotide sequence to neoepitope, & 
    -pretty sure I need the sliding window algorithm for extracting the neoepitopes
    -Q: How do I retrieve the start & end translation positions for the gene
        -IDEA 1: just retrieve the coding sequence (CDS)
        -IDEA 2: compare the exon to CDS, align the starting of the CDS to all the exons & 
        -IDEA 3: need another array that records genomic position that corresponds to mRNA -> relative position 
        -IDEA 4: Just append the end of the nucleotide 
        -IDEA 5: make up backup system to evaluate possibility for NSD - Need exon, genomic position
    -Q: should I include raw amino acid change (e.g. V/A where original_aa/mut_aa) for double checking?
    -Q: Should I make an Isoform instance with SimpleNeoepitopeIsoform()? Perhaps just make it separate...
"""
class SimpleNeoepitopeIsoform():
    def __init__( self, hash_isoform_alt, path_genomeidx ):
        """
        -Args:
            -hash_isoform_alt = hash retrieved from SimpleNeoepitope's function retrieve_isoform_alt(), where the format for hash_isoform_alt can be found
        """
        self.gene_symbol = hash_isoform_alt["gene_symbol"]
        self.isoform_id = hash_isoform_alt["isoform_id"].split(".")[0] if "." in hash_isoform_alt["isoform_id"] else hash_isoform_alt["isoform_id"]
        self.gene_strand = hash_isoform_alt["gene_strand"]
        self.allele_string = hash_isoform_alt["allele_string"]
        self.variant_class = hash_isoform_alt["variant_class"]
        self.nuc_orig = hash_isoform_alt["nuc_orig"]
        self.nuc_alt = hash_isoform_alt["nuc_alt"]

        #CONJ: for cds_start & cds_end, I think I need to subtract by 1 since VEP returns relative position of mRNA starting at 1 while the array index starts at 0 --> do the -1 in the methods below
        self.chrom = hash_isoform_alt["chrom"]
        self.genome_start = hash_isoform_alt["genome_start"]
        self.genome_end = hash_isoform_alt["genome_end"]
        self.cds_start = hash_isoform_alt["cds_start"] if not str(hash_isoform_alt["cds_start"]).isdigit() else int( hash_isoform_alt["cds_start"] )
        self.cds_end = hash_isoform_alt["cds_end"] if not str(hash_isoform_alt["cds_end"]).isdigit() else int( hash_isoform_alt["cds_end"] )
        self.codon_change = hash_isoform_alt["codon_change"]

        self.aa_change = hash_isoform_alt["aa_change"]      #this is the amino acid change in the format "ORG/ALT" where "ORG" is the original amino acid & "ALT" is the altered amino acid
        self.aa_start = hash_isoform_alt["aa_start"] if not str(hash_isoform_alt["aa_start"]).isdigit() else int( hash_isoform_alt["aa_start"] )
        self.aa_end = hash_isoform_alt["aa_end"] if not str(hash_isoform_alt["aa_end"]).isdigit() else int( hash_isoform_alt["aa_end"] )

        self.consequence = hash_isoform_alt["consequence"]      #this is the consequence of the alteration (e.g missense_variant, synonymous_variant)
        self.hgvsc = hash_isoform_alt['hgvsc']      #format: NM_002880.3:c.935T>C
        self.hgvsp = hash_isoform_alt['hgvsp']      #format: NP_002871.1:p.Val312Ala
        self.path_genomeidx = path_genomeidx        #FIX THIS: it is inefficient to save string this for each isoform
        
        ##MAY DELETE LATER BECAUSE OF METHOD 2: check if amino acid change has occurred!!
        # #METHOD 1: determine if the genomic alteration changes the amino acid
        # try:
        #     ##TEST::
        #     print "\t^^^^^^SN_ISO TEST:: list_aa_change = ", list_aa_change, " & list_codon = ", list_codon
        #     #check amino acid
        #     list_aa_change = self.aa_change.split('/')
        #     #check codon
        #     list_codon = self.codon_change.split('/')
        #     check_aa_orig = str( Seq( list_codon[0] ).translate( to_stop = False ) )
        #     check_aa_alt = str( Seq( list_codon[1] ).translate( to_stop = False ) )


        #     ##TEST::
        #     print "\t^^^^^^ SN_ISO TEST:: list_aa_change = ", list_aa_change, " & list_codon = ", list_codon

        #     #determine if amino acid alteration occurs
        #     self.bool_alteration = False
        #     if list_aa_change[0] != list_aa_change[1]:
        #         self.bool_alteration = True
        #     elif check_aa_orig != check_aa_alt:     #this may be a bit excessive, not sure
        #         self.bool_alteration = True
        # except:
        #     self.bool_alteration = False

        #METHOD 2: check if amino acid change has occurred (if no amino acid change has occurred with genomic alteration, then what is the point of making an "altered mRNA")
        if not '/' in self.aa_change:
            self.bool_alteration = False
        else:
            list_aa_change = self.aa_change.split('/')
            if list_aa_change[0] != list_aa_change[1]:
                self.bool_alteration = True
            else:
                self.bool_alteration = False



        ##DELETE - JUST FOR REFERENCE
        # hash_isoform_alt[isoform_id] = {
        #     "gene_symbol": gene_symbol, 
        #     "gene_id": gene_id,
        #     "isoform_id": isoform_id,
        #     "gene_strand": hash_transcript['strand'],
        #     "nuc_orig": list_allele[0],
        #     "nuc_alt": list_allele[1],
        #     "codon_change": codon_change,
        #     "cds_start": cds_start,
        #     "cds_end": cds_end,
        #     "aa_change": aa_change,
        #     "aa_start": aa_start,
        #     "aa_end": aa_end,
        #     "hgvsc": hash_transcript['hgvsc'],
        #     "hgvsp": hash_transcript['hgvsp']
        #     }

        #create the mRNA sequence - original & mutation
        # self.list_mRNA = self.create_mRNA()
        # #NEED TO MAKE MUTATION TO mRNA using HGVSc!! - need to see how insertion & deletions work as well!!
        # self.list_mRNA_alt = self.create_mRNA_alt()

        self.mRNA = self.create_mRNA()
        bool_insert = False if self.genome_start <= self.genome_end else True      #determine whether only insert (& no deletion) should occur or not

        #METHOD 1: retrieve the altered mRNA based on mutation
        # self.mRNA_alt = [] if not self.mRNA else SimpleNeoepitopeIsoform.create_alt_mRNA( self.mRNA, self.nuc_alt, self.cds_start, self.cds_end, bool_insert )
        
        #METHOD 2: retrieve the altered mRNA based on mutation
        if not self.mRNA or not self.bool_alteration:
            self.mRNA_alt = []
        else:
            self.mRNA_alt = SimpleNeoepitopeIsoform.create_alt_mRNA( self.mRNA, self.nuc_alt, self.cds_start, self.cds_end, bool_insert )
        

    def __str__( self ):
        str_info = self.gene_symbol + " - " + self.isoform_id
        str_info+= " | cds_start = " + str( self.cds_start ) + " & cds_end = " + str( self.cds_end ) + " | nuc_orig = " + str( self.nuc_orig ) + " & nuc_alt = " + str( self.nuc_alt )
        str_info+= " | codon_change = " + str( self.codon_change )
        str_info+= " | aa_change = " + str( self.aa_change ) + " & aa_start = " + str( self.aa_start ) + " & aa_end = " + str( self.aa_end )
        str_info+= " | hgvsc = " + str( self.hgvsc ) + " & hgvsp = " + str( self.hgvsp )
        str_info+= " | variant_class = " + str( self.variant_class )
        str_info+= " | bool alteration = " + str( self.bool_alteration )
        return str_info

    #INCOMPLETE - MAY DELETE THIS
    # @staticmethod
    # def split_HGVSc( str_HGVSc ):
    #     list_1 = str_HGVSc.split(":")
    #     list_2 = list_1[1].split(".")

    def create_mRNA( self ):
        """
        NOTE: this only records one mRNA sequence for the isoform. I think I've seen cruzdb retrieve multiple sequences for a specific isoform, but I will just record the first isoform
        """
        if "ENST" in self.isoform_id:
            info_isoform = Isoform.obj_cruzdb.ensGene.filter_by( name = self.isoform_id ).all()
        else:
            info_isoform = Isoform.obj_cruzdb.refGene.filter_by( name = self.isoform_id ).all()

        convert_ss = { '+' : 1, '-' : -1 }
        list_mRNA_seq = []
        for each_iso in info_isoform:
            if each_iso.name != self.isoform_id:
                continue

            ##TEST::print ">>>>>>SN_ISOFORM.create_mRNA 3: each_iso.name = ", each_iso.name, " & self.isoform_id = ", self.isoform_id

            iso_strand = convert_ss[each_iso.strand]
            arr_pos = each_iso.cds      #other possible array values - each_iso.exons, each_iso.introns
            str_mRNA_seq = SimpleNeoepitopeIsoform.retrieve_seq( each_iso.chrom, arr_pos, self.path_genomeidx, iso_strand, True )

            ##TEST::
            # print ">>>>>>SN_ISOFORM.create_mRNA 3: str_mRNA_seq = ", arr_pos
            # print ">>>>>>SN_ISOFORM.create_mRNA 4: str_mRNA_seq = ", str_mRNA_seq

            #split string into an array so I can make mutations to it later -> strings are immutable whereas arrays are not
            list_mRNA_seq = list( str_mRNA_seq )

            ##TEST::print ">>>>>>SN_ISOFORM.create_mRNA 5: list_mRNA_seq = ", list_mRNA_seq

            break

        return list_mRNA_seq

    def extract_changed_codon_self( self ):
        """
        Applies the function extract_changed_codon() for this specific instance of the alteration
        """
        #determine if an insertion is the type of alteration
        bool_insert = False if self.genome_start <= self.genome_end else True      #determine whether only insert (& no deletion) should occur or not
        return SimpleNeoepitopeIsoform.extract_changed_codon( self.mRNA, self.nuc_alt, self.cds_start, self.cds_end, bool_insert )

    def extract_changed_aa_self( self ):
        """
        Applies the function extract_changed_codon() for this specific instance of the alteration
        """
        #determine if an insertion is the type of alteration
        bool_insert = False if self.genome_start <= self.genome_end else True      #determine whether only insert (& no deletion) should occur or not
        return SimpleNeoepitopeIsoform.extract_changed_aa( self.mRNA, self.nuc_alt, self.cds_start, self.cds_end, bool_insert )

    @staticmethod
    def extract_changed_aa( mrna_orig, nuc_alt, cds_start, cds_end, bool_insert = False ):
        """
        determines the amino acid sequence based on the 
        """
        [seq_orig, seq_alt] = SimpleNeoepitopeIsoform.extract_changed_codon( mrna_orig, nuc_alt, cds_start, cds_end, bool_insert )
        if not seq_orig or not seq_alt:
            return [None, None]

        aa_orig = str( Seq( seq_orig.upper() ).translate( to_stop = False ) )
        aa_alt = str( Seq( seq_alt.upper() ).translate( to_stop = False ) )

        return [aa_orig, aa_alt]

    @staticmethod
    def extract_changed_codon( mrna_orig, nuc_alt, cds_start, cds_end, bool_insert = False ):
        """
        Retrieves the codon that is affected by the alteration, before (original) & after (altered) the alteration --> this displays the original & altered codon in more of a Ensembl VEP style! (i.e. lowercase letters are unchanged nucleotides, whereas the capital letters are changed)
        Args:
            -mrna_orig = array of nucleotide sequences. It is assumed these have been corrected for gene strand.
            -nuc_alt = string that is the new nucleotide alteration that will be inputed in cds_start.
                -if it is a deletion, then this string will be '-'
            -cds_start & cds_end = integers that are the relative mRNA position, where cds_start <= cds_end uness there is an insertion that does not remove a base.
                -NOTE: the cds_start is of base-1 whereas an array is base-0
            -bool_insert = boolean that handles if there is an insertion that doesn't replace any bases, this happens when end_pos = start_pos - 1 (e.g. "11:67046908-67046907")
                -True = only true when end_pos is 1 less than start pos (e.g. "11:67046908-67046907"), will 
        PROTOCOL:
            -get the mRNA seq (there may be multiple)
            -mutation: str_mRNA[pos] = new_base
            -insertion: somehow need to add a gap position so I can add new nucleotides
                -if an insertion is replace nucleotides, then need to "remove bases" & then add nucleotides
            -deletion: use start & end replace nucleotides with '-', and when need to translate just remove these '-'
            -MAKE SURE TO PRESERVE mRNA string length just in case there are multiple alterations --> this may be hard with insertions UNLESS I make it into an array
        """
        if not str(cds_start).isdigit() or not str(cds_end).isdigit():
            return [None, None]
        else:
            cds_start = int( cds_start )
            cds_end = int( cds_end )

        #copy the entire array that is the mRNA sequence (each element is a nucleotide sequence)
        mRNA_orig = copy.copy( mrna_orig )
        mRNA_orig = list( ''.join( mRNA_orig ).lower() )
        mRNA_alt = copy.copy( mRNA_orig )
        # mRNA_alt = self.mRNA[:]

        #retrieve the position of a complete reading window
        [rf_start, rf_end] = SimpleNeoepitopeAll.calculate_full_rf_pos( cds_start, cds_end, len(mrna_orig), True )
        if not rf_start or not rf_end:
            return [None, None]

        if bool_insert:
            mRNA_alt[ cds_start - 1 ] += nuc_alt.upper()
        else:
            for i in range( cds_start - 1, cds_end ):
                ##TEST:: print "\tSN_Iso.extract_changed_codon - BEFORE change nucleotide = ", mRNA_alt[i-1:i+2], " | pos = ", i

                mRNA_orig[i] = mRNA_orig[i].upper()
                mRNA_alt[i] = "-"

                ##TEST:: print "\tSN_Iso.extract_changed_codon - AFTER change nucleotide = ", mRNA_alt[i-1:i+2], " | pos = ", i

            mRNA_alt[ cds_start - 1 ] = nuc_alt.upper()

        ##TEST:: print "\tSN_Iso.extract_changed_codon - ADD = ", mRNA_alt[cds_start-1:cds_end+2], " & bool_insert = ", bool_insert

        #extract the nucleotide sequence for the full window 
        seq_orig = ''.join( [mRNA_orig[x] for x in range( rf_start, rf_end + 1 ) ] ).replace('-', '')
        seq_alt = ''.join( [mRNA_alt[x] for x in range( rf_start, rf_end + 1 ) ] ).replace('-', '')


        ##TEST::
        # print "\tSN_Iso.extract_changed_codon - FINAL VERSION = ", seq_orig
        # print "\tSN_Iso.extract_changed_codon - FINAL VERSION = ", seq_alt

        return [seq_orig, seq_alt]

    @staticmethod
    def create_alt_mRNA( mrna_orig, nuc_alt, cds_start, cds_end, bool_insert = False ):
        """
        Create the alteration in the mRNA sequence
        NEED TO FIX THIS FOR INSERTIONS, DELETIONS - I think I did fix this for insertions/deletions
        Args:
            -mrna_orig = array of nucleotide sequences. It is assumed these have been corrected for gene strand.
            -nuc_alt = string that is the new nucleotide alteration that will be inputed in cds_start.
                -if it is a deletion, then this string will be '-'
            -cds_start & cds_end = integers that are the relative mRNA position, where cds_start <= cds_end uness there is an insertion that does not remove a base.
                -NOTE: the cds_start is of base-1 whereas an array is base-0
            -bool_insert = boolean that handles if there is an insertion that doesn't replace any bases, this happens when end_pos = start_pos - 1 (e.g. "11:67046908-67046907")
                -True = only true when end_pos is 1 less than start pos (e.g. "11:67046908-67046907"), will 
        PROTOCOL:
            -get the mRNA seq (there may be multiple)
            -mutation: str_mRNA[pos] = new_base
            -insertion: somehow need to add a gap position so I can add new nucleotides
                -if an insertion is replace nucleotides, then need to "remove bases" & then add nucleotides
            -deletion: use start & end replace nucleotides with '-', and when need to translate just remove these '-'
            -MAKE SURE TO PRESERVE mRNA string length just in case there are multiple alterations --> this may be hard with insertions UNLESS I make it into an array
        """
        if not str(cds_start).isdigit() or not str(cds_end).isdigit():
            return []
        else:
            cds_start = int( cds_start )
            cds_end = int( cds_end )

        mRNA_alt = copy.copy( mrna_orig )
        # mRNA_alt = self.mRNA[:]

        if bool_insert:
            mRNA_alt[ cds_start - 1 ] += nuc_alt
        else:
            for i in range( cds_start - 1, cds_end ):
                ##TEST:: print "\tSN_Iso.create_mRNA_alt - BEFORE change nucleotide = ", mRNA_alt[i-1:i+2], " | pos = ", i
                mRNA_alt[i] = "-"
                ##TEST:: print "\tSN_Iso.create_mRNA_alt - AFTER change nucleotide = ", mRNA_alt[i-1:i+2], " | pos = ", i

            mRNA_alt[ cds_start - 1 ] = nuc_alt

        ##TEST:: print "\tSN_Iso.create_mRNA_alt - ADD = ", mRNA_alt[cds_start-1:cds_end+2], " & bool_insert = ", bool_insert

        #this will remove all '-' and return an array 
        # return list( ''.join( [x for x in mRNA_alt if x != '-'] ) )
        mRNA_alt = list( ''.join( [x for x in mRNA_alt if x != '-'] ) )

        ##TEST::
        print "\tSN_Iso.create_mRNA_ORG - FINAL VERSION = ", mrna_orig[cds_start-4:cds_end+5]
        print "\tSN_Iso.create_mRNA_alt - FINAL VERSION = ", mRNA_alt[cds_start-4:cds_end+5]

        return mRNA_alt

    def compare_codon_orig_alt( self, return_str = False ):
        """
        Compares the codon before & after alteration
        Args:
            -return_str = boolean where - True = will return codons in string form; False = returns codon in array form (this is the default)
        -NOTE: how numbers change between CDS, array index, & reading frame
            -cds_position = 1 2 3 4 5 6 7 8 9
            -array_index  = 0 1 2 3 4 5 6 7 8
            -read_frame   = 0 1 2 0 1 2 0 1 2
        """
        ##TEST:: print "\t~~~~~SNIso.compare_codon_orig_alt 1 -> self.cds_start = ", self.cds_start, " | self.cds_end = ", self.cds_end

        if not str(self.cds_start).isdigit() or not str(self.cds_end).isdigit():
            return [None, None]

        #retrieve the start & end 
        [rf_start, rf_end] = SimpleNeoepitopeAll.calculate_full_rf_pos( self.cds_start, self.cds_end, len( self.mRNA ), True )
        if not rf_start or not rf_end:
            return [None, None]
        #get the original codon sequence
        codon_orig = self.mRNA[rf_start:rf_end + 1] if not return_str else ''.join( self.mRNA[rf_start:rf_end + 1] )
        #get the mutated codon sequence
        if rf_start > len( self.mRNA_alt ) or rf_end > len( self.mRNA_alt ):
            codon_alt = None
        else:
            codon_alt = self.mRNA_alt[rf_start:rf_end + 1] if not return_str else ''.join( self.mRNA_alt[rf_start:rf_end + 1] )

        ##TEST:: print "\t~~~~SNIso.compare_codon_orig_alt 8 -> self.cds_start = ", self.cds_start, " | self.cds_end = ", self.cds_end, " | rf_start = ", rf_start, " | rf_end = ", rf_end, " | len( self.mRNA ) = ", len( self.mRNA )

        #MAY DELETE THIS
        # #double check to see if range is not out of bounds
        # if rf_start < 0 or rf_end < 0:
        #     return [None, None]
        # #get the original codon
        # if rf_start > len( self.mRNA ) or rf_end > len( self.mRNA ):
        #     codon_orig = None
        # else:
        #     codon_orig = self.mRNA[rf_start:rf_end + 1] if not return_str else ''.join( self.mRNA[rf_start:rf_end + 1] )
        # #get the mutated codon
        # if rf_start > len( self.mRNA_alt ) or rf_end > len( self.mRNA_alt ):
        #     codon_alt = None
        # else:
        #     codon_alt = self.mRNA_alt[rf_start:rf_end + 1] if not return_str else ''.join( self.mRNA_alt[rf_start:rf_end + 1] )

        return [codon_orig, codon_alt]

    def compare_aa_orig_alt( self ):
        """
        Compares the amino acid sequence before & after alteration
        -NOTE: how numbers change between CDS, array index, & reading frame
            -cds_position = 1 2 3 4 5 6 7 8 9
            -array_index  = 0 1 2 3 4 5 6 7 8
            -read_frame   = 0 1 2 0 1 2 0 1 2
        """
        [codon_orig, codon_alt] = self.compare_codon_orig_alt( True )

        if codon_orig:
            aa_orig = str( Seq( codon_orig ).translate(to_stop = False) )
        else:
            aa_orig = None

        if codon_alt:
            aa_alt = str( Seq( codon_alt ).translate(to_stop = False) )
        else:
            aa_alt = None

        return [aa_orig, aa_alt]



    #MAY DELETE because of def retrieve_comparative_neoeps()
    # def retrieve_comparative_neoeps_BACKUP( self, mrna_orig, mrna_alt, aa_len ):
    #     """
    #     calculate neoepitopes
    #     PROTOCOL:
    #         -split into codons
    #         -combine the codons & see if the nucleotide sequences match between original & mutation
    #         -
    #     """

    #     codons_orig = [mrna_orig[i:i+3] for i in range(0, len( )) if i % 3 == 0]
    #     codons_mut = [mrna_alt[i:i+3] for i in range(0, len( )) if i % 3 == 0]

    #     for i in range(0, len( codons_mut )):
    #         seq_orig = ''.join( codons_orig[i:i + aa_len] )
    #         seq_mut = ''.join( codons_mut[i:i + aa_len] )
    #         #if sequence does not match, then record both AAs
    #         if seq_orig != seq_mut:
    #             pass
    #     pass


    @staticmethod
    def retrieve_mRNA_neoep_subseq( mrna_orig, mrna_alt, aa_len ):
        """
        compares the nucleotide sequence between "mrna_orig" & "mrna_alt" and returns the subsequence that is different between them
        Args:
            -mrna_orig & mrna_alt = array of nucleotide sequences. It is assumed these have been corrected for gene strand.
            -aa_len = integer that is the number of codons (aa) to retrieve before nucleotide difference between mrna_orig & mrna_alt
        Returns: an array of 2 elements where each is the subsequence "aa_len" positions before the genomic alteration, where [0] = array that is original mRNA sequence (each element is a nucleotide base), [1] = same as [0] but for altered mRNA sequence
        """
        #find the first position where 
        i_diff = -1
        for i in range( 0, len(mrna_alt) ):
            if mrna_orig[i] != mrna_alt[i]:
                i_diff = i
                break

        #TEST:: print "SimpleNeoepitope.retrieve_mRNA_neoep_subseq 1: i_diff = ", i_diff

        if i_diff < 0:
            return [[], []]

        start_i = 0 if i_diff < (aa_len * 3) else i_diff - (aa_len * 3)
        adj_start_i = start_i + 1       #need to add +1 because array starts at index 0, not 1. For example, at index 2 of array it is actually nucleotide 3 (where rf = 2, therefore need to +1 to make sure algorithm considers this as the end of a codon

        ##TEST:: print "SimpleNeoepitope.retrieve_mRNA_neoep_subseq 2a: BEFORE - start_i = ", start_i, " & start_i % 3 = ", (start_i % 3)

        #need to find complete open reading frame, therefore need to start at rf = 0
        if adj_start_i % 3 != 0:
            adj_start_i = adj_start_i - (adj_start_i % 3)
        start_i = adj_start_i
        
        ##TEST:: print "SimpleNeoepitope.retrieve_mRNA_neoep_subseq 2b: AFTER - start_i = ", start_i, " & start_i % 3 = ", (start_i % 3)

        #TEST::  print "SimpleNeoepitope.retrieve_mRNA_neoep_subseq: start_i = ", start_i, " | mrna_orig = ", mrna_orig[start_i], " | mrna_alt = ", mrna_alt[start_i], " |||| i_diff = ", i_diff, " --> mrna_orig = ", mrna_orig[i_diff], " | mrna_alt = ", mrna_alt[i_diff]

        #retrieve the subsequence that is coding a different
        # end_i = i + (aa_len * 3) if alt_type == 0 else -1
        end_i = -1
        subseq_orig = mrna_orig[ start_i:end_i ]
        subseq_alt = mrna_alt[ start_i:end_i ]

        return [subseq_orig, subseq_alt]


    ##QUES: why am I working on "retrieve_mRNA_neoep_subseq_v3()" again? - CONJ: look at playAtlas.py post "##17.10.13 - be able to retrieve amino acid sequence". It's suppose to be a way to find differences between mrna_orig & mrna_alt by using the AA sequence instead
    @staticmethod
    def retrieve_mRNA_neoep_subseq_v3( mrna_orig, mrna_alt, aa_len ):
        """
        Returns the subsequence of the original & altered nucleotide sequence that is altered
        Args:
            -mrna_orig & mrna_alt = array of nucleotide sequences, where each element is a nucleotide base. Assumed that this is the sequence based on the gene strand
            -aa_len = integer that is the number of codons (aa) to retrieve before nucleotide difference between mrna_orig & mrna_alt
        PROTOCOL:
            -get mRNA, original & altered
            -translate sequence, and calculate the minimum # of AAs needed to find difference & same sequences
            -
        """
        aa_orig = str( Seq( mrna_orig ).translate( to_stop = False ) )
        aa_alt = str( Seq( mrna_alt ).translate( to_stop = False ) )
        len_aa = max( [len(aa_orig), len(aa_alt)] )

        #find start of dissimilarity between both sequences
        pass



    @staticmethod
    def retrieve_mRNA_neoep_subseq_v2( mrna_orig, mrna_alt, aa_len ):
        """
        compares the nucleotide sequence between "mrna_orig" & "mrna_alt" and returns the subsequence that is different between them
        Args:
            -mrna_orig & mrna_alt = array of nucleotide sequences, where each element is a nucleotide base
            -aa_len = integer that is the number of codons (aa) to retrieve before nucleotide difference between mrna_orig & mrna_alt
        Returns: an array of 2 elements where each is the subsequence "aa_len" positions before the genomic alteration, where [0] = array that is original mRNA sequence (each element is a nucleotide base), [1] = same as [0] but for altered mRNA sequence
        """
        #find the first position where 
        i_diff_start = -1
        num_bases = aa_len * 3      #converts the number of amino acids (number of codons) to number of nucleotide bases
        for i in range( 0, len(mrna_alt) ):
            # sub_start_i = i if i < num_bases else i - num_bases
            sub_end_i = len(mrna_alt) if len(mrna_alt) < (i + num_bases) else i + num_bases

            # if mrna_orig[sub_start_i:sub_end_i] != mrna_alt[sub_start_i:sub_end_i]:
            if mrna_orig[i:sub_end_i] != mrna_alt[i:sub_end_i]:
                i_diff_start = i
                break

        #if no difference found, then return None
        if i_diff_start < 0:
            return [[], []]

        #Find where mrna_orig & mrna_alt stop being different from each other
        i_diff_end = -1
        for i in range( i_diff_start, len(mrna_alt) ):
            # sub_start_i = i if i < num_bases else i - num_bases
            sub_end_i = len(mrna_alt) if len(mrna_alt) < (i + num_bases) else i + num_bases

            # if mrna_orig[sub_start_i:sub_end_i] != mrna_alt[sub_start_i:sub_end_i]:
            if mrna_orig[i:sub_end_i] == mrna_alt[i:sub_end_i]:
                i_diff_end = sub_end_i
                break

        ##QUES: Should I do this?? - ACTUALLY, NO, I DON'T THINK I SHOULD DO THIS
        # [rf_start, rf_end] = SimpleNeoepitopeIsoform.calculate_full_rf_pos( cds_start, cds_end, len_seq, True )

        #calculate the start position
        # start_i = 0 if i_diff_start < num_bases else i_diff_start - num_bases
        adj_start_i = i_diff_start + 1       #need to add +1 because array starts at index 0, not 1. For example, at index 2 of array it is actually nucleotide 3 (where rf = 2, therefore need to +1 to make sure algorithm considers this as the end of a codon
        #retrieve correct reading frame for start_i
        if adj_start_i % 3 != 0:
            adj_start_i = adj_start_i - (adj_start_i % 3)
        start_i = adj_start_i

        #calculate the end position - NOTE: I did not calculate this to end exactly at reading frame = 2
        # end_i = -1 if len(mrna_alt) < i_diff_end + num_bases else i_diff_end + num_bases
        end_i = -1 if len(mrna_alt) < i_diff_end else i_diff_end
        

        ##TEST:: print "\t\tSN_Iso.retrieve_mRNA_neoep_subseq_v2(): start_i = ", start_i, " & end_i = ", end_i, " || i_diff_start = ", i_diff_start, " & i_diff_end = ", i_diff_end, " & len(mrna_alt) = ", len(mrna_alt)


        #retrieve subsequences for original & altered mRNA and return
        subseq_orig = mrna_orig[ start_i:end_i ]
        subseq_alt = mrna_alt[ start_i:end_i ]

        return [subseq_orig, subseq_alt]


    @staticmethod
    def retrieve_comparative_neoeps( mrna_orig, mrna_alt, aa_len ):
        """
        calculate neoepitopes between the original & altered mRNA sequence
        Args:
            -mrna_orig & mrna_alt = array that are RNA sequences (each element is a nucleotide), where "orig" is the original mRNA sequence & "mut" is the mutated mRNA sequence. These sequences should be corrected based on strand (+ or -), should only contain nucleotide bases (no "-")
            -aa_len = integer that is the number of amino acids to retrieve before the genomic alteration
            -alt_type = integer that is a certain type of change (mutation, insertion, deletion)
                -0 = mutation
                -1 = insertion
                -2 = deletion
                -3 = aberrant splicing
        PROTOCOL:
            -find the starting position of nucleotide difference
            -find the ending position where nucleotide
            -combine the codons & see if the nucleotide sequences match between original & mutation
            -
        """
        [subseq_orig, subseq_mut] = SimpleNeoepitopeIsoform.retrieve_mRNA_neoep_subseq_v2( mrna_orig, mrna_alt, aa_len )
        subseq_orig = ''.join( subseq_orig )
        subseq_mut = ''.join( subseq_mut )

        #translate sequence
        aa_orig = str( Seq( subseq_orig ).translate( to_stop = False ) )
        aa_alt = str( Seq( subseq_mut ).translate( to_stop = False ) )

        #retrieve the sequence before the termination codon "*"
        aa_alt_2 = aa_alt.split('*')[0]
        aa_orig_2 = aa_orig[0:len(aa_alt_2)]

        # return [aa_orig_2, aa_alt_2]
        """
        For the test output:
        -aa_orig = from mrna_orig, the full amino acid sequence that is the original, non-mutated sequence
        -aa_alt = from mrna_alt, the full alterated amino acid sequence
        -aa_orig_2 = same as aa_orig, but only the subsequence before the first stop signal identified in "aa_alt_2"
        -aa_alt_2 = same as aa_alt, but only the subsequence before the first stop signal identified in "aa_alt_2"
        """
        return [aa_orig, aa_alt, aa_orig_2, aa_alt_2]        #THIS IS JUST A TEST


    @staticmethod
    def retrieve_seq( chrom, arr_pos, path_genomeidx, strand_sign = 1, rev_comp = True ):
        """
        Retrieves the nucleotide sequence based on the genomic positions in array "arr_pos". BE CAREFUL: sequence.sequence() method is a cruzdb method that makes a request to the UCSC database for sequences, so this is not meant for high-throughput
        Args:
            -chrom = string that is the chromosome (format: chrNum, "chr9")
            -arr_pos = array of genomic positions
            -path_genomeidx = string that is the path to the samtools faidx genome file. NOTE: Make sure path_genomeidx & "build" match
            -strand_sign = integer that is the strand sign
            -rev_comp = reverse complement, only applies when strand_sign == -1. a boolean where True = returns the reverse complement of the sequence, False = returns just the complement
        """
        compile_seq = ""        #compile DNA sequence
        for c, genome_range in enumerate( arr_pos ):
            #retreive sequence
            exon_range = chrom + ":" + str(genome_range[0] + 1) + "-" + str(genome_range[1])
            compile_seq += SimpleNeoepitopeIsoform.sequence_dna( exon_range, path_genomeidx )

            ##TEST:: print "exon ", c, " - strand_sign = ", strand_sign, " & type( strand ) = ", type(strand_sign), " exon_range = ", exon_range, " & seq = ", SimpleNeoepitopeIsoform.sequence_dna( exon_range, path_genomeidx )

        #create BioSeq sequence object
        if strand_sign < 0:     #if strand sign is negative, then get reverse complement
            # seq_transcript = seq_transcript[::-1]     #reverse the string of nucleotides (assuming the exons entered are from right to left aka higher to lower)
            return str( Seq(compile_seq.upper(), IUPAC.unambiguous_dna).reverse_complement() ) if rev_comp else str( Seq(compile_seq.upper(), IUPAC.unambiguous_dna ).complement() )     ##QUES: what is IUPAC.unambiguous_dna?? WHAT DOES IT DO?
        else:       #else strand is positive
            return str( Seq(compile_seq.upper(), IUPAC.unambiguous_dna) )


    @staticmethod
    def sequence_dna( str_pos, path_genomeidx ):
        """ 
        Args:
            -str_pos = string that is the position of interest, in the format 'chrom:posStart-posEnd'
            -path_genomeidx = string that is the path to the samtools faidx genome file
        Function:
            retrieves the nucleotide sequence from the samtools index file in the path path_genome
        """
        #retrieve the genomic sequence
        exon_seq = subprocess.check_output( ["samtools", "faidx", path_genomeidx, str_pos] )

        #manipulate output to retrieve nucleotide sequence
        exon_seq = re.sub( r"\r\n?", "\n", exon_seq )     #replace any different carriage returns with the unix carriage return
        arr_exon_seq = exon_seq.split("\n")
        del arr_exon_seq[0]            #remove the first element, as this contains ">chr#:posStart-posEnd"
        exon_seq = "".join(arr_exon_seq)

        ##TEST:: print "str_pos = ", str_pos, ": ", Seq(exon_seq).reverse_complement().upper()

        return exon_seq



class SimpleNeoepitopeAll():
    def __init__( self, genome_range, orig, alt, path_genomeidx, strand = 1, opt_param = "variant_class=1&hgvs=1&refseq=1", build_hg38 = False, spec_isoform_id = None ):
        """
        -need to retrieve transcript ID (for RefSeq = "NM_", for Ensembl = "ENST")
        Args:
            -genome_pos = string that is genomic position of genomic alteration (format: chrom:start-end, BUT exclude 'chr' (e.g. 9:25-40, NOT chr9:25-40))
            -orig = string that is either None or '-', where '-' means it is an insertion. If it is any other alteration (mutation, deletion), then this should be "None"
                -None 
                -the param "orig" is suppose to be the reference nucleotide base, but I do not need the original nucleotide base for Ensembl's VEP
            -alt = string that is the genomic string change. Can be either nucleotide string (A,C,G,T) or '-' (meaning a deletion)
            -path_genomeidx = string that is the path to the samtools-indexed genome
            -opt_param = the optional parameteres for EnsemblVEP class (the EnsemblVEP class uses "vep/human/region" method for finding )
            -build_hg38 = string that retrieves
            -(MAYBE DELETE) spec_isoform_id (specific_isoform_id) = if I'm interested in only a specific isoform_id, specify it here. Make sure the isoform_id is of the same form as the database being queried (e.g. if Ensembl, then "ENST", else if Refseq, then "NM_")
        """
        self.obj_vep = EnsemblVEP( genome_range, orig, alt, strand, opt_param, build_hg38 )
        self.path_genomeidx = path_genomeidx
        self.spec_isoform_id = spec_isoform_id
        self.list_isoforms = []    #records a list of SimpleNeoepitopeIsoform() 
        self.r_stat, self.list_isoforms = self.define_isoform_alts()     #records a list of SimpleNeoepitopeIsoform()


    @staticmethod
    def calculate_full_rf_pos( cds_start, cds_end, len_seq, zero_base = True ):
        """
        this calculates a window starting around "cds_start" and ends around "cds_end", but will start at reading frame 0 & end at reading frame 2
        Args:
            -cds_start & cds_end = integers that refer to relative CDS positions on "self.mRNA", where cds_start <= cds_end.
            -len_seq = integer that is the length of the string (or array) that is the full sequence of interest, usually the full mRNA sequence.
            -zero_base = boolean where:
                -True = cds_start & cds_end are based on 0-based genomic positioning
                -False = cds_start & cds_end are based on 1-based genomic positioning
        """
        #retrieve the start & end 
        start_i = cds_start - 1 if zero_base else cds_start            #CONJ: I think -1 is for correcting for 0-base
        rf_start = start_i - (start_i % 3)      #need to find rf = 0
        end_i = cds_end - 1 if zero_base else cds_end 
        rf_end = end_i - (end_i % 3) + 2        #need to find rf = 2

        if rf_start < 0 or rf_end < 0:
            return [None, None]
        elif rf_start >= len_seq or rf_end >= len_seq:
            return [None, None]
        else:
            return [rf_start, rf_end]
    

    @staticmethod
    def determine_alteration( allele_string, alt_strand, gene_strand ):
        """
        Returns the original & altered nucleotide sequence
        Args:
            -allele_string = string in the format "orig_nuc/alt_nuc" (e.g. G/A, T/C, what about insertions & deletions???)
            -alt_strand = integer where the strand the alteration takes place (-1 = minus strand, 1 = plus strand)
            -gene_strand = integer that is the gene strand (-1 = minus strand, 1 = plus strand)
        Returns: a list where [0] = original nucleotide & [1] = altered nucleotide
        """
        #determine mutation
        str_allele = allele_string.replace('-', '')
        list_allele = str_allele.split('/')
        if alt_strand != gene_strand:
            list_allele = [ str( Seq(x.upper(), IUPAC.unambiguous_dna ).reverse_complement() ) for x in list_allele ]

        return list_allele

    def retrieve_isoform_alt( self ):
        """
        returns a hash of isoforms & HGVS information about the mutation
        """
        [r_stat, list_info] = self.obj_vep.get_all_alt_info()

        ##TEST::
        print "SN_ALL.RIA 0: r_stat = ", r_stat, " & list_info = ", list_info

        if r_stat < 1:
            return [r_stat, {}]

        hash_isoform_alt = {}      #k = isoform ID, v = hash where k2 = 'hgvsc' & 'hgvsp' & v2 = value for coding nucleotide & protein change, respectively
        for (hash_gen, hash_transcript) in list_info:
            ##TEST::
            print "\nSN_ALL.RIA: hash_gen = ", hash_gen.keys()
            print "\nSN_ALL.RIA: hash_transcript = ", hash_transcript.keys()


            for k_t, v_t in hash_transcript.iteritems():        #k_t = isoform_id, v_t = hash that contains information about transcript consequences for each isoform
                gene_symbol = v_t['gene_symbol']
                gene_id = v_t['gene_id']
                isoform_id = v_t['transcript_id']
                
                #determine relative mRNA change
                codon_change = '-' if not 'codons' in v_t else v_t['codons']
                cds_start = '-' if not 'cds_start' in v_t else v_t['cds_start']
                cds_end = '-' if not 'cds_end' in v_t else v_t['cds_end']

                list_allele = SimpleNeoepitopeAll.determine_alteration( hash_gen['allele_string'], self.obj_vep.strand, v_t['strand'] )

                #determine alteration
                # str_allele = hash_gen['allele_string'].replace('-', '')
                # list_allele = str_allele.split('/')
                # if obj_vep.strand != v_t['strand']:
                #     list_allele = [ str( Seq(x.upper(), IUPAC.unambiguous_dna ).reverse_complement() ) for x in list_allele ]
                
                #amino acid change
                aa_change = '-' if not 'amino_acids' in v_t else v_t['amino_acids']
                aa_start = '-' if not 'protein_start' in v_t else v_t['protein_start']
                aa_end = '-' if not 'protein_end' in v_t else v_t['protein_end']

                hash_isoform_alt[k_t] = {
                "gene_symbol": str( gene_symbol ), 
                "gene_id": str( gene_id ),
                "isoform_id": str( isoform_id ),
                "gene_strand": v_t['strand'],
                "allele_string": str( hash_gen['allele_string'] ),
                "variant_class": None if not "variant_class" in hash_gen else str( hash_gen['variant_class'] ),
                "nuc_orig": str( list_allele[0] ),
                "nuc_alt": str( list_allele[1] ), 
                "codon_change": str( codon_change ),
                "chrom": str( hash_gen['seq_region_name'] ),
                "genome_start": int( hash_gen['start'] ),
                "genome_end": int( hash_gen['end'] ),
                "cds_start": cds_start,
                "cds_end": cds_end,
                "aa_change": str( aa_change ),
                "aa_start": aa_start,
                "aa_end": aa_end,
                "consequence": hash_gen['most_severe_consequence'], 
                "hgvsc": None if not "hgvsc" in v_t else str( v_t['hgvsc'] ),
                "hgvsp": None if not "hgvsp" in v_t else str( v_t['hgvsp'] )
                }

        return [r_stat, hash_isoform_alt]


    # def retrieve_isoform_hgvs( self ):
    #     """
    #     returns a hash of isoforms & HGVS information about the mutation
    #     """
    #     [r_stat, list_info] = self.obj_vep.get_all_alt_info()
    #     if r_stat < 1:
    #         return {}

    #     hash_isoform_hgvs = {}      #k = isoform ID, v = hash where k2 = 'hgvsc' & 'hgvsp' & v2 = value for coding nucleotide & protein change, respectively
    #     for (hash_gen, hash_transcript) in list_info:
    #         isoform_id = hash_transcript['transcript_id']
    #         HGVSc = hash_transcript['hgvsc']        #format: NM_002880.3:c.935T>C
    #         HGVSp = hash_transcript['hgvsp']        #format: NP_002871.1:p.Val312Ala
    #         hash_isoform_hgvs[isoform_id] = {"HGVSc": HGVSc, "HGVSp": HGVSp}

    #     return hash_isoform_hgvs

    def define_isoform_alts( self ):
        """
        defines each instance of SimpleNeoepitopeIsoform by using the isoform ID & the HGVS associated information
        Returns: an array where [0] = the request status (see class EnsemblVEP def get_vep_request()) & [1] = an array of SimpleNeoepitopeIsoform instances
        """
        [r_stat, hash_isoform_alt] = self.retrieve_isoform_alt()

        list_isoforms = []
        for k,v in hash_isoform_alt.iteritems():       #k = isoform ID, v = hash where k2 = contains information about gene including gene symbol, isoform ID, cds_start & cds_end of genomic alteration

            if self.spec_isoform_id:
                if k != self.spec_isoform_id:
                    continue

            ##TEST:: print "SN_ALL.define_isoform_alts 3 loop:: k = ", k, " & v = ", v

            list_isoforms.append( SimpleNeoepitopeIsoform( v, self.path_genomeidx ) )

        return [r_stat, list_isoforms]

    