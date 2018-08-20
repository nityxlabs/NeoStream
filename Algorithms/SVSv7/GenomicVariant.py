#/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Isoform import Isoform
from IsoformSJ import IsoformSJ
from MultiIsoform import MultiIsoform
from TranscribeTranslate_V5 import TranscribeTranscript, TranslateTranscript

"""
Needs for TTV5/GenomicVariant
-need a way to see which isoforms contain the position of interest (e.g. in the exons) - CONJ: I think GenomicVariant
-do I need this function to find overlapping insertions/deletions
    def is_overlapping(x1,x2,y1,y2):
        return max(x1,y1) <= min(x2,y2)
Questions for KM:
    -should I just pick the first isoform? I've noticed sometimes 
    -Should I consider multiple mutations for a single neoepitope? Right now I only consider 1 alteration at a time (well let's just tell KM that), but for the neoepitopes with gene reactivation, maybe we should consider multiple mutations
"""

"""
Class: Think of this as GeneSNV part 2, however this works
"""
class GenomicVariant():
    #class variables
    hash_change_type = {'SNP': 0, 'DNP': 0, 'ONP': 0, 'INS': 1, 'DEL': 2}       #SNP = Singe Nucleotide Polymorphism (basically treated as a mutation aka SNV), DNP = Dinucleotide Polymorphism (2 base changes), ONP = Oligoucleotide Polymorphism (more than 2 nucleotide changes)

    def __init__( self, chrom, start, end, strand, base_orig, base_alt, alt_type, path_genomeidx, gene_sym = None, isoform_id = None ):
        """
        Args:
            -chrom = string that is the chromosome (e.g. chr12, chr4)
            -base_orig = character that is the original nucleotide base. Should be any of the following: A, T, G, C
            -base_alt = character that is the mutation nucleotide base. Should be any of the following: A, T, G, C
            -alt_type = string that describes the genomic alteration type
            -path_genomeidx = string that is the path to the samtools-indexed genome, will be used with TranscribeTranscript class
            -gene_sym = string that, if defined, will only record isoforms that are assigned to this gene symbol
            -isoform_id = optional parameter, if this is provided then will only record one isoform. If gene_sym is also defined, then isoform_id needs to be an isoform of gene_sym else it will not be created
        Before Using this Class
            -make sure to use Isoform.set_cruzdb() to set obj_cruzdb         
        """
        #assign cruzdb instance to global variable in Isoform
        # Isoform.set_cruzdb( obj_cruzdb )      
        self.gene_sym = gene_sym
        self.isoform_id = isoform_id
        self.snv_chrom = chrom if 'chr' in chrom else 'chr' + str(chrom)
        self.snv_start = int( start )
        self.snv_end = int( end )
        self.snv_strand = -1 if strand == '-' else 1        #+1 if plus strand & -1 if minus strand
        self.base_orig = base_orig
        self.base_alt = base_alt
        self.alt_type = GenomicVariant.hash_change_type[alt_type]
        self.path_genomeidx = path_genomeidx
        self.obj_mi = MultiIsoform( self.snv_chrom, self.snv_start, self.snv_end, self.gene_sym, self.isoform_id )
        

        ##TEST::
        # print "GV info:"
        # print "chrom = ", self.snv_chrom
        # print "snv_start = ", self.snv_start, " & type = ", type( self.snv_start )
        # print "snv_end = ", self.snv_end, " & type = ", type( self.snv_end )
        # print "snv_strand = ", self.snv_strand
        print "base_orig = ", self.base_orig
        print "base_alt = ", self.base_alt
        # print "path_genomeidx = ", self.path_genomeidx
        # print "gene_sym = ", self.gene_sym
        # print "isoform_id = ", self.isoform_id

    def get_genomic_range( self ):
        """
        Returns string format of the genomic alteration
        """
        return self.snv_chrom + ':' + str( self.snv_start ) + '-' + str( self.snv_end )

    def check_orig_base_isoform( self, isoform_id ):
        """
        Checks if the nucleotide from instance of TranslateTranscript & self.base_orig are the same

        Args:
            isoform_id = string that is the isoform ID, will be used to recreate transcript
        Output:
            returns boolean, where True = nucleotide from instance of TranslateTranscript & self.base_orig match, else False if they do not match
        """
        #create transcript isoform based on isoform_id
        obj_tt = self.create_transcript_instances( isoform_id )

        #double check to see if the original nucleotide reported (self.base_orig) & nucleotide in hg19.fa
        i_genome_pos = obj_tt.arr_genome_pos.index( self.snv_start )
        check_orig_nuc = obj_tt.arr_nuc_seq[ i_genome_pos ]
        if self.snv_strand != obj_tt.iso_sj.strand:
            check_orig_nuc = str( Seq(check_orig_nuc).reverse_complement().upper() )
        #if True, means the original base are the same, else if they are different then I don't why this would happen
        status_orig_base = True if check_orig_nuc.upper() == self.base_orig.upper() else False

        return status_orig_base

    def determine_aa_change( self ):
        """
        Determine the effect of the variant on the AA sequence

        PROTOCOL: retrieve the transcript & position (from class TranscribeTranscript ) -> make changes to positions of interest
        """
        for k,v in self.obj_mi.hash_isoforms.iteritems():        #k = string that is isoform_id, v = Isoform instance
            obj_tt = self.create_transcript_instances( k )

            #METHOD 1: get the original codon & mutated codon
            # orig_codon = obj_tt.retrieve_containing_codon( self.snv_start, self.snv_strand )
            # i_genome_pos = obj_tt.arr_genome_pos.index( self.snv_start )
            # obj_tt.arr_nuc_seq[ i_genome_pos ] = self.base_alt
            # mut_codon = obj_tt.retrieve_containing_codon( self.snv_start, self.snv_strand )


            #METHOD 2: get the mutated codon
            full_pos = self.snv_chrom + ':' + str( self.snv_start ) + '-' + str( self.snv_end )
            hash_codon_info = obj_tt.get_mutated_codon( self.base_orig, self.base_alt, full_pos, self.snv_strand, True )      #output is hash -> {'codon_orig': codon_orig, 'codon_mut': codon_mut, 'aa_orig': aa_orig, 'aa_mut': aa_mut}



            ##TEST:: show the AA change based on mutation
            # print "hash_codon_info: "
            # print hash_codon_info

            # print "gene strand & snv strand: ", obj_tt.iso_sj.strand, " & ", self.snv_strand
            # print "original base > mutated base: ", self.base_orig, " > ", self.base_alt
            # print "original codon > mutated codon: ", hash_codon_info['codon_orig'], " > ", hash_codon_info['codon_mut']
            # print "original AA > mutated AA: ", hash_codon_info['aa_orig'], " > ", hash_codon_info['aa_mut']


            ##TEST:: determine consequence
            print "GV_DAAC 1: "
            obj_tt.alteration_consequence( self.base_alt, self.get_genomic_range(), self.snv_strand, self.alt_type )
            

            ##TEST METHOD - SEE WHAT STEPS I NEED TO PERFORM
            #TEST:: retrieve the original base & the mutated base
            # i_genome_pos = obj_tt.arr_genome_pos.index( self.snv_start )
            # orig_base = obj_tt.arr_nuc_seq[ i_genome_pos ]
            # print "k = ", k, " & i_genome_pos = ", i_genome_pos, " | orig_base = ", orig_base, " & double_check = ", self.base_orig, " & iso_sj.strand = ", obj_tt.iso_sj.strand, " & mut strand = ", self.snv_strand
            # hash_orig_codon = obj_tt.find_containing_codon( self.snv_start )
            # print "hash_orig = ", hash_orig_codon
            # get_orig_codon = obj_tt.arr_nuc_seq[ hash_orig_codon['i_genome_start']:hash_orig_codon['i_genome_end'] + 1 ]
            # str_orig_codon = ''.join( get_orig_codon ) if obj_tt.iso_sj.strand > 0 else ''.join( get_orig_codon[::-1] )
            # print "seq_orig = ", str_orig_codon, " & type = ", type( get_orig_codon ), " & rf = ", obj_tt.arr_rf[ hash_orig_codon['i_genome_start']:hash_orig_codon['i_genome_end'] + 1 ], " & list_orig_codon = ", get_orig_codon

            # ##TEST:: make mutation
            # obj_tt.arr_nuc_seq[ i_genome_pos ] = self.base_alt
            # hash_mut_codon = obj_tt.find_containing_codon( self.snv_start )
            # print "hash_muts = ", hash_mut_codon
            # get_mut_codon = obj_tt.arr_nuc_seq[ hash_mut_codon['i_genome_start']:hash_mut_codon['i_genome_end'] + 1 ]
            # str_mut_codon = ''.join( get_mut_codon ) if obj_tt.iso_sj.strand > 0 else ''.join( get_mut_codon[::-1] )
            # print "seq_muts = ", str_mut_codon, " & type = ", type( get_mut_codon ), " & rf = ", obj_tt.arr_rf[ hash_mut_codon['i_genome_start']:hash_mut_codon['i_genome_end'] + 1 ], " & list_mut_codon = ", get_mut_codon 

            # ##TEST:: retrieve 
            # print "AA: from ", Seq( str_orig_codon ).translate( to_stop = False ), ">", Seq( str_mut_codon ).translate( to_stop = False )

            # try:
            #     i_genome_pos = obj_tt.arr_genome_pos.index( self.snv_start )
            #     orig_base = obj_tt.arr_nuc_seq[ i_genome_pos ]
            #     print "k = ", k, " & i_genome_pos = ", i_genome_pos, " | orig_base = ", orig_base, " & double_check = ", self.base_orig, " & iso_sj.strand = ", obj_tt.iso_sj.strand, " & mut strand = ", self.snv_strand
            #     hash_orig_codon = obj_tt.find_containing_codon( self.snv_start )
            #     print "hash_orig = ", hash_orig_codon
            #     get_orig_codon = obj_tt.arr_nuc_seq[ hash_orig_codon['i_genome_start']:hash_orig_codon['i_genome_end'] ]
            #     print "seq_orig = ", get_orig_codon

            #     ##TEST:: make mutation
            #     obj_tt.arr_nuc_seq[ i_genome_pos ] = self.base_alt
            #     hash_mut_codon = obj_tt.find_containing_codon( self.snv_start )
            #     print "hash_muts = ", hash_mut_codon
            #     get_mut_codon = obj_tt.arr_nuc_seq[ hash_mut_codon['i_genome_start']:hash_mut_codon['i_genome_end'] ]
            #     print "seq_muts = ", get_mut_codon 

            #     ##TEST:: retrieve 
            #     print "AA: from ", Seq( orig_codon ).translate( to_stop = False ), ">", Seq( mut_codon ).translate( to_stop = False )
            # except:
            #     print "ERROR:: for ", k, ", position does not exist: ", self.snv_start
            #     continue

            print "////////////////////\n"


    def create_isosj_instance( self, isoform_id ):
        """
        Creates an instance of the class IsoformSJ. Can be used for creating an instance of TranscribeTranscript or TranslateTranscript
        """
        hash_pos = { 'chrom': self.snv_chrom, 'pos_oi': self.snv_start }      #used to find the closest position for an isoform
        bool_simulant_sj = False
        group_sj = 0        #this means splicing events will NOT be grouped into 5' or 3' competitive splicing
        iso_sj = IsoformSJ( isoform_id, [], -10, hash_pos, bool_simulant_sj, group_sj )
        
        return iso_sj

    def create_transcript_instances( self, isoform_id ):
        """
        Args:
            -isoform_id = string that is the isoform ID - this isoform_id should be a key in the MultiIsoform instance 'self.obj_mi' (though I don't check that in this function do I...)
        Function: creates an instance of classes IsoformSJ & TranscribeTranscript
        """
        iso_sj = self.create_isosj_instance( isoform_id )
        canon_transcript = iso_sj.create_canon_transcript()
        #reconstruct canonical transcript
        obj_tt = TranslateTranscript( canon_transcript, iso_sj, self.path_genomeidx, {} )

        ##TEST::
        # print "GeneSNV CTI - SHOW EXONS: "
        # for i, each_exon in enumerate(obj_tt.list_exons):
        #     print "exon ", i, " - ", each_exon
        # print "GeneSNV CTI - self.snv_start = ", self.snv_start

        return obj_tt

    def define_variant_type( self ):
        """
        Determine if the variant type is synonymous, nonsynonymous, insertion, or deletion
        """
        pass

