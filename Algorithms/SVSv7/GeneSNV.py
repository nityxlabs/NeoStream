#/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Isoform import Isoform
from IsoformSJ import IsoformSJ
from MultiIsoform import MultiIsoform
from TranscribeTranslate_V4 import TranscribeTranscript


"""
Class: Looks at Single Nucleotide Variants
"""
class GeneSNV():
    def __init__( self, chrom, start, end, base_alt, obj_cruzdb, path_genomeidx, gene_sym = None, isoform_id = None ):
        """
        Args:
            chrom = string that is the chromosome (e.g. chr12, chr4)
            base_alt = character that is the mutation nucleotide base. Should be any of the following: A, T, G, C
            obj_cruzdb = instance of CruzDB, will be used by MultiIsoform class
            path_genomeidx = string that is the path to the samtools-indexed genome, will be used with TranscribeTranscript class
        """
        #assign cruzdb instance to global variable in Isoform
        Isoform.set_cruzdb( obj_cruzdb )
        self.snv_chrom = chrom
        self.snv_start = start
        self.snv_end = end
        self.base_alt = base_alt
        self.path_genomeidx = path_genomeidx
        self.obj_mi = MultiIsoform( chrom, start, start, gene_sym, isoform_id )


    def pos_in_cds( self, pos, bool_contained ):
        """
        Args:
            -pos = integer that is genomic position (assumed the chromosome is the same as the isoform for self object). USE = find the genomic position that contains mutation
            -bool_contained = boolean where
                -True = returns hash where isoforms that contain position in exon will be return.
                -False = returns hash of all isoforms in MultiIsoform
        Function: determines if position is in CDS of all isoforms, and returns hash where k = isoform id (string), v = Exon object or None, depends on if position is found in cds
        """
        if bool_contained:
            hash_isoforms = self.obj_mi.in_cds_all( pos )      #returns a hash where k = isoform id (string), v = Exon object or None, depends on if position is found in cds

            ##TEST::
            print "GeneSNV pos_in_cds TEST:"
            for k,v in hash_isoforms.iteritems():
                print "k = ", k, " & v = ", v

            ##TEST::
            filter_hash = {str(k):v for k,v in hash_isoforms.iteritems() if v}
            print "GeneSNV pos_in_cds TEST - filter = ", filter_hash

            # return {k:v for k,v in hash_isoforms.iteritems() if not v}
            return filter_hash
        else:
            return self.obj_mi.in_cds_all( pos )      #returns a hash where k = isoform id (string), v = Exon object or None, depends on if position is found in in_cds_all

    def pos_in_exon( self, pos ):
        """
        Args:
            -pos = integer that is genomic position (assumed the chromosome is the same as the isoform for self object). USE = find the genomic position that contains mutation
            -bool_contained = boolean where
                -True = returns hash where isoforms that contain position in exon will be return.
                -False = returns hash of all isoforms in MultiIsoform
        Function: determines if position is in CDS of all isoforms, and returns hash where k = isoform id (string), v = Exon object or None, depends on if position is found in exon
        """
        if bool_contained:
            hash_isoforms = self.obj_mi.in_exon_all( pos )      #returns a hash where k = isoform id (string), v = Exon object or None, depends on if position is found in exon
            return {str(k):v for k,v in hash_isoform.iteritems() if not v}
        else:
            return self.obj_mi.in_exon_all( pos )      #returns a hash where k = isoform id (string), v = Exon object or None, depends on if position is found in exon

    def pos_in_intron( self, pos ):
        """
        Args:
            -pos = integer that is genomic position (assumed the chromosome is the same as the isoform for self object). USE = find the genomic position that contains mutation
            -bool_contained = boolean where
                -True = returns hash where isoforms that contain position in exon will be return.
                -False = returns hash of all isoforms in MultiIsoform
        Function: determines if position is in CDS of all isoforms
        """
        if bool_contained:
            hash_isoforms = self.obj_mi.in_intron_all( pos )      #returns a hash where k = isoform id (string), v = Exon object or None, depends on if position is found in intron
            return {str(k):v for k,v in hash_isoforms.iteritems() if not v}
        else:
            return self.obj_mi.in_intron_all( pos )      #returns a hash where k = isoform id (string), v = Exon object or None, depends on if position is found in intron

    """
    Function: determine SNV type (synonymous, non-synonymous)
    """

    def determine_snv_type( self, pos, mut_nuc, list_isoform_id ):
        """
        Args:
            -pos = integer that is genomic position (assumed the chromosome is the same as the isoform for self object). USE = find the genomic position that contains mutation
            -mut_nuc = character that is the mutated nucleotide. Should be any of the following (in upper case) - A, T, C, G
            -list_isoform_id = an array of isoform IDs of interest. This should come from pos_in_cds()
        Function: determines if SNV is synonymous or non-synonymous
        """
        #find the codon of the SNV - 1st create canonical transcript of isoform, 2nd 
        #create IsoformSJ instance so I can retrieve the canonical transcript for the isoform
        isoform_snv_type = {}      #k = isoform id, v = hash where {'codon_orig': str(codon_orig), 'codon_mut': str(codon_mut), 'aa_orig': aa_orig, 'aa_mut': aa_mut, 'snv_type': None}
        list_isoform_id = [str(x) for x in list_isoform_id]     #make sure all elements are string
        for iso_id in list_isoform_id:
            #skip any non-coding isoform IDs (i.e. contain 'NR_', these are RefSeq IDs)
            if 'NR_' in iso_id:
                isoform_snv_type[iso_id] = {}
                continue

            iso_sj, obj_tt = self.create_transcript_instances( iso_id )

            #get codon associated with position
            hash_codon = self.mutate_codon_from_pos( pos, mut_nuc, obj_tt )
            if not hash_codon:
                isoform_snv_type[iso_id] = {}
                continue


            ##TEST:: print "GeneSNV dSNVType 1 - hash_codon['str_codon'] = ", hash_codon['str_codon'], " & type = ", type( hash_codon['str_codon'] ), " & hash_codon['str_codon_mut'] = ", hash_codon['str_codon_mut'], " & type = ", type( hash_codon['str_codon_mut'] )

            codon_orig = Seq( hash_codon['str_codon'], IUPAC.unambiguous_dna )
            codon_mut = Seq( hash_codon['str_codon_mut'], IUPAC.unambiguous_dna )
            aa_orig = codon_orig.translate()
            aa_mut = codon_mut.translate()

            ##TEST:: print "GeneSNV dSNVType 2 - str_codon = ", codon_orig, " & mut_codon = ", codon_mut, " & aa_orig = ", aa_orig, " & aa_mut = ", aa_mut

            get_strand = None if not (iso_id in self.obj_mi.hash_isoforms) else self.obj_mi.hash_isoforms[iso_id].strand
            isoform_snv_type[iso_id] = {'hash_pos_nuc': hash_codon['hash_pos_nuc'], 'hash_pos_nuc_mut': hash_codon['hash_pos_nuc_mut'], 'strand': get_strand, 'codon_orig': str(codon_orig), 'codon_mut': str(codon_mut), 'aa_orig': str(aa_orig), 'aa_mut': str(aa_mut), 'snv_type': None}

            #need to translate codon before & after mutation to see if changes amino acid
            if aa_orig != aa_mut:
                isoform_snv_type[iso_id]['snv_type'] = 'non-synonymous'
            else:
                isoform_snv_type[iso_id]['snv_type'] = 'synonymous'

            ##TEST:: print "GeneSNV dSNVType 3 - isoform_snv_type = ", isoform_snv_type

        return isoform_snv_type

    
    #PUT THIS FUNCTION IN TranscrbieTranslate_V4.py
    @staticmethod
    def mutate_codon_from_pos( pos, mut_nuc, obj_tt ):
        """
        Args:
           -pos = integer that is the genomic position in 
           -mut_nuc = character that is the mutated nucleotide. Should be any of the following (in upper case) - A, T, C, G
           -obj_tt = instance of TranscribeTranscript class, will be used to retrieve the codon
        Function: creates point mutation based on position 'pos' & nucleotide at 'nuc'. Returns result from function 'def get_codon_from_pos()' but also adds a key where it contains the position and nucleotide.
        """
        hash_codon = GeneSNV.get_codon_from_pos( pos, obj_tt )
        if not hash_codon:
            return {}

        ##TEST:: print "GeneSNV: mutate_codon_from_pos -> hash_codon['hash_pos_nuc'] = ", hash_codon['hash_pos_nuc'], " & pos = ", pos, " & mut_nuc = ", mut_nuc

        #create mutated version of codon
        codon_mut = hash_codon['hash_pos_nuc'].copy()
        codon_mut[pos] = mut_nuc
        arr_codon_mut = codon_mut.values()[::-1] if obj_tt.iso_sj.strand < 0 else codon_mut.values()
        str_codon_mut = ''.join( arr_codon_mut )

        #record mutation information to hash codon
        hash_codon['hash_pos_nuc_mut'] = codon_mut
        hash_codon['str_codon_mut'] = str_codon_mut

        return hash_codon

    @staticmethod
    def get_codon_from_pos( pos, obj_tt ):
        """
        Args:
            -pos = integer that is the genomic position (do not need 'chrom')
            -obj_tt = instance of TranscribeTranscript class, will be used to retrieve the codon
        Function: get the full codon that resides at a position
        """
        try:
            i_pos = obj_tt.arr_genome_pos.index( pos )
        except ValueError:      #this means the position is not found in arr_genome_pos (contains all exon positions)
            print "GeneSNV get_codon_from_pos() Error: position ", pos, " is not in obj_tt.arr_genome_pos"
            return {}

        rf = obj_tt.arr_rf[i_pos]
        #if no reading frame associated with position, return None. rf = -1 means there is no reading frame present
        if rf < 0:
            print "GeneSNV get_codon_from_pos() Error: position ", pos, " does not have reading frame. RF = ", rf
            return {}
        
        #check if the codon has the possibility of being out of range of the arrays obj_tt.arr_genome_pos (and effectively obj_tt.arr_rf & obj_tt.arr_nuc_seq)
        if i_pos < 2 or len( obj_tt.arr_genome_pos ) - i_pos <= 2 - rf:
            print "GeneSNV get_codon_from_pos() Error: position ", pos, " out of range. i_pos = ", i_pos, " & rf = ", rf, " & len( obj_tt.arr_genome_pos ) = ", len( obj_tt.arr_genome_pos )
            return {}

        #retrieve the starting & ending position for each codon
        if obj_tt.iso_sj.strand < 0:
            i_pos_start = i_pos + (rf - 2)
            # i_pos_end = i_pos + rf
            i_pos_end = i_pos_start + 2
        else:
            i_pos_start = i_pos - rf
            # i_pos_end = i_pos + (2 - rf)
            i_pos_end = i_pos_start + 2

        #get codon for position
        arr_codon = obj_tt.arr_nuc_seq[i_pos_start:i_pos_end + 1][::-1] if obj_tt.iso_sj.strand < 0 else obj_tt.arr_nuc_seq[i_pos_start:i_pos_end + 1]      #need to add +1 else will only return nucleotides instead of 3 (need 3 for codon)
        str_codon = ''.join( arr_codon )

        #create hash where key = genomic position (integer only) & value = nucleotide
        hash_pos_nuc = {}       #key = genomic position (integer only) & value = nucleotide
        hash_TEST_rf = {}
        for i in range( i_pos_start, i_pos_end + 1 ):
            actual_pos = obj_tt.arr_genome_pos[i]
            hash_pos_nuc[actual_pos] = obj_tt.arr_nuc_seq[i]
            hash_TEST_rf[actual_pos] = obj_tt.arr_rf[i]


        ##TEST::
        print "GeneSNV - Get Codon: "
        print "pos = ", pos, " & rf = ", rf
        print "hash_pos_nuc = ", hash_pos_nuc
        print "hash_TEST_rf = ", hash_TEST_rf
        print ">>>>>>>>>>>>>>>>>>>--------->>>>>>>>>>>"

        return {'str_codon': str_codon, 'hash_pos_nuc': hash_pos_nuc}


    def create_transcript_instances( self, isoform_id ):
        """
        Args:
            -isoform_id = string that is the isoform ID - this isoform_id should be a key in the MultiIsoform instance 'self.obj_mi' (though I don't check that in this function do I...)
        Function: creates an instance of classes IsoformSJ & TranscribeTranscript
        """
        hash_pos = { 'chrom': self.snv_chrom, 'pos_oi': self.snv_start }      #used to find the closest position for an isoform
        bool_simulant_sj = False
        group_sj = 0
        iso_sj = IsoformSJ( isoform_id, [], -10, hash_pos, bool_simulant_sj, group_sj )
        canon_transcript = iso_sj.create_canon_transcript()


        ##TEST::
        # print "GeneSNV canon_transcript = ", canon_transcript,  " & iso_sj: strand = ", iso_sj.strand, " & isoform_id = ", iso_sj.isoform_id, " & gene_sym = ", iso_sj.gene_sym, " & is_coding = ", iso_sj.is_coding, " & boundary = ", iso_sj.boundary, " & coding_boundary = ", iso_sj.coding_boundary
        # for i, each in enumerate( canon_transcript ):
        #     print "SJ ", i, "  - ", each


        obj_tt = TranscribeTranscript( canon_transcript, iso_sj, self.path_genomeidx, {} )

        ##TEST::
        print "GeneSNV CTI - SHOW EXONS: "
        for i, each_exon in enumerate(obj_tt.list_exons):
            print "exon ", i, " - ", each_exon
        print "GeneSNV CTI - self.snv_start = ", self.snv_start


        return [iso_sj, obj_tt]

    # def check_in_gene():
    #     obj_mi = MultiIsoform( chrom, start, end, gene_sym = None, isoform_id = None )

    #     obj_mi.in_exon_all( position )
    #     obj_mi.in_cds_all( position )

    #     obj_tt = TranscribeTranscript( transcript_sj, iso_sj, DIR_GENOME, {} )

    #     hash_pos = { 'chrom': each_gene.chrom, 'pos_oi': each_gene.introns[0][0] }      #used to find the closest position for an isoform
    #     bool_simulant_sj = False
    #     group_sj = 0
    #     iso_sj = IsoformSJ( isoform_id, [], -10, hash_pos, bool_simulant_sj, group_sj )
    #     canon_transcript = iso_sj.create_canon_transcript()
    #     obj_tt = TranscribeTranscript( canon_transcript, iso_sj, DIR_GENOME, {} )