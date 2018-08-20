#/usr/bin/python

from Bio.SeqFeature import SeqFeature, FeatureLocation
# from Bio.SeqFeature import SeqFeature, FeatureLocation, BeforePosition, AfterPosition

#NOTE: using ( object ) is the new style of identifying a class
class Exon( object ):
    def __init__( self, hash_exon_info ):
        """
        Args:
            hash_exon_info = hash with keys & the following information
                -posStart = integer that is the position, posStart < posEnd
                -posEnd = integer 
                -strand_sign = either 1 (plus) or -1 (minus)
                -exon_frame = integer that is 0, 1, or 2. Refers to the reading frame of the 5' end of exon (lower position for plus genes & higher position for minus genes). It is -1 if the exon is non-coding & None if reading frame is unknown
                -exon_num = exon number in gene. For plus genes, exon number increases from lower to higher position. For minus genes, exon number decreases from lower to higher position. If None, means unknown what the exon number is (could be non-canonical exon number)
                -chrom = chromosome number, format "chr#" (e.g chr9)
                -canonical = boolean that, if true, means the exon is a canonical exon in the gen
                -isoformID (optional) = gene isoform
                -feat_type = string that conveys the type of feature the element is, either 'exon' or 'intron'
        -CAUTION:
            -Need to be careful of 0-based positions of exons because the first position of the exon isn't the start of the exon but instead the last nucleotide of the previous exon. Therefore, especially when retrieving nucleotide sequence for the exon, need to add +1 to the start of the exon.
        """
        self.chrom = hash_exon_info['chrom']
        self.canonical = hash_exon_info['canonical']

        self.readFrame = hash_exon_info['exon_frame']      #reading frame for the 5' end of the exon (lower position for + genes & higher position for - genes)
        self.exonNum = int( hash_exon_info['exon_num'] ) if hash_exon_info['exon_num'] else str( hash_exon_info['exon_num'] )      #exon number. If None, then exon number is not known

        self.constitutive = None        #boolean that, if True, means it is present in all isoforms of a particular gene, else if False, means it is not constitutive element -> assigned by def set_constitutive()

        #NOTE: for FeatureLocation, it is 0-based (not 1-based), and the range is [start, end), therefore anything within start to end-1 will be within this FeatureLocation
        self.exonPos = SeqFeature( FeatureLocation( int(hash_exon_info["pos_start"]), int(hash_exon_info["pos_end"]) ), strand = hash_exon_info["strand_sign"], type = hash_exon_info["feat_type"] )
        self.cdsPos = None      #SeqFeature object similar to self.exonPos, but records CDS range (0-based - position flanking nucleotides)
        self.cds_base_1 = None      #same as cdsPos, but is 1-based, meaning position is on nucleotides (as opposed to flanking nucleotides in 0-based)
        
        #append isoformID if it exists
        self.arrAllIsoforms = []
        if hash_exon_info[ "isoform_id" ]:
        # if "isoform_id" in hash_exon_info:
            self.arrAllIsoforms.append( str( hash_exon_info[ "isoform_id" ] ) )

        ##QUESTION: what about exons that are modified? Should I use a boolean value to say true or false? Also, with modified exons, there could be an isoform with the same exon number for a normal exon & a modified exon
        # self.boolModExon = hash_exon_info["boolModExon"]

    def get_exon_info( self ):
        """
        Function: creates & returns hash that is information about this exon. Adopted from function Isoform.get_feat_info()
        """
        #prepare hash that will record elements into exon list
        exon_info = {
        "pos_start": self.exonPos.location.start,
        "pos_end": self.exonPos.location.end,
        "strand_sign": self.exonPos.strand,
        "chrom": self.chrom,
        "isoform_id": self.arrAllIsoforms,
        "exon_frame": self.readFrame,
        "exon_num": self.exonNum,
        "canonical": self.canonical,
        "feat_type": self.exonPos.type }

        return exon_info

    def __str__( self ):
        #retrieve the 5' end for the exon
        rf_pos = None
        if isinstance( self.exonPos.strand, int ):
            rf_pos = self.exonPos.location.end if self.exonPos.strand < 0  else self.exonPos.location.start
        return self.exonPos.type + " " + str( self.exonNum ) + " - " + self.chrom + ":" + str( self.exonPos.location.start ) + "-" + str( self.exonPos.location.end ) + ", canonical = " + str( self.canonical ) + ", constitutive = " + str( self.constitutive ) + ", isoforms = " + ','.join( self.arrAllIsoforms ) + ", strand = " + str( self.exonPos.strand ) + ", & at 5' end = " + str( rf_pos ) + ", the reading frame is " + str( self.readFrame )

    def __eq__( self, objExon ):
        """ Function: checks if this is the same exon with comparison '==' """
        if not objExon:     #this means the comparison is between Exon object & a None type object
            return False

        if ( self.chrom == objExon.chrom ) and ( self.exonPos.location.start == objExon.exonPos.location.start ) and ( self.exonPos.location.end == objExon.exonPos.location.end ):
            return True
        else:
            return False

    def str_genomic_pos( self, str_only = True, add_1 = False ):
        """
        Args:
            str_only = boolean that:
                -True = will return in string format (e.g chrom:start-end)
                -False = returns in array format, where [chrom, start, end], where chrom = string, start & end = int
            -add_1 = boolean that is used to add +1 to the start of the exon. I do this because the position reported for the exon range is usually 0-based. Therefore if an exon position is reported, it may report not from the start of the exon to the end, but from -1 of the start exon (the end of the intron) to the end of the exon. This problematic, especially when retrieving nucleotide sequences using genomic position.
                -0-based is where the gap between nucleotides is enumerated, whereas 1-based means the actual nucleotides are enumerated
                -True = will add +1 to the start of the exon
                -False = will not add +1 to the start of the exon
        Function: returns position of gene in string form
        """
        start_pos = self.exonPos.location.start + 1 if add_1 else self.exonPos.location.start

        if str_only:
            return str( self.chrom ) + ":" + str( start_pos ) + "-" + str( self.exonPos.location.end )
        else:
            return {'chrom': str( self.chrom ), 'start': int( start_pos ), 'end': int(self.exonPos.location.end) }
            # return [str( self.chrom ), int( self.exonPos.location.start ), int(self.exonPos.location.end)]

    def add_cds( self, pos_start, pos_end ):
        """ 
        Args:
            pos_start & pos_end = integers that are the start & end genomic positions of the CDS, respectively. These need to be integers as SeqFeature.FeatureLocation only accepts integers
        Function: adds start & end position of CDS 
        """
        pos_start = int( pos_start )
        pos_end = int( pos_end )
        if self.in_exon( pos_start ) and self.in_exon( pos_end ):
            self.cdsPos = SeqFeature( FeatureLocation( int(pos_start), int(pos_end) ), strand = self.exonPos.strand, type = "cds" )
            self.cds_base_1 = SeqFeature( FeatureLocation( int(pos_start + 1), int(pos_end) ), strand = self.exonPos.strand, type = "cds_base_1" )

    def set_constitutive( self, stat_constitutive ):
        """ 
        Function: sets constitutive attribute to either True or False
        """
        self.constitutive = stat_constitutive

    """
    Functions: see if position is located in exon/cds
    """
    def in_exon( self, position ):
        """ Function: checks if position is within exon """
        if position == self.exonPos.location.start:
            return "exonLeft"
        elif position == self.exonPos.location.end:
            return "exonRight"
        elif position in self.exonPos:
            return "withinElem"
        else:
            return None

    def in_cds( self, position ):
        """ Function: checks if position is within cds """
        if not self.cdsPos:
            return None

        if position == self.cdsPos.location.start:
            return "exonLeft"
        elif position == self.cdsPos.location.end:
            return "exonRight"
        elif position in self.cdsPos:
            return "withinElem"
        else:
            return None

    def in_cds_base_1( self, position ):
        """ Function: checks if position is within cds_base_1 """
        if not self.cds_base_1:
            return None

        if position == self.cds_base_1.location.start:
            return "exonLeft"
        elif position == self.cds_base_1.location.end:
            return "exonRight"
        elif position in self.cds_base_1:
            return "withinElem"
        else:
            return None

    """
    Functions: calculate reading frame
    """
    def calc_reading_frame( self, position ):
        """ Function: if position within exon, then check the reading frame at position of interest """
        #if reading frame is -1, that means there is no reading frame
        if self.readFrame == -1:
            return -1

        ##TEST:: print "Exon.calc_reading_frame(): ", position, " & in cds = ", self.in_cds_base_1( position )

        #make sure position is within exon
        # if not position in self.cds_base_1.location:
        if not self.in_cds_base_1( position ):
            return None

        #calculate the distance of the position from the 5' end of exon
        if self.exonPos.location.strand == -1:      #minus genes, 5' end on higher numerical position of exon (right side of exon)
            dist = self.cds_base_1.location.end - position
        else:       #plus genes: 5' end on lower numerical position of exon (left side of exon)
            dist = position - self.cds_base_1.location.start

        return (dist + self.readFrame) % 3

    ##DELETE? I think this is obsolete now
    def calc_reading_frame_backup( self, position ):
        """ Function: if position within exon, then check the reading frame at position of interest """
        #if reading frame is -1, that means there is no reading frame
        if self.readFrame == -1:
            return -1

        #make sure position is within exon
        # if not position in self.cdsPos.location:
        if not self.in_cds( position ):
            return None

        #calculate the distance of the position from the 5' end of exon
        if self.exonPos.location.strand == -1:      #minus genes, 5' end on higher numerical position of exon (right side of exon)
            dist = self.cdsPos.location.end - position
        else:       #plus genes: 5' end on lower numerical position of exon (left side of exon)
            dist = position - self.cdsPos.location.start

        return (dist + self.readFrame) % 3

    """
    Function: calculate exon expression
    """
    @staticmethod
    def quant_rpkm( bam_reader, genomic_range, unique_reads = False ):
        """
        Args:
            bam_reader = HTSeq.BAM_Reader instance, i.e. HTSeq.BAM_Reader( bam_path )
            genomic_range = string with the format "chr_num:start-end"
            unique_reads = boolean
                -True = only consider uniquely mapped reads that map to genomic_range (i.e. a.optional_field( "NH" ) == 1)
                -False = consider all reads that map to genomic_range (unique + multimapped)
        Function: quantify the number of single reads mapping to region "genomic region"
        """
        #determine if only unique reads should be considered (unique_reads = True) or all reads (unique_reads = False)
        if unique_reads:
            unique_count = 0
            for i, a in enumerate( bam_reader.fetch( region = genomic_range ) ):
                if a.optional_field( "NH" ) == 1:
                    unique_count += 1
            return unique_count
        else:
            return sum( 1 for x in bam_reader.fetch( region = genomic_range ) )

    @staticmethod
    def quant_fpkm( bam_reader, genomic_range, unique_reads = False ):
        """
        Args:
            bam_reader = HTSeq.BAM_Reader instance, i.e. HTSeq.BAM_Reader( bam_path )
            genomic_range = string with the format "chr_num:start-end"
            unique_reads = boolean
                -True = only consider uniquely mapped reads that map to genomic_range (i.e. a.optional_field( "NH" ) == 1)
                -False = consider all reads that map to genomic_range (unique + multimapped)
        Function: same as def quant_genes_rpkm(), but uses pair-end reads (FPKM) instead of single-end reads (RPKM)
        """
        #retrieve the number of reads that mape to genomic region
        if unique_reads:
            unique_count = 0
            pair_reader = HTSeq.pair_SAM_alignments( bam_reader.fetch( region = genomic_range ) )
            for i, (a, b) in enumerate( pair_reader ): 
                if a and b:     #make sure a & b are not NONE
                    if a.optional_field( "NH" ) == 1 and b.optional_field( "NH" ) == 1:
                        unique_count += 1
            
            return unique_count
        else:
            pair_reader = HTSeq.pair_SAM_alignments( bam_reader.fetch( region = genomic_range ) )
            return sum( 1 for a,b in pair_reader if a and b )       #make sure both reads in pair are not NONE

    def quant_mapped_reads( self, bam_reader, bool_rpkm = True, unique_reads = False ):
        """
        Args:
            bam_reader = HTSeq.BAM_Reader instance, i.e. HTSeq.BAM_Reader( bam_path )
            genomic_range = string with the format "chr_num:start-end"
            bool_rpkm = boolean
                -True = quantify the # of single-end reads that map to exon
                -False = quantify the # of pair-end reads that map to exon
            unique_reads = boolean
                -True = only consider uniquely mapped reads that map to genomic_range (i.e. a.optional_field( "NH" ) == 1)
                -False = consider all reads that map to genomic_range (unique + multimapped)
        Function: quantify the number of reads that maps 
        """
        genomic_range = self.str_genomic_pos( True )
        if bool_rpkm:
            return Exon.quant_rpkm( bam_reader, genomic_range, unique_reads )
        else:
            return Exon.quant_fpkm( bam_reader, genomic_range, unique_reads )

    def calc_read_density( self, bam_reader, library_size, bool_rpkm = True, unique_reads = False ):
        """
        Args:
            -bam_reader = HTSeq.BAM_Reader instance, i.e. HTSeq.BAM_Reader( bam_path )
            -count = number of reads mapped to region, calculated by quant_exons()
            -genomic_range = string in the format "chrom:start-end"
            -library_size = integer that is the total number of reads in library. This is calculated by Isoform.total_mapped_reads( bam_path )
            -bool_rpkm = boolean
                -True = quantify the # of single-end reads that map to exon
                -False = quantify the # of pair-end reads that map to exon
            -unique_reads = boolean
                -True = only consider uniquely mapped reads that map to genomic_range (i.e. a.optional_field( "NH" ) == 1)
                -False = consider all reads that map to genomic_range (unique + multimapped)
        Function: calculates RPKM or FPKM, depending on how parameter 'count' was generated (def quant_genes_rpkm() or def quant_genes_fpkm())
        """
        #calculate RPKM = 10^9 * num_of_reads / (region_length * library_size )
        count_reads = self.quant_mapped_reads( bam_reader, bool_rpkm, unique_reads )
        #split genomic range 
        exon_len = self.exonPos.location.end - self.exonPos.location.start
        return ( 10**9 * float( count_reads ) ) / ( library_size * exon_len )

