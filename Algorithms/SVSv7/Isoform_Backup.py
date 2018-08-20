#/usr/bin/python
import tabix

from Exon import Exon

exon_base = 1           #if 0, then exons are 0-based, if 1, then exons are 1-based

class Isoform( object ):
    obj_cruzdb = None
    #make exon range string ( format: chrom:start-end )
    make_key_range = staticmethod( lambda chrom, exon_start, exon_end : str(chrom) + ':' + str(exon_start) + '-' + str(exon_end) )

    def __init__( self, isoform_id ):
        # isoform_info = Isoform.get_isoform( isoform_id )
        info_isoform = Isoform.obj_cruzdb.refGene.filter_by( name = isoform_id ).first()

        self.isoform_id = str( info_isoform.name )
        self.gene_sym = str( info_isoform.name2 )
        self.chrom = str( info_isoform.chrom )

        convert_ss = { '+' : 1, '-' : -1 }
        self.strand = convert_ss.get( info_isoform.strand, 0 )
        self.boundary = ( info_isoform.txStart, info_isoform.txEnd )       #tuple that records the gene boundary

        # self.hashExonList = self.organize_exons( info_isoform )     #key = string that is exon range (chrom:start-end), value = Exon object

        self.hashExonList = self.organize_features( info_isoform, True )     #key = string that is exon range (chrom:start-end), value = Exon object
        self.hashIntronList = self.organize_features( info_isoform, False )     #key = string that is exon range (chrom:start-end), value = Exon object

        # self.last_exon_num = len( self.hashExonList.keys() ) - 1      #for 0-based exons, -1 because exon starts at 0, not 1
        self.last_exon_num = len( self.hashExonList.keys() )     #-1 because exon starts at 0, not 1

        self.donor_acceptor_sites = self.get_donor_acceptor_sites()

    
    def __str__( self ):
        pos = self.chrom + ':' + str( self.boundary[0] ) + '-' + str( self.boundary[1] )
        strand_sign = '+' if self.strand > 0 else '-'
        return self.isoform_id + " (" + self.gene_sym + "): " + pos + ", (" + strand_sign + ")"

    def display_exon_list( self ):
        """
        Function: displays the exons in the isoform, returns string that gives exon number & the genomic range (e.g. chrom:start-end)
        """
        str_exons = ""
        for k,v in self.hashExonList.iteritems():       #key = string that is exon range (chrom:start-end), value = Exon object
            str_exons += "exon " + str( v.exonNum ) + ": " + k + "\n"

        return str_exons

    def json_isoform_info( self ):
        """
        Function: returns a hash that contains information about the isoform - gene name, isoform id, genomic position, exons & introns, etc. 
        """
        info = {}
        info['isoform_id'] = self.isoform_id
        info['gene_sym'] = self.gene_sym
        info['chrom'] = self.chrom
        info['pos_start'] = self.boundary[0]
        info['pos_end'] = self.boundary[1]
        info['strand_sign'] = self.strand
        info['last_exon_num'] = self.last_exon_num

        info['features'] = {}

        #get all the exons in the isoform
        for i in range(exon_base, self.last_exon_num):
            obj_feat_num = self.get_exon_num( i )
            obj_feat = {
            'type': 'exon',
            'chrom': obj_feat_num.chrom,
            'start': obj_feat_num.exonPos.location.start,
            'end': obj_feat_num.exonPos.location.end,
            'strand': obj_feat_num.exonPos.strand,
            'featureNum': obj_feat_num.exonNum
            }

            #add object to hash of features
            feat_name = 'exon' + str( i )
            info['features'][feat_name] = obj_feat

        #get all the introns in the isoform
        for i in range(exon_base, self.last_exon_num - 1):
            obj_feat_num = self.get_intron_num( i )
            obj_feat = {
            'type': 'intron',
            'chrom': obj_feat_num.chrom,
            'start': obj_feat_num.exonPos.location.start,
            'end': obj_feat_num.exonPos.location.end,
            'strand': obj_feat_num.exonPos.strand,
            'featureNum': obj_feat_num.exonNum
            }

            #add object to hash of features
            feat_name = 'intron' + str( i )
            info['features'][feat_name] = obj_feat

        return info


    """
    Function: organizing genomic information
    """
    def organize_features( self, isoform, bool_exon = True ):
        """
        Args:
            isoform = cruzdb object that contains information for a specific isoform
        Function: this will organize the feats, cds, & reading frame for a specific gene
        """
        #hash_feats = will contain feats sorted by start position (key = index, value = feat)
        hash_feats = {}       #key = string that is feat range (chrom:start-end), value = feat object

        #get reading frames
        

        ##TEST:: print "organize_feats: isoform = ", isoform.feats
        if bool_exon:
            list_features = isoform.exons
            feat_frames = map( int, [x for x in isoform.exonFrames.split(',') if x] )
        else:
            list_features = isoform.introns

        for i, feat in enumerate( list_features ):      #i = feat number
            #check to see if the key_range exists 
            key_range = Isoform.make_key_range( self.chrom, feat[0], feat[1] )

            #calculate the feat number
            # feat_num = i if self.strand == 1 else ( len(list_features) - i - 1 )        #for 0-based exons, use this
            # feat_num = i + 1 if self.strand == 1 else ( len(list_features) - i )        #for 1-based exons, use this
            feat_num = i + exon_base if self.strand == 1 else ( len(list_features) - i - ( exon_base - 1 ) )        #this handles both 0-based & 1-based

            #NOTE: UCSC has 0-based genome, meaning the first position of the feat is actually the last position in the previous intron, that is why I add '+1'
            if bool_exon:
                feat_info = self.get_feat_info( feat[0], feat[1], feat_frames[i], feat_num, True, 'exon' )
            else:
                feat_info = self.get_feat_info( feat[0], feat[1], None, feat_num, True, 'intron' )
            hash_feats[ key_range ] = Exon( feat_info )

        #go through all cds, assign CDS to each feat
        if bool_exon:
            for each_cds in isoform.cds:     #each_cd = tuple where [0] = start position & [1] = end position
                hash_feats = self.organize_exons_cds( each_cds[0], each_cds[1], hash_feats )

        return hash_feats


    
    ##MAY DELETE THIS BECAUSE NOW I HAVE - organize_features
    # def organize_exons( self, isoform ):
    #     """
    #     Args:
    #         isoform = cruzdb object that contains information for a specific isoform
    #     Function: this will organize the exons, cds, & reading frame for a specific gene
    #     """
    #     #hashExonList = will contain exons sorted by start position (key = index, value = Exon)
    #     hashExonList = {}       #key = string that is exon range (chrom:start-end), value = Exon object

    #     #get reading frames
    #     exon_frames = map( int, [x for x in isoform.exonFrames.split(',') if x] )

    #     ##TEST:: print "organize_exons: isoform = ", isoform.exons

    #     for i, exon in enumerate( isoform.exons ):      #i = exon number
    #         #check to see if the key_range exists 
    #         key_range = Isoform.make_key_range( self.chrom, exon[0], exon[1] )

    #         #calculate the exon number
    #         exon_num = i if self.strand == 1 else ( len(isoform.exons) - i - 1 )

    #         #NOTE: UCSC has 0-based genome, meaning the first position of the exon is actually the last position in the previous intron, that is why I add '+1'
    #         exon_info = self.get_feat_info( exon[0], exon[1], exon_frames[i], exon_num, True )
    #         hashExonList[ key_range ] = Exon( exon_info )

    #     #go through all cds, assign CDS to each exon
    #     for each_cds in isoform.cds:     #each_cd = tuple where [0] = start position & [1] = end position
    #         hashExonList = self.organize_exons_cds( each_cds[0], each_cds[1], hashExonList )

    #     return hashExonList


    def organize_exons_cds( self, start_cds, end_cds, hashExonList ):
        """
        Args:
            start_cds & end_cds = the starting (lower numerical position) & ending CDS positions (higher numerical position) that is fused by the splice junction.
            tuple_cds = tuple where [0] = start position of cds & [1] = 
            hashExonList = hash that records the exon for the gene, where key = string that is exon range (chrom:start-end), value = Exon object
        Function: maps the cds exons to the original exons in hashExonList
        """
        for k,v in hashExonList.iteritems():        #key = exon range (chrom:start-end), value = Exon object
            v.add_cds( start_cds, end_cds )

        return hashExonList

    """
    Function: get exons based on exon number
    """
    def get_exon_num( self, exon_num ):
        """ 
        Args:
            exon_num = integer that is the 
        Function: returns Exon object with exon number
        """
        for k,v in self.hashExonList.iteritems():       #k = string (chrom:start-end), v = Exon object
            if v.exonNum == exon_num:
                return v

        #if no exon found with that exon number, then return None
        return None

    def get_intron_num( self, intron_num ):
        """ 
        Args:
            intron_num = integer that is the 
        Function: returns intron object with intron number
        """
        for k,v in self.hashIntronList.iteritems():       #k = string (chrom:start-end), v = Exon object
            if v.exonNum == intron_num:
                return v

        #if no intron found with that intron number, then return None
        return None

    """
    Function Start: finding CDS that contain position
    """
    def get_element( self, position, bool_exon = True, direction = 0, bool_strand = True, rel_pos = None ):
        """
        Args:
            position = integer that is the genomic position of interest
            bool_exon = boolean
                -True = retrieve exon features
                -False = retrieve cds features
            direction = integer that is the direction to find the next exon
                - 0 = direction does not matter
                - 1 = find next exon to the right (highest position - higher exon number)
                - -1 = find next exon to the left (lowest position - lower exon number)
            bool_strand = boolean (Note: direction = 0 will not pick a specific )
                -True = will consider strand sign to find exon
                    -if direction = 1, then will find the higher exon number (higher position for + & lower position for -)
                    -if direction = -1, then will find lower exon number (lower position for + & higher position for -)
                -False = will consider finding the higher-positioned exon regardless of exon number & strand sign
                    -if direction = 1, then will just find exon with highest position
                    -if direction 
            rel_pos = string that is relative position of the element. Can be the following values:
                -None = position can be anywhere within exon
                -"exonLeft" = left side of exon 
                -"exonRight" = right side of exon
        Function: retrieve exons from all isoforms, returns the first Exon object that contains integer 'position' 
        Assumption: assumes that exons are not overlapping (therefore only 1 exon, if any, will contain position)
        """
        elems = self.get_exon( position, rel_pos ) if bool_exon else self.get_cds( position, rel_pos )

        #if considering gene strand, then configure direction to take into account the strand sign
        if bool_strand:
            direction *= self.strand

        #determine which exon to return based on if exons are found & the direction
        if not elems:
            elem = None
        elif direction > 0:      #find the elem with the highest numerical position
            elem = max( elems, key = lambda x: x.exonPos.location.start )
        elif direction < 0:     #find the elem with the lowest numerical position
            elem = min( elems, key = lambda x: x.exonPos.location.end )
        else:
            elem = elems[0]

        ##TEST:: print "Isoform.get_element: gene = ", self.gene_sym," & position = ", position, " & elems = ", elems

        return elem

    def get_exon( self, position, rel_pos = None ):
        """
        Args:
            position = integer that is the genomic position of interest
            rel_pos = string that is relative position of the element. Can be the following values:
                -None = position can be anywhere within exon
                -"exonLeft" = left side of exon 
                -"exonRight" = right side of exon
        Function: retrieve exons from all isoforms, returns the first Exon object that contains integer 'position' 
        Assumption: assumes that exons are not overlapping (therefore only 1 exon, if any, will contain position)
        """
        if not rel_pos:
            exons = [self.hashExonList[k] for k in self.hashExonList if self.hashExonList[k].in_exon( position )]
        else:
            exons = [self.hashExonList[k] for k in self.hashExonList if self.hashExonList[k].in_exon( position ) == rel_pos]

        #Assumption: that the exons are not overlapping, therefore only 1 exon will be retrieved
        # if len( exon ) > 1:
        #     print "!!Isoform.get_exon Warning: More than 1 exon found. Returning only first occurrence!!"

        return exons

    ##MAY NEED TO DELETE
    # def get_exon_backup( self, position, rel_pos = None ):
    #     """
    #     Args:
    #         position = integer that is the genomic position of interest
    #         rel_pos = string that is relative position of the element. Can be the following values:
    #             -None = position can be anywhere within exon
    #             -"exonLeft" = left side of exon 
    #             -"exonRight" = right side of exon
    #     Function: retrieve exons from all isoforms, returns the first Exon object that contains integer 'position' 
    #     Assumption: assumes that exons are not overlapping (therefore only 1 exon, if any, will contain position)
    #     """
    #     if not rel_pos:
    #         exon = [self.hashExonList[k] for k in self.hashExonList if self.hashExonList[k].in_exon( position )]
    #     else:
    #         exon = [self.hashExonList[k] for k in self.hashExonList if self.hashExonList[k].in_exon( position ) == rel_pos]

    #     #Assumption: that the exons are not overlapping, therefore only 1 exon will be retrieved
    #     # if len( exon ) > 1:
    #     #     print "!!Isoform.get_exon Warning: More than 1 exon found. Returning only first occurrence!!"

    #     return exon[0] if exon else None

    def in_exon( self, position ):
        """ Function: returns the relative position within the exon ('exonLeft', 'exonRight', 'withinExon') """
        rel_pos = None      #records the position relative in the exon ('exonLeft', 'exonRight', 'withinExon')
        for k,v in self.hashExonList.iteritems():       #k = string (chrom:start-end), v = Exon object
            #if position is found, then return relative position
            rel_pos = v.in_exon( position )
            if rel_pos:
                break

        return rel_pos

    def get_cds( self, position, rel_pos = None ):
        """
        Args:
            position = integer that is the genomic position of interest
            rel_pos = string that is relative position of the element. Can be the following values:
                -None = position can be anywhere within exon
                -"exonLeft" = left side of exon 
                -"exonRight" = right side of exon
        Function: retrieve exons from all isoforms, returns an array of Exon objects that contain integer 'position'
        Assumption: assumes that exons are not overlapping (therefore only 1 exon, if any, will contain position)
        """

        #NOTE: cdss is just the plural of CDS (Coding Sequence)
        if not rel_pos:
            cdss = [self.hashExonList[k] for k in self.hashExonList if self.hashExonList[k].in_cds( position )]
        else:
            cdss = [self.hashExonList[k] for k in self.hashExonList if self.hashExonList[k].in_cds( position ) == rel_pos]
        
        # #Assumption: that the exons are not overlapping, therefore only 1 exon will be retrieved
        # if len( cds ) > 1:
        #     print "!!Isoform.get_cds Warning: More than 1 cds found. Returning only first occurrence!!"

        return cdss


    #MAYBE DELETE THIS AS THIS DOESN'T USE "self.hashIntronList"
    # def get_intron( self, position ):
    #     """ Function: determines if position is within intron """
    #     #as the donor-acceptor sites span the intron, use it to see if position is within introns
    #     da_sites = self.get_donor_acceptor_sites()
    #     intron_num = None
    #     for i, intron in enumerate( da_sites ):
    #         if intron[0] < position <= intron[1]:
    #             intron_num = i
    #             break

    #     return intron_num

    def get_intron( self, position, rel_pos = None ):
        #NOTE: cdss is just the plural of CDS (Coding Sequence)
        if not rel_pos:
            obj_intron = [self.hashIntronList[k] for k in self.hashIntronList if self.hashIntronList[k].in_exon( position )]
        else:
            obj_intron = [self.hashIntronList[k] for k in self.hashIntronList if self.hashIntronList[k].in_exon( position ) == rel_pos]


        ##TEST::
        print "Isoform.get_intron: gene = ", self.gene_sym," & position = ", position, " & obj_intron = ", obj_intron
        
        # #Assumption: that the exons are not overlapping, therefore only 1 exon will be retrieved
        # if len( cds ) > 1:
        #     print "!!Isoform.get_cds Warning: More than 1 cds found. Returning only first occurrence!!"

        return obj_intron[0] if obj_intron else None

    def in_intron( self, position ):
        """ Function: similar to def in_exon(), but does this for introns. returns the relative position within the exon ('exonLeft', 'exonRight', 'withinExon') """
        rel_pos = None      #records the position relative in the exon ('exonLeft', 'exonRight', 'withinExon')
        for k,v in self.hashIntronList.iteritems():       #k = string (chrom:start-end), v = Exon object
            #if position is found, then return relative position
            rel_pos = v.in_exon( position )
            if rel_pos:
                break

        return rel_pos

    """
    Function Start: calculate reading frame
    """
    def calc_reading_frame( self, position, direction ):
        """
        Args:
            position = integer that is the genomic position of interest
            direction = integer that is the direction to find the next exon
                - 0 = direction does not matter
                - 1 = find next exon to the right (highest position - higher exon number)
                - -1 = find next exon to the left (lowest position - lower exon number)
        Function: return reading frame for isoform. Returns integer 0, 1, or 2. If position not in any exons, then return None, else if no reading frame for exon then return -1
        """
        #parameters for get_element( self, position, bool_exon = True, direction = 0, bool_strand = True, rel_pos = None ) --> don't consider strand sign & retrieves CDS, not exon
        bool_exon = False
        bool_strand = False
        rel_pos = None
        cds = self.get_element( position, bool_exon, direction, bool_strand, rel_pos )    #retrieves Exon object

        if cds:
            return cds.calc_reading_frame( position )
        else:
            return None

    ##DELETE?? - 
    # def calc_reading_frame( self, position, isoform ):
    #     """ Function: return reading frame for isoform. Returns integer 0, 1, or 2 """
    #     #get exons that contain position
    #     cds_exon = self.get_cds( position, isoform )       #k = exon number, v = Exon Model that contains position & is from specific isoform

    #     #if no exons contain position, then return
    #     if not cds_exon:
    #         return None

    #     #if -1, this means there is no reading frame
    #     if cds_exon.read_frame == -1:
    #         return -1

    #     zero_base = position - cds_exon.cdsPos.location.start
    #     if self.strand == -1:      #minus genes, 5' end on higher numerical position of exon (right side of exon)
    #         dist = cds_exon.cdsPos.location.end - position
    #         if zero_base == 0:       #similar to dist = exons[cds_num].exonPos.location.end - ( position + 1 )
    #             dist -= 1       #CORRECTION (ONLY A TEMPORARY FIX): Since UCSC (and as a result cruzdb) is 0-based, the first nucleotide position for each exon is actually the last position nucleotide of the previous intron, so need to add +1 in order to be at the correct position

    #         ##TEST:: print "CRF minus gene --> dist = ", dist, " & position = ", position, " & end_pos = ", exons[cds_num].exonPos.location.end, " & cds_num = ", cds_num
    #     else:       #plus genes: 5' end on lower numerical position of exon (left side of exon)
    #         dist = position - cds_exon.cdsPos.location.start
    #         if zero_base > 0:       #similar to dist = position - ( exons[cds_num].exonPos.location.start + 1 )
    #             dist -= 1       #CORRECTION (ONLY A TEMPORARY FIX): Since UCSC (and as a result cruzdb) is 0-based, the first nucleotide position for each exon is actually the last position nucleotide of the previous intron, so need to add +1 in order to be at the correct position

    #     return (dist + cds_exon.read_frame) % 3

    """
    Function Start: determine if reading frame is preserved
    """
    ##MAY NEED TO DELETE THIS FUNCTION
    # def sj_position_reading_frame( self, position, bool_start = True ):
    #     """
    #     Args:
    #         position = integer that is the genomic position of interest - will be used to find the reading frame at that position
    #         bool_start = is this the starting or ending position of a splice junction 
    #     Helpful when: this function is helping when considering the position of a splice junction
    #     """

    #     if bool_start:      #if this is the starting position
    #         if self.strand > 0:     #if plus strand
    #             start_direction = -1        #retrieve the lowest-positioned cds (left-most cds)
    #         else:                   #if minus strand
    #             start_direction = 1     #retrieve the highest-positioned cds (right-most cds)
    #     else:           #else if this is the ending position
    #         if self.strand > 0:     #if plus strand
    #             end_direction = 1           #retrieve the highest-positioned cds (right-most cds)
    #         else:                   #if minus strand
    #             end_direction = -1      #retrieve the lowest-positioned cds (left-most cds)

    #     return self.calc_reading_frame( position, direction )


    def frame_preserved( self, pos_a, pos_b ):
        """ 
        Args:
            pos_a & pos_b = integer the is the genomic position of interest, where pos_a < pos_b
            isoform_id: string that isoform id (RefSeq ID)
        Function: determines if the reading frame is preserved. Returns True if reading frame is preserved, else returns false
        """
        #need to take strand sign into account (if plus then reading frame goes left to right, else if minus then goes right to left)
        if self.strand > 0:     #if plus strand
            start = pos_a
            end = pos_b
            start_direction = -1        #retrieve the lowest-positioned cds (left-most cds)
            end_direction = 1           #retrieve the highest-positioned cds (right-most cds)
        else:                   #if minus strand
            start = pos_b
            end = pos_a
            start_direction = 1     #retrieve the highest-positioned cds (right-most cds)
            end_direction = -1      #retrieve the lowest-positioned cds (left-most cds)

        #calculate the reading frame for position A & B
        frame_a = self.calc_reading_frame( start, start_direction )
        frame_b = self.calc_reading_frame( end, end_direction )

        #check if frame found for both positions
        if frame_a == None or frame_b == None:
            return { 'rf_start': frame_a, 'rf_end': frame_b, 'rf_preserved': None }

        #calculate next frame
        next_frame = ( frame_a + 1 ) % 3
        if next_frame == frame_b:
            stat = True
        else:
            stat = False

        return { 'rf_start': frame_a, 'rf_end': frame_b, 'rf_preserved': stat }

    """
    Function: check if exons skipped in isoform
    """

    def get_donor_acceptor_sites( self ):
        """ Function: returns a tuple of exon ends 
        Notes: for each tuple: for + genes, [0] = donor site (3' end) & [1] = acceptor site (5' end). For - gene, [0] = acceptor site (5' end) & [1] = donor site (5' end)
        """
        sj_sites = []       #array of tuples, where for each tuple contains ends of exon
        exons = sorted( self.hashExonList.values(), key = lambda k: k.exonPos.location.start, reverse = False )
        for i in range( 0 , len( exons ) - 1 ):
            sj_sites.append( ( exons[i].exonPos.location.end, exons[i + 1].exonPos.location.start )  )

        return sj_sites

    def is_exon_skip( self, start, end ):
        """
        Args:
            start & end = genomic positions of start & end position of SJ, where start < end
        Function: returns an array of integers, where each integer represents the exons skipped. Blank array means no exons were skipped 
        NOTE: the exon numbers are zero-based, so exon 0 = 1st exon, exon 1 = 2nd exon
        """
        #retrieve Exon objects that contain position
        start_num = self.get_element( start, True, -1, False, None )    #retrieve the left-most exon (exon with lowest numerical position)
        end_num = self.get_element( end, True, 1, False, None )     #retrieve the right-most exon (exon with highest numerical position)
        if start_num == None or end_num == None:
            return None

        #determine the lower & higher exon number
        skip_range = range( start_num.exonNum + 1, end_num.exonNum ) if start_num.exonNum < end_num.exonNum else range( end_num.exonNum + 1, start_num.exonNum )

        return [i for i in skip_range]


    """
    Function: determine the effect of a splice junction in the isoform
    """
    def is_canon_sj( self, start, end ):
        """ Function: checks if positions of splice sites land in canonical donor-acceptor sites - returns True if position exists array of donor-acceptor sites, else returns False """
        #check start & end position of splice junction to see if it is canonical
        sj_pos = (start, end)
        if sj_pos in self.donor_acceptor_sites:
            return True
        else:
            return False

    def is_aberrant_sj( self, start, end, intronic = False ):
        """
        Args:
            start & end = integers that are the start & end genomic position 
            intronic = boolean
                -True = return an isoform ID even if it does land in the intronic space
                -False = do not return isoform ID if only lands in intronic space
        Function: retrieve information about start & end splicing to see if how the splice junction is aberrant, including whether splicing is in a exon (donor site, acceptor site, in exon) or intronic, any exons skipped, and any frameshifts """
        #see where position of start & end positions is located
        exon_start = self.in_exon( start )
        exon_start = 'exon_canon' if exon_start == 'exonRight' else exon_start      #'exon_canon' means splicing occurring at canonical position
        exon_end = self.in_exon( end )
        exon_end = 'exon_canon' if exon_end == 'exonLeft' else exon_end         #'exon_canon' means splicing occurring at canonical position

        #if intronic == False, then check if SJ is intronic - if so then discard SJ
        if not intronic:
            #if no exon (or intron) is returned, return None
            if exon_start == None or exon_end == None:
                return None
        
        exon_skips = self.is_exon_skip( start, end )        #check for exons skipped
        check_frame = self.frame_preserved( start, end )    #check to see if reading frame is preserved

        #exon_skips = array of integers referring to exons skipped
        #frame_preserved = True (SJ preserves frame), False (SJ does not preserve frame), or None (no frame, probably b/c SJ is intronic)
        #intronic (when I do implement it) = True (means SJ is intronic), False = (SJ is not intronic)
        calc_rf = ( check_frame['rf_start'], check_frame['rf_end'] )
        aberrations = { 'exon_start': exon_start, 'exon_end': exon_end, 'exon_skip': exon_skips, 'frame_preserved': check_frame['rf_preserved'], 'calc_rf': calc_rf }
        aberrations.update( self.score_aberrations( aberrations ) )

        return aberrations

    def score_aberrations( self, aberrations ):
        """ 
        Args:
            aberrations = hash of all aberrations, keys are 'exon_start', 'exon_end', 'exon_skips', 'frame_preserved'
        Function: scores the aberrations associated with the abnormalities associated with splice junction
        Note on scoring: 0 = no aberrations, 1 = exon skip, 100 = exonic (not on exon ends), 200 = frameshift, 500 = intronic
        """
        score = 0
        effects = []        #record the effects
        #quantify number of exon skips - make sure not aberrations['exon_skip'] is not 'NoneType'
        score += len( aberrations['exon_skip'] ) if aberrations['exon_skip'] is list else 0
        if score > 0:
            effects.append( 'exon_skip: ' + ','.join( aberrations['exon_skip'] ) )

        #check if splice junction is intronic (one of the positions lands in the intron) or exonic (within exon body, not at ends of exon)
        if aberrations['exon_start'] == None or aberrations['exon_end'] == None:        #splicing in intron
            score += 500
            effects.append( 'intronic' )
        elif aberrations['exon_start'] != 'exon_canon' or aberrations['exon_end'] != 'exon_canon':      #splicing in exon (not at donor-acceptor sites)
            score += 100
            effects.append( 'exonic' )
        #check if reading frame is preserved
        if not aberrations['frame_preserved']:
            score += 200
            effects.append( 'frameshift' )

        return { 'score': score, 'effects': effects }
    
    """
    Functions: Retrieve information for gene & isoform, including start & end position 
    """

    ##DELETE? - MAY DELETE THIS LATER AS I THINK I CAN RETRIEVE THE ISOFORM JUST BY USING CRUZDB
    @classmethod
    def get_isoform( cls_obj, isoform_id ):
        """ retrieves all isoforms that have the same gene id 'isoform_id' """
        #retrieve information from cruzdb
        #NOTE: I am using "refGene.filter_by" because for some reason more isoform IDs retrieved by it than by .bin_query
        # cruzdbInfo = ExonList.objCruzDB.bin_query( "refGene", chrom, pos_start, pos_end ).all()
        info = Isoform.obj_cruzdb.refGene.filter_by( name = isoform_id ).all()

        #save all isoforms with geneSym into array
        isoform_info = None
        for x in info:
            #Method 1: If I only want to consider mRNA protein-coding isoforms
            # if x.name2 == isoform_id and "NM_" in x.name:        #if want to consider only mRNA isoforms
            #     isoform_info = x
            
            #Method 2: if I want all isoforms
            """
            IMPORTANT: the reason I decided this is because there could a splice junction could correctly splice a non-mRNA transcript, therefore reducing erroneously stating an SJ is aberrant when indeed it maps to a non-mRNA transcript. So need to consider reads that map to non-mRNA transcripts as well.
            """
            if x.name == isoform_id:      #if want to consider all isoforms (mRNA, ncRNA, etc.)
                isoform_info = x
                break

        return isoform_info

    @classmethod
    def set_cruzdb( cls_obj, obj_cruzdb ):
        """ Function: assigns cruzdb object to class variable """
        cls_obj.obj_cruzdb = obj_cruzdb

    def get_feat_info( self, pos_start, pos_end, exon_frame, exon_num, canonical, feat_type ):
        """
        Args:
            pos_start & pos_end = integer that is the lower & higher nucleotide position for the exon, respectively
            exon_frame = integer that is 0, 1, or 2. Refers to the reading frame of the 5' end of exon (lower position for plus genes & higher position for minus genes). It is -1 if the exon is non-coding & None if reading frame is unknown
            exon_num = exon number in gene. For plus genes, exon number increases from lower to higher position. For minus genes, exon number decreases from lower to higher position. If None, means unknown what the exon number is (could be non-canonical exon number)
            canonical = boolean that, if True, means this exons is canonical (known to exist) in the gene, else if False, then is not a canonical exon in the gene
            feat_type = string that conveys the type of feature the element is, either 'exon' or 'intron'
        Function: creates & returns hash to create Exon
        """
        #prepare hash that will record elements into exon list
        exon_info = {
        "pos_start": pos_start,
        "pos_end": pos_end,
        "strand_sign": self.strand,
        "chrom": self.chrom,
        "isoform_id": self.isoform_id,
        "exon_frame": exon_frame,
        "exon_num": exon_num,
        "canonical": canonical,
        "feat_type": feat_type }

        return exon_info