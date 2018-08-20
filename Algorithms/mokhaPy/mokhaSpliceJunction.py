#/usr/bin/python

import os
import sys
import subprocess
import time
import math
import re

"""
from inspect import currentframe, getframeinfo  --> use this to find the current filename & line number. 
frameinfo = getframeinfo(currentframe()) --> this retrieves information about the current line
print frameinfo.filename," ",frameinfo.lineno --> filename prints the file name & lineno prints the line number
"""
from inspect import currentframe, getframeinfo 

import Bio.Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, BeforePosition, AfterPosition
import mygene
from cruzdb import Genome, sequence         #Genome = retrieves genomic information, sequence = can retrieve nucleotide sequence based on chromosome # & position

sys.path.insert(0,"/home/mokha/Documents/Krauthammer_Lab/PythonClasses")
import mokhaPy
from mokhaGeneModel import GeneModel


"""
QUESTIONS:
-Python: for class methods, why do I need an object as the parameter (e.g. cls_obj): CONJ: becuase when the class method is called, it is preceded by the class name, therefore referencing the class that contains the class method
"""

class SpliceJunction(GeneModel):
    #class variables
    cls_sjID=0          #numerical ID for each SJ event
    objCruzDB=None      #will record the cruzdb object
    path_genomeIndex=None       #string that is the path to the samtools indexed genome
    temp_arrGeneList=None       #this will record the list of genes that program should be looking for (e.g. list of kinases)

    @classmethod
    def setCruzDB(cls_obj,cruzdb_genome):
        """
        Args:
            cruzdb_genome = a cruzDB object that is assigned to the database where genomic information is being queried (e.g. cruzdb_genome=cruzdb.Genome("hg19"), or for temporary databases - cruzdb_genome=cruzdb.Genome("sqlite:////tmp/hg19.db")
        Function:
            this function will define the class variable objCruzDB, cruzdb object that will query UCSC resources via cruzdb class
        """
        cls_obj.objCruzDB=cruzdb_genome
        # SpliceJunction.objCruzDB=cruzdb_genome
        GeneModel.objCruzDB=cruzdb_genome

    @classmethod
    def setPath_genomeIndex(cls_obj,strPath):
        """
        Function:
            this function will define the class variable path path_genomeIndex. This should be to a samtools indexed genome file (e.g. samtools faidx hg19.fa)
        """
        cls_obj.path_genomeIndex=strPath
        # SpliceJunction.path_genomeIndex=strPath
        GeneModel.strPath=strPath

    @staticmethod
    def checkSJ_assocGene(hashInfo):
        """
        Args:
            hashInfo is a hash (python dictionary) with the following properties:
                chrNum: the chromosome ID, which will be "chr#" (e.g. chr2, chr9, chr13)
                posStart:
                posEnd:
                strandSign: should be 1 (+ strand), -1 (- strand), or 0 (unknown)
                getType: the type of entity this object is, e.g. CDS (coding sequence), exon, intron, splice_site
                readCount: SHOULD I ADD??
        Function:
            if a gene is found associated with the splice junction, then it will not return "None", else it will return "None"
        """
        #VERSION 1:
        return SpliceJunction.objCruzDB.bin_query("refGene",hashInfo["chrNum"],hashInfo["posStart"],hashInfo["posEnd"]).first()

        #VERSION 2: this will only return true (if gene exists) or false 
        # geneInfo = SpliceJunction.objCruzDB.bin_query("refGene",hashInfo["chrNum"],hashInfo["posStart"],hashInfo["posEnd"]).first()

        # if geneInfo:
        #   return True
        # else:
        #   return False

    @classmethod
    def tempfunc_loadGeneList(cls_obj,arrGeneList):
        """
        Args:
            arrGeneList = array of all genes
        Function:
            this will load a list of genes that the splice variant should be looking for
        """
        cls_obj.temp_arrGeneList = arrGeneList
        # SpliceJunction.temp_arrGeneList = arrGeneList


    @classmethod
    def tempfunc_checkGene(cls_obj,geneName):
        """
        Arg:
            geneName = string that is the gene name
        Function:
            this function determines if the gene of interest is a kinase
        """
        #STEP: see if geneName is present in list of genes 
        if geneName.upper() in (eachGene.upper() for eachGene in cls_obj.temp_arrGeneList):     #this capitalizes all letters as to make search case-insensitive
        # if geneName.upper() in (eachGene.upper() for eachGene in SpliceJunction.temp_arrGeneList):        #this capitalizes all letters as to make search case-insensitive
            return True
        else:
            return False

    @classmethod    
    def tempfunc_checkBeforeCreate(cls_obj,hashInfo):
        geneInfo=cls_obj.checkSJ_assocGene(hashInfo)
        if not geneInfo:        #if no gene information is found, then return nothing
            print "Mokha ERROR: No associated RefSeq gene found for this splicing event"
            return False
        else:       #if the gene symbol is not located within the list of genes SpliceJunction.temp_arrGeneList
            print "gene Sym = ",geneInfo.name2
            if not cls_obj.tempfunc_checkGene(geneInfo.name2):
                print "Mokha HALT: Gene is not in list"
                return False
            else:       #this means there is a gene symbol & it is within the list of genes SpliceJunction.temp_arrGeneList
                return True 

    def __init__(self,hashInfo):
        """
        Args:
            hashInfo is a hash (python dictionary) with the following properties:
                chrNum: the chromosome ID, which will be "chr#" (e.g. chr2, chr9, chr13)
                posStart:
                posEnd:
                strandSign: should be 1 (+ strand), -1 (- strand), or 0 (unknown)
                getType: the type of entity this object is, e.g. CDS (coding sequence), exon, intron, splice_site
                readCount: SHOULD I ADD??
            cruzdb_genome = a cruzDB object that is assigned to the database where genomic information is being queried (e.g. cruzdb_genome=cruzdb.Genome("hg19"), or for temporary databases - cruzdb_genome=cruzdb.Genome("sqlite:////tmp/hg19.db")
        Properties of the SpliceJunction:
            posStart
            posEnd
            associated_gene
            exonSkip
            frameShift
            readCount
        """
        #STEP: check to see if there are any annotations associated with this splice junction
        if not SpliceJunction.tempfunc_checkBeforeCreate(hashInfo):
            print "Mokha - SJ does not map to a gene of interest"
            return None

        #STEP: retrieve information for gene associated with this splice junction event
        GeneModel.__init__(self,hashInfo)

        print "kinase name = ",self.geneSym

        #STEP: recored properties of this object: id, chrNum, posStart, posEnd, strandSign, gene, gene isoform, exon
        #count the number of events
        self.sjID=SpliceJunction.cls_sjID
        SpliceJunction.cls_sjID+=1      #increment counter 
        #record genomic position
        self.posInfo=SeqFeature(FeatureLocation(hashInfo["posStart"],hashInfo["posEnd"]),strand=hashInfo["strandSign"],type=hashInfo["type"])
        self.chrNum=hashInfo["chrNum"]          #should be in the form "chr#" (e.g. chr3, chr16, chrX)
        self.readCount=hashInfo["readCount"]

        #STEP: find the nearest feature. Returns an array with 2 values: exonNum:0/1, & distance from exonNum:0/1 
        self.start_nearestExon=self.aid_findNearestFeatureV3(self.posInfo.location.start,self.arrUniqExons)
        self.end_nearestExon=self.aid_findNearestFeatureV3(self.posInfo.location.end,self.arrUniqExons)

        #STEP: get feature information for starting position
        self.start_EI=self.posExonIntron(self.posInfo.location.start)
        self.end_EI=self.posExonIntron(self.posInfo.location.end)

        if(self.start_nearestExon[1]>0):
            self.modExon_start=self.modifiedExonSJ(True)
        else:
            self.modExon_start=None

        if(self.end_nearestExon[1]>0):
            self.modExon_end=self.modifiedExonSJ(True)
        else:
            self.modExon_end=None


    def posExonIntron(self,posOI):
        """
        Function:
            will return if the position is within the intron or exon of the gene of interest
        Output:
            will return a hash with 3 elements, isEI, featStr, & distFeat: 
                -isEI = is Exon Intron, is either "exon" or "intron"
                -featNum = the feature string in the format [featNum:0/1, distance from end position], where featNum = exon/intron number and 0 = on left side of feature (numerically lower nucleotide position) & 1 = on right side of feature (numerically higher nucleotide position)
                -distFeat = the numerical distance from the feature in featStr
        """

        #STEP: go through each exon & intron to see where the position is located
        hashPosExon=self.aid_posInFeat(posOI,self.arrUniqExons)
        hashPosIntron=self.aid_posInFeat(posOI,self.arrUniqIntrons)

        #STEP: if the position is found in the exon, then return it with information about which exon & how far it is from the closest end. ELSE if not found in the exon, then look into introns
        if hashPosExon["boolPosPresent"]:       #this means that posOI is present within exon
            for num,eachExon in enumerate(self.arrUniqExons):           #eachExon = each element is a tuple where [0] = start position & [1] = end position
                if eachExon[0] <= posOI <= eachExon[1]:
                    #STEP: retrieve the relative position within the exon, the "exon number" (based on variable "num") & whether in exon or not
                    #STEP: find the relative position from the ends of the exon & see if posOI is on the exon end
                    distLeft=posOI-eachExon[0]
                    distRight=posOI-eachExon[1]
                    #if posOI lands on the exon ends, the boolOnExonEnd=True, else it is False
                    if distLeft==0 or distRight==0:
                        boolOnExonEnd=True
                    else:
                        boolOnExonEnd=False
                    relPos_left=str(num)+":"+str(distLeft)
                    relPos_right=str(num)+":"+str(distRight)
                    hashEI={"isEI":"exon", "index_arrExon":num ,"relPos_left":relPos_left, "relPos_right":relPos_right, "boolOnExonEnd":boolOnExonEnd}
                    break
        elif hashPosIntron["boolPosPresent"]:       #this means that posOI is present within intron
            for num,eachIntron in enumerate(self.arrUniqIntrons):           #eachIntron = each element is a tuple where [0] = start position & [1] = end position
                if eachIntron[0] <= posOI <= eachIntron[1]:
                    #STEP: retrieve the relative position within the intron, the "intron number" (based on variable "num") & whether in intron or not
                    relPos_left=str(num)+":"+str(posOI-eachIntron[0])
                    relPos_right=str(num)+":"+str(posOI-eachIntron[1])
                    hashEI={"isEI":"intron", "index_arrIntron":num , "relPos_left":relPos_left, "relPos_right":relPos_right, "boolOnExonEnd":False}
                    break
        else:       #this means the position is outside the range of the gene
            hashEI={"isEI":None, "relPos_left":None, "relPos_right":None, "boolOnExonEnd":None}


        return hashEI


    def modifiedExonSJ(self, boolSJStart):
        """
        Args:
            arrSJ (array splice junctions) = array of splice junctions
            arrGIE (array Gene Isoform Exons) = array of all unique exons within a gene. These are all exons from all isoforms for a given gene
            rangeGeneBoundary = range that is the start & end position of the gene
        Function:
            this function will output all modified exons based on position of the splice junction. 
        Output:
            outputs an array of all modified unique exons based on the collection of splice junctions. By unique, there will be no duplicate exon 

        arrSJ: class SpliceJunctions
            -retrieve splice junctions within an range of positions
        arrGIE: class GeneModel
            -record all the exons associated with 
            -identify exons that are constitutive in all isoforms
            -identify exons that are unique to specific isoforms
            -return the range (lowest & highest position)
        NOTE:
            -SJ start position should be on 3' exon end, & SJ end positino should be on 5' exon end (regardless of strand sign of gene)
            -for "+" gene, for each exon the 5' end < 3' end (numerical nucleotide position). The opposite is true for "-" genes (5' end > 3' end)
        """

        #STEP: retrieve the position of interest
        posOI=self.posInfo.location.start if boolSJStart else self.posInfo.location.end

        #STEP: go through each exon & intron to see where the position is located
        hashPosExon=self.aid_posInFeat(posOI,self.arrUniqExons)

        print "TEST mESJ 1 - hashPosExon = ",hashPosExon
    
        #STEP: if the position is found in the exon, then return it with information about which exon & how far it is from the closest end. ELSE if not found in the exon, then look into introns
        hashExonRange_mod={}
        if hashPosExon["boolPosPresent"]:       #this means that posOI is present within exon
            ##REDUNDANT: Don't need to search through all exons because exon index found from sefl.aid_posInFeat()
            # for num,eachExon in enumerate(self.arrUniqExons):     #go through each exon 
            #   #STEP: see if it is within range, then see if it is on one of the acceptor end (acceptor/donor site). If yes, then return unmodified exon, else if within exon then return modified exon
            #   if eachExon[0] <= posOI <= eachExon[1]:
            #       hashExonRange_mod=self.aid_relExonPos_withinExon(posOI,eachExon[0],eachExon[1],boolSJStart) #record the modified exon range
            #       hashExonRange_mod["index_arrExon"]=num      #add the index where the tuple range was found in this array
            exonOI=self.arrUniqExons[hashPosExon["arrIndex"]]
            hashExonRange_mod=self.aid_relExonPos_withinExon(posOI,exonOI[0],exonOI[1],boolSJStart) #record the modified exon range
            hashExonRange_mod["index_arrExon"]=hashPosExon["arrIndex"]
        else:           #this means position not found in exon
            for num,eachIntron in enumerate(self.arrUniqIntrons):       #go through each intron
                #STEP: see if it is within range, then see if it is on one of the acceptor end (acceptor/donor site). If yes, then return unmodified exon, else if within exon then return modified exon
                if eachIntron[0] <= posOI <= eachIntron[1]:
                    # exonRange_mod=self.aid_relExonPos_withinIntron(posOI,eachIntron[0],eachIntron[1],boolSJStart) #record the modified exon range         #NOTE: i'm saving this because what if I need the ends of the intron
                    hashExonRange_mod=self.aid_relExonPos_withinIntron(posOI,boolSJStart)   #record the modified exon range
                    hashExonRange_mod["index_arrIntron"]=num        #add the index where the tuple range was found in this array

        ##TEST:: print "TEST mESJ 2 - hashExonRange_mod = ",hashExonRange_mod 

        #OUTPUT: hashExonRange_mod is a hash that contains multiple elements of information - see functions  
        return hashExonRange_mod

    def aid_relExonPos_withinExon(self,posOI,exonStart,exonEnd,boolSJStart):
        """
        Args:
            posOI = integer that is the numerical genomic position of interest
            exonStart = integer that is the lower numerical genomic position of the exon (+ gene: 5' end, - gene: 3' end)
            exonEnd = integer that is the higher numerical genomic position of the exon (+ gene: 3' end, - gene: 5' end)
            boolSJStart = boolean that, if true, means posOI is the start position of the splice junction (the lower nucleotide position)
        Function:
            this function will find where a position lies relative to one of the ends of the exon, depending on boolSJStart & if the gene is a "+" gene or "-" gene
        Output:
            this function should output a hash with 3 elements: description of exon start, numerical genomic position of exon position of interest, & the relative distance from that exon position
        NOTE:

        """
        #STEP: check if this is the start or the end of the splice junction 
        if boolSJStart:     #means this is the start of the splice junction position
            if self.posInfo.location.strand>0:      #means this is a plus strand
                descripExonBound="+ gene:exon 3' site:SJstart"      #description of the exon boundary of interest
                exonBoundary=exonStart              #this is one of the ends of the new exon boundary
                relPos=posOI-exonStart          #this is the relative position with respect to the appropriate exon boundary
                tupleRangeOI=(exonStart,posOI)      #this is the new modified exon boundary
            else:           #means this is a minus strand
                descripExonBound="- gene:exon 5' site:SJstart"      #description of the exon boundary of interest
                exonBoundary=exonEnd                #this is one of the ends of the new exon boundary
                relPos=posOI-exonEnd            #this is the relative position with respect to the appropriate exon boundary
                tupleRangeOI=(posOI,exonEnd)        #this is the new modified exon boundary

                ##TEST:
                print "SJ Start Minus Gene - posOI =",posOI," && exonEnd =",exonEnd
        else:       #means this is the end of the splice junction position
            if self.posInfo.location.strand>0:      #means this is a plus strand
                descripExonBound="+ gene:exon 5' site:SJend"        #description of the exon boundary of interest
                exonBoundary=exonEnd                #this is one of the ends of the new exon boundary
                relPos=posOI-exonEnd            #this is the relative position with respect to the appropriate exon boundary
                tupleRangeOI=(posOI,exonEnd)        #this is the new modified exon boundary
            else:           #means this is a minus strand
                descripExonBound="- gene:exon 3' site:SJend"        #description of the exon boundary of interest
                exonBoundary=exonStart              #this is one of the ends of the new exon boundary
                relPos=posOI-exonStart          #this is the relative position with respect to the appropriate exon boundary
                tupleRangeOI=(exonStart,posOI)      #this is the new modified exon boundary

                ##TEST:
                print "SJ End Minus Gene - posOI =",posOI," && exonStart =",exonStart

        # return {"relPos":relPos, "tupleRangeOI":tupleRangeOI}
        return {"relExonPos":relPos, "tupleRangeOI":tupleRangeOI}


    def aid_relExonPos_withinIntron(self,posOI,boolSJStart):
        """
        Args:
            posOI = integer that is the numerical genomic position of interest
            exonStart = integer that is the lower numerical genomic position of the exon (+ gene: 5' end, - gene: 3' end)
            exonEnd = integer that is the higher numerical genomic position of the exon (+ gene: 3' end, - gene: 5' end)
            boolSJStart = boolean that, if true, means posOI is the start position of the splice junction (the lower nucleotide position)
        Function:
            this function will find where a position lies relative to one of the ends of the exon, depending on boolSJStart & if the gene is a "+" gene or "-" gene
        Output:
            this function should output a hash with 3 elements: description of exon start, numerical genomic position of exon position of interest, & the relative distance from that exon position
        NOTE:

        """
        #STEP: determine which range should be selected based on scoring scheme below
        #(1) = left nucleotide position of left exon to posOI
        #(2) = posOI to right nucleotide position of right exon
        #if boolSJStart true & + gene: (1)
        #if boolSJStart false & - gene: (1)
        #if boolSJStart false & + gene: (2)
        #if boolSJStart true & - gene: (2)
        #scoreRC = score Range Category. If abs(scoreRC)==2, then its option (1), else if abs(scoreRC)==0, then option (2)
        scoreRC=0
        scoreRC+=1 if boolSJStart else -1
        scoreRC+=1 if self.posInfo.location.strand>0 else -1

        #STEP: retrieve the range based on the value of "scoreRC"
        #VAR: expandExonRange will record the new modified exon range due to splicing within the intron
        if abs(scoreRC)==2:
            nearestExon=self.aid_findNearestFeatureV2(posOI,self.arrUniqExons,False)        #retrieve the exon before the splice position
            expandExonRange=(nearestExon[0],posOI)
            relPos=posOI-nearestExon[1]     #the position from SJ position to the right end of the exon
        else:   #else if abs(scoreRC)==0
            nearestExon=self.aid_findNearestFeatureV2(posOI,self.arrUniqExons,True)     #retrieve the exon after the splice position
            expandExonRange=(posOI,nearestExon[1])
            relPos=posOI-nearestExon[0]     #the position from SJ position to the left end of the exon

        # return {"relPos":relPos, "expandExonRange":expandExonRange}
        return {"relIntronPos":relPos, "expandExonRange":expandExonRange}


    ##AID FUNCTIONS##
    def aid_posInFeat(self,posOI,arrTupleRange):
        """
        Args:
            posOI (position of interest) = this integer that will be check if it is present
            arrTupleRange = array of features where each element in arrTupleRange is a tuple with 2 elements, [0] = start position & [1] = end position. arrTupleRange[0] < arrTupleRange[1]
        Function:
            this function will return True if posOI is within 
        """
        boolPosPresent=False            #boolPosPresent = False if posOI is not present in array of ranges
        arrIndex=None
        for index,eachFeat in enumerate(arrTupleRange):
            if eachFeat[0] <= posOI <= eachFeat[1]:
                boolPosPresent=True
                arrIndex=index
                break

        return {"boolPosPresent":boolPosPresent, "arrIndex":arrIndex}

    def aid_findNearestFeatureV2(self,posOI,arrTupleRange,boolPosGreater):
        """
        Function:
            this function will find the nearest genomic feature (i.e. exon) to the splice junction within a specific direction
        """
        #STEP: add all start feature positions into an array, & insert posOI into that arrays
        arrPosExon=[x[0] for x in arrTupleRange]        #retrieve the first position (starting p)
        arrPosExon.append(posOI)

        #STEP: sort from least to greatest & find the index of the value
        arrPosExon.sort(reverse=False)
        index_posOI=arrPosExon.index(posOI)

        #STEP: retrieve the nearest
        if boolPosGreater:
            nearestPos=arrPosExon[index_posOI+1] if index_posOI!=len(arrPosExon) else None
        else:
            nearestPos=arrPosExon[index_posOI-1] if index_posOI!=0 else None

        #STEP: retrieve the feature range by retrieving the array element that has the position of interest, & return that feature range
        if nearestPos:
            for featTuple in self.arrUniqExons:     #featTuple = feature tuple (i.e. exon tuple)
                if nearestPos in featTuple:
                    nearestFeat=featTuple
        else:
            nearestFeat=None

        return nearestFeat


    def aid_findNearestFeatureV3(self,posOI,arrTupleRange):
        """
        Args:
            posOI (position of interest): integer that is the position of interest. No need to worry about chromosome number since posOI & the positions in arrExons
            arrTupleRange: an array that contains the position ranges, where the first number is the left position & the second number is the right position. Note that for "+" genes the left position is the 5' end and the right end is the 3' end, whereas for "-" gene the contrary is true
        Function:
            this will find the nearest genomic feature to the splice junction based on the 
        Output:
            returns an array where the first element = featNum:0/1 and the second value = distance from either either the left-end (featNum:0) or right-end (featNum:1). featNum usually is the exon/intron number and 0 is the left-end whereas 1 is the right-end
        """

        #STEP: go through each element in the array "arrTupleRange", calculating the difference between both positions
        hashRelativePos={}      #this will record all the distances between posOI & each number in the tuple, & find the smallest value
        for num,eachTuple in enumerate(arrTupleRange):
            #here, :0 signifies the numerically lower nucleotide position
            keyExon_start=str(num)+":0"
            hashRelativePos[keyExon_start]=abs(int(posOI)-int(eachTuple[0]))
            #here, :1 signifies the numerically higher nucleotide position
            keyExon_end=str(num)+":1"
            hashRelativePos[keyExon_end]=abs(int(posOI)-int(eachTuple[1]))

        #STEP: find the feature that is the closest to posOI
        nearestFeature_key=min(hashRelativePos,key=hashRelativePos.get)     #this finds the minimum value in the hash & the associated key

        #OUTPUT:
        #nearestFeature_key = string that is the following format: index_arrTupleRange:index_tuple, where "index_arrTupleRange" can be used to access the element in arrTupleRange & "index_tuple" accesses the value within the element in arrTuple Range
        #hashRelativePos[nearestFeature_key] = this is the distance between posOI & the exon end "nearestFeature_key"
        return [nearestFeature_key, hashRelativePos[nearestFeature_key]]


    def aid_findNearestFeature(self,posOI,arrFeatIndex,arrFeatRange):
        """
        Args:
            posOI (position of interest): integer that is the position of interest. No need to worry about chromosome number since posOI & the positions in arrExons
            arrFeatIndex: an array that contains the indices that will access the numerical range in the array "arrFeatRange". This will usually a list of exon/intron indices that come from aid_accessExons() or aid_accessIntrons()
            arrFeatRange: an array that contains the position ranges for each exon/intron feature, where the first number is the left position & the second number is the right position. Note that for "+" genes the left position is the 5' end and the right end is the 3' end, whereas for "-" gene the contrary is true
        Function:
            this will find the nearest genomic feature to the splice junction, regardless of exon or intron
        Output:
            returns an array where the first element = featNum:0/1 and the second value = distance from either either the left-end (featNum:0) or right-end (featNum:1). featNum usually is the exon/intron number and 0 is the left-end whereas 1 is the right-end
        """

        #STEP: go through each element in the hash, calculating the difference between both positions
        hashRelativePos={}
        for featNum in arrFeatIndex:
            if self.cruzdb_isoformInfo.strand=="+":     #for "+" strand genes, 0 will refer to the numerically lower nucleotide position & 1 will refer to the numerically higher nucleotide position
                #here, :0 signifies 5' end AND is the numerically lower nucleotide position
                keyExon_start=str(featNum)+":0"
                hashRelativePos[keyExon_start]=abs(int(posOI)-int(arrFeatRange[featNum][0]))
                #here, :1 signifies 3' end AND is the numerically higher nucleotide position
                keyExon_end=str(featNum)+":1"
                hashRelativePos[keyExon_end]=abs(int(posOI)-int(arrFeatRange[featNum][1]))
            else:               #for "-" strand genes, 0 will refer to the numerically higher nucleotide position & 1 will refer to the numerically lower nucleotide position
                #here, :0 signifies 5' end BUT is the numerically higher nucleotide position (opposite of "+" strand)
                keyExon_start=str(featNum)+":0"
                hashRelativePos[keyExon_start]=abs(int(posOI)-int(arrFeatRange[featNum][1]))
                #here, :1 signifies 3' end BUT is the numerically lower nucleotide position  (opposite of "+" strand)
                keyExon_end=str(featNum)+":1"
                hashRelativePos[keyExon_end]=abs(int(posOI)-int(arrFeatRange[featNum][0]))

        nearestFeature_key=min(hashRelativePos,key=hashRelativePos.get)     #this finds the minimum value in the hash & the associated key

        #OUTPUT:
        #nearestFeature_key = the exon number & exon start/end code (e.g. 14:0 means the 5' end (0) of exon 14, 7:1 means 3' end (1) of exon 7)
        #hashRelativePos[nearestFeature_key] = this is the distance between posOI & the exon end "nearestFeature_key"
        return [nearestFeature_key, hashRelativePos[nearestFeature_key]]

    def displayObj(self):
        """
        Function:
            this will display all the properties and their values in the object
        """
        for property, value in vars(self).iteritems():
            print property, ": ", value


