#/usr/bin/python

import math
import time
import subprocess
import re


"""
To Kill a Python script on Linux, because Ctrl+c doesn't work:
	press Ctrl+Shift+\
	pressing Ctrl+z only pauses the process apparently
"""

# %c char single character
# %d (%i) int signed integer
# %e (%E) float or double exponential format
# %f float or double signed decimal
# %g (%G) float or double use %f or %e as required
# %o int unsigned octal value
# %p pointer address stored in pointer
# %s array of char sequence of characters
# %u int unsigned decimal
# %x (%X) int unsigned hex value


##CLASSES - manipulate matrix hash
class manipMatrixHash:

	def __init__(self,getMatrix):
		"""
		Function:
			initializes the manipMatrix class 
		"""
		self.rowSize=len(getMatrix)
		self.colSize=len(getMatrix[0])
		self.matrix={}		#matrix will actually be a hash (i.e. dictionary)
		
		#STEP: go through each element
		for i,row in enumerate(getMatrix):
			for j,col in enumerate(row): 
				strIndex=i+","+j
				self.matrix[strIndex]=col 		#OR can use self.matrix[strIndex]=getMatrix[i][j]
	
	def printMatrix(self):
		"""
		Function:
			prints the matrix
		"""
		for k,v in self.matrix.iteritems():
			print k,"=",v


class manipMatrix:

	def __init__(self,getMatrix):
		"""
		Function:
			initializes the manipMatrix class 
		"""
		self.matrix=getMatrix
		self.rowSize=len(getMatrix)
		self.colSize=len(getMatrix[0])

	def createArr_column(self,colOI):
		"""
		Args:
			colOI = integer that is the column of interest
		Function:
			this function will return an array with all the values from a specific column from self.matrix
		"""
		arrColOI=[]
		for i,row in enumerate(self.matrix):
			arrColOI.append(row[colOI])

		return arrColOI

	def createHash_twoColumns(self,colKey,colVal):
		"""
		Args
			colKey & colVal: integers that refer to a column in the matrix "self.matrix"
		Function:
			This will correspond column 1 "colKey" and column 2 "colVal" from the matrix
		Errors:
			if colKey or colVal are out of range
		"""
		hashCorrespond={}		#hashCorrespond = hash (i.e. dictionary) that will record colKey as the key & colVal as the value

		#STEP: go through each element in for both columns
		for i,row in enumerate(self.matrix):
			hashCorrespond[row[colKey]]=row[colVal]

		return hashCorrespond 

	def printMatrix(self):
		"""
		Function:
			prints the matrix
		"""
		#STEP: go through each row of the 2D array
		for i,row in enumerate(self.matrix):
			for j,col in enumerate(row):
				if j==0:
					print "[",i,"]: ",
				print self.matrix[i][j],"\t",		#OR can do print col,"\t"
			print "\n"



###FUNCTIONS###

##START: Time functions##
def timeElapse_convertToHMS(getSeconds):
	"""
	Function:
		Will display the time in hours, minutes, and seconds
	Note:
		This will require "start_time=time.time()" in the beginning of the code, "elapse_time=time.time() - start_time" at the end of the code. Note that the variable names are arbitrary
	"""
	calcHours = int(getSeconds/3600)
	calcMin = int(getSeconds/60 - calcHours*60)
	calcSeconds = int(getSeconds-(calcMin*60 + calcHours*3600))
	# calcMilliSec = (((getSeconds-(calcMin*60 + calcHours*3600)))-calcSeconds)*1000
	
	print "Hours = %d: Min = %d: Seconds = %d" % (calcHours, calcMin, calcSeconds)

##END: Time functions##

##START: Data Measure functions##

def dataMeasure_countLineByLine(filePath):
	"""
	Function:
		this function will quantify the number of lines a file has and return that number
	"""
	return int(subprocess.check_output(['wc','-l',filePath]).strip().split(" ")[0])

def dataMeasure_countAllLines(filePath):
    """     
    Function:
        quantifies the number of lines a file has and return that number
    Note:
        this is the same as the "def dataMeasure_countLineByLine(filePath)"

    """
    with open(filePath) as f:
        for i, r in enumerate(f):
            pass

    return i + 1

def dataMeasure_rowsPeriodicUpdate( totalNumRows, currRow, updateStat, updateRate = 0.1):
	"""
	Args:
		totalNumRows (integer) = total # of rows in dataset
		currRow (integer) = current row position in dataset, updateStat (decimal #) = the threshold to send an update of percentage number of completion
		updateStat = a float that will be updated by the rate of "updateRate"
		updateRate = refers to percent rate increase, e.g. 0.1 = 10%
	Function:
		considers the # of lines in a dataset, and periodically returns update of completion of analysis
	"""
	#STEP: retrieve the # of significant figures based on the updateRate - this makes the assumption that updateRate < 1
	roundNum=int(math.fabs(math.log10(updateRate))+1)

	if(currRow>=(totalNumRows*updateStat)):
		calcPercentage=(float(currRow)/float(totalNumRows))*100
		calcPercentage=round(calcPercentage,roundNum)
		print "%u of %u analyzed - %f percent" % (currRow,totalNumRows,calcPercentage)
		updateStat+=updateRate;

	return updateStat

##END: Data Measure functions##


##START: File Manipulation functions##
def fileManip_handleTextFile(filePath):
	"""
	Function:
		this function will replace all new line characters from Mac OS (\r) and Windows (\r\n)
	NOTE: for indicating new lines, each operating system does something differently:
		\r = use as new line character in Mac OS
		\n = used as new line character in Unix
		\r\n = used as new line character in Windows
	"""

	#STEP: read in file and convert each return carrier to the Unix return carrier (\n)
	getFile=open(filePath,"r")
	fileStr=re.sub(r"\r\n?","\n",getFile.read())

	return fileStr


def fileManip_createMatrix(filePath,skipRowNum=0,minRowLen=1):
	"""
	Args:
		filePath: string that is the path to the file
		skipRowNum: integer that is the number of rows to skip in the file before recording to the matrix "matrixFile"
		minRowLen: integer that is the minimum length of each row, if it is not this length then it will not be recorded into the matrix "matrixFile"
	Function:
		this function will read in a tab-delimited file and split it into a 2D array (i.e. matrix)
	NOTE: for indicating new lines, each operating system does something differently:
		\r = use as new line character in Mac OS
		\n = used as new line character in Unix
		\r\n = used as new line character in Windows
	"""

	#STEP: read in the file
	fileRead=fileManip_handleTextFile(filePath)		#this will convert all new line carriers (e.g. \r & \r\n) to \n
	arrFile=fileRead.split("\n")
	matrixFile=[]

	#STEP: go through each line of the file
	i2=0			#this is an auxiliary counter to add elements to the matrix "matrixFile"
	for i,row in enumerate(arrFile):
		if i>=skipRowNum:
			if len(row.split("\t"))>=minRowLen:
				matrixFile.append([[]])
				matrixFile[i2]=row.split("\t")
				i2+=1

	objMM=manipMatrix(matrixFile)		#objMM = object Manip Matrix

	return objMM

##END: File Manipulation functions##

##START: Math functions##
def mathTools_roundNearestNum(valToRound, baseNum, boolRoundUp = True):
	"""
	Args:
		valToRound =  integer/float value to round to the nearest base value "baseNum"
		baseNum = integer/float value to round to.
		boolRoundUp = if true, then will round up, else will round down
	Function:
		this will round the number "valToRound" to the nearest value baseNum
	"""
	if boolRoundUp:
		return int((valToRound + baseNum) / baseNum) * baseNum
	else:
		return int(valToRound / baseNum) * baseNum

##END: Math functions


##START: Array functions##
def arrManip_convertStrToInt(getArr):
	"""
	Function:
		will convert an array of strings to an array of integers
	"""
	return map(int, getArr)		#another way to do this: getArr = [int(i) for i in getArr]

def arrManip_convertIntToStr(getArr):
	"""
	Function:
		will convert an array of strings to an array of integers
	"""
	return map(str, getArr)		#another way to do this: getArr = [int(i) for i in getArr]

def arrManip_convertStrToFloat(getArr):
	"""
	Function:
		will convert an array of strings to an array of integers
	"""
	return map(float, getArr)		#another way to do this: getArr = [int(i) for i in getArr]

def arrManip_calcSum(getArr):
	getArr_float = arrManip_convertStrToFloat(getArr)
	return sum(getArr_float)

def arrManip_calcAvg(getArr):
	getArr_float = arrManip_convertStrToFloat(getArr)
	return sum(getArr_float) / float(len(getArr))

def arrManip_maxVal(getArr,toNum=False):
	"""
	Args:
		getArr = array of values
		toNum = boolean that, if true, will convert the values to float, else if false will keep the values 
	Function:
		will return the max value in the list "getArr", either for string or number
	"""
	if(toNum):
		getArr2=arrManip_convertStrToFloat(getArr)
		return max(getArr2)
	else:
		return max(getArr)

def arrManip_minVal(getArr,toNum=False):
	"""
	Args:
		getArr = array of values
		toNum = boolean that, if true, will convert the values to float, else if false will keep the values 
	Function:
		will return the max value in the list "getArr", either for string or number
	"""
	if(toNum):
		getArr2=arrManip_convertStrToFloat(getArr)
		return min(getArr2)
	else:
		return min(getArr)

	
def arrManip_findStrInElem(listStr,strOI,boolCaseInsen=False):
	"""
	Args:
		listStr = array where each element is a string, looking to see if the string "strOI" is found
		strOI (string of interest) = search for the string of interest in each element in the array "listStr"
		boolCaseInsen = boolean that, if true, will mean to disregard the case sensitivity, else if false will still consider case sensitivity
	Function:
		this function will find and return the element in a list (i.e. array) that contains the string of interest
	Returns:
		returns a dictionary (i.e. hash) with 2 elements, list index & the string found
	"""

	#STEP: go through each element in the list (i.e. array) until the string is found
	for c,x in enumerate(listStr):
		#STEP: first see if case insensitivity is true or not
		if(boolCaseInsen):
			x2=x.lower()
			strOI_2=strOI.lower()
		else:
			x2=x
			strOI_2=strOI

		#STEP: if string is found, then break out of loop
		if(x2.find(strOI_2)>=0):
			listIndex=c
			strFound=x2.strip()		#remove the whitespaces from the left & right side of the string
			break

	return {"listIndex":listIndex,"strFound":strFound}

def arrManip_removeDuplicates(getArr):
	"""
	Function:
		this function will remove duplicates from the array "getArr"
	"""
	return list(set(getArr))

def arrManip_removeBlanks(getArr):
	return filter(None, getArr)

def arrManip_removeVal(getArr,valOI):
	newArr=[x for x in getArr if x!=valOI]
	return newArr
##END: Array functions##

##START: Hash functions##
def hashManip_sortByValThenKey(hash,boolGtoL=False):
	"""
	Args:
		hash = hash to be sorted
		boolGtoL (boolean Greatest to Least) = if true, will return hash with keys sorted from greatest to least, else will return keys from least to greatest
	Function:
		will sort hash based on values then keys
	Output:
		will return an array of tuples, where the 1st element is the key & the 2nd element is the value
	"""
	return sorted(hash.items(), key=lambda k:(k[1],k[0]), reverse=boolGtoL)

##START: String functions##
def strFunc_isStrPresent(fullStr,findStr,boolCaseInsen=True):
	"""
	NOTE:
		re.search(findString,fullStr) will return an object called a "match object". Can extract information from the match object using the following:
			-group(): Return the string matched by the RE (e.g. re.search(findString,fullStr).group())
			-start(): Return the starting position of the match (e.g. re.search(findString,fullStr).group())
			-end(): Return the ending position of the match (e.g. re.search(findString,fullStr).group())
			-span(): Return a tuple containing the (start, end) positions of the match (e.g. re.search(findString,fullStr).group())
	NOTE: r"%s" % findStr gives a different answer than re.escape(findStr)
	Args:
		fullStr: the full string to search into
		findStr: the string to find in "fullStr"
		boolCaseInsen: boolean that, if True, will ignore letter casing. Else will be case-sensitive
	Function:
		will check if string "findStr" is present in "fullStr"
	"""
	#STEP: check if string "findStr" is present in "fullStr"
	#NOTE: r"%s" % findStr gives a different answer than re.escape(findStr)
	if boolCaseInsen:
		checkStrPresent=re.search(r"%s" % findStr,fullStr,re.IGNORECASE)		#can also do: re.search(r""+findStr,fullStr,re.IGNORECASE)
	else:
		checkStrPresent=re.search(r"%s" % findStr,fullStr)		#can also do: re.search(r""+findStr,fullStr)

	if checkStrPresent is None:
		return False
	else:
		return True

def strFunc_strWBPresent(fullStr,findStr,boolCaseInsen=True):
	"""
	strWBPresent = string Word Boundary Present
	Args:
		fullStr: the full string to search into
		findStr: the string to find in "fullStr"
		boolCaseInsen: boolean that, if True, will ignore letter casing. Else will be case-sensitive
	Function: 
		this function will see if only the word (and not part of the word) is present by using the regular expression word boundaries "\b"
	"""
	#STEP: check if string "findStr" is present in "fullStr"
	if boolCaseInsen:
		checkStrPresent=re.search(r'\b%s\b' % findStr,fullStr,re.IGNORECASE)
	else:
		checkStrPresent=re.search(r'\b%s\b' % findStr,fullStr)

	if checkStrPresent is None:
		return False
	else:
		return True


def strFunc_countStr(strOI,fullStr):
	"""
	Args:
		strOI = string of interest that will be quantified in string "fullStr"
		fullStr = longer string that may or may not contain strOI
	Function:
		this function will count the frequency of strOI in fullStr
	"""
	return fullStr.count(strOI)

def strFunc_replaceStr(oldStr,newStr,fullStr):
	"""
	Args:
		fullStr: the full string to search into
		findStr: the string to find in "fullStr"
		boolCaseInsen: boolean that, if True, will ignore letter casing. Else will be case-sensitive
	Function: 
		this function will see if only the word (and not part of the word) is present by using the regular expression word boundaries "\b"
	"""
	return re.sub(r"%s" % oldStr,r"%s" % newStr,fullStr)

##END: String functions##


##START: directory functions##
def dirFunc_dirInfo(dirPath):
	"""
	Args:
		dirPath = string that is path to the directory of interest
	Function:
		will retrieve the directory path, name of directories, and name of files
	CREDIT: http://stackoverflow.com/questions/3207219/how-to-list-all-files-of-a-directory-in-python
	"""
	for (dirpath, dirnames, filenames) in os.walk(dirPath):
		break

	return [dirpath, dirnames, filenames]
##START: directory functions##

##START: Print functions##

def printList(getList,boolReverse=False):
	"""
	Function:
		displays content of an array aka python list
	"""
	if not boolReverse:		#equivalent to "if !boolReverse"
		for c,elem in enumerate(getList):
			print c,": ",elem
	else:
		for i in range(len(getList)*-1,1):
			print i,": ",getList[i]

def printDict(getDict):
	"""
	Function:
		prints the hash aka python dictionary
	"""
	#STEP: go through each row of the 2D array
	for k in getDict:
		print k,": ",getDict[k]

##END: Print functions##


##REFERENCE FUNCTIONS - These are functions just for understanding how things work, don't really have to use these##


##BIOLOGY FUNCTIONS##

##START: Samtools##

def samTools_quantRegion(path_bamFile,chrNum,posStart,posEnd):
	"""
	Args:
		path_bamFile = string that is the path to the .bam file
		chrNum = string that is the chromosome number, in the format "chr#" (e.g. chr9,chr12)
		posStart = integer that is the start nucleotide position
		posEnd = integer that is the posEnd nucleotide position
	Function:
		will quantify the number reads within a specific region and return an integer of the number of reads within that region
	"""
	#STEP: create the genomic range & then quantify the number of reads within this range
	strGenomicRange=str(chrNum)+":"+str(posStart)+"-"+str(posEnd)

	return int(subprocess.check_output(['samtools','view','-c',path_bamFile,strGenomicRange]))

def samTools_totalReadCount(path_bamFile):
	"""
	Args:
		path_bamFile = string that is the path to the .bam file
	Function:
		This will sum the total number of reads in a sample
	NOTE: the output of "samtools idxstats file.bam" --> The output is TAB-delimited with each line consisting of reference sequence name (column 0, e.g. chr3), sequence length (column 1), # mapped reads (column 2) and # unmapped reads (column 3). It is written to stdout.
	"""
	colChrNum=0
	colChrNumLength=1
	colMappedReads=2
	colUnmappedReads=3
	#STEP: create table for read count information 
	tableStats=subprocess.check_output(['samtools','idxstats',path_bamFile]).split("\n")

	#STEP: sum the total number of reads
	totalReadCount=0		#totalReadCount = will keep track of the total read count
	for count,row in enumerate(tableStats):
		arrRowCol=row.split("\t")
		if len(arrRowCol)==4:		#NOTE: this is 4 because each row should have 4 columns. Look at the notes in this function
			totalReadCount+=int(arrRowCol[colMappedReads])

	return totalReadCount

def samTools_mappedReadCount(filePath):
    """
    
    Function:
        this function will quantify the number of mapped reads in a sample. Actually, check function samTools_totalReadCount() as this usually produces the same number at is faster than this command
    Note:
        -0x04 = smatools flag for unmapped read
        -"-F" means to look for the opposite, so -F 0x04 means to look for mapped reads. "-f 0x04" means to look for unmapped reads
        -"-c" means count 
        -samtools view -F 0x04 -c [file.bam] = count the number of mapped reads
    """
    return int( subprocess.check_output(["samtools", "view", "-F", "0x04", "-c", filePath]) )