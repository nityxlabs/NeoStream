# NeoStream
Mutanome &amp; Transcriptome-Derived Neoepitope Discovery System

## Getting Started

The following

### Prequisites 

In order to run NeoStream, the following libraries will also need to be installed:
* cruzdb (https://github.com/brentp/cruzdb)
* HTSeq (https://htseq.readthedocs.io/en/release_0.10.0/)
* Tabix (https://pypi.org/project/pytabix/)
* requests
* pandas
* numpy
* scipy

### Installing

To begin using NeoStream, clone this Github repository to your computer.

```
git clone https://github.com/nityxlabs/NeoStream.git
```

The directory structure should be as follows:

```
/Algorithm
/Data - contains file for neoepitope comparison 180106_Proteasome_TAP_Chart.txt
/Genome - contains hg19.fa & hg19.fa.fai (samtools-indexed) genome build file
/Results
```

The hg19.fa build and hg19.fa.fai (samtools-indexed) files should be located in the /Genome directory. To create the samtools

Next, run the following script to set up local genomic databases
```
python setup_cruzdb_databases.py
```

## Deployment

To execute NeoStream, perform the following procedure:

```py
# python script "180820_Pipeline_NeoStream_V2.py [filename] [prefix_label]"
# filename = the file in the ../Data path that is being analyzed (example: PDAC.final.somatic.maf)
# prefix_label = string that will be appended to each output file in the /Results directory (example: 180901)
python 180820_Pipeline_NeoStream_V2.py file.maf 180901
```

## Output

The algorithm will output information into the /Results directory:

```
* [prefix_label]_Thres0_[filename]_ProcessGenomeAlt_V1.txt - generates a list of potential neoepitopes from a given mutation
* [prefix_label]_Thres0_[filename]_NeoepCompare_V1.txt - contains the processing efficiency & MHC affinity for each neoepitope derived from a given mutation
* [prefix_label]_Thres0_[filename]_NeoepCompare_V1_Extended.txt - further identifies information associated with neoepitopes, including assignment of high MHC affinity allele
```

---------
```
NOTE: note to self - finish documentation looking at following link: https://gist.github.com/PurpleBooth/109311bb0361f32d87a2
```
