#/usr/bin/python
#Script: setup_cruzdb_databases.py
from cruzdb import Genome

arrTables_hg19 = ["refGene","knownGene","ensGene","ccdsKgMap","knownGeneMrna","kgProtAlias","knownToEnsembl","knownToRefSeq", "wgEncodeGencodeBasicV19"]
#NOTE: "ensGene" is not present in hg38

#NOTE: "ensGene"
arrTables_hg38 = ["refGene","knownGene","ccdsKgMap","knownGeneMrna","kgProtAlias","knownToEnsembl","knownToRefSeq"]


# db_hg19_path = "sqlite:////tmp/hg19.db"
# db_hg38_path = "sqlite:////tmp/hg38.db"
db_hg19_path = "sqlite:////tmp/hg19_v2.db"
db_hg38_path = "sqlite:////tmp/hg38_v2.db"
Genome( db = "hg19" ).mirror( arrTables_hg19, db_hg19_path )
Genome( db = "hg38" ).mirror( arrTables_hg38, db_hg38_path )
