#/usr/bin/python
import sys

from cruzdb import Genome

sys.path.insert( 0, "/home/mokha/Documents/Krauthammer_Lab/PythonClasses" )
from SVSv5 import Exon, Isoform, MultiIsoform, SpliceJunction, IsoformSJ

print "------------ Algorithm: 160919_Isoform_1.py ------------"
""" Reconstruct transcripts based on Splice Junctions """


#assign all splice junctions to specific gene: go through cruzdb & find end points for each gene --> assign 
# g = Genome( 'sqlite:////tmp/hg19.db' )
g = Genome( 'sqlite:////tmp/hg19_v2.db' )
Isoform.set_cruzdb( g )

#retrieve gene & print information on it based on 
gene = g.refGene.filter_by( name2 = 'BRAF' ).all()
# all_genes = g.refGene.filter_by( name2 = 'TTN' ).first()
# all_genes = g.refGene.filter_by( name2 = 'AGRN' ).all()
# all_genes = g.refGene.filter_by( name2 = 'AGRN' ).first()
# all_genes = g.refGene.filter_by( name2 = 'DIXDC1' ).all()


for each_isoform in gene:
    obj_iso = Isoform( each_isoform.name )

    #print name
    print obj_iso.isoform_id, ":", obj_iso.gene_sym

    print "obj_iso = ", obj_iso

    # print "exons: "
    # for i, (k,v) in enumerate( obj_iso.hashExonList.iteritems() ):        #k = feature name, v = Exon object
    #     print i, ": ", k, " -> ", v

    # print "introns: "
    # for i, (k,v) in enumerate( obj_iso.hashIntronList.iteritems() ):        #k = feature name, v = Exon object
    #     print i, ": ", k, " -> ", v

    print "exons: "
    for x in range( 0, obj_iso.last_exon_num ):
        get_exon = obj_iso.get_exon_num( x )
        print x, ": ", get_exon

    print "introns: "
    for x in range( 0, obj_iso.last_exon_num - 1 ):
        get_intron = obj_iso.get_intron_num( x )
        print x, ": ", get_intron

print "------------ Algorithm Completed: 160919_Isoform_1.py ------------"