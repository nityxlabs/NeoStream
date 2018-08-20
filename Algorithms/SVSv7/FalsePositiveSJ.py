#/usr/bin/python
import tabix

class FalsePositiveSJ( SpliceJunction ):
    """ retrieves all aberrant splice junctions & quantifies how many samples """
    def __init__( self, path_sj, control_sj ):
        #retrieve all splice junction files -> get all splice junctions for specific gene -> find splice junctions that are aberrant -> determine how prevalent each SJ is across all samples
        self.path_sj = path_sj
        self.control_sj = control_sj

    def sj_prevalence( self, gene_sym ):
        """ Function: determines prevalence of aberrant splice junctions in specific gene across all samples """
        pass

        gene = obj_cruzdb.refGene.filter_by( name2 = gene_sym ).all()

        for isoform in gene:
            if not isoform.is_coding:
                continue


    def sj_prevalence_all( self ):
        """ Function: quantifies how prevalent a splice junction is across all samples & control samples """
        all_genes = obj_cruzdb.refGene.filter_by().all()

        list_genes = []     #records all
        for i, each_gene in enumerate( all_genes ):
            gene_sym = each_gene.name2

            #make sure gene is protein coding
            if not each_gene.is_coding:
                continue

            #make sure not gene name duplicate, and record gene name as to not repeat
            if gene_sym in list_genes:
                continue
            else:
                list_genes.append(gene_sym)