#!/usr/bin/env python3
import function_parser as fct
import taxonomy_fct as tax


def getPolyproteinFaa(gb_file):
    record, polyprot, matpep, poly_di = fct.extractPolyMatpep(gb_file)

    organism = record.annotations['organism']

    file_name = "polyprotein_{}.faa".format(organism.replace(' ', "_"))

    print(record)
    fct.write_prot(polyprot, file_name, organism)


if __name__ == '__main__':

    taxon = 'Viruses'
    # taxon = "Nidovirales"
    gbff_iter = tax.getAllRefseqFromTaxon(taxon)

    # print(len(list(gbff_iter)))
    # gb_file = next(gbff_iter)
    types = set()
    for i, gb_file in enumerate(gbff_iter):
        # print(gb_file)
        record, polyprot, matpep, poly_di = fct.extractPolyMatpep(gb_file, types)
        if i%100 == 0:
            print(i)
        
        # print('pol', polyprot)
    print(types)

    # types = set()
    # gb_file = "/mirror/ncbi/current/genomes/refseq/viral/Lactate_dehydrogenase-elevating_virus/latest_assembly_versions/GCF_000850185.1_ViralProj14702/GCF_000850185.1_ViralProj14702_genomic.gbff.gz"
    # record, polyprot, matpep, poly_di = fct.extractPolyMatpep(gb_file, types)
