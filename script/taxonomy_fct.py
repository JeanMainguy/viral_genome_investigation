import os
import logging
import re
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
from Bio.Alphabet import generic_protein


logging.basicConfig(filename='refseq_genome_path.log')
def see_objet(obj):
    print('SEE OBJET ', type(obj))

    for attr in dir(obj):
        if attr[0] != '_' and attr != "seq":
            print(attr, " ", getattr(obj, attr))


def getGbffFile(virus_name, ncbi_refseq_db="/mirror/ncbi/current/genomes/refseq/viral/"):
    # ncbi_refseq_db = "/mirror/ncbi/current/genomes/refseq/viral/"
    # virus_name = "Ageratum_yellow_vein_China_virus"

    virus_path = os.path.join(ncbi_refseq_db, virus_name)

    assembly_path = os.path.join(virus_path, 'latest_assembly_versions')
    gb_files = []
    # Checking if the path exist..
    if not os.path.isdir(virus_path):
        logging.warning('Cannot find the virus folder: '+virus_path)
        # raise Exception(virus_path, 'the virus folder does not exist..')
        return []

    if not os.path.isdir(assembly_path):
        logging.warning('No latest_assembly_versions folder for {}'.format(virus_name))
        # raise Exception(assembly_path, 'the virus folder does not have "latest_assembly_versions" directory...')
        return []
    if len(os.listdir(assembly_path)) == 0:
        logging.warning('the path {} is empty'.format(assembly_path))

    for assembly in os.listdir(assembly_path):
        gb_file = os.path.join(assembly_path, assembly, assembly+'_genomic.gbff.gz')
        # print(gb_file)
        if os.path.exists(gb_file):
            gb_files.append(gb_file)
        else:
            logging.warning('No gbff file have been found in the directory: '+gb_file)
            # raise Exception(assembly_path, 'No gbff file have been found in the directory')
    return gb_files

def getTaxonomy(gb_file):
    polyprot = []
    matpep = []
    poly_di = {} # dict with poly as key and list of matpep as value

    # print('='*20)
    taxonomy_set = set()
    taxon_set = set()
    with gzip.open(gb_file, "rt") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            taxonomy = record.annotations['taxonomy']
            taxonomy_set.add(";".join(taxonomy))
            organism = record.annotations['organism']
            # print(record.annotations)
            for f in record.features:
                if f.type == 'source':
                    taxon = f.qualifiers['db_xref'][0].replace('taxon:', '')
                    # taxon_set.add(taxon)
            # see_objet(record)
        if len(taxonomy_set) == 1:
            return taxonomy_set.pop(), organism, taxon
        else:
             logging.warning('Taxonomy is divergente in '+gb_file)
             print(taxonomy_set)

def createTaxonomyFile(output_dir, ncbi_refseq_db="/mirror/ncbi/current/genomes/refseq/viral/"):
    """
    create csv file from the genome refseq db with info on avaible genome
    """
    taxonomy_file = os.path.join(output_dir, 'taxonomy_virus.txt')
    with open(taxonomy_file, 'w') as handle:
        for i, virus in enumerate(os.listdir(ncbi_refseq_db)):
            gb_files = getGbffFile(virus, ncbi_refseq_db)
            # print(gb_files)
            for gb in gb_files:
                taxonomy, organism, taxon = getTaxonomy(gb)
                handle.write( "{}\t{}\t{}\t{}\n".format(taxon, organism, taxonomy, gb))

            if i%100 == 0:
                print(i)
        # print i
    return taxonomy_file

def getAllRefseqFromTaxon(wanted_taxonomy, taxonomy_file="/home/user/mainguy/Documents/Data_Analysis/data/taxonomy/taxonomy_virus.txt"):
    with open(taxonomy_file, 'r') as taxfl:

        for l in taxfl:
            # print(l)
            (tax_id, organism, taxonomy, gbff) = l.split("\t")
            # print(taxonomy)
            if wanted_taxonomy in taxonomy.split(';'):
                # print('taxonomy')
                yield gbff.rstrip()



if __name__ == '__main__':

    # taxonomy_file = createTaxonomyFile("/home/user/mainguy/Documents/Data_Analysis/data/taxonomy/")
    taxonomy_file = "/home/user/mainguy/Documents/Data_Analysis/data/taxonomy/taxonomy_virus.txt"
    taxon = "Nidovirales"
    gb_file_iter = getAllRefseqFromTaxon(taxon, taxonomy_file)

    for g in gb_file_iter:
        print(g)


# print(len(list(gb_file_iter)))
