import taxonomy as tax
import object_analysis as obj

import os, gzip, logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys, csv, re
from Bio.Alphabet import generic_protein
from Bio.SeqFeature import SeqFeature, FeatureLocation
from operator import attrgetter



def extractProteins(gb_file, handle_poly, handle_prot, genetic_code):


    genome = obj.Genome( gb_file)

    with gzip.open(gb_file, "rt") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "genbank")):


            segment = obj.Segment(record, gb_file)
            genome.segments.append(segment)

            segment.getMatpeptidesAndPolyproteins()
            if not segment.peptides:
                continue

            segment.checkPeptideRedundancy() #remove the redundant peptide
            segment.checkSubPeptides()
            segment.associatePepWithProt()

            segment.checkForSlippage()
            segment.identifySubProtein()

            segment.getCleavageSites()



    # genome.visualisation(1, genetic_code)
    print(segment.taxon_id)
    for i, segment in enumerate(genome.segments):
        header_info = segment.taxon_id

        for p, cds in enumerate(segment.cds):
            # for pep in cds.peptides:
            #     print(pep)
            key = '{}_{}'.format(i+1, p+1) # first segment and then the polyprotein number

            # header_info += '|annotation' if .segments else '|no_annotation'

            header = '|'.join([header_info, cds.protein_id, key])

            seq_to_write = cds.getSequenceRecord(segment.organism, header, segment.record, genetic_code)

            SeqIO.write(seq_to_write, handle_prot,"fasta")

            if cds.peptides:
                SeqIO.write(seq_to_write, handle_poly,"fasta")


if __name__ == '__main__':

    # logging.basicConfig(filename='log/genbankparser.log',level=logging.INFO)
    try:
        taxon = sys.argv[1]
    except IndexError:
        taxon = "Viruses"

    taxonomy_file="data/taxonomy/taxonomy_virus.txt"

    output_dir = 'data/polyprotein_sequence'

    polyprotein_db = '{}_polyprotein_db.faa'.format(taxon)
    protein_db = "{}_protein_db.faa".format(taxon)

    handle_poly = open(os.path.join(output_dir,polyprotein_db), "w")
    handle_prot = open(os.path.join(output_dir,protein_db), "w")

    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)

    for i, gb_dico in enumerate(gbff_iter):
        # print(gb_dico)
        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']
        # print(gb_file)
        extractProteins(gb_file, handle_poly, handle_prot, genetic_code)
        if i%100 == 0:
            print(i)
    print('Polyproteins extraction completed for taxon', taxon)

    handle_poly.close()
    handle_prot.close()
