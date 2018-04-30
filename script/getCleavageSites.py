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



def extractCleavageSites(gb_file,file_handle, genetic_code, gff_file, taxon_expectation, polyprotein_number):


    genome = obj.Genome( gb_file)

    with gzip.open(gb_file, "rt") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "genbank")):


            segment = obj.Segment(record, gb_file)
            genome.segments.append(segment)

            segment.getMatpeptidesAndPolyproteins()
            segment.checkPeptideRedundancy() #remove the redundant peptide
            segment.checkSubPeptides()
            segment.associatePepWithProt()

            segment.checkForSlippage()
            segment.identifySubProtein()

            segment.getCleavageSites()

            # segment.writeCleavageSite(file_handle, genetic_code, window_step)

    genome.getTaxonExpectation(taxon_expectation)

    genome.identifyExpectedElement()

    genome.getMatchObject(gff_file)
    genome.associateMatchWithPolyprotein()
    genome.visualisation(1, genetic_code)

    for segment in genome.segments:
        for cds in segment.cds:
            if cds.polyprotein_number == polyprotein_number:
                seq_to_write = cds.getSequenceRecord(segment.organism, segment.taxon_id, segment.record, genetic_code)
                SeqIO.write(seq_to_write, file_handle,"fasta")



if __name__ == '__main__':

    # logging.basicConfig(filename='log/genbankparser.log',level=logging.INFO)

    polyprotein_number = 2.0

    taxonomy_file="data/taxonomy/taxonomy_virus.txt"
    expected_file =  "polyprotein_expectation_by_taxon.csv"

    output_dir = '/home/user/mainguy/Documents/Data_Analysis/data/Cleavage_site_sequences'
    gff_file = 'data/interproscan_result/annotated_polyprotein_by_taxon/polyprotein_All.gff3'


    taxon = 'ssRNA viruses'
    taxon = "Togaviridae"
    taxon = 'Rubivirus'
    taxon = "Marafivirus"
    taxon= "Togaviridae"
    taxon='ssRNA positive-strand viruses, no DNA stage'
    taxon='Alphavirus'
    # taxon="11036"

    taxon_expectation = tax.expectedPeptide(expected_file)

    file_handle = open(os.path.join(output_dir,'{}_polyprotein_{}.faa'.format(taxon,polyprotein_number  )), "w")

    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)


    for i, gb_dico in enumerate(gbff_iter):
        # print(gb_dico)
        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']

        if i%50 == 0:
            print(i)

        extractCleavageSites(gb_file,file_handle, genetic_code, gff_file, taxon_expectation, polyprotein_number)

    print(i+1, 'Genome anlalysed from taxon', taxon)
