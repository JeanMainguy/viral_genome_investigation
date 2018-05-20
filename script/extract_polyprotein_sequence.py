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



def extractCleavageSites(gb_file,dico_file, genetic_code, gff_file, taxon_expectation, polyprotein_number):


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

    # genome.visualisation(1, genetic_code)
    if not genome.polyprotein_expectation:
        print("NO POLYPROTEIN NUMBER FOR: ",genome.expectation_node+'|'+genome.organism)
        input()
    header_info = genome.expectation_node+'|'+genome.taxon_id
    print(header_info)
    for i, segment in enumerate(genome.segments):
        for cds in segment.cds:

            if cds.polyprotein_number: # == polyprotein_number and i+1 == wanted_segment:
                key = '{}_{}'.format(i+1, cds.polyprotein_number) # first segment and then the polyprotein number
                if key not in dico_file:
                    file_name = '{}_{}.faa'.format(dico_file["prefix_file"],key)
                    print(file_name)

                    dico_file[key] = open(file_name, "w")
                # header_info += '|annotation' if .segments else '|no_annotation'

                seq_to_write = cds.getSequenceRecord(segment.organism, header_info, segment.record, genetic_code)
                SeqIO.write(seq_to_write, dico_file[key] ,"fasta")



def polyproteinsFileExist(taxon, output_dir):
    for file in os.listdir(output_dir):
        if file.startswith(taxon):
            return True
    return False

if __name__ == '__main__':

    # logging.basicConfig(filename='log/genbankparser.log',level=logging.INFO)
    polyprotein_number = 1.0

    taxonomy_file="data/taxonomy/taxonomy_virus.txt"
    expected_file =  "data/taxonomy/polyprotein_expectation_by_taxon.csv"

    output_dir = '/home/user/mainguy/Documents/Data_Analysis/data/polyprotein'
    gff_file = 'data/interproscan_result/annotated_polyprotein_by_taxon/polyprotein_All.gff3'






    taxon_expectation = tax.expectedPeptide(expected_file)

    for taxon in taxon_expectation:

        if not taxon_expectation[taxon]['polyprotein']:
            logging.info('Information about {} are not complete in the expectation file: {}. Polyprotein are then not extracted'.format(taxon, expected_file))
            continue
        if polyproteinsFileExist(taxon, output_dir):
            logging.warning('polyproteins from {} have been already extracted'.format(taxon))
            # continue
        taxon= "12139"
        prefix_file = os.path.join(output_dir,'{}_polyprotein'.format(taxon))

        dico_file = {'prefix_file':prefix_file}

        gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)


        for i, gb_dico in enumerate(gbff_iter):
            # print(gb_dico)
            gb_file = gb_dico['gb_file']
            genetic_code = gb_dico['genetic_code']

            extractCleavageSites(gb_file,dico_file, genetic_code, gff_file, taxon_expectation, polyprotein_number)


        print('Polyproteins extraction completed for taxon', taxon)


        [ v.close() for k, v in dico_file.items() if k != 'prefix_file']
        break
