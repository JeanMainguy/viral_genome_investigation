import taxonomy as tax
import object_analysis as obj

import os, gzip, logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys, csv
from Bio.Alphabet import generic_protein
from Bio.SeqFeature import SeqFeature, FeatureLocation
from operator import attrgetter


def extractSpecificPolyprotein(gb_file,file_handle, genetic_code, gff_file, taxon_expectation, window_step, nb_line):

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
    # genome.getTaxonExpectation(taxon_expectation)
    # # print(genome.expectation_node)
    # print(genome.expectation_node)
    # if not genome.hasEnoughPeptide():
    genome.getMatchObject(gff_file)
    genome.associateMatchWithPolyprotein()
    genome.visualisation(nb_line)

    for segment in genome.segments:
        for cds in segment.cds:
            if cds.position_number == 1.0:
                seq_to_write = polyprotein.getSequenceRecord(segment.organism, segment.taxon_id, segment.record, genetic_code)
                SeqIO.write(seq_to_write, file_handle,"fasta")

        # print(genome.matchs)




def writeAnnotatedProtein(gb_file, file_handle, genetic_code, taxon_expectation):
        genome = obj.Genome( gb_file)



        with gzip.open(gb_file, "rt") as handle:
            for i, record in enumerate(SeqIO.parse(handle, "genbank")):


                segment = obj.Segment(record, gb_file)


                segment.getMatpeptidesAndPolyproteins()
                segment.checkPeptideRedundancy() #remove the redundant peptide
                segment.checkSubPeptides()
                segment.associatePepWithProt()

                segment.checkForSlippage()
                segment.identifySubPolyprotein()
                genome.segments.append(segment)

        genome.getTaxonExpectation(taxon_expectation)
        # print(genome.expectation_node)
        for segment in genome.segments:
            segment.writeAnnotatedProteins(file_handle, genetic_code)
            # file_handle = files_dico[genome.expectation_node]
            # #we don't want to write polyprotein that come from taxon where we know that there is not normally polyprotein
            # #But we still write the polyprotein where the taon hasn't been found in taxon expectation
            # if segment.mat_peptides and genome.expectation_node != "Not_Found" and taxon_expectation[genome.expectation_node] != 0:
            #     if isinstance(files_dico[genome.expectation_node ], str): #if file handle is a string then we create the file
            #         files_dico[genome.expectation_node ] = open(files_dico[genome.expectation_node ], "w")
            #     segment.writeAnnotatedProteins(files_dico[genome.expectation_node ], genetic_code)



if __name__ == '__main__':

    # logging.basicConfig(filename='log/genbankparser.log',level=logging.INFO)

    taxonomy_file="data/taxonomy/taxonomy_virus.txt"
    expected_file =  "viruses_w_polyproteins.txt"

    output_dir = "/scratch/polyproteins_interpro/annotated_genomes_by_taxon"

    output_file = 'data/protein_test.faa'

    taxon = 'ssRNA viruses'
    taxon = "Togaviridae"
    taxon = 'Rubivirus'
    taxon = "Alphavirus"
    taxon= "Retroviridae"
    taxon='Coronavirinae'
    taxon="Viruses"
    # taxon='Blueberry virus A'


    taxon_expectation = tax.expectedPeptide(expected_file)

    file_handle = open(os.path.join(output_dir,'polyprotein_All.faa'), "w")


    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)

    # file_handle = open(output_file, "w")
    genome = set()

    for i, gb_dico in enumerate(gbff_iter):
        # print(gb_dico)
        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']

        if i == 50:
            print(i)
            # break

        # writeAnnotatedProtein(gb_file,file_handle, genetic_code, taxon_expectation)
        # print(gb_file)
        # investigateGenomes(gb_file, taxon_expectation)
        extractSpecificPolyprotein(gb_file,file_handle, genetic_code, gff_file, taxon_expectation, window_step, nb_line)
    print(len(genome))
    print(i+1, 'Genome prossessed from taxon', taxon)
