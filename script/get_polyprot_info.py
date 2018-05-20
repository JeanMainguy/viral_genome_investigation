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


def writeInfo(gb_file, genetic_code, gff_file, taxon_expectation, domain_header, writer_domain):

    genome = obj.Genome( gb_file)
    with gzip.open(gb_file, "rt") as handle:
    # with open("/home/user/mainguy/Documents/Data_Analysis/GCF_000885175.1_ViralMultiSegProj39867_genomic_MODIFIED.gbff", "rt") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "genbank")):


            segment = obj.Segment(record, gb_file)
            genome.segments.append(segment)

        genome.getTaxonExpectation(taxon_expectation)
        if not genome.peptide_expectation:
            return
        for segment in genome.segments:
            segment.getMatpeptidesAndPolyproteins()
            segment.checkPeptideRedundancy() #remove the redundant peptide
            segment.checkSubPeptides()
            segment.associatePepWithProt()

            segment.checkForSlippage()
            segment.identifySubProtein()

            segment.getCleavageSites()

            segment.identifyExpectedPolyprotein(genome.polyprotein_expectation)  # genome.identifyExpectedElement()

    genome.getTaxonExpectation(taxon_expectation)
    genome.getMatchObject(gff_file)

    genome.associateMatchWithPolyprotein()

    # print(taxon_expectation)

    dico_genome = genome.getGenomeInfo()

    for m in genome.matchs:
        d = m.get_csv_dico(domain_header)
        d.update(dico_genome)
        writer_domain.writerow(d)
        # if genome.hasEnoughPeptide() or genome.numberOf('peptides'):
        #     genome.visualisation(1, genetic_code)



if __name__ == '__main__':

    # logging.basicConfig(filename='log/genbankparser.log',level=logging.INFO)


    taxonomy_file="data/taxonomy/taxonomy_virus.txt"
    expected_file =  "data/taxonomy/polyprotein_expectation_by_taxon.csv"

    output_dir = '/home/user/mainguy/Documents/Data_Analysis/data/Cleavage_site_sequences'
    gff_file = 'data/interproscan_result/annotated_polyprotein_by_taxon/polyprotein_All.gff3'

    try:
        taxon = sys.argv[1]
    except IndexError:
        taxon = "Viruses"
    taxon_expectation = tax.expectedPeptide(expected_file)

    # file_handle = open(os.path.join(output_dir,'cleavage_site_{}_window_{}.faa'.format(taxon, window_step_clavage_site*2)), "w")

    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)



    genome_header = ['taxon_id', 'organism', "taxon_of_expectation",
                    "expected_number_of_peptide", "peptide", "number_of_final_peptide",
                    "expected_number_reached", "polyprotein_fully_covered", "unannotated_parts", "relevant_annotation"] # genome header


    domain_header = ["matchID", "method" ,"score" ,"Dbxref" , "name","signature_desc" ,
                    "start_in_prot", "end_in_prot", "start", "end", "duplicated",
                    "overlapping", "OverlappingDistance",
                    "left_overlaps_peptide", "right_overlaps_peptide"]

    header = genome_header + domain_header

    fl_domain = open('result/polyproteins_interpro_domains2.csv', "w")

    writer_domain = csv.DictWriter(fl_domain, fieldnames=header, delimiter='\t')
    writer_domain.writeheader()

    for i, gb_dico in enumerate(gbff_iter):
        # print(gb_dico)

        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']
        # print(gb_file)

        if i%200== 0:
            # continue
            print(i)
        # print('genetic code', genetic_code)
        # print(gb_file)
        writeInfo(gb_file, genetic_code, gff_file, taxon_expectation, domain_header, writer_domain)

    print(i+1, 'Genome analysed from taxon', taxon)
