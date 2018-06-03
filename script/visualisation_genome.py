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



def visualisation(gb_file, genetic_code, gff_file, taxon_expectation, window_step, nb_line, minimum_nb_peptide):

    genome = obj.Genome( gb_file)

    with gzip.open(gb_file, "rt") as handle:
    # with open("test_genome.gbff", "rt") as handle:
    # with open("/home/user/mainguy/Documents/Data_Analysis/GCF_000885175.1_ViralMultiSegProj39867_genomic_MODIFIED.gbff", "rt") as handle:
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

    genome.getTaxonExpectation(taxon_expectation)
    genome.identifyExpectedElement()
    genome.getMatchObject(gff_file)
    genome.associateMatchWithPolyprotein()



    if genome.numberOf("peptides") >= minimum_nb_peptide:
        genome.visualisation(nb_line, genetic_code)

    for cds in segment.cds:
        if cds.number == 1:
            print(cds.number, cds.protein_id)
            print(' sub prot:',[s.number for s in cds.sub_prot])
            print(' peptide:',[(p.number, p.location) for p in cds.peptides])
            print(' unannotated_region:',[(u.number, str(u.bp_obj.location)) for u in cds.unannotated_region])
            for c in cds.cleavage_sites:
                if c.number == 3:
                    print(c)
                    for p in c.peptides:
                        print(p.number)
                    for p in c.proteins:
                        print(p.number)
                    print(c.get_position_prot_relative(cds))

if __name__ == '__main__':
    # logging.basicConfig(filename='log/genbankparser.log',level=logging.INFO)
    try:
        nb_line = sys.argv[2]
    except IndexError:
        nb_line = 1

    try:
        minimum_nb_peptide = int(sys.argv[3])
    except IndexError:
        minimum_nb_peptide = 0

    window_step_clavage_site = 10

    taxonomy_file="data/taxonomy/taxonomy_virus.txt"
    expected_file =  "data/taxonomy/polyprotein_expectation_by_taxon.csv"

    output_dir = '/home/user/mainguy/Documents/Data_Analysis/data/Cleavage_site_sequences'
    gff_file = '/scratch/polyproteins_interpro/interproscan_result/All_polyprotein.gff3'

    try:
        taxon = sys.argv[1]
    except IndexError:
        taxon = 'ssRNA viruses'
        taxon = "Togaviridae"
        taxon = 'Rubivirus'
        taxon = "Marafivirus"
        # taxon= "Togaviridae"
        # taxon='ssRNA positive-strand viruses, no DNA stage'
        # taxon='Alphavirus'
        # taxon="11036"

    taxon_expectation = tax.expectedPeptide(expected_file)

    # file_handle = open(os.path.join(output_dir,'cleavage_site_{}_window_{}.faa'.format(taxon, window_step_clavage_site*2)), "w")

    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)

    print('*-'*100)
    print('VISUALISATION OF THE RefSEq GENOME FROM THE TAXON {} THAT HAVE AT LEAST {} ANNOTATED PEPTIDE'.format(taxon,minimum_nb_peptide ))
    print('*-'*100)
    for i, gb_dico in enumerate(gbff_iter):
        # print(gb_dico)

        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']
        # print(gb_file)

        if i%200== 0:
            # continue
            print(i)
        # print('genetic code', genetic_code)
        print(gb_file)
        visualisation(gb_file, genetic_code, gff_file, taxon_expectation, window_step_clavage_site, nb_line, minimum_nb_peptide)

    print(i+1, 'Genome analysed from taxon', taxon)
