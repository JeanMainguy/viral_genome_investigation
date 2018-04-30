import function_parser as fct
import taxonomy as tax
import gb_parser
import os, re, subprocess, gzip, logging, csv
from operator import attrgetter

logging.basicConfig(filename='log/interproscan_parser.log')

def colorString(category, string):
    COLOR = {"peptide0":'\033[94m',"peptide1":'\033[34m', "peptide2":'\033[35m', "no_peptide_region":'\033[93m', "border":'\033[31m', "match2":'\033[35m',
                "match0":'\033[31m', "match4":'\033[33m', "match1":'\033[32m',"CDD":'\033[34m',"ProDom":'\033[35m' }

    # OKBLUE = '\033[94m'
    # OKGREEN = '\033[92m'
    # WARNING = '\033[93m'
    # FAIL = '\033[91m'
    ENDC = '\033[0m'
    # BOLD = '\033[1m'
    # UNDERLINE = '\033[4m'

    if category in COLOR:
        return COLOR[category]+string+ENDC
    elif category == "peptideborder":
        return COLOR["border"]+string[:3]+ COLOR["peptide"]+string[3:-3]+COLOR["border"]+string[-3:]+ENDC
    return string

def see_objet(obj):
    print('SEE OBJET ', type(obj))

    for attr in dir(obj):
        if attr[0] != '_' and attr != "seq":
            print(attr, " ", getattr(obj, attr))


def polyproteinVisualisation(record, pol, mat_peptides, matchs):
    #Pol is a gb feature of a polyprotein
    #mat_peptides is a sorted list of mat peptide from pol
    organism = record.annotations['organism']
    pep_draw = ""
    segmented_seq = []

    polyprot_seq = pol.location.extract(record).seq.translate()
    pep_draw = ' '*len(polyprot_seq)


    pep_draws = [pep_draw]

    # for unannotated in pol.qualifiers['unannotated_region']:


    #Build PEP draw string
    for pep in mat_peptides:
        draw_indice = 0
        # print("==============")
        pep_start = int(pep.qualifiers['start_aa'])
        pep_end = int(pep.qualifiers['end_aa'])
        pep_seq = pep.location.extract(record).seq.translate()
        if not pep_seq == polyprot_seq[pep_start-1:pep_end]:
            # print('NEIN')
            pep_start = pep.qualifiers['alternative_start_aa']
            pep_end = pep.qualifiers['alternative_end_aa']
            if not pep_seq == polyprot_seq[pep_start-1:pep_end]:
                logging.warning('peptide sequence and polyprotein sequence are not similar in : "{}"'.format(organism))
                print("PROBELEM WITH SEQUENCE")
                print(pep_seq)
                print(polyprot_seq[pep_start-1:pep_end])
                input()

                continue

        # pep_seq = colorString('peptideborder', pep_seq)
        # print('draw indice', draw_indice)
        while pep_draws[draw_indice][pep_start-5:pep_end+5] != " "*len(pep_draws[draw_indice][pep_start-5:pep_end+5]):
            # print('--draw indice', draw_indice)
            if len(pep_draws) == draw_indice+1:
                pep_draws.append(' '*len(polyprot_seq))
                draw_indice += 1
                break
            draw_indice += 1

        pep_draw = pep_draws[draw_indice]
        pep_draw = pep_draw[:pep_start-1] + pep_seq + pep_draw[pep_end:]
        pep_draws[draw_indice] = pep_draw

    #Build Match draw string
    symbol = "="
    match_draw = ' '*len(polyprot_seq)
    match_draws = [match_draw]

    for i, m in enumerate(matchs):
        draw_indice = 0
        m_len = m["end"] - m['start'] +1
        start = m['start']-1
        end = m['end']-1
        while symbol in match_draws[draw_indice][start:end+1]:
            if len(match_draws) == draw_indice+1:
                match_draws.append(' '*len(polyprot_seq))
                draw_indice += 1
                break
            draw_indice += 1
        match_draw = match_draws[draw_indice]

        string = '{} {}:{} '.format(m['matchID'],m['app'],m['Name'])
        symbol_len = m_len - len(string)
        if symbol_len < 0:
            string = string[:symbol_len] # remove character to the string because it is too large
        elif symbol_len%2 == 0:
            string = int(symbol_len/2)*symbol + string + int(symbol_len/2)*symbol

        else: # odd length
            string = int(symbol_len/2)*symbol + string + (int(symbol_len/2)+1)*symbol # fct int rounds down so we had 1 to get the proper len


        match_draw = match_draw[:start] + string + match_draw[end+1:]
        match_draws[draw_indice] = match_draw
        # print(colorString('peptideborder', pep_seq))
    step = 120
    # print("len draws", len(pep_draws))

    print(display_prot_info(pol, mat_peptides, matchs))
    # return
    for i in range(0, len(polyprot_seq), step):

        print ('.'*step)
        print()
        for k, match_draw in enumerate(match_draws):
            if symbol in match_draw[i:i+step]:
                print(colorString('match'+str(k%3), match_draw[i:i+step]))

        print(polyprot_seq[i:i+step])

        for j, pep_draw in enumerate(pep_draws):
            # if not peptides
            print(colorString('peptide'+str(j%3), pep_draw[i:i+step]))

        print('\033[0m\n')


def display_prot_info(pol, mat_peptides, matchs):
    covered_peptide = pol.qualifiers['covered_peptide']
    uncovered_peptide = pol.qualifiers['uncovered_peptide']
    overflowed_peptide = pol.qualifiers['overflowed_peptide']

    only_covered = covered_peptide - overflowed_peptide
    only_overflowed = overflowed_peptide - covered_peptide
    convered_overflowed = covered_peptide & overflowed_peptide



    info_string = "" #"="*200+"\n"
    info_string += "{} | {} | {}\n".format(pol.qualifiers['protein_id'][0], pol.qualifiers.setdefault('product', [''])[0], pol.qualifiers.setdefault('note', [''])[0])
    info_string = info_string if len(info_string)<205 else info_string[:200]+"...\n"
    # print(pol)
    # print('covered pep', len(covered_peptide),"\nuncovered_peptide",len(uncovered_peptide) ,"\noverflowed_peptide",len(overflowed_peptide) )

    info_string += "length:{}aa peptides:{} domains:{}\n".format(len(pol)/3, len(pol.qualifiers["mat_peptide"]), len(pol.qualifiers["domains"]) )
    info_string += "Peptide:{}  only covered:{} convered and overflowed {} | only overflowed:{} | uncovered {}\n".format(len(pol.qualifiers["mat_peptide"]), len(only_covered), len(convered_overflowed), len(only_overflowed), len(uncovered_peptide)  )

    info_string += "Percentage of coverage:{} | unannotated regions:{} \n".format(pol.qualifiers['peptide_coverage']/(len(pol)/3 -1)*100, len(pol.qualifiers['unannotated_region']))


    return info_string

def getStatisticFromInterpro(taxid, d_match, counter, writer_match, expected_number, taxon, writer_prot):
    """

    """

    gb_file = next(tax.getAllRefseqFromTaxon(taxid))
    records = gb_parser.extractPolyMatpep(gb_file)
    poly_di = {}
    translationflag = False

    #OUTPUT of match result
    # match_fl = open("match_result.csv", 'w')
    # header =['start', 'end', 'app', 'seqid', 'taxid', 'matchID', 'signature_desc', 'Name', 'Dbxref', 'polyprotein', 'organism', 'included_by', 'overflow_pep', 'overflow_info',  "included", 'organism_id', 'overflow', 'overflow_distance']
    # writer_match = csv.DictWriter(match_fl, fieldnames=header, delimiter='\t')
    # writer_match.writeheader()


    for d_record in records:
        # poly_di.update(d_record['poly_di'])
        poly_di = d_record['poly_di']

        ##VISUALIsation
        for polyprot, mat_peptides in poly_di.items():
            polyprot.qualifiers["Expected peptides number"] = True if len(mat_peptides) >= expected_number else False
            protid = polyprot.qualifiers["protein_id"][0]
            prot_match = [] if protid not in d_match else d_match[protid]
            fct.checkForOverlaping(prot_match, poly_di[polyprot], polyprot, counter)

            polyproteinVisualisation(d_record['record'], polyprot, mat_peptides, prot_match)
            # print(display_prot_info(polyprot, mat_peptides, prot_match))
            # writeMatchInfo(prot_match, writer_match, taxon)
            # writePolyproteinInfo(polyprot, mat_peptides, writer_prot, taxon, counter)
            # print(prot_match[0].keys())

    # fct.transformePolySeq(records, polyprot, poly_di[polyprot])
def writeMatchInfo(matchs, writer, taxon):
    for match in matchs:
        match['taxon'] = taxon
        match.pop('included_by', None)
        match.pop('overflow_pep', None)
        match.pop('overflow_info', None)
        # match = {key:str(value) for key, value in match.items()}
        writer.writerow(match)

def writePolyproteinInfo(pol, mat_peptides, writer, taxon, counter):
    dico_csv = {}
    covered_peptide = pol.qualifiers['covered_peptide']
    uncovered_peptide = pol.qualifiers['uncovered_peptide']
    overflowed_peptide = pol.qualifiers['overflowed_peptide']

    dico_csv['total_peptide'] = len(mat_peptides)
    dico_csv['only_covered'] = len(covered_peptide - overflowed_peptide)
    dico_csv['only_overflowed'] = len(overflowed_peptide - covered_peptide)
    dico_csv['convered_and_overflowed'] = len(covered_peptide & overflowed_peptide)
    dico_csv['uncovered'] = len(uncovered_peptide)

    # counter = {k:v+dico_csv[k] for k,v in counter.items()}
    for k, v in counter.items():
        counter[k] += dico_csv[k]


    dico_csv['taxon'] = taxon
    dico_csv["organism"] = pol.qualifiers['organism'][0]
    dico_csv["organism_id"] = pol.qualifiers['taxid'][0]
    dico_csv["Expected_peptides_number"] =  pol.qualifiers["Expected peptides number"]
    dico_csv["Polyprotein_Fully_covered"] = True if pol.qualifiers['peptide_coverage'] == len(pol)/3 - 1 else False
    dico_csv['protein_id'] = pol.qualifiers["protein_id"][0]

    # print(dico_csv.keys())
    writer.writerow(dico_csv)


def writeCounter(counter_dico, taxon, fl_result, headers):
    string = taxon

    for k in headers:
        if isinstance(counter_dico[k], dict):
            string += '\t'
            for app, nb in counter_dico[k].items():
                string += "{} {} |".format(app, nb)
        else:
            string += '\t' + str(counter_dico[k])
    string += '\n'
    fl_result.write(string)

def expectedPeptide(expected_pep_file= "viruses_w_polyproteins.txt"):
        expected_pep_file = "viruses_w_polyproteins.txt"
        taxon_expectation = {}
        nb_correct_poly = 0
        with open(expected_pep_file, "r") as fl:
            for l in fl:
                (taxon, expected_nb) = l.split("\t")
                expected_nb = int(expected_nb.rstrip())
                if expected_nb > 0:
                    taxon_expectation[taxon] = expected_nb
        return taxon_expectation


if __name__ == '__main__':
    gff_outdir = "/home/user/mainguy/Documents/Data_Analysis/data/interproscan_result/annotated_genomes_by_taxon"
    # gff_outdir = "/home/user/mainguy/Documents/Data_Analysis/data/interproscan_result/test_Nidovirales"

    taxon_expectation = expectedPeptide()


    # header = ["total_polyprotein","total_peptide",  "covered_peptide", "uncovered_peptide", 'overflow_match_peptide', "included_method", "overflow_method", "total_method"]
    # poly_result_fl.write('Taxon'+'\t'+'\t'.join(header)+'\n')
    #"Peptide:{}  only covered:{} convered and overflowed {} | only overflowed:{} | uncovered {}\n"
    counter = {"total_peptide":0,'only_covered':0, "convered_and_overflowed":0, "only_overflowed":0, "uncovered":0}


    patern_taxon = re.compile("polyprotein_([\d\w]*).gff3")

    input_files = os.listdir(gff_outdir)

    ##MATCH OUTPUT
    match_fl = open("match_result.csv", 'w')
    header_match =['start', 'end', 'app', 'seqid','taxon',  'taxid', 'matchID', "Polyprotein_Fully_covered", 'signature_desc', 'Name', 'Dbxref',  'organism', 'included_by', "Expected_peptides_number",   "included", 'organism_id', 'overflow', 'overflow_distance']
    writer_match = csv.DictWriter(match_fl, fieldnames=header_match, delimiter='\t')
    writer_match.writeheader()

    ## Protein output
    prot_fl = open("polyprotein_result.csv", 'w')
    header_prot =['total_peptide', 'only_covered', 'only_overflowed', 'convered_and_overflowed', 'uncovered', 'taxon', 'organism', 'organism_id', 'Expected_peptides_number', 'Polyprotein_Fully_covered', 'protein_id']
    writer_prot = csv.DictWriter(prot_fl, fieldnames=header_prot, delimiter='\t')
    writer_prot.writeheader()


    for gff_file in input_files[:1]: #[:1]:
        if not gff_file.endswith('.gff3'):
            continue
        rematch = patern_taxon.search(gff_file)
        if rematch:
            print(rematch[1])
            taxon = rematch[1]
        else:
            taxon = gff_file
        expected_number = taxon_expectation[taxon.replace('_', ' ')]
        gff_file = os.path.join(gff_outdir, gff_file)
        match_by_taxid = fct.gffParser(gff_file)

        # counter = dict(original_counter)
        for taxid, d_match in match_by_taxid.items():

            getStatisticFromInterpro(taxid, d_match, counter, writer_match, expected_number, taxon, writer_prot)


        # writeCounter(counter, taxon, poly_result_fl, header)
        #
        #     # print(counter)
        # for k, item in counter.items():
        #     if isinstance(item, dict):
        #         for k2, item2 in item.items():
        #             global_counter[k].setdefault(k2, 0)
        #             global_counter[k][k2] += item2
        #     else:
        #         global_counter[k] += item
        #
        # print(global_counter['total_peptide'])
        # writeCounter(global_counter, 'Total', poly_result_fl, header)

    print(counter)
    # print(global_counter)
    with open('Global_counter.txt', 'w') as fl:
        for k, v in counter.items():
            fl.write(k+':'+str(v)+'\n')
    prot_fl.close()
    match_fl.close()
