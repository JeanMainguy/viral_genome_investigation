import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
from Bio.Alphabet import generic_protein

def extractPolyMatpep(gb_file, types):
    records = []
    # print('Extracting Polyproteins and mature peptides:')
    with gzip.open(gb_file, "rt") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "genbank")):
            polyprot = []
            matpep = []
            poly_di = {} # dict with poly as key and list of matpep as value
            record_dict = {}

            # if i > 0:
            #     print(">>>")
            #     print(gb_file)
            #     print(record)
            #     input("more than one record")
            organism = record.annotations['organism']
            # print(type(record))

            for feat in record.features:
                # print("feat")
                if feat.type == "mat_peptide" or feat.type == "sig_peptide":
                    # print("mat")
                    matpep.append(feat)
                types.add(feat.type)

                if feat.type == "CDS":
                    for k in feat.qualifiers:
                        if "polyprotein" in feat.qualifiers[k][0]:
                            polyprot.append(feat)
                            break
                    #print("cds")
                    #if 'note' in feat.qualifiers:
                        #if "polyprotein" in feat.qualifiers['note'][0]:
                            #polyprot.append(feat)
                    #if "polyprotein" in feat.qualifiers['product'][0]:


            #Grouping polyprotein with their mat pep
            for poly in polyprot:
                poly_di[poly] = []
                # print("POLY****"*9)
                # see_objet(poly.location)

                for pep in matpep:

                    # print("PEP=="*5)
                    # see_objet(pep.location)
                    if pepBelongsToPoly(pep, poly, record):
                        poly_di[poly].append(pep)


            record_dict["record"] = record
            record_dict["poly_di"] = poly_di
            records.append(record_dict)
        return records, polyprot, matpep, poly_di


def pepBelongsToPoly(pep, poly, record):
    try:
        if pep.qualifiers['gene'] == poly.qualifiers['gene']:
            # print(pep.qualifiers['gene'])
            return True
        else:
            return False
    except KeyError:
        # Start < end always. The coordinate are given according the strand of the gene
        if not (pep.location.start in poly.location and pep.location.end in poly.location and pep.location.strand ==  poly.location.strand):

        # if not (poly.location.start <= pep.location.start < poly.location.end and poly.location.start < pep.location.end <= poly.location.end and pep.location.strand ==  poly.location.strand):
            # print("pep is not in poly..")
            return False


        start_flag = False
        for sub_location in poly.location.parts:
            if pep.location.start in sub_location and pep.location.start%3 == sub_location.start%3:
                # print("start seems to be correct")
                start_flag = True
            if start_flag and pep.location.end in sub_location and pep.location.end%3 == sub_location.end%3:
                # print('PEPTIDE BELONG TO POLYP')
                # print('==record=='*10)
                # print(record)
                # print('==pep=='*20)
                # see_objet(pep)
                # print('==poly=='*20)
                # see_objet(poly)
                # input()
                return True

        return False


    # elif not (poly.location.start => pep.location.start > poly.location.end and poly.location.start > pep.location.end >= poly.location.end)

    # for locprot in poly.location.parts:
    #     for locpep in pep.location.parts:
    #         if locpep.start


#     try:
#         if pep.qualifiers['gene'] == poly.qualifiers['gene']:
#             poly_di[poly].append(pep)
#     except KeyError:
#
# # if hasattr(pep.qualifiers, 'parts') and hasattr(poly.qualifiers, 'parts'):
# for loc in pep.qualifiers.parts:



def write_prot(proteins, file_name, organism):
    ##proteins: list of feature with their
    protseqs = []
    for prot in proteins:
        seq = prot.qualifiers["translation"][0]
        protseqs.append(SeqRecord(Seq(seq,generic_protein) , id=str(prot.qualifiers["protein_id"][0]), description="{} {}".format( prot.qualifiers["product"][0], organism)))

        SeqIO.write(protseqs, file_name,"fasta")



def see_objet(obj):
    print('SEE OBJET ', type(obj))

    for attr in dir(obj):
        if attr[0] != '_' and attr != "seq":
            print(attr, " ", getattr(obj, attr))


def transformePolySeq(record, pol, mat_peptides):
    #Pol is a gb feature of a polyprotein
    #mat_peptides is a sorted list of mat peptide from pol
    current_position = pol.location.start
    pep_draw = ""
    union_sequence = ''
    segmented_seq = []

    # peptide need to be sorted by their start
    # normally they are coming sorted in the gb file
    for pep in mat_peptides:
        print("**"*20)

        if pep.location.start >= current_position:
            if pep.location.start > current_position:
                # add the sequence to join the 2 peptides
                poly_adjust = record[current_position+1:pep.location.start]
                poly_adjust = poly_adjust.seq.translate()
                union_sequence += poly_adjust.upper()
                pep_draw += '{}'.format('-'*(len(poly_adjust)))
                segmented_seq.append(poly_adjust)
                print('ADJUST WITH PROT')

            pep_seq = pep.location.extract(record)
            pep_seq = pep_seq.seq.translate()

            pep_seq = pep_seq[:4].lower() + pep_seq[4:-4].upper() + pep_seq[-4:].lower()
            pep_draw += '<{}>'.format('='*(len(pep_seq)-2))
            union_sequence += pep_seq
            current_position = pep.location.end
            segmented_seq.append(pep_seq)
        else:
            print(pep.location)
            print('start pep', pep.location.start, 'current pos', current_position)
            print(pep.location.start)
            raise Exception('Pep Overlap.. need to take that into account')

    return union_sequence, pep_draw

def gffParser(gff_file):
    sequences_match = {}
    with open(gff_file, "r") as fl:
        for l in fl:
            if l.startswith("#"):
                continue
            elif l[0] == ">":
                break

            (gff_elementid, method, feature_type, begin, end, score, strand, something, comment) = l.split("\t")
            if splited_l[2] == "protein_match" and splited_l[5] != ".":

                match = {"start":int(splited_l[3]), "end":int(splited_l[4]), "app":splited_l[1] }
                sequences_match.setdefault(splited_l[0], [])
                sequences_match[splited_l[0]].append(match)
        return sequences_match

def checkForOverlaping(seq_matchs, peptides):
    for i, p in enumerate(peptides):
        pstart = p.location.start
        pend = p.location.end

        for j, m in enumerate(seq_matchs):
            mstart = m["start"]
            mend = m['end']
            if mstart <  pstart <  mend or mstart <  pend <  mend:
                print('Match {} from {} to {} overflows peptide {}'.format(j, mstart, mend, i+1))
