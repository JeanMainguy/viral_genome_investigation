#!/usr/bin/env python3

import taxonomy as tax
import os, gzip, logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys, csv, re, collections
from Bio.Alphabet import generic_protein
from Bio.SeqFeature import SeqFeature, FeatureLocation
from operator import attrgetter

SCREEN_SIZE = 200


class Genome:
    def __init__(self, gb_file):
        self.gb_file = gb_file
        self.segments = []
        self.matchs = []
        self.taxonomy = None
        self.organism = None
        self.taxon_id = None
        self.expectation_node = None
        self.peptide_expectation = None
        self.polyprotein_expectation = None
        self.variable_polyprotein_expectation = None
        self.expectation_info = {}
        Protein.COUNTER = 0

    def __len__(self):
        return sum([len(s) for s in self.segments])


    def numberOf(self, attribute):
        return sum([len(getattr(s, attribute)) for  s in self.segments])


    def extractObjectFromGbff(self, gb_file):

        with gzip.open(gb_file, "rt") as handle:
            for i, record in enumerate(SeqIO.parse(handle, "genbank")):
                segment = Segment(record)
                genome.segments.append(segment)


    def hasEnoughPeptide(self):

        total_final_peptides = sum([len(segment.peptides)-len(segment.parent_peptides) for segment in self.segments])
        return True if self.peptide_expectation <=  total_final_peptides else False


    def getTaxonExpectation(self, taxon_expectation):
        """
        We presume that the taxonomy is homogenous in the differents segments of a genome
        For one genome it is not the case but we ignore this case
        """
        # print('TAXON ....')
        if not self.segments:
            logging.warning( 'ERROR no segment stored in genome')
            return 'ERROR no segment stored in genome'
        segment = self.segments[0]
        record = segment.record


        taxonomy = record.annotations['taxonomy']
        self.taxonomy = taxonomy
        self.organism = record.annotations['organism']
        self.taxon_id = segment.taxon_id
        # print("TAXONOMY")
        # print(taxonomy)
        self.expectation_node='Not_Found'
        # self.peptide_expectation = 1
        for node in taxonomy:

            if node in taxon_expectation:
                self.expectation_node = node
                self.peptide_expectation = taxon_expectation[node]['peptide']
                self.polyprotein_expectation = taxon_expectation[node]['polyprotein']
                self.variable_polyprotein_expectation = taxon_expectation[node]['variable_polyprotein_expectation']
                self.expectation_info = taxon_expectation[node]
        # print(len(self.segments)," == ",self.expectation_info['segment'])
        if not  self.expectation_info: # if the genome taxonmy has no match with any taxon where we expect peptides then we stop here
            return
        if len(self.segments) > 0 and len(self.segments) == self.expectation_info['segment']:
            for i in range(self.expectation_info["segment"]):
                # BIG approximation here! the segment in the genbank file need to be in the same order all the time
                # WARNING

                nb_polyprotein_key = 'segment{}_polyprotein'.format(i+1)
                nb_peptide_key = 'segment{}_peptide'.format(i+1)
                # print(nb_polyprotein_key,  self.expectation_info[nb_peptide_key])
                self.segments[i].polyprotein_expectation = self.expectation_info[nb_polyprotein_key]
                self.segments[i].peptide_expectation = self.expectation_info[nb_peptide_key]
                # print( self.segments[i].polyprotein_expectation)
        # print('expectation node', self.expectation_node)
        # print(taxonomy)


    def getGenomeInfo(self):
        """
        genome_header = ['taxon_id', 'organism', "taxon_of_expectation",
                    "expected_number_of_peptide", "peptide", number_of_final_peptide",
                    "expected_number_reached", "polyprotein_fully_covered", "unannotated_parts", "relevant_annotation"] # genome header
        """

        csv_dico = {}

        csv_dico['taxon_id'] = self.taxon_id
        csv_dico['organism'] = self.organism
        csv_dico['taxon_of_expectation'] = self.expectation_node

        csv_dico['expected_number_of_peptide'] = self.peptide_expectation

        csv_dico['peptide'] = self.numberOf('peptides')
        csv_dico['number_of_final_peptide'] = csv_dico['peptide'] - self.numberOf('parent_peptides')


        csv_dico['expected_number_reached'] = True if csv_dico['number_of_final_peptide'] >= self.peptide_expectation else False
        csv_dico['unannotated_parts'] = self.numberOf('unannotated_region')
        csv_dico['polyprotein_fully_covered'] = False if csv_dico['unannotated_parts'] > 0 else True

        csv_dico['relevant_annotation'] = self.isGenomeAnnotationRelevant()



        return csv_dico


    def visualisation(self, nb_line=1, genetic_code=None):
        '''
        Visualisation of the genome in a friendly way
        Example:
        Segment 1:
        ------------------------------------------------------
        ======>   ======================>               ====>
                  ==========<-1====================>
        Segment 2:
        ----------------------------------------
        ======>   ======================>
        '''
        print('GENOME of {}|{}'.format(self.organism, self.taxon_id))
        longuest_segment_len = max([len(s) for s in self.segments])

        conversion = longuest_segment_len/SCREEN_SIZE
        for i, segment in enumerate(self.segments):
            print('SEGMENT {}/{}: {} cds | {} polyprotein | {}/{} final mat_peptides'.format(i+1,
                len(self.segments),
                len(segment.cds),
                len(segment.polyproteins),
                len(segment.peptides)-len(segment.parent_peptides),
                len(segment.peptides)))
            print('taxonomy:{} \nnode giving the nb of peptides {}'.format(self.taxonomy, self.expectation_node))
            print('Expected number of mature peptide {} '.format(self.peptide_expectation))
            # for p in segment.cds:
            #     print('Is annotation relevant? ',p.polyprotein_number, p.isAnnotationRelevant())
            segment.visualisation(nb_line, genetic_code)
            # print(features_strg)


    def getMatchObject(self, gff_file):
        #retreive all the match that belong to the genome
        flag = False
        with open(gff_file) as fl:
            for l in fl:

                if l.startswith(">") or l.startswith("##FASTA"):
                    break
                if l.startswith('##sequence-region'):
                    flag = True if l.startswith('##sequence-region {}|{}'.format(self.expectation_node,self.taxon_id)) else False
                    continue
                if flag:
                    (taxseqid, method, feature_type, start, end, score, strand, phase, attributes) = l.split("\t")
                    #match = {"start":int(begin), "end":int(end), "app":method, "seqid":seqid, "taxid":taxid}
                    if feature_type == 'protein_match':
                        (expectation_node, taxid, seqid, polyprotein_number) = taxseqid.split("|")
                        re_result = re.search("ID=match\$([\d]+)_[\d]+_[\d]+", attributes)
                        matchID = "" if not re_result else int(re_result.group(1))

                        # signature_desc
                        re_result = re.search("signature_desc=([^;]+)", attributes)
                        signature_desc= "" if not re_result else re_result.group(1)

                        # Domain name
                        re_result = re.search("Name=([\d\w]+)", attributes)
                        Name = "" if not re_result else re_result.group(1)

                        # Dbxref
                        re_result = re.search('Dbxref="([^"]+)"', attributes)
                        Dbxref= "" if not re_result else re_result.group(1)

                        self.matchs.append(Match(taxid, seqid,  method, start, end, score, Dbxref, Name, matchID, signature_desc))


    def associateMatchWithPolyprotein(self):
        for segment in self.segments:
            for poly in segment.polyproteins:
                # poly.matchs = [match for match in self.matchs if match.seqid == poly.protein_id]
                for match in self.matchs:
                    if match.seqid == poly.protein_id:
                        poly.matchs.append(match)
                        match.protein = poly
                        match.getGenomicPositions(poly.start)
                        segment.matchs.add(match)
            segment.getDomainOverlappingInfo()
            segment.identifyDuplicatedMatch()


    def identifyExpectedElement(self):
        #The question is when their is more than 1 segment with more than 1 polyprotein...
        for segment in self.segments:
            segment.identifyExpectedPolyprotein(self.polyprotein_expectation, self.variable_polyprotein_expectation)


    def isGenomeAnnotationRelevant(self):
        return False if any([s for s in self.segments if not s.isSegmentAnnotationRelevant() ]) else True


class Segment:

    POSSIBLE_TYPE = set()

    def __init__(self, record, gb_file):
        self.record = record
        self.gb_file = gb_file
        self.source = None
        self.polyproteins = set() # polyproteins list
        self.parental_polyproteins = []
        self.peptides = set() # mat_peptides list
        self.matchs = set()
        self.parent_peptides = set()
        self.sub_peptides = set()
        self.cds = set()
        self.unannotated_region = set()
        self.taxon_id = None
        self.organism = record.annotations['organism']
        self.cleavage_sites = [] # list of cleavage_site object
        self.relevant_annotation = False

        self.polyprotein_expectation = None



    def __len__(self):
        return len(self.source)


    def identifyExpectedPolyprotein(self, polyprotein_expectation, variable_polyprotein_expectation):
        # print("self.polyprotein_expectation", self.polyprotein_expectation)
        if self.polyprotein_expectation is not None:
            # print("self.polyprotein_expectation", self.polyprotein_expectation)
            polyprotein_expectation = self.polyprotein_expectation
            # print("polyprotein_expectation", polyprotein_expectation)
        if not polyprotein_expectation:
            # print(self.polyprotein_expectation)

            return
        # if the segment has it own polyprotein_expectation then we take them as reference


        # print(polyprotein_expectation, 'Expected polyprotein in this genome')
        parental_prot =  [ prot for prot in self.cds if not prot.parental_prot] # if this prot is not included in another one
        if variable_polyprotein_expectation and len(parental_prot) > 1:
            parental_prot = sorted(parental_prot, key=lambda x: len(x), reverse=False)
            length_difference = [len(parental_prot[i]) - len(parental_prot[i-1]) for i in range(len(parental_prot))[1:]]
            index_boundary  = length_difference.index(max(length_difference))
            self.parental_polyproteins = parental_prot[index_boundary+1:]


        elif len(parental_prot) < polyprotein_expectation:
            logging.warning('Not enough parental protein found {}/{}'.format(len(parental_prot), polyprotein_expectation))
            self.parental_polyproteins = parental_prot
        elif len(parental_prot) == polyprotein_expectation:
            self.parental_polyproteins = parental_prot
        elif len(parental_prot) > polyprotein_expectation:# rule: take the biggest protein until the number match the expectation and consider them as polyprotein
            self.parental_polyproteins = sorted(parental_prot, key=lambda x: len(x), reverse=True)[:polyprotein_expectation]

        self.parental_polyproteins = sorted(self.parental_polyproteins, key=lambda x: x.start, reverse=False)
        #giving a positional number to the polyprotein
        #if the polyprotein is included it get the positional number of the parent + a decimal number between 0.1


        # self.polyproteins = []
        for i, p_poly in enumerate(self.parental_polyproteins):
            p_poly.polyprotein_number = i+1.0
            self.polyproteins.add(p_poly)
            # print(p_poly)
            for j, sub_poly in enumerate(sorted(p_poly.sub_prot, key=lambda x: len(x), reverse=True)):

                sub_poly.polyprotein_number = (i+1)+(j+1)/10
                self.polyproteins.add(sub_poly)
                # print(sub_poly)
                if j > 9: logging.error('polyprotein_number is not working') # subpolyprotein > 9


    def getSubProteins(self):
        return [p for p in self.cds if p.parental_prot ] # return protein that have parental prot


    def getPolyWithFeature(self, feature):
        return [p for p in self.polyproteins if getattr(p, feature)]


    def visualisation(self, nb_line, genetic_code):

        try:
            display_len = int(nb_line)*SCREEN_SIZE
            conversion = int(len(self.record)) / display_len
        except ValueError:
            display_len = int(len(self.record)/3) if nb_line == 'aa' else int(len(self.record))
            conversion = 3 if nb_line == 'aa' else 1



        compatible_dico = collections.OrderedDict()
        compatible_dico['match'] = self.buildCompatibleGroup(self.matchs)[::-1]
        compatible_dico['cds'] = self.buildCompatibleGroup(self.cds)
        compatible_dico['pep'] = self.buildCompatibleGroup(self.peptides)
        compatible_dico['unannotated_region'] = self.buildCompatibleGroup(self.unannotated_region)

        strings = []
        seq_type_dico = {
            'match': ('[', '#', ']' , ' '),
            'cds': ('=', '=', '>', '-') ,
            'pep':('|', '+', '|' , ' '),
            'unannotated_region':('|', '~', '|' , ' ')
            }
        color_end = '\033[0m'
        overlap_col = '\033[91m'

        # color = '$'
        # color_end= '$'
        for seq_type, compatible_group in compatible_dico.items():
            left, central, right, neutral = seq_type_dico[seq_type] # different symbol if it protein or peptide
            for group in compatible_group:
                # print('GROUPPP')
                string = neutral*display_len
                # print(string)

                for seq in group:
                    i_start, i_end = seq.getStringIndices(conversion)

                    string = string[:i_start] + left+central*(i_end - i_start-1)+ right + string[i_end+1:]

                    if seq.__class__.__name__ ==  'Protein':

                        # print(seq.getSequenceRecord(self.organism, self.taxon_id, self.record, genetic_code).seq)
                        if  nb_line == 'aa':
                            seqaa = seq.getSequenceAA(self.record, genetic_code)

                            string = string[:i_start] + seqaa + string[i_start+len(seqaa):]
                            continue
                        if  nb_line == 'nt':
                            seqaa = seq.getSequenceAA(self.record, genetic_code)
                            seqaa = ''.join([a*3 for a in seqaa])
                            string = string[:i_start] + seqaa + string[i_start+len(seqaa):]
                            continue

                        #display positional number
                        # if seq.polyprotein_number:
                        #     nb = str(seq.polyprotein_number)
                        #     string = string[:i_start+2] +nb+ string[i_start+2+len(nb):]
                        #
                        #Display protein number
                        nb = str(seq.number)
                        string = string[:i_start+2] +nb+ string[i_start+2+len(nb):]

                        if seq.ribosomal_slippage:
                            for position, shift in seq.ribosomal_slippage.items():
                                string_postion = round(position/conversion)
                                string_postion = string_postion if string_postion<i_end else string_postion-1
                                shift_string = '+' if shift > 0 else "-"
                                shift_string += str(abs(shift))

                                string = string[:string_postion-len(shift_string)+1] + shift_string + string[string_postion+1:]

                    else:
                        name = seq.getNameToDisplay()
                        string = string[:i_start+2] +name+ string[i_start+2+len(name):]
                # if seq_type == "match" and nb_line == 1:
                    # print(string)

                    # string = re.sub(r"(\[#<[A-Za-z0-9]+#+\]?)", r"{}\1{}".format(color, color_end), string)
                # print(string)
                strings.append(string)




        ##REGEX
        # overlap_match re.compile(r"\[#<"g)
        # match = re.compile(r"\[#[^<]"g)
        # sub_string = re.sub(r"(\[#<[A-Za-z0-9]+#+\]?)", r"{}\1{}".format(color, color_end), sub_string)

        unfinished_col_list = ["" for s in strings]
        for i in range(0, display_len, SCREEN_SIZE):
            print('/'*SCREEN_SIZE)
            print()
            for i_line, string in enumerate(strings):
                sub_string = string[i:i+SCREEN_SIZE]
                if unfinished_col_list[i_line]:
                    sub_string = re.sub(r"^([^ \[]+)", r"{}\1{}".format(unfinished_col_list[i_line], color_end), sub_string)
                sub_string = re.sub(r"(\[?#?<[A-Za-z0-9]+#*\]?)", r"{}\1{}".format(overlap_col, color_end), sub_string)

                ansi_list = re.findall(r"(\x1B\[[0-?]*[ -/]*[@-~])" ,sub_string[:-1]) # -1 not not catch the end color character of the end of the line

                unfinished_col_list[i_line] = '' if not ansi_list or ansi_list[-1] == color_end else ansi_list[-1]
                # if ansi_list and ansi_list[-1] != color_end:
                #     unfinished_col_dict[i_line] = ansi_list[-1]
                print(sub_string)

        # protein = '\n'.join(strings)

        print(color_end)
        # print(protein)
        # return 'protein'


    def buildCompatibleGroup(self, set_of_sequence):
            sequences = set_of_sequence.copy() #sorted(list(set_of_sequence), key=lambda x: x.bp_obj.location.start, reverse=False)
            #Build the set of non overlaping protein

            for i, seq in enumerate(sequences):

                for seq_next in list(sequences)[i+1:]:

                    if not seq.overlap(seq_next):
                        seq.non_overlapping_prot.add(seq_next)
                        seq_next.non_overlapping_prot.add(seq)

            #Build list of sequences that don't overlap
            compatible_groupes = []
            used_seq = []
            while len(sequences) > 0:
                seq = sequences.pop()
                group = [seq]
                # finding the sequence that are compatible with seq and within each other
                group_compatible = seq.non_overlapping_prot # set of sequence that are compatible with the sequence stored in the list group
                while group_compatible:
                    try:
                        compatible_seq = (group_compatible & sequences).pop() #extarct one seq that is compatibe with group
                    except KeyError:
                        break
                    group_compatible = group_compatible & compatible_seq.non_overlapping_prot & sequences # build new set 'group_compatible' that store compatible set

                    group.append(compatible_seq)
                    sequences.remove(compatible_seq)

                compatible_groupes.append(sorted(group, key=lambda x: x.end, reverse=True))
            return compatible_groupes


    def getMatpeptidesAndPolyproteins(self):
        for feat in self.record.features:
            Segment.POSSIBLE_TYPE.add(feat.type)

            if feat.type == 'source':
                self.taxon_id = feat.qualifiers['db_xref'][0].replace('taxon:', '')
                self.source = feat

            elif feat.type == "mat_peptide" or feat.type == "sig_peptide" or feat.type == "proprotein":
                self.peptides.add(Peptide(feat))
                # self.peptides.append(feat)
            elif feat.type == 'CDS':
                prot_obj = Protein(feat)
                self.cds.add(prot_obj)
                # Check if polyprotein appears in qualifiers maybe a bit risky..
                # if any(k for k in feat.qualifiers if "polyprotein" in feat.qualifiers[k][0]):
                #     self.polyproteins.add(prot_obj)
                #     prot_obj.polyprotein = True
                #     prot_obj.reasons.add('"polyprotein" appears in qualifiers')


    def associatePepWithProt(self):
        #Should be changed at the end to catch polyprot that are not annotated
        #Here we want to find the polyprot that have a mat peptide so we are pretty sure that thz protein is a polyprotein
        # EDIT: The protein that has some peptide is no longer consider as polyprotein
        # To be consider as polyprotein the cds undergoes the identifyExpectedPolyprotein() function
        # No Anymore the case

        for cds in self.cds:
            # print('CDS', cds)
            for pep in self.peptides:
                # print('  PEP', pep)
                if pep in cds:
                    cds.peptides.add(pep)
                    pep.polyproteins.add(cds)

            if cds.peptides:
                cds.status = 'Polyprotein'
                cds.polyproteinCoverage()
                self.unannotated_region.update(cds.unannotated_region)
                # cds.getCleavageSites()

        ##SMALL check up to be sure that every peptides have been assigned to at least one protein
        for pep in self.peptides:
            if not pep.polyproteins:
                logging.warning("The peptide {} {} was not assigned to any cds in {}".format(pep.number, str(pep.bp_obj.location), self.taxon_id))


    def checkSubPeptides(self):
        #Check if peptide is included in a bigger peptide
        for i, pep in enumerate(self.peptides):
            for pep_next in self.peptides:
                #try to determine if pep is included in pep_next
                if pep == pep_next:
                    continue

                if pep_next.start <= pep.start <= pep_next.end and pep_next.start <= pep.end <= pep_next.end:

                    if not(pep_next.start%3 == pep.start%3 and pep.end%3 == pep_next.end%3):
                        # print("PROBLEM/////")
                        # print('\n', pep.bp_obj,'\nin\n', pep_next.bp_obj, '\n' )
                        logging.warning('Peptide {} is included in peptide {} but does not share the same strand'.format(pep.number, pep_next.number))
                        # print(' \ \ \ \ \ ')
                        continue

                     # and pep_next.start%3 == pep.start%3 and pep.end%3 == pep_next.end%3:
                    # print('\n', pep.bp_obj,'\nin\n', pep_next.bp_obj, '\n' )
                    self.parent_peptides.add(pep_next)

                    pep_next.parent_peptide = True
                    self.sub_peptides.add(pep)
                    # print(pep_next)
                    # print(pep)


    def checkPeptideRedundancy(self):
        #Remove Peptide that have similat start and end...
        peptides = list(self.peptides)

        for i, pep in enumerate(peptides):

            for pep_next in peptides[i+1:]:

                if pep.location.start ==  pep_next.location.start and pep.location.end ==  pep_next.location.end:

                    self.peptides.remove(pep_next)
                    pep.redundant_pep.append(pep_next)


    def checkForSlippage(self):
        """
        Find the protein that are included in another protein due to ribosomal_slippage are alternative start
        ribosomal_slippage attr is a dico with as a key the position of the end of the part of the protein and as a key the shift
        """

        for poly in self.cds:
            # print(poly)
            # print(poly.bp_obj.location)
            if len(poly.bp_obj.location.parts) == 1:
                # print('PROTEIN HAS ONLY ONE PART')
                continue
            first_part = poly.bp_obj.location.parts[0]
            # print(first_part)
            for part_next in poly.bp_obj.location.parts[1:]:
                # print(part_next)
                # print('CHECK NEXT PART')
                shift = part_next.start - first_part.end
                # print('SHIFT of ',shift)

                if abs(shift) < 3:

                    poly.ribosomal_slippage[first_part.end] = shift
                elif shift == 3:
                    poly.readthrough[first_part.end] = shift
                    logging.info('Shift of 3.. Should be a readthrough but better to check {}'.format(self.gb_file))
                else:

                    logging.info('shift of {}, probably a spliced gene  {}'.format(shift, self.gb_file))

                first_part = part_next
            # if len(poly.bp_obj.location.parts) > 2:
            #     print('PROT WITH MORE THAN 2 PART')




            # self.ribosomal_slippage = True if len(poly.bp_obj.location.parts)>1 else False



            # for poly_next in poly_list[i+1:]:

                # if poly.bp_obj.location.end == poly_next.bp_obj.location.end and poly.bp_obj.location.start != poly_next.bp_obj.location.start:
                    # self.alternative_start = True


    def identifySubProtein(self):
        list_prot = list(self.cds)

        for i, prot in enumerate(list(list_prot)):
            for prot_next in list(list_prot)[i+1:]:

                if prot.isIncludedIn(prot_next):
                    prot.parental_prot.append(prot_next)
                    prot_next.sub_prot.append(prot)
                    prot.checkforAlternativeStart(prot_next)

                elif prot_next.isIncludedIn(prot):
                    prot_next.parental_prot.append(prot)
                    prot.sub_prot.append(prot_next)
                    prot_next.checkforAlternativeStart(prot)


    def writeAnnotatedProteins(self, file_handle, genetic_code):
        for polyprotein in self.polyproteins:
            if polyprotein.peptides: # we write every protein that have at least one peptide annotation
                seq_to_write = polyprotein.getSequenceRecord(self.organism, self.taxon_id, self.record, genetic_code)
                SeqIO.write(seq_to_write, file_handle,"fasta")


    def writeIdentifiedPolyproteins(self, file_handle, genetic_code):
        for polyprotein in self.polyproteins:
            seq_to_write = polyprotein.getSequenceRecord(self.organism, self.taxon_id, self.record, genetic_code)
            SeqIO.write(seq_to_write, file_handle,"fasta")


    def writeCleavageSite(self, file_handle, genetic_code, window_step):
        for cleavage_site in self.cleavage_sites:
            site_seq = cleavage_site.extractSeq(self.organism, self.taxon_id, self.record, genetic_code, window_step) #extractSeq(self, organism, taxon, record, genetic_code, window_step)
            # print(site_seq)
            SeqIO.write(site_seq, file_handle, "fasta")


    def getCleavageSites(self):

        start_cleavage_sites = {}
        for pep in sorted(list(self.peptides | self.unannotated_region), key=lambda x: x.start, reverse=False): #(self.peptides | self.unannotated_region):  # sorted(list(self.peptides | self.unannotated_region), key=lambda x: x.start, reverse=False):

            for type, border in [("start", pep.start), ('end', pep.end)]:

                start = border - 3
                end = border + 3
                # if start == 6077:
                #     print('START and END', start, end)
                #     print(type, border)
                polyproteins = {poly for poly in pep.polyproteins if poly.start+9 < start < poly.end-9} #polyprot that are compatible with the cleavage site
                if polyproteins:
                    if start in start_cleavage_sites:
                        site = start_cleavage_sites[start]

                        start_cleavage_sites[start].update(type, pep, polyproteins) # this cleavage site has been already treated
                        continue


                    cleavage_site = CleavageSite(start, end, polyproteins, pep)
                    # print("Cleavage site:", type, border,' start end', start, end, '\t',  cleavage_site.start_aa(poly), cleavage_site.end_aa(poly))
                    start_cleavage_sites[start] =  cleavage_site
                    self.cleavage_sites.append(cleavage_site)


    def getDomainOverlappingInfo(self):

        for poly in self.polyproteins:
            for pep in list(poly.peptides) + poly.unannotated_region:

                if getattr(pep, 'parent_peptide', False):
                    continue

                pep_start = pep.start_aa(poly)
                pep_end   = pep.end_aa(poly)

                for m in poly.matchs:

                    if pep_end < m.start_in_prot  or m.end_in_prot < pep_start: # the match is not on peptide seq
                        continue

                    if pep_start <= m.start_in_prot and m.end_in_prot <= pep_end: # The match fall into the peptide
                        # print(m.name)
                        m.including_peptides.append(pep)
                        continue

                    m.overlapping = True
                    # here 3 possible scenario
                    # The match overlap on the left, on the right or overlap oon the left and on the right

                    # We find the length of the sequence that the domain annotation share with pep
                    # to later determine which peptide share the biggest sequence with annotatuon
                    max_start = max(pep_start, m.start_in_prot)
                    min_end = min(pep_end, m.end_in_prot)

                    m.partially_included_in[pep] = min_end - max_start +1
                    # if m.name == "PF13608":
                    #
                    #     print("PEPTIDE ", pep.number, pep_start, pep_end)
                    #     print("match", m.start_in_prot, m.end_in_prot)
                    #     print("max satrt ", max_start,'min end', min_end)
                    #     print('distance', min_end - max_start +1)
                    #     print('*******')

            for m in poly.matchs:
                #For the match that overlap 2 or more peptide
                if not m.overlapping:
                    continue
                # to which peptide the match belong to? based on the length
                max_length = 0
                for pep, length in m.partially_included_in.items():

                    if length > max_length:
                        max_length = length
                        peptide = pep
                m.including_peptides.append(peptide)
                pep_start = peptide.start_aa(poly)
                pep_end   = peptide.end_aa(poly)

                if  pep_end < m.end_in_prot: # the annotation overlap on the right
                    m.right_overlaps_peptide =  m.end_in_prot - pep_end

                if  m.start_in_prot < pep_start: # it overlaps on the left
                    m.left_overlaps_peptide = pep_start - m.start_in_prot


    def identifyDuplicatedMatch(self):
        # Due to ribosomal_slippage 2 proteins have been given to interpro but they share a similar part and then may share same somain
        # To not count twice this domain we identify the domain that have the same start and end in the genome and have same name

        # The structure is not that efficient because the duplicated match check the other match and it is not usefull they should be remove
        matchs = list(self.matchs)
        for i, match in enumerate(matchs):
            for match_next in matchs[i+1:]:
                if match.isEqualTo(match_next):
                    match_next.duplicated = True


    def isSegmentAnnotationRelevant(self):
        return False if any([p for p in self.polyproteins if not p.isAnnotationRelevant() ]) else True


class Sequence:

    def __str__(self):

        string = 'Sequence: from {} to {}'.format(self.start,self.end)
        # string += '\n'.join([str(p) for p in self.peptides])
        # print('ddd')
        return string


    def __len__(self):
        return len(self.bp_obj)
        # return self.end - self.start +1


    def overlap(self, seq):
        if self.start <= seq.end <= self.end or seq.start <= self.end <= seq.end:
            return True
        return False


    def getStringIndices(self, conversion):
        # print(conversion,' start', self.start/conversion, 'end', (self.end-1)/conversion)
        # print(conversion,' start', int(self.start/conversion), 'end', int((self.end-1)/conversion))
        return (int(self.start/conversion), int((self.end-1)/conversion))
        # return round(self.bp_obj.location.start/conversion), round(self.bp_obj.location.end/conversion)


    def getGenomicPositions(self, prot_start):
        #for uncovered region and match to get the genomic position and not the prot relative position
        self.start = prot_start -1 + self.start_in_prot*3 -2 -1
        self.end = prot_start + self.end_in_prot*3 -1


    def getNameToDisplay(self):
        return 'None'


class Protein(Sequence):

    COUNTER=0

    def __init__(self, biopyth_obj):
        # self.start = start
        # self.end = end

        self.bp_obj = biopyth_obj
        self.start, self.end = (biopyth_obj.location.start, biopyth_obj.location.end)

        self.gene_type = biopyth_obj.type
        self.protein_id = "Unknown" if 'protein_id' not in biopyth_obj.qualifiers else biopyth_obj.qualifiers['protein_id'][0]
        self.product =  "Unknown" if 'product' not in biopyth_obj.qualifiers else biopyth_obj.qualifiers['product'][0]
        self.peptides = set()
        self.matchs = []
        self.unannotated_len = len(biopyth_obj)/3
        self.polyprotein = True
        self.unannotated_region = []
        self.non_overlapping_prot = set() #for the visualisation
        self.ribosomal_slippage = {}
        self.readthrough = {}
        self.parental_prot = []
        self.sub_prot = []
        self.alternative_start = False
        self.cleavage_sites = []
        self.polyprotein_number = 0.0
        self.sequence = ''

        Protein.COUNTER+=1
        self.number = Protein.COUNTER


    def __len__(self):
        return len(self.bp_obj)


    def __str__(self):

        # status += str(self.polyprotein_number)
        string = 'Protein: {} {}aa   | {} \n'.format(self.number, len(self)/3, self.bp_obj.location)
        # for p in self.peptides:
        #     string += str(p) +"\n"
        # string += "Position_number {}\n".format(self.polyprotein_number)
        # string += "{}:{}\n".format(self.gene_type, ' and '.join(self.reasons))


        # string += '{} annotated peptides | {} covered by peptide\n'.format(len(self.peptides), 100 - (self.unannotated_len/(len(self)/3))*100) if self.polyprotein else ''
        # string += '{} domain annotations\n'.format(len(self.matchs)) if self.polyprotein else ''
        # string += '\n'.join([str(p) for p in sorted(list(self.peptides), key=lambda x: x.start, reverse=False)])
        # print('ddd')
        return string


    def __contains__(self, seq):
        ##test if a peptide is in polyprotein
        # the qualifiers gene is no more used to know if a peptide belong to a polyprotein
        # Only the position of the peptide tel us if the peptide belongs to a poly
        # Start < end always. The coordinate are given according the strand of the gene
        # pep = seq.bp_obj # extract biopython object of peptide objet
        # seq_location = seq if not seq.__class__.__name__ ==  'Peptide' else seq.bp_obj.location
        if not seq.__class__.__name__ ==  'Peptide':
            logging.warning("the magic method from Porotein class is used by {}".format(seq.__class__.__name__))

        poly = self.bp_obj
        if not (seq.start in poly.location and seq.end in poly.location and seq.location.strand ==  poly.location.strand):
            # print("direct F")
            return False


        pep_parts = iter(seq.location.parts)
        pep_part = next(pep_parts)
        case2 = False
        len_prot = 0
        # print('---------------------------')
        # print('PEPTIDE: ', seq.number, seq.location)
        #
        # print('Protein: ', self.number, poly.location)

        for sub_location in poly.location.parts:
            len_prot += len(sub_location) # len prot: length of the protein from the first part to the i_part
            modulo_check = sub_location.end - len_prot%3 # len_prot%3: correction to be always in the correct frame. modulo_check is always in a

            # print(sub_location)
            # print('LEN',len(sub_location),"MODULO", len(sub_location)%3)
            # print('\n==Protein Sub location==')
            # print('pap part', pep_part)
            # print('prot part', sub_location)
            if case2:
                # print('\n-WE are in case 2', )

                if pep_part.start == sub_location.start: # check if the 2 parts start are the same. So if the join() is similar
                    # print('--the start are the same')
                    # print("--pep parts", pep_parts)
                    try:
                        pep_part = next(pep_parts)
                    except StopIteration:
                        if seq.end in sub_location and seq.end%3 == len_prot - modulo_check%3: # No more part for the peptide and it end is included in the subprot prat
                            # print('---end is included in the prot part')
                            # print('TRUUE')
                            return True
                    else:
                        # print('--We go to next part')

                        continue
                else: #The peptide in somehow not folowing the prot part
                    # print('-The peptide in somehow not folowing the prot part')
                    # print('FALSE')
                    return False
            # Case 1 : The peptide is included in one part of the protein
            if seq.start in sub_location and seq.end in sub_location and  len(seq.location.parts) == 1 and seq.end%3 == modulo_prot%3:
                # print("-the peptide is included in the subprot part")
                return True

            # Case 2 Peptide is overlapping this part with the next one
            if seq.start in sub_location and seq.end not in sub_location and pep_part.end == sub_location.end:
                pep_part = next(pep_parts)
                case2 = True
                # print('-the peptide is overlaping the prot part')
        return False


    def getSequenceAA(self, record, genetic_code):

        if self.sequence:
            return self.sequence

        elif 'translation' in self.bp_obj.qualifiers:
            seq = self.bp_obj.qualifiers["translation"][0]
        else:
            seq = self.bp_obj.location.extract(record).seq.translate(table=genetic_code)
            seq = str(seq)
        seq = seq.upper()
        #WRITE CLEAVAGE SITE IN lower CASE
        # print(self.cleavage_sites)
        # if self.isAnnotationRelevant():
        print(seq)
        for site in self.cleavage_sites:
            site_position = site.start_aa(self)-1 # -1 to be in base 0

            site_extraction = site.bp_obj.location.extract(record).seq.translate(table=genetic_code)
            site_seq = seq[site_position:site_position+2]
            if site_extraction.upper() != site_seq.upper():
                logging.warning('cleavage site is not similar in extraction from genome sequence and the protein sequence: extraction:{} and prot_seq:{}'.format(site_extraction, site_seq))
            # # print(site.peptide.bp_obj.location.extract(record).seq.translate(table=genetic_code))
            # print("site position",site_position)
            seq = seq[:site_position] + seq[site_position:site_position+2].lower() + seq[site_position+2:]


            print('cleavage site', site.number, site_position, seq[site_position-10:site_position+11])
            print('Extraction:', site.bp_obj.location.extract(record).seq.translate(table=genetic_code))
            for pep in site.peptides:
                print('peptide : ',pep.number)
                print("pep sequence", pep.bp_obj.location.extract(record).seq.translate(table=genetic_code))

        self.sequence = seq

        return seq


    def extractCleavageSites(self, record, genetic_code, window):
        cleavage_sites_seq = []
        if self.isAnnotationRelevant():
            for site in self.cleavage_sites:
                site_position = site.peptide.end_aa(self)-1 # -1 to be in base 0
                print(site.peptide.bp_obj.location.extract(record).seq.translate(table=genetic_code))

                cleavage_sites_seq.append(seq[site_position-window:site_position+window])

                print(site_position, seq[site_position-10:site_position+11])


    def isAnnotationRelevant(self):
        #Simple rule to be sure that the annotation is relevant
        #  the number of uncvered region < number of peptide
        #To not take into acount this kind of annotation
        # ---==1.0============================================================================================================================================================================>------------------------
        #           |+5|         |+6+++|     |+7++++|
        #    |~8~~~~|  |~9~~~~~~~|     |~10~~|      |~11~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

        if len(self.unannotated_region) >= len([p for p in self.peptides if not p.parent_peptide]):
             return False
        else:
            return True


    def getSequenceRecord(self,organism, header, record, genetic_code, subPosition=False):


        seq = self.getSequenceAA(record, genetic_code)

        # print(seq)
        if subPosition:
            seq = seq[subPosition[0]:subPosition[1]]
            header += '|{}:{}'.format(subPosition[0],subPosition[1])

        return  SeqRecord(Seq(seq,generic_protein) , id=header, description="{}|{}|{}".format( self.protein_id, self.product, organism))
        # print("PROTTTT")
        # print(seq_to_write)
        # SeqIO.write(seq_to_write, file_handle,"fasta")


    def polyproteinCoverage(self):
        ##Check if the mat peptide associated with the protein are covering all the prot or not
        #If not create unannotated_region to fill the gaps
        unannotated_position = []

        self.unannotated_len = 0
        iter_parts = iter(self.bp_obj.location.parts)
        protpart = next(iter_parts)
        current_po = protpart.start
        partial_location = None
        # for pep in sorted(list(self.peptides), key=lambda x: x.start, reverse=False) :
        peptides = sorted(list(self.peptides), key=lambda x: x.start, reverse=False)
        pep_i = 0
        while pep_i < len(peptides):
            pep = peptides[pep_i]
            if current_po < pep.start: # if pep start is after the current positon then a unnatotade region need to be created

                if protpart.end < pep.start: #that means the next peptide start is located in the next part of the protein
                    location = FeatureLocation(current_po, protpart.end)
                    partial_location = partial_location + location if partial_location else location

                    protpart = next(iter_parts)
                    current_po = protpart.start
                    continue

                else:
                    #an unannotated position is created.
                    location = FeatureLocation(current_po, pep.start)
                    if partial_location:
                        location = partial_location + location
                        partial_location = None
                    unannotated_seq = UnannotatedRegion(location, self)
                    current_po = pep.end
                    self.unannotated_region.append(unannotated_seq)

            elif current_po < pep.end:
                if protpart.end < pep.end: #the peptide is overlapping two part of the prot
                    protpart = next(iter_parts) # we take the next part

                current_po = pep.end
            pep_i += 1

            # else: #Another case where curent_po > pep.end but then we don't change the current po and we move to the next pep
            #     current_po =  current_po


        if current_po < protpart.end: # we didn't not arrive at the end of the
            location = FeatureLocation(current_po, protpart.end)
            final_location = partial_location + location if partial_location else location
            current_po = protpart.end
            for protpart in iter_parts:
                # print(protpart)
                final_location += FeatureLocation(current_po, protpart.end)
                current_po = protpart.end

            if len(final_location) > 3:
                unannotated_seq = UnannotatedRegion(final_location, self)

                self.unannotated_region.append(unannotated_seq)
            # unannotated_seq = UnannotatedRegion(current_po, int(len(self.bp_obj)/3 -1), self) #SeqFeature(FeatureLocation(current_po, int(len(self.bp_obj)/3 -1)), type="unannotated_region", qualifiers={'note':'Position given in amino acid', "start_aa":current_po-1, 'end_aa':int(len(self.bp_obj)/3 -1-1)})
            # unannotated_seq = Peptide(unannotated_seq_feature)

            # self.unannotated_len += len(self.bp_obj)/3 - current_po
        # poly.qualifiers['unannotated_region'] = unannotated_region
        # poly.qualifiers['peptide_coverage'] = len(poly)/3 - unannotated_len -1 # -1 because we dont count the stop codon
        # print(unannotated_position)


    def isIncludedIn(self, poly_next):

        next_location = poly_next.bp_obj.location
        start, end = self.bp_obj.location.start, self.bp_obj.location.end
        if start in next_location and end in next_location and start%3 == next_location.start%3:
            return True
        return False


    def checkforAlternativeStart(self, poly):

        if self.end ==  poly.end and  self.start !=  poly.start and self.start%3 ==  poly.start%3:
             self.alternative_start = True


    # def getCleavageSites(self):
    #
    #     pep_border  = set()
    #     print('Polyprot:',self.start, self.end)
    #     # print( len(self.peptides))
    #     for pep in sorted(list(self.peptides)+self.unannotated_region, key=lambda x: x.start, reverse=False):
    #         if pep.__class__.__name__ ==  'Peptide':
    #             print("Peptide", pep.start, pep.end, '\t', pep.start_aa(self), pep.end_aa(self))
    #         else:
    #             print("Unannot", pep.start, pep.end, '\t', pep.start_aa(self), pep.end_aa(self))
    #         print(pep.end)
    #     print(pep_border)


class Peptide(Sequence):

    COUNTER = 0
    def __init__(self, bp_obj):
        #by default the start and end is in amino acid
        self.bp_obj = bp_obj
        self.location = bp_obj.location
        self.start = bp_obj.location.start
        self.end =  bp_obj.location.end
        self.redundant_pep = []
        self.non_overlapping_prot = set() #for the visualisation
        self.parent_peptide = False

        self.polyproteins = set()

        self.position_prot_relative = {}
        Peptide.COUNTER +=1
        self.number = Peptide.COUNTER
        # self.genome_start = None
        # self.genome_end = None


    def __str__(self):

        string = 'Peptide {}: from {} to {} | belongs to {}\n'.format(self.number, self.start,self.end, [prot.number for prot in self.polyproteins])
        string += '{} nt| location :{}  \n'.format(len(self), self.bp_obj.location)
        string += '  Position in protein\n'
        for prot in self.polyproteins:
            string += '    protein {} from {} to {}\n'.format(prot.number, self.start_aa(prot), self.end_aa(prot))

        return string


    def start_aa(self, prot):
        if prot.protein_id not in self.position_prot_relative:
            self.getProteinPosition(prot)

        return self.position_prot_relative[prot.protein_id][0]


    def end_aa(self, prot):

        if prot.protein_id not in self.position_prot_relative:
            self.getProteinPosition(prot)

        return self.position_prot_relative[prot.protein_id][1]


    def get_position_prot_relative(self, prot):
        if prot.protein_id not in self.position_prot_relative:
            self.getProteinPosition(prot)

        return self.position_prot_relative[prot.protein_id]


    def getProteinPosition(self, prot):
        print("==================getProteinPosition======")
        len_previous_part = 0
        bp_prot = prot.bp_obj
        ## Searching the peptides position protein relative
        ## Due to ribosomal slippage the conversion is not trivial
        for subprotpart in bp_prot.location.parts:
            print("subprotpart ", subprotpart)
            print(self.start, self.start%3)
            print(subprotpart.start, subprotpart.start%3 )

            if subprotpart.start <=  self.start <= subprotpart.end:

                pstart = len_previous_part +  self.start-subprotpart.start +1
                p_end = pstart + len(self.bp_obj) -1

                if p_end <= bp_prot.location.end:

                    self.position_prot_relative[prot.protein_id] =  (int((pstart-1)/3+1), int((p_end-2-1)/3+1))
                    break
            # In some genome the cleavage site is located at the location of the ribosomal slippage
            # Then the first codon of the site is in frame of the next subprotpart but is not included in this part
            # and consequently if above don't catch it
            # I test here the end codon if it is included in the part and is in frame
            elif subprotpart.start <=  self.end <= subprotpart.end:
                pstart = len_previous_part +  self.start-subprotpart.start +1
                p_end = pstart + len(self.bp_obj) -1

                if p_end <= bp_prot.location.end:
                    self.position_prot_relative[prot.protein_id] =  (int((pstart-1)/3+1), int((p_end-2-1)/3+1))
                    break

            len_previous_part += len(subprotpart)


    def getNameToDisplay(self):
        return str(self.number)


class UnannotatedRegion(Peptide):

    def __init__(self, location, protein):
        self.start =  location.start
        self.end   = location.end
        self.polyproteins = [protein] # this region belongs to this protein

        self.non_overlapping_prot = set()
        self.position_prot_relative = {}

        Peptide.COUNTER +=1
        self.number = Peptide.COUNTER

        self.bp_obj = SeqFeature(location,
            type="unannotated_region",
            qualifiers={})


class CleavageSite(Peptide):
    COUNTER = 0
    def __init__(self, start, end, proteins, peptide):
        self.start = start
        self.end = end
        self.proteins = proteins # Set of protein
        self.peptides = {peptide}
        self.position_prot_relative = {}
        self.bp_obj = SeqFeature(FeatureLocation(self.start, self.end),
            type="cleavage site",
            qualifiers={})
        #NOTIFY PROTEINS THAT THEY hAVE A CLEAVAGE SITE
        for poly in proteins:
            poly.cleavage_sites.append(self)

        CleavageSite.COUNTER +=1
        self.number = CleavageSite.COUNTER

    def __eq__(self, other):
        return True if self.position == other.position else False

    def __str__(self):
        return 'Cleavage site {} from {} to {} ({}nt) belongs to {}, And have been made from peptide {}'.format(self.number, self.start, self.end, len(self), [p.number for p in self.proteins], [p.number for p in self.peptides] )

    def update(self, type, pep, polyproteins):

        polyproteins = polyproteins - self.proteins # new proteins that haven't been notify of the cleavage site

        for poly in polyproteins:
            poly.cleavage_sites.append(self)
        self.peptides.add(pep)
        self.proteins |= polyproteins

    # def extractSeq(self, organism, taxon, record, genetic_code, window_step):
    #     #window in aa to extract
    #     sites = []
    #     # proteins = sorted(list(self.proteins), key=lambda x: x.polyprotein_number, reverse=True)
    #
    #     for protein in proteins:
    #         # print(protein.polyprotein_number)
    #         position_prot_relative = int((self.position - protein.start)/3) #posiiton starting at 0
    #         if (self.position - protein.start)%3 != 0:
    #             logging.warning('cleavage_site seems to not point the first nt of a codon '+str(self))
    #
    #         # print(position_prot_relative)
    #         header = taxon+'|protein_{}'.format(protein.polyprotein_number)
    #         seq =  protein.getSequenceRecord('', header, record, genetic_code, (position_prot_relative-window_step,position_prot_relative+window_step +1))
    #         # print(dir(seq))
    #         seq.description = ''
    #         sites.append(seq)
    #         # seq.id = header #seq.id +'prot'+str(protein.number)
    #     for i, seq in enumerate(sites):
    #         for seq_next in sites[i+1:]:
    #             if seq_next.seq != seq.seq:
    #                 logging.warning("A same cleavage site give different sequences")
    #
    #     return seq
    #         # print(position_prot_relative/3)


class Match(Sequence):
    def __init__(self, taxid, seqid,  method, start, end, score, Dbxref, Name, matchID, signature_desc):
        self.taxid =taxid
        self.seqid=seqid
        self.method = method
        self.start_in_prot = int(start)
        self.end_in_prot = int(end)
        self.score=score
        self.Dbxref=Dbxref
        self.name =Name
        self.matchID = matchID
        self.signature_desc = signature_desc
        self.protein = None

        self.duplicated = False

        self.overlapping = False
        self.partially_included_in = {} # all the peptides that have a part of the domains annotation in their seq
        self.including_peptides = []
        self.right_overlaps_peptide = 0 # distance overlap on the right with the peptide
        self.left_overlaps_peptide = 0 # distance on the left

        self.non_overlapping_prot = set()

    def __len__(self):
        return self.end - self.start +1

    def __str__(self):
        string = '\n==Domain {} from {} and {}\n'.format(self.name, self.start_in_prot, self.end_in_prot)
        # string += 'Is included in: ' + str([p.number for p in self.including_peptides]) + '\n'
        # for p in self.including_peptides:
        #     string += '  pep:'+str(p.number) +':   '+ str(p.get_position_prot_relative( self.protein) ) + '\n'
        # string += '\nIs overlaping on the right: ' + str({p.number:dist for p, dist in self.right_overlaps_peptides.items()}) + '\n'
        # for p, dist in self.right_overlaps_peptides.items():
        #     string +='  pep:'+str(p.number) +'  overlap of '  + str(dist) +':   '+str(p.get_position_prot_relative( self.protein) ) + '\n'
        #     string += "  "+str(p)+'\n'
        # string += '\nIs overlaping on the left: ' + str({p.number:dist for p, dist in self.left_overlaps_peptides.items()}) + '\n'
        # for p, dist in self.left_overlaps_peptides.items():
        #
        #     string +='  pep:'+str(p.number) +'  overlap of '  + str(dist) +':   '+str(p.get_position_prot_relative( self.protein) ) + '\n'
        #     string += "  "+str(p)+'\n'
        return string

    def isEqualTo(self, other_match):
        # If 2 match has same genomic start and end and have teh same interprot id they are duplicated
        # they would be only different in the seqid attr and potentially in the protein coordinate
        if self.start == other_match.start and  self.end == other_match.end and self.Dbxref and other_match.Dbxref:
            return True

        return False

    def getNameToDisplay(self):

        if self.duplicated:
            name = 'd'*len(self.name)
        else:
            name = self.name
        return name if not self.overlapping else '<' + name

    def get_csv_dico(self, match_header):
        dico = {}
        for attribute_name in match_header:
            if hasattr(self,attribute_name):
                attribute = getattr(self, attribute_name)

                if callable(attribute):
                    dico[attribute_name] = attribute()
                else:
                    dico[attribute_name] = attribute

        return dico
    def OverlappingDistance(self):
        return self.left_overlaps_peptide + self.right_overlaps_peptide
