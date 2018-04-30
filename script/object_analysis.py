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

SCREEN_SIZE = 250


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

        # self.taxon_id = None
        # self.organism = None


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
        # print(taxonomy)²²²
        self.expectation_node='Not_Found'
        # self.peptide_expectation = 1
        for node in taxonomy:

            if node in taxon_expectation:
                self.expectation_node = node
                self.peptide_expectation = taxon_expectation[node]['peptide']
                self.polyprotein_expectation = taxon_expectation[node]['polyprotein']

        # print('expectation node', self.expectation_node)
        # print(taxonomy)


    def getTaxonomyInformation(self, dico, node_level):
        """
        We presume that the taxonomy is homogenous in the differents segments of a genome
        For one genome it is not the case but we ignore this case
        """
        taxon_expectation = TAXON_EXPECTATION

        if not self.segments:
            return 'ERROR no segment stored in genome'
        segment = self.segments[0]
        record = segment.record


        taxonomy = record.annotations['taxonomy']
        dico['taxonomy'] = taxonomy
        dico['organism'] = record.annotations['organism']
        dico['taxon_id'] = segment.taxon_id

        for node in taxonomy:
            if node.endswith('virales') and ' ' not in node:
                dico["order"] = node
            if node.endswith('viridae') and ' ' not in node:
                dico["family"] = node
            if node.endswith('virinae') and ' ' not in node:
                dico["subfamily"] = node
            if node.endswith('virus') and ' ' not in node:
                dico["genus"] = node

            if node in taxon_expectation:
                dico["expectation_node"] = node
                dico["peptide_expectation"] = taxon_expectation[node]

        last_node = taxonomy[0]
        last_first_node =taxonomy[-1]

        for n in range(node_level):
            try:
                dico[n+1] = taxonomy[n]
                last_node = taxonomy[n]
            except IndexError:
                dico[n+1] = last_node

            try:
                dico[(n+1)*-1] = taxonomy[(n+1)*-1]
                last_first_node = taxonomy[(n+1)*-1]
            except IndexError:
                dico[(n+1)*-1] = last_first_node


    def writecsvfile(self, writer_genome,writer_prot,  writer_poly):
        #['taxon_id', 'organism', 'taxonomy', 1, 2, 3, 4, -1, -2, -3, -4, 'segment',
        #'genome_length', 'cds', 'polyprotein', 'subpolyprotein', 'peptide', 'peptide_child',
        # 'ribosomal_slippage', 'alternative_start']

        #Take one segement

        node_level = max([k for k in writer_genome.fieldnames if isinstance(k, int)])
        csv_dico = {}
        self.getTaxonomyInformation(csv_dico, node_level)

        csv_dico['gb_file'] = self.gb_file

        csv_poly_dico = dict(csv_dico)
        for s in self.segments:

            s.writecsvProt(csv_poly_dico,writer_prot, writer_poly)
            # print(writer_poly.fieldnames)



        csv_dico['segment'] = len(self.segments)
        csv_dico['genome_length'] = len(self)
        csv_dico['cds'] = self.numberOf('cds')

        csv_dico['polyprotein'] = self.numberOf('polyproteins')
        csv_dico['subpolyprotein'] =  sum([len(s.getSubPolyprotein()) for  s in self.segments])
        csv_dico['peptide'] = self.numberOf('peptides')
        csv_dico['peptide_child'] = self.numberOf('child_peptides')


        csv_dico['ribosomal_slippage'] = sum([len(s.getPolyWithFeature( 'ribosomal_slippage')) for  s in self.segments])
        csv_dico['readthrough'] =sum([len(s.getPolyWithFeature( 'readthrough')) for  s in self.segments])
        csv_dico['alternative_start'] =sum([len(s.getPolyWithFeature( 'alternative_start')) for  s in self.segments])


        writer_genome.writerow(csv_dico)


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
            features_strg = segment.visualisation(nb_line, genetic_code)
            # print(features_strg)


    def getMatchObject(self, gff_file):
        #retreive all the match that belong to the genome
        flag = False
        with open(gff_file) as fl:
            for l in fl:
                if l.startswith(">"):
                    break
                if l.startswith('##sequence-region'):
                    flag = True if l.startswith('##sequence-region {}'.format(self.taxon_id)) else False
                    continue
                if flag:
                    (taxseqid, method, feature_type, start, end, score, strand, phase, attributes) = l.split("\t")
                    #match = {"start":int(begin), "end":int(end), "app":method, "seqid":seqid, "taxid":taxid}
                    if feature_type == 'protein_match':
                        (taxid, seqid) = taxseqid.split("|")
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
                poly.matchs = [match for match in self.matchs if match.seqid == poly.protein_id]
                for match in self.matchs:
                    if match.seqid == poly.protein_id:
                        poly.matchs.append(match)
                        match.protein = poly
                        match.getGenomicPositions(poly.start)
                        segment.matchs.add(match)


    def identifyExpectedElement(self):
        #The question is when their is more than 1 segment with more than 1 polyprotein...
        for segment in self.segments:
            segment.identifyExpectedPolyprotein(self.polyprotein_expectation)


    def checkDomainOverlapping(self):
        for segment in self.segments:
            segment.getDomainOverlappingInfo()
        for m in self.matchs:
            print(m)


class Segment:

    POSSIBLE_TYPE = set()

    def __init__(self, record, gb_file):
        self.record = record
        self.gb_file = gb_file
        self.source = None
        self.polyproteins =set() # polyproteins list
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


    def __len__(self):
        return len(self.source)


    def identifyExpectedPolyprotein(self, polyprotein_expectation):

        if not polyprotein_expectation:
            return
        # print(polyprotein_expectation, 'Expected polyprotein in this genome')
        parental_prot =  [ prot for prot in self.cds if not prot.parental_prot] # if this prot is not included in another one

        if len(parental_prot) < polyprotein_expectation:
            logging.warning('Not enough parental protein found {}/{}'.format(len(parental_prot), polyprotein_expectation))
            self.parental_polyproteins = parental_prot
        elif len(parental_prot) == polyprotein_expectation:
            self.parental_polyproteins = parental_prot
        elif len(parental_prot) > polyprotein_expectation:# rule: take the biggest protein until the number match the expectation and consider them as polyprotein
            self.parental_polyproteins = sorted(parental_prot, key=lambda x: len(x), reverse=True)[:polyprotein_expectation]

        self.parental_polyproteins = sorted(self.parental_polyproteins, key=lambda x: x.start, reverse=False)
        #giving a positional number to the polyprotein
        #if the polyprotein is included it get the positional number of the parent + a decimal number between 0.1

        for i, p_poly in enumerate(self.parental_polyproteins):
            p_poly.polyprotein_number = i+1.0
            # print(p_poly)
            for j, sub_poly in enumerate(sorted(p_poly.sub_prot, key=lambda x: x.start, reverse=False)):

                sub_poly.polyprotein_number = (i+1)+(j+1)/10
                # print(sub_poly)
                if j > 9: logging.error('polyprotein_number is not working')


    def writecsvProt(self, dico, writer_prot, writer_poly, genetic_code):
        '''
        Write csv file relative to Protein
        ['taxon_id', 'organism', 'taxonomy', 'order', 'family', 'subfamily', 'genus', 1, 2, 3, 4, -1, -2, -3, -4,
         'protein_id', 'protein_length', 'polyprotein', 'reason', 'coverage', 'peptide', 'peptide_child', 'unannotated_region'
        'ribosomal_slippage', 'alternative_start', 'readthrough', 'slippage1', 'slippage2', 'gb_file']
        '''
        for cds in self.cds:
            dico['protein_id'] = cds.protein_id
            dico['protein_length'] = len(cds)
            dico['polyprotein'] = 'TRUE' if cds.reasons else 'FALSE'
            dico['reason'] = ' and '.join(cds.reasons)
            dico['coverage'] =  round(((len(cds)/3 - cds.unannotated_len)/(len(cds)/3))*100, 2)
            dico['peptide'] = len(cds.peptides)
            dico['peptide_child'] = len([p for p in cds.peptides if not p.parent_peptide])
            dico['unannotated_region'] = len(cds.unannotated_region)

            dico['ribosomal_slippage'] = len(cds.ribosomal_slippage)
            dico['alternative_start'] = cds.alternative_start
            dico['readthrough'] =  len(cds.readthrough)
            shifts = list(cds.ribosomal_slippage.values())
            for i in range(2):
                dico['slippage'+str(i+1)] = None if len(shifts)<i+1 else shifts[i]
            # for i, shift in enumerate(cds.ribosomal_slippage.values()):
            #     dico['slippage'+str(i+1)] = shift


            if len(cds.ribosomal_slippage) > 2:
                print(cds.ribosomal_slippage)
                input()

            writer_prot.writerow(dico)
            if  cds.reasons:
                writer_poly.writerow(dico)


    def getSubPolyprotein(self):
        return [p for p in self.polyproteins if p.parental_prot ] # return protein that have parental prot


    def getPolyWithFeature(self, feature):
        return [p for p in self.polyproteins if getattr(p, feature)]


    def visualisation(self, nb_line, genetic_code):

        try:
            display_len = int(nb_line)*SCREEN_SIZE
        except ValueError:
            display_len = int(len(self.record)/3) if nb_line == 'aa' else int(len(self.record))

        conversion = int(len(self.record)) / display_len

        compatible_dico = collections.OrderedDict()
        compatible_dico['match'] = self.buildCompatibleGroup(self.matchs)[::-1]
        compatible_dico['cds'] = self.buildCompatibleGroup(self.cds)
        compatible_dico['pep'] = self.buildCompatibleGroup(self.peptides)
        compatible_dico['unannotated_region'] = self.buildCompatibleGroup(self.unannotated_region)

        strings = []
        seq_type_dico = {
            'match': ('[', '#', ']' , ' '),
            'cds': ('=', '=', '>', '-' ) ,
            'pep':('|', '+', '|' , ' '),
            'unannotated_region':('|', '~', '|' , ' ')
            }

        for seq_type, compatible_group in compatible_dico.items():
            left, central, right, neutral = seq_type_dico[seq_type]#different symbol if it protein or peptide
            for group in compatible_group:
                # print('GROUPPP')
                string = neutral*display_len
                # print(string)

                for seq in group:
                    i_start, i_end = seq.getStringIndices(conversion)
                    # print(i_start, i_end)
                    # print(len(seq),i_end - i_start +1 )
                    # print(i_start, i_end)
                    #match_draw[:start] + m_len*str(i%10) + match_draw[end+1:]
                    string = string[:i_start] +left+central*(i_end - i_start-1)+right + string[i_end+1:]

                    if seq.__class__.__name__ ==  'Protein':

                        # print(seq.getSequenceRecord(self.organism, self.taxon_id, self.record, genetic_code).seq)
                        if  nb_line == 'aa':
                            seqaa = seq.getSequenceAA(self.record, genetic_code)
                            # print('longueur aa*3', len(seqaa), 'longueur symb', len(left+central*(i_end - i_start-1)+right),'longueur seq', len(seq), 'longueur apres conv',i_end - i_start +1 )
                            # print(len(seqaa))
                            # print(len(seq))
                            # print(i_end+1, 'istart + seqaa len',i_start+len(seqaa) )
                            string = string[:i_start] + seqaa + string[i_start+len(seqaa):]
                            continue
                        if  nb_line == 'nt':
                            seqaa = seq.getSequenceAA(self.record, genetic_code)
                            seqaa = ''.join([a*3 for a in seqaa])
                            string = string[:i_start] + seqaa + string[i_start+len(seqaa):]
                            continue

                        #display positional number
                        if seq.polyprotein_number:
                            nb = str(seq.polyprotein_number)
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

                strings.append(string)

        for i in range(0, display_len, SCREEN_SIZE):
            print('/'*SCREEN_SIZE)
            print()
            for string in strings:
                print(string[i:i+SCREEN_SIZE])

        protein = '\n'.join(strings)

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
                if any(k for k in feat.qualifiers if "polyprotein" in feat.qualifiers[k][0]):
                    self.polyproteins.add(prot_obj)
                    prot_obj.polyprotein = True
                    prot_obj.reasons.add('"polyprotein" appears in qualifiers')


    def associatePepWithProt(self):
        #Should be changed at the end to catch polyprot that are not annotated
        #Here we want to find the polyprot that have a mat peptide so we are pretty sure that thz protein is a polyprotein
        for cds in self.cds:

            for pep in self.peptides:
                if pep in cds:
                    self.polyproteins.add(cds)
                    cds.peptides.add(pep)
                    pep.polyproteins.add(cds)


                    # print('cds', cds.bp_obj.location)
                    # print(len(cds.bp_obj))
                    # print("poly", pep.bp_obj)
            if cds.peptides:
                cds.reasons.add('has mature peptide annotation')
                cds.status = 'Polyprotein'
                cds.polyproteinCoverage()
                self.unannotated_region.update(cds.unannotated_region)


    def checkSubPeptides(self):
        for i, pep in enumerate(self.peptides):
            for pep_next in self.peptides:
                #try to determine if pep is included in pep_next
                if pep == pep_next:
                    continue

                if pep_next.start <= pep.start <= pep_next.end and pep_next.start <= pep.end <= pep_next.end:

                    if not(pep_next.start%3 == pep.start%3 and pep.end%3 == pep_next.end%3):
                        # print("PROBLEM/////")
                        # print('\n', pep.bp_obj,'\nin\n', pep_next.bp_obj, '\n' )
                        logging.warning('One peptide is included in another one but does not share the same strand: {}')
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
        # print("LEN PEP before removal",len(self.peptides))
        for i, pep in enumerate(peptides):
            # print("PEPTIDE:", pep)
            for pep_next in peptides[i+1:]:
                # pep.qualifiers['redundant_pep'] = []
                # print(pep_next)
                if pep.location.start ==  pep_next.location.start and pep.location.end ==  pep_next.location.end:
                    # print("REMove peptide")
                    self.peptides.remove(pep_next)
                    pep.redundant_pep.append(pep_next)

        # print("LEN PEP after removal",len(self.peptides))

    def checkForSlippage(self):
        """
        Find the protein that are included in another protein due to ribosomal_slippage are alternative start
        ribosomal_slippage attr is a dico with as a key the position of the end of the part of the protein and as a key the shift
        """

        for poly in self.polyproteins:
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
                    #Probably exon....
                    # print(poly.bp_obj)
                    logging.info('shift of {}, probably a spliced gene  {}'.format(shift, self.gb_file))
                    # input()
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


    def writeCleavageSite(self, file_handle, genetic_code, window_step):
        for cleavage_site in self.cleavage_sites:
            site_seq = cleavage_site.extractSeq(self.organism, self.taxon_id, self.record, genetic_code, window_step) #extractSeq(self, organism, taxon, record, genetic_code, window_step)
            # print(site_seq)
            SeqIO.write(site_seq, file_handle, "fasta")


    def getCleavageSites(self):
        # self.cleavage_sites
        threshold_distance = 5*3 #threshold distance between 2 cleavage sites to get a # WARNING:  logging

        treated_positions = []
        for pep in self.peptides | self.unannotated_region:  # sorted(list(self.peptides | self.unannotated_region), key=lambda x: x.start, reverse=False):
            # print(pep)
            position = pep.end

            # print('CURRENT PO:', position)
            polyproteins = {poly for poly in pep.polyproteins if poly.start+9 < position < poly.end-6} #polyprot that are compatible with the cleavage site
            if polyproteins:
                # updated_site = {site.proteins.update(polyproteins) for site in self.cleavage_sites if position == site.position}
                # if not updated_site:
                    ## CHECKING ///
                adjacent_sites = [site for site in self.cleavage_sites if position-threshold_distance < site.position < position+threshold_distance]
                if adjacent_sites:
                    logging.warning("A cleavage site has been found very close to another one... {} with {}".format(position, [p.position for p in adjacent_sites]))
                ## \\\\

                cleavage_site = CleavageSite(position, polyproteins, pep)
                for poly in polyproteins: poly.cleavage_sites.append(cleavage_site)
                self.cleavage_sites.append(cleavage_site)

    def getDomainOverlappingInfo(self):

        for poly in self.polyproteins:

            for pep in poly.peptides:
                pep_start = pep.start_aa(poly)
                pep_end   = pep.end_aa(poly)

                for m in poly.matchs:
                    # start_in_prot
                    # end_in_prot

                    if pep_end < m.start_in_prot  or m.end_in_prot < pep_start: # the match is not on peptide seq
                        continue
                    if pep_start <= m.start_in_prot and m.end_in_prot <= pep_end: # The match fall into the peptide
                        m.including_peptides.append(pep)
                        continue
                    if  pep_start <= m.start_in_prot <= pep_end and pep_end < m.end_in_prot: # the annotation overlap on the right
                        m.right_overlaps_peptides[pep] = m.end_in_prot - pep_end

                    if  pep_start <= m.end_in_prot <= pep_end and m.start_in_prot < pep_start: # it overlaps on the left
                        m.left_overlaps_peptides[pep] = pep_start - m.start_in_prot




class Sequence:
    def __str__(self):

        string = 'Sequence: from {} to {}'.format(self.start,self.end)
        # string += '\n'.join([str(p) for p in self.peptides])
        # print('ddd')
        return string

    def __len__(self):
        return self.end - self.start +1

    def overlap(self, seq):
        if self.start <= seq.end <= self.end or seq.start <= self.end <= seq.end:
            return True
        return False

    def getStringIndices(self, conversion):
        return (int(self.start/conversion), int((self.end-1)/conversion))
        # return round(self.bp_obj.location.start/conversion), round(self.bp_obj.location.end/conversion)

    def getGenomicPositions(self, prot_start):
        #for uncovered region and match to get the genomic position and not the prot relative position
        self.start = prot_start + self.start_in_prot*3 -2 -1
        self.end = prot_start + self.end_in_prot*3 #-1
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
        self.peptides = set()
        self.matchs = []
        self.unannotated_len = len(biopyth_obj)/3
        self.polyprotein = True
        self.reasons = set()
        self.unannotated_region = []
        self.non_overlapping_prot = set() #for the visualisation
        self.ribosomal_slippage = {}
        self.readthrough = {}
        self.parental_prot = []
        self.sub_prot = []
        self.alternative_start = False
        self.cleavage_sites = []
        self.polyprotein_number = 0.0

        Protein.COUNTER+=1
        self.number = Protein.COUNTER


    def __len__(self):
        return len(self.bp_obj)


    def __str__(self):

        # status += str(self.polyprotein_number)
        string = 'Protein: {} {}aa   | from {} to {} \n'.format(self.protein_id, len(self)/3, self.bp_obj.location.start,self.bp_obj.location.end)
        string += "Position_number {}\n".format(self.polyprotein_number)
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
        seq_location = seq if not seq.__class__.__name__ ==  'Peptide' else seq.bp_obj.location

        poly = self.bp_obj
        if not (seq_location.start in poly.location and seq_location.end in poly.location and seq_location.strand ==  poly.location.strand):
            return False


        start_flag = False
        for sub_location in poly.location.parts:
            if seq_location.start in sub_location and seq_location.start%3 == sub_location.start%3:
                # print("start seems to be correct")
                start_flag = True
            if start_flag and seq_location.end in sub_location and seq_location.end%3 == sub_location.end%3:
                # print('PEPTIDE BELONG TO POLYP')

                return True

        return False

    def getSequenceAA(self, record, genetic_code):

        if 'translation' in self.bp_obj.qualifiers:
            seq = self.bp_obj.qualifiers["translation"][0]
        else:
            seq = self.bp_obj.location.extract(record).seq.translate(table=genetic_code)
            seq = str(seq)
        seq = seq.upper()
        #WRITE CLEAVAGE SITE IN CAPITAL LETTER
        # print(self.cleavage_sites)
        for site in self.cleavage_sites:
            site_position = site.peptide.end_aa(self)-1 # -1 to be in base 0
            print(site.peptide.bp_obj.location.extract(record).seq.translate(table=genetic_code))

            seq = seq[:site_position] + seq[site_position:site_position+2].lower() + seq[site_position+2:]
            print(site_position, seq[site_position-10:site_position+11])
        return seq

    def getSequenceRecord(self,organism, taxon, record, genetic_code, subPosition=False):

        header = "{}|{}|{}".format(taxon, self.bp_obj.qualifiers["protein_id"][0], self.polyprotein_number)
        seq = self.getSequenceAA(record, genetic_code)

        # print(seq)
        if subPosition:
            seq = seq[subPosition[0]:subPosition[1]]
            header += '|{}:{}'.format(subPosition[0],subPosition[1])

        return  SeqRecord(Seq(seq,generic_protein) , id=header, description="{}|{}".format(self.bp_obj.qualifiers["product"][0], organism))
        # print("PROTTTT")
        # print(seq_to_write)
        # SeqIO.write(seq_to_write, file_handle,"fasta")


    def polyproteinCoverage(self):
        ##Check if the mat peptide associated with the protein are covering all the prot or not
        #If not create unannotated_region to fill the gaps
        unannotated_region = []
        current_po = 1
        self.unannotated_len = 0
        # pep_position = [(pep.qualifiers['start_aa'], pep.qualifiers['end_aa']) for pep in peptides] # if "parent_peptide" not in pep.qualifiers].sort()
        for pep in sorted(list(self.peptides), key=lambda x: x.bp_obj.location.start, reverse=False) :
            start, end = pep.get_position_prot_relative(self)
            if start > current_po: # if start is beyond the current positon
                unannotated_seq = UnannotatedRegion(current_po, start-1, self) #SeqFeature(FeatureLocation(current_po-1, start-1-1), type="unannotated_region", qualifiers={'note':'Position given in amino acid', "start_aa":current_po-1, 'end_aa':start-1-1})
                # unannotated_seq = Peptide(unannotated_seq_feature)
                self.unannotated_region.append(unannotated_seq)
                self.unannotated_len += start - current_po # nex peptide start -1 - current_po +1
            current_po = end+1 if end +1 > current_po else current_po


        if current_po < len(self.bp_obj)/3 -1: # the stop codon is take into account in the len of the prot so -1
            unannotated_seq = UnannotatedRegion(current_po, int(len(self.bp_obj)/3 -1), self) #SeqFeature(FeatureLocation(current_po, int(len(self.bp_obj)/3 -1)), type="unannotated_region", qualifiers={'note':'Position given in amino acid', "start_aa":current_po-1, 'end_aa':int(len(self.bp_obj)/3 -1-1)})
            # unannotated_seq = Peptide(unannotated_seq_feature)
            self.unannotated_region.append(unannotated_seq)
            self.unannotated_len += len(self.bp_obj)/3 - current_po
        # poly.qualifiers['unannotated_region'] = unannotated_region
        # poly.qualifiers['peptide_coverage'] = len(poly)/3 - unannotated_len -1 # -1 because we dont count the stop codon


    def isIncludedIn(self, poly_next):

        next_location = poly_next.bp_obj.location
        start, end = self.bp_obj.location.start, self.bp_obj.location.end
        if start in next_location and end in next_location and start%3 == next_location.start%3:
            return True
        return False


    def checkforAlternativeStart(self, poly):

        if self.end ==  poly.end and  self.start !=  poly.start and self.start%3 ==  poly.start%3:
             self.alternative_start = True


    def write(self, file_handler):
        pass


    def getProportionOfCoverage(self):
        pass

    def getProteinSequence(self):
        pass


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

        string = 'Peptide: from {} to {} | belongs to {}'.format(self.bp_obj.location.start,self.bp_obj.location.end, [prot.number for prot in self.polyproteins])

        # string += '\n'.join([str(p) for p in self.peptides])
        # print('ddd')
        return string


    def start_aa(self, prot):
        if self.__class__.__name__ ==  'UnannotatedRegion':
            return self.start_in_prot
        if prot.protein_id not in self.position_prot_relative:
            self.getProteinPosition(prot)

        return self.position_prot_relative[prot.protein_id][0]

    def end_aa(self, prot):
        if self.__class__.__name__ ==  'UnannotatedRegion':
            return self.end_in_prot

        if prot.protein_id not in self.position_prot_relative:
            self.getProteinPosition(prot)

        return self.position_prot_relative[prot.protein_id][1]

    def get_position_prot_relative(self, prot):
        if prot.protein_id not in self.position_prot_relative:
            self.getProteinPosition(prot)

        return self.position_prot_relative[prot.protein_id]

    def getProteinPosition(self, prot):
        len_previous_part = 0
        bp_prot = prot.bp_obj
        ## Searching the peptides position protein relative
        ## Due to ribosomal slippage the conversion is not trivial
        for subprotpart in bp_prot.location.parts:

            if subprotpart.start <=  self.bp_obj.location.start <= subprotpart.end and  self.bp_obj.location.start%3 == subprotpart.start%3 :
                pstart = len_previous_part +  self.bp_obj.location.start-subprotpart.start +1
                p_end = pstart + len( self.bp_obj) -1
                if p_end <= bp_prot.location.end:

                    self.position_prot_relative[prot.protein_id] =  (int((pstart-1)/3+1), int((p_end-2-1)/3+1))

                # pstart /= 3 # amino acid
            len_previous_part += len(subprotpart)

    def getNameToDisplay(self):
        return str(self.number)

class UnannotatedRegion(Peptide):

    def __init__(self, start_in_prot, end_in_prot, protein):
        self.start_in_prot =  start_in_prot
        self.end_in_prot   = end_in_prot
        self.polyproteins = [protein] # this region belongs to this protein
        self.getGenomicPositions(protein.start) #give start and end attr
        self.non_overlapping_prot = set()
        self.position_prot_relative = {}
        Peptide.COUNTER +=1
        self.number = Peptide.COUNTER

        self.bp_obj = SeqFeature(FeatureLocation(self.start, self.end),
            type="unannotated_region",
            qualifiers={})


class CleavageSite():
    def __init__(self, position, proteins, peptide):
        self.position = position
        self.proteins = proteins #Set of protein
        self.peptide = peptide

    def __eq__(self, other):
        return True if self.position == other.position else False

    def __str__(self):
        return 'Cleavage site {} belonging to proteins {}'.format(self.position, [p.number for p in self.proteins])

    def extractSeq(self, organism, taxon, record, genetic_code, window_step):
        #window in aa to extract
        sites = []
        proteins = sorted(list(self.proteins), key=lambda x: x.polyprotein_number, reverse=True)

        for protein in proteins:
            # print(protein.polyprotein_number)
            position_prot_relative = int((self.position - protein.start)/3) #posiiton starting at 0
            if (self.position - protein.start)%3 != 0:
                logging.warning('cleavage_site seems to not point the first nt of a codon '+str(self))

            # print(position_prot_relative)
            header = taxon+'|protein_{}'.format(protein.polyprotein_number)
            seq =  protein.getSequenceRecord('', header, record, genetic_code, (position_prot_relative-window_step,position_prot_relative+window_step +1))
            # print(dir(seq))
            seq.description = ''
            sites.append(seq)
            # seq.id = header #seq.id +'prot'+str(protein.number)
        for i, seq in enumerate(sites):
            for seq_next in sites[i+1:]:
                if seq_next.seq != seq.seq:
                    logging.warning("A same cleavage site give different sequences")

        return seq
            # print(position_prot_relative/3)


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
        self.including_peptides = []
        self.right_overlaps_peptides = {} # match overlap on the right with this peptide
        self.left_overlaps_peptides = {} # match overlap on the left with this peptide

        self.non_overlapping_prot = set()

    def __len__(self):
        return self.end - self.start +1

    def __str__(self):
        string = '\n==Domain {} from {} and {}\n'.format(self.name, self.start_in_prot, self.end_in_prot)
        string += 'Is included in: ' + str([p.number for p in self.including_peptides]) + '\n'
        for p in self.including_peptides:
            string += '  pep:'+str(p.number) +':   '+ str(p.get_position_prot_relative( self.protein) ) + '\n'
        string += '\nIs overlaping on the right: ' + str({p.number:dist for p, dist in self.right_overlaps_peptides.items()}) + '\n'
        for p, dist in self.right_overlaps_peptides.items():
            string +='  pep:'+str(p.number) +'  overlap of '  + str(dist) +':   '+str(p.get_position_prot_relative( self.protein) ) + '\n'
            string += "  "+str(p)+'\n'
        string += '\nIs overlaping on the left: ' + str({p.number:dist for p, dist in self.left_overlaps_peptides.items()}) + '\n'
        for p, dist in self.left_overlaps_peptides.items():

            string +='  pep:'+str(p.number) +'  overlap of '  + str(dist) +':   '+str(p.get_position_prot_relative( self.protein) ) + '\n'
            string += "  "+str(p)+'\n'
        return string

    def getNameToDisplay(self):
        return self.name


if __name__ == '__main__':

    logging.basicConfig(filename='log/genbankparser.log',level=logging.INFO)

    node_level = 4 # number of node to write in the output file 4 first one and 4 last one
    nodes_header = list(range(1,node_level+1))
    nodes_header += [n*-1 for n in nodes_header]
    taxonomy_file="data/taxonomy/taxonomy_virus.txt"
    expected_file =  "viruses_w_polyproteins.txt"

    # taxon = 'ssRNA viruses'
    # taxon = "390845"
    taxon="Viruses"
    # taxon='Blueberry virus A'


    #OUTPUT FILE GENOME: 1 line/genome
    fl_genome = open('RefSEQ_csv_stat_genomes.csv', 'w')
    header_genome = ['taxon_id', 'organism', 'taxonomy', 'order', 'family', 'subfamily', 'genus', 'expectation_node', "peptide_expectation"] +  nodes_header
    header_genome += ["segment",'genome_length', 'cds', 'polyprotein', 'subpolyprotein', 'peptide', 'peptide_child']
    header_genome += ['ribosomal_slippage','alternative_start', 'readthrough', 'gb_file']

    writer_genome = csv.DictWriter(fl_genome, fieldnames=header_genome, delimiter='\t')
    writer_genome.writeheader()

    #OUTPUT FILE Polypritein: 1 line/polyprot
    fl_poly = open('RefSEQ_csv_stat_polyproteins.csv', 'w')
    fl_prot = open('RefSEQ_csv_stat_proteins.csv', 'w')
    header_poly = ['taxon_id', 'organism', 'taxonomy', 'order', 'family', 'subfamily', 'genus', 'expectation_node', "peptide_expectation"] +  nodes_header


    header_poly += ['protein_id', 'protein_length', 'polyprotein', 'reason','coverage', 'peptide', 'peptide_child', 'unannotated_region']
    header_poly += ['ribosomal_slippage','alternative_start', 'readthrough', 'slippage1', 'slippage2', 'gb_file']

    writer_poly = csv.DictWriter(fl_poly, fieldnames=header_poly, delimiter='\t')
    writer_poly.writeheader()

    writer_prot = csv.DictWriter(fl_prot, fieldnames=header_poly, delimiter='\t')
    writer_prot.writeheader()


    taxon_expectation = tax.expectedPeptide(expected_file)



    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)


    for i, gb_dico in enumerate(gbff_iter):
        # print(gb_dico)
        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']
        # print(gb_file)
        # print('**'*20)
        # print(gb_file)
        if i%500 == 5:
            print(i)
            # break
            # break
        print(gbff)
        extractObjectFromGbff(gb_file, writer_genome, writer_prot, writer_poly, genetic_code, taxon_expectation)
    print(Segment.POSSIBLE_TYPE)
    fl_poly.close()
