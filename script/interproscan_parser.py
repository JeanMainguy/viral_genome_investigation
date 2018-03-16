import function_parser as fct
                
    
if __name__ == '__main__':


    db_path = "/mirror/ncbi/current/genomes/refseq/viral/"
    gb_file = "/mirror/ncbi/current/genomes/refseq/viral/Abisko_virus/latest_assembly_versions/GCF_002270725.1_ViralProj399942/GCF_002270725.1_ViralProj399942_genomic.gbff.gz"
    #gb_file = sys.argv[1]
    file_path = "Feline_astrovirus_2/latest_assembly_versions/GCF_000910975.1_ViralProj218014/GCF_000910975.1_ViralProj218014_genomic.gbff.gz"
    
    file_path = "Betacoronavirus_1/latest_assembly_versions/GCF_000862505.1_ViralProj15385/GCF_000862505.1_ViralProj15385_genomic.gbff.gz"

    gb_file = db_path+file_path
    
    record, polyprot, matpep, poly_di = fct.extractPolyMatpep(gb_file)
    
    gff_file = "polyprotein_Bovine_coronavirus.gff3"
    #gff_file = "polyprotein_Middle_East_respiratory_syndrome-related_coronavirus.faa.html.tar.gz"
    
    sequences_match = fct.gffParser(gff_file)
    #print(sequences_match) 
    
    
    
    for pol in poly_di:
        union_sequence, pep_draw = fct.transformePolySeq(record, pol, poly_di[pol])
        protein_id = pol.qualifiers["protein_id"][0]
        
        match_draw = ' '*len(union_sequence)
        match_draws = [match_draw]
        draw_indice = 0
        symbol = '%'
        for i, m in enumerate(sequences_match[protein_id]):
            draw_indice = 0
            m_len = m["end"] - m['start'] +1
            start = m['start']-1
            end = m['end']-1
            while symbol in match_draws[draw_indice][start:end+1]:
                if len(match_draws) == draw_indice+1:
                    match_draws.append(' '*len(union_sequence))
                    draw_indice += 1 
                    break
                draw_indice += 1 
            match_draw = match_draws[draw_indice]
            
            match_draw = match_draw[:start] + m_len*str(i%10) + match_draw[end+1:]
            match_draw = match_draw[:start] + m_len*symbol + match_draw[end+1:]
            match_draws[draw_indice] = match_draw 
            
        #print(pep_draw)
        #print(union_sequence)
        #print(match_draw)
        step = 180
        
        for i in range(0, len(union_sequence), step):
            print(union_sequence[i:i+step])
            print(pep_draw[i:i+step])
            for md in match_draws:
                if '%' in md[i:i+step]:
                    print(md[i:i+step])
            print('\n')
    
            
        fct.checkForOverlaping(sequences_match[protein_id], poly_di[pol])
