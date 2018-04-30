
out_file=data/taxonomy/taxonomy_virus.txt
ncbi_refseq_db=/mirror/ncbi/current/genomes/refseq/viral/
taxon_id_gencode_file=data/taxonomy/taxon_id_gencode.txt
final_output=data/taxonomy/taxonomy_virus_gencode.txt
python3 script/taxonomy.py ${out_file}  ${ncbi_refseq_db} ${taxon_id_gencode_file}



# join -1 1 -2 1 -t$'\t'  <(sort -n -k 1  ${out_file}) <(sort -n -k 1 ${taxon_id_gencode_file}) > ${final_output}
#
#
#
# wc -l ${out_file}
#
# wc -l ${final_output}
