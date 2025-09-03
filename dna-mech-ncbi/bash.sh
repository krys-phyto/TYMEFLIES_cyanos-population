### Generate reference genes from ncbi geneIDs
# list currently in notebook


### conda download of ncbi bash interface programs
# get sequences from ncbi
datasets download gene gene-id --inputfile 2025-03-03_geneIDs.txt --include gene,cds,protein

### blastdb
conda install -y blast

cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/dna-mech-ncbi/ncbi_dataset/data

makeblastdb -in gene.fna -dbtype nucl

blastn -query /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/cyano_diamond_76.fna -db gene.fna -out cyano-diamond_dna-mech_blast.txt

blastn -query /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/cyano_diamond_76.fna -db gene.fna -out cyano-diamond_dna-mech_blast.txt \
   -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

### fix up blast output
header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"
awk "BEGIN {print \"$header\"} {print}" cyano-diamond_dna-mech_blast.txt > header.tmp
mv header.tmp cyano-diamond_dna-mech_blast.txt

### create blastn_seqIDs.txt
grep '>' gene.fna > blastn_seqIDs.txt

sed -i 's/>//g' blastn_seqIDs.txt
header="gene_seq"
awk "BEGIN {print \"$header\"} {print}" blastn_seqIDs.txt > header.tmp
mv header.tmp blastn_seqIDs.txt
