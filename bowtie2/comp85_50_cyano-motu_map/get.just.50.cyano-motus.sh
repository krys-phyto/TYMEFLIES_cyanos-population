### For mapping need a concatenated just 50 cyano-motus
# diamond set cyanos w/ above 85% completeness

cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/comp85_50_cyano-motu_map


for line in $(cat cyano_motus.tsv)
do
 cat /home/glbrc.org/trina.mcmahon/TYMEFLIES/data/cyano_diamond_76_fastas/$line.fna >> cyano_diamond_50_fasta.fna
 echo "Line: $line"
done


### create stb file for inStrain
for line in $(cat /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/comp85_50_cyano-motu_map/cyano_motus.tsv)
do
 #cat /home/glbrc.org/trina.mcmahon/TYMEFLIES/data/cyano_diamond_76_fastas/$line.fna >> cyano_diamond_50_fasta.fna
 echo "Line: $line"
 grep $line cyano_diamond_76.stb >> cyano_diamond_50.stb
done
