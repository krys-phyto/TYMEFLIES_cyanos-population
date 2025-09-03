### How to combine instrain output files
for f in *.IS
do
  cp $f/output/${f}_genome_info.tsv /mnt/bigdata/linuxhome/kjkibler/cyanoSTRONG/krys/analysis/instrain/cyano_diamond_50/instrain-95ANI/output/diamond_50/
  cp $f/output/${f}_linkage.tsv /mnt/bigdata/linuxhome/kjkibler/cyanoSTRONG/krys/analysis/instrain/cyano_diamond_50/instrain-95ANI/output/diamond_50/
  cp $f/raw_data/detailed_mapping_info.csv.gz /mnt/bigdata/linuxhome/kjkibler/cyanoSTRONG/krys/analysis/instrain/cyano_diamond_50/instrain-95ANI/output/diamond_50/${f}_detailed_mapping_info.csv.gz
done

find . -name '*.IS_*.tsv' -type f -empty -delete

for f in *detailed_mapping_info.csv.gz
do
gunzip $f
done

mkdir files
for f in *genome*
do
  name=$(basename -s .IS_genome_info.tsv "$f")
  mv $f files/
  mv $name.IS_detailed_mapping_info.csv files/
  mv $name.IS_linkage.tsv files/
done

# genome instrain
for i in *.IS_genome_info.tsv
do
    name=${i%.IS_genome_info.tsv}
    name=${name#cyano_50-vs-}
    awk -v name="$name" 'BEGIN { FS=OFS="\t" } { print $0, name }' $i > edit_$i
done

echo -e "genome\tcoverage\tbreadth\tnucl_diversity\tlength\ttrue_scaffolds\tdetected_scaffolds\tcoverage_median\tcoverage_std\tcoverage_SEM\tbreadth_minCov\tbreadth_expected\tnucl_diversity_rarefied\tconANI_reference\tpopANI_reference\tiRep\tiRep_GC_corrected\tlinked_SNV_count\tSNV_distance_mean\tr2_mean\td_prime_mean\tconsensus_divergent_sites\tpopulation_divergent_sites\tSNS_count\tSNV_count\tfiltered_read_pair_count\treads_unfiltered_pairs\treads_mean_PID\tdivergent_site_count\treads_unfiltered_reads\tsampleID" > cyano_50_readANI95_genome.tsv

awk 'FNR>1 || NR==1' edit*.tsv >> cyano_50_readANI95_genome.tsv

# linkage instrain
for i in *.IS_linkage.tsv
do
    name=${i%.IS_linkage.tsv}
    name=${name#cyano_50-vs-}
    awk -v name="$name" 'BEGIN { FS=OFS="\t" } { print $0, name }' $i > edit_$i
done

echo -e "scaffold\tposition_A\tposition_B\tdistance\tr2\td_prime\tr2_normalized\td_prime_normalized\tallele_A\tallele_a\tallele_B\tallele_b\tcountab\tcountAb\tcountaB\tcountAB\ttotal" > cyano_50_readANI95_linkage.tsv

awk 'FNR>1 || NR==1' edit*.tsv >> cyano_50_readANI95_linkage.tsv

# mapping details

for i in *.IS_linkage.tsv
do
    name=${i%.IS_linkage.tsv}
    name=${name#cyano_50-vs-}
    awk -v name="$name" 'BEGIN { FS=OFS="\t" } { print $0, name }' $i > edit_$i
done

echo -e "scaffold\tposition_A\tposition_B\tdistance\tr2\td_prime\tr2_normalized\td_prime_normalized\tallele_A\tallele_a\tallele_B\tallele_b\tcountab\tcountAb\tcountaB\tcountAB\ttotal" > cyano_50_readANI95_linkage.tsv

awk 'FNR>1 || NR==1' edit*.tsv >> cyano_50_readANI95_linkage.tsv
