#! /bin/bash


#################################################################
###                  Krys Kibler 2024-07-22                   ###
###                    Condor BASH SCRIPT                     ###
### Purpose: obtain readani in figure 3 of instrain           ###
### https://github.com/MrOlm/inStrain/issues/131              ###
#################################################################

source /home/GLBRCORG/kjkibler/miniconda3/etc/profile.d/conda.sh
export PATH=//home/GLBRCORG/kjkibler/miniconda3/bin:$PATH
unset PYTHONPATH

PYTHONPATH=""
PERL5LIB=""

# arguments
instrain_out=$1
sampleID=$2

# csv output
csv=`echo /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/output/readANI_$sampleID.csv`

cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/

# run python script
python instrain_readANI.py --instrain_out=$instrain_out --csv=$csv



#code to create instrain_output.tsv for queue
#for f in *IS/
#do
#  instrain_out=/home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/output/$f
#  sampleID=${f%.IS/}
#  echo -e $instrain_out '\t' $sampleID >> /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/instrain_output.tsv
#done

# example inputs
#instrain_out=/home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/output/diamond_76/cyano_76-vs-rr0692.IS
#sampleID=cyano_76-vs-rr0692
#csv=/home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/output/readANI_cyano_76-vs-rr0692

# for loop for running this python script since condor is a pain
for f in *.IS
do
  instrain_out=$f
  csv=`echo /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/output/readANI_${f%.IS/}.csv`
  echo "working on $f"
  python /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/instrain_readANI.py --instrain_out=$instrain_out --csv=$csv
done


# concatenating output files
awk '
    FNR == 1 {
        ID = FILENAME
        gsub(/.IS.csv$/,"",ID)
        next
    }
    { print ID "," $0 }
' *.IS.csv > all_readANI.csv

echo -e "sampleID\tscaffold\tgene\tgene_length\tcoverage\tbreadth\tbreadth_minCov\tnucl_diversity\tstart\tend\tdirection\tpartial\tdNdS_substitutions\tpNpS_variants\tSNV_count\tSNV_S_count\tSNV_N_count\tSNS_count\tSNS_S_count\tSNS_N_count\tdivergent_site_count" > all_gene_info.tsv

for f in *.IS
do
  echo "working on $f"
  cd $f/output
  awk '
  FNR == 1 {
      ID = FILENAME
      gsub(/.tsv$/,"",ID)
      next
  }
  { print ID "\t" $0 }
' *.IS_gene_info.tsv >> /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/output/diamond_76/all_gene_info.tsv
 cd ../..
done

echo -e "sampleID\tgenome\tcoverage\tbreadth\tnucl_diversity\tlength\ttrue_scaffolds\tdetected_scaffolds\tcoverage_median\tcoverage_std\tcoverage_SEM\tbreadth_minCov\tbreadth_expected\tnucl_diversity_rarefied\tconANI_reference\tpopANI_reference\tiRep\tiRep_GC_corrected\tlinked_SNV_count\tSNV_distance_mean\tr2_mean\td_prime_mean\tconsensus_divergent_sites\tpopulation_divergent_sites\tSNS_count\tSNV_count\tfiltered_read_pair_count\treads_unfiltered_pairs\treads_mean_PID\treads_unfiltered_reads\tdivergent_site_count" > all_genome_info.tsv

for f in *.IS
do
  echo "working on $f"
  cd $f/output
  awk '
  FNR == 1 {
      ID = FILENAME
      gsub(/.tsv$/,"",ID)
      next
  }
  { print ID "\t" $0 }
' *.IS_genome_info.tsv >> /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/output/diamond_76/all_genome_info.tsv
 cd ../..
done

echo -e "sampleID\tscaffold\tlength\tcoverage\tbreadth\tnucl_diversity\tcoverage_median\tcoverage_std\tcoverage_SEM\tbreadth_minCov\tbreadth_expected\tnucl_diversity_median\tnucl_diversity_rarefied\tnucl_diversity_rarefied_median\tbreadth_rarefied\tconANI_reference\tpopANI_reference\tSNS_count\tSNV_count\tdivergent_site_count\tconsensus_divergent_sites\tpopulation_divergent_sites" > all_scaffold_info.tsv

for f in *.IS
do
  echo "working on $f"
  cd $f/output
  awk '
  FNR == 1 {
      ID = FILENAME
      gsub(/.tsv$/,"",ID)
      next
  }
  { print ID "\t" $0 }
' *.IS_scaffold_info.tsv >> /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/output/diamond_76/all_scaffold_info.tsv
 cd ../..
done
