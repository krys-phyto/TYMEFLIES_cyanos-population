#! /bin/bash

#################################################################
###                  Krys Kibler 2024-07-08                   ###
###           CONDOR SUBMIT EXECUTABLE BASH SCRIPT            ###
### Purpose: microdiversity in the cyano_diamond_76           ###
### https://instrain.readthedocs.io/en/latest/                ###
#################################################################

# load up instrain
source /home/GLBRCORG/kjkibler/miniconda3/etc/profile.d/conda.sh

INSTRAIN=/home/glbrc.org/kjkibler/instrain/bin
#version 1.9.0

export PATH=//home/GLBRCORG/kjkibler/miniconda3/bin:$PATH
unset PYTHONPATH

PYTHONPATH=""
PERL5LIB=""

cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/cyano_diamond_50

# generate an stb file
# python /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/parse_stb.py --reverse -f /mnt/bigdata/linuxhome/trina.mcmahon/TYMEFLIES/data/cyano_diamond_76_fastas/*.fna -o cyano_diamond_76.stb

# prodigal version 2.6.3
# prodigal -i /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/comp85_50_cyano-motu_map/cyano_diamond_50_fasta.fna -d cyano_diamond_50_fasta.fna.genes.fna -a cyano_diamond_50_fasta.fna.genes.faa


# prepare qeueu file for sam files
#for f in *sam
#do
#  map=/home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/comp85_50_cyano-motu_map/$f
#  sampleID=${f%.sam}
#  echo -e $map '\t' $sampleID >> /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/cyano_diamond_50/mapping.files.tsv
#done

#sed -i 's#bt2#cyano_diamond_76/bt2#g' mapping.files.tsv


# arguments
map=$1
sampleID=$2

resultName=cyano_50-vs-$sampleID.IS


# profile command

$INSTRAIN/inStrain profile --pairing_filter paired_only $map \
  /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/comp85_50_cyano-motu_map/cyano_diamond_50_fasta.fna \
  -g cyano_diamond_50_fasta.fna.genes.fna \
  -o /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/cyano_diamond_50/instrain-95ANI/output/diamond_50/$resultName \
  -p 5 -s cyano_diamond_50.stb \
  --min_genome_coverage 10 \
  --min_read_ani 0.95 \
  --min_cov 5 \
  --detailed_mapping_info
