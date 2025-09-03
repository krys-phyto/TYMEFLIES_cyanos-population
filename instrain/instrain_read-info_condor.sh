#! /bin/bash

#################################################################
###                  Krys Kibler 2025-05-29                   ###
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

cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/cyano_diamond_50/instrain-92ANI/output/diamond_50

# arguments
map=$1



# profile command

$INSTRAIN/inStrain filter_reads --pairing_filter paired_only $map \
  /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/comp85_50_cyano-motu_map/cyano_diamond_50_fasta.fna \
  -o /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/cyano_diamond_50/instrain-92ANI/output/diamond_50/ \
  -p 6 \
  --min_read_ani 0.92 \
  --detailed_mapping_info

inStrain filter_reads --pairing_filter paired_only /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/comp85_50_cyano-motu_map/bt2/rr0906.sorted.bam \
  /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/comp85_50_cyano-motu_map/cyano_diamond_50_fasta.fna \
  -o /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/cyano_diamond_50/instrain-92ANI/output/diamond_50/rr0906/ \
  -p 6 \
  --min_read_ani 0.92 \
  --detailed_mapping_info
