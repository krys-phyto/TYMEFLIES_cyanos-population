#! /bin/bash

#################################################################
###                  Krys Kibler 2024-07-15                   ###
###           CONDOR SUBMIT EXECUTABLE BASH SCRIPT            ###
### Purpose: running antismash on the diamond_76              ###
### https://antismash.secondarymetabolites.org/#!/start       ###
#################################################################

source /home/GLBRCORG/kjkibler/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
PERL5LIB=""

conda activate antismash


cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/antismash


# code for antismash
antismash --genefinding-tool prodigal  \
    --output-dir output \
    --fullhmmer \
    --tfbs \
    --cb-subclusters \
    --cb-knownclusters \
    --rre \
    --pfam2go \
    /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/cyano_diamond_76.fna
