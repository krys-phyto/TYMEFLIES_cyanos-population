#! /bin/bash

#################################################################
###                  Krys Kibler 2025-03-07                   ###
###           CONDOR SUBMIT EXECUTABLE BASH SCRIPT            ###
### Purpose: mapping the tymeflies APHAN-134                  ###
### https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml   ###
#################################################################


BOWTIE2=/mnt/bigdata/bifxapps/bowtie2-2.4.5
#version 2.4.5


cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/STRONG/VULCA_20/



# arguments
reads=$1
sampleID=$2


# bowtie2 code
# make an index
#mkdir bt2
#$BOWTIE2/bowtie2-build /home/glbrc.org/trina.mcmahon/TYMEFLIES/data/cyano_diamond_76_fastas/VULCA_20.fna bt2/VULCA_20.fa


# do the mapping
$BOWTIE2/bowtie2 -p 4 -x bt2/VULCA_20.fa --sensitive --interleaved \
$reads > /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/STRONG/VULCA_20/bt2/$sampleID.sam

# --sensitive -D 15 -R 2 -N 0 -L 22 -i S,1,1.15
