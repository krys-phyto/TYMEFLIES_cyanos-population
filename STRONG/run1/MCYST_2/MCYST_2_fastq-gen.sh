#! /bin/bash

#################################################################
###                  Krys Kibler 2025-02-06                   ###
###           CONDOR SUBMIT EXECUTABLE BASH SCRIPT            ###
### Purpose: creating fastq files from sorted trimmed bam     ###
### https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml   ###
#################################################################


SAMTOOLS=/home/glbrc.org/kjkibler/miniconda3/bin/samtools
#version 1.20

cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/STRONG/MCYST_2/

# arguments
maps=$1
sampleID=${maps%.sam}

#convert from sam to bam
$SAMTOOLS view -bS bt2/$maps > bt2/$sampleID.bam

# remove reads that did not map
$SAMTOOLS view -b -F 4 bt2/$sampleID.bam > bt2/$sampleID.onlymap.bam

# sorting mapping files
$SAMTOOLS sort -n bt2/$sampleID.onlymap.bam -o  bt2/$sampleID.sorted.onlymap.bam


mkdir /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/STRONG/MCYST_2/$sampleID


$SAMTOOLS fastq -@ 8 bt2/$sampleID.sorted.onlymap.bam \
    -1 /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/STRONG/MCYST_2/${sampleID}/${sampleID}.1.fastq.gz \
    -2 /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/STRONG/MCYST_2/${sampleID}/${sampleID}.2.fastq.gz \
    -0 /dev/null -s /dev/null -n
