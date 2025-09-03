#! /bin/bash

#################################################################
###                  Krys Kibler 2025-06-11                   ###
###           CONDOR SUBMIT EXECUTABLE BASH SCRIPT            ###
### Purpose: creating fastq files from sorted trimmed bam     ###
### https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml   ###
#################################################################


SAMTOOLS=/home/glbrc.org/kjkibler/miniconda3/bin/samtools
#version 1.20

cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/STRONG/run2/MCYST_31/

# arguments
maps=$1
sampleID=${maps%.sorted.bam}
sampleID=${sampleID#/mnt/bigdata/linuxhome/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/comp85_50_cyano-motu_map/bt2/}


# code
mkdir $sampleID

$SAMTOOLS fastq -@ 8 $maps \
    -1 /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/STRONG/run2/MCYST_31/${sampleID}/${sampleID}.1.fastq.gz \
    -2 /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/STRONG/run2/MCYST_31/${sampleID}/${sampleID}.2.fastq.gz \
    -0 /dev/null -s /dev/null -n

/mnt/bigdata/bifxapps/bbmap-38.32/filterbyname.sh in=$sampleID/${sampleID}.1.fastq.gz \
                                                  in2=$sampleID/${sampleID}.2.fastq.gz \
                                                  out=$sampleID/${sampleID}_filtered.1.fastq.gz \
                                                  out2=$sampleID/${sampleID}_filtered.2.fastq.gz include=t ths=t\
                                                  names=/home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/instrain/cyano_diamond_50/instrain-92ANI/output/diamond_50/files/STRONG_reads/MCYST_31/${sampleID}_filtered.reads
