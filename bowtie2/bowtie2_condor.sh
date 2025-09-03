#! /bin/bash

#################################################################
###                  Krys Kibler 2024-07-08                   ###
###           CONDOR SUBMIT EXECUTABLE BASH SCRIPT            ###
### Purpose: mapping the tymeflies mags                       ###
### https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml   ###
#################################################################


BOWTIE2=/mnt/bigdata/bifxapps/bowtie2-2.4.5
#version 2.4.5
SAMTOOLS=/opt/bifxapps/samtools-1.9/bin
#version 1.9

#cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2
cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/comp85_50_cyano-motu_map


# arguments
reads=$1
sampleID=$2


# bowtie2 code
# make an index
#mkdir bt2
#$BOWTIE2/bowtie2-build /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/comp85_50_cyano-motu_map/cyano_diamond_50_fasta.fna bt2/cyano_diamond_50_fasta.fa


# do the mapping
$BOWTIE2/bowtie2 -p 16 -x bt2/cyano_diamond_50_fasta.fa --sensitive --interleaved \
$reads > /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/comp85_50_cyano-motu_map/bt2/$sampleID.sam

# --sensitive -D 15 -R 2 -N 0 -L 22 -i S,1,1.15


# remove reads that did not map
$SAMTOOLS/samtools view -F 4 bt2/$sampleID.sam -o bt2/$sampleID.onlymap.sam

#convert sam to bam
$SAMTOOLS/samtools view -S -b bt2/$sampleID.onlymap.sam -o bt2/$sampleID.bam

# sorting mapping files
$SAMTOOLS/samtools  sort  bt2/$sampleID.bam -o  bt2/$sampleID.sorted.bam

#rm -f bt2/$sampleID.sam
