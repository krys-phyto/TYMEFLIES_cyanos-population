#! /bin/bash

#################################################################
###                  Krys Kibler 2025-07-24                   ###
###           CONDOR SUBMIT EXECUTABLE BASH SCRIPT            ###
### Purpose: filter bam files to only include 95%ANI reads    ###
#################################################################

# code is from pre-mcorr.sh and adpated for a condor run
# qeueu file is here analysis/mcorr/APHAN_134/linkage.maps

mOTU=$1
sampleID=$2
maps=$3

### Generating a file of fiiltered reads
#a. detailed mapping info for each read in bam files is produced from instrain

# itterate over detailed mapping files from instrain
#cyano_50-vs-rr0283.IS_detailed_mapping_info.csv
detailed_mapping_file=/mnt/bigdata/linuxhome/kjkibler/cyanoSTRONG/krys/analysis/instrain/cyano_diamond_50/instrain-92ANI/output/diamond_50/files
detailed_mapping=cyano_50-vs-${sampleID}.IS_detailed_mapping_info.csv
destination=/home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/mcorr

cd $destination
mkdir $mOTU/

#line,read_pair,scaffold,mm,insert_dist,mapq,length,reads,start,stop
#b. grep tax_abbr_name from file - gets all reads mapped to it
grep "$mOTU" $detailed_mapping_file/$detailed_mapping > $mOTU/filtered_$detailed_mapping


#c. filter for just paired reads
awk -F, '$8==2 {print}' $mOTU/filtered_$detailed_mapping > $mOTU/temp_file && mv -f $mOTU/temp_file $mOTU/filtered_$detailed_mapping


#d. filter for reads with less than 95% ANI
# read length * 95% = 7.5 different bps = < 8 mm
# remove reads with less than 8mm
awk -F, '$4<8 {print}' $mOTU/filtered_$detailed_mapping > $mOTU/temp_file && mv -f $mOTU/temp_file $mOTU/filtered_$detailed_mapping


#e. filter for min/max insert size, min = => 50, max =< 3* mean insert size
mean=`awk -F, '{sum += $5} END {if (NR > 0) print sum / NR}' $mOTU/filtered_$detailed_mapping`
max=$(echo "$mean * 3" | bc)
awk -v var="$max" -F, '$5<var {print}' $mOTU/filtered_$detailed_mapping > $mOTU/temp_file && mv -f $mOTU/temp_file $mOTU/filtered_$detailed_mapping
awk -F, '$5>49 {print}' $mOTU/filtered_$detailed_mapping > $mOTU/temp_file && mv -f $mOTU/temp_file $mOTU/filtered_$detailed_mapping


#f. print out list of reads that belong in metag to screen in bam/fastq files

cut -d "," -f2 $mOTU/filtered_$detailed_mapping > $mOTU/${mOTU}_${sampleID}_filtered.reads

#g. final output is a list of reads that passed instrain filter for an mOTU for each metagenome


### Generating filtered bam files from filtered reads file
SAMTOOLS=/home/glbrc.org/kjkibler/miniconda3/bin/samtools
#version 1.20

mkdir $mOTU/bt2/

# filter bam file based on generated ^^^
$SAMTOOLS view -b $maps -N $mOTU/${mOTU}_${sampleID}_filtered.reads > $mOTU/bt2/${mOTU}_${sampleID}_filtered.bam
