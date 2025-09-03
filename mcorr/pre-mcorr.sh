#! /bin/bash

#################################################################
###                  Krys Kibler 2025-07-24                   ###
###           CONDOR SUBMIT EXECUTABLE BASH SCRIPT            ###
### Purpose: filter bam files to only include 95%ANI reads    ###
#################################################################

# I think I can only run this pipeline on one genome at a time
# 1. produce a GFF3 file for each mOTU in the snp linkage analysis
# 2. filter bam files for only 95% read ani for each mOTU (instrain filter process)
# 3. run mcorr on the mOTUs individually

# bam files
# /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/comp85_50_cyano-motu_map/bt2

# mOTUs in snp linkage analysis
# APHAN_134, MCYST_2, PSEUDA_13, NODOS_105, CYANO_106, DOLIS_187, MCYST_56, MCYST_62, MCYST_31
# CYBIM_90, CYBIM_157, CYBIM_104, CYBIM_89, CYBIM_101, CYBIM_63, CYBIM_190,
# CYBIM_200, CYBIM_119_1, CYBIM_199, CYBIM_73, VULCA_20, VULCA_28


### 1. produce a GFF3 file for each mOTU in the snp linkage analysis
# prodigal version 2.6.3
# prodigal -i <input_file.fasta> -f gff -o <output_file.gff> -p meta
prodigal -i /mnt/bigdata/linuxhome/kjkibler/cyanoSTRONG/krys/analysis/final_50_tree/cyano_mOTUs/APHAN_134.fna \
         -f gff -o /mnt/bigdata/linuxhome/kjkibler/cyanoSTRONG/krys/analysis/mcorr/APHAN_134/APHAN_134.gff -p meta


### 2. filter bam files for only 95% read ani for each mOTU (instrain filter process)
# process to filter bam files

#a. detailed mapping info for each read in bam files is produced from instrain
      # example: /mnt/bigdata/linuxhome/kjkibler/cyanoSTRONG/krys/analysis/instrain/cyano_diamond_50/instrain-92ANI/output/diamond_50/files/cyano_50-vs-rr0283.IS_detailed_mapping_info.csv
      #line,read_pair,scaffold,mm,insert_dist,mapq,length,reads,start,stop
#b. grep tax_abbr_name from file - gets all reads mapped to it
#c. filter for just paired reads
#d. filter for reads with less than 95% ANI
#e. filter for min/max insert size, min = => 50, max =< 3* mean insert size
#f. print out list of reads that belong in metag to screen in bam/fastq files
#g. final output is a filtered metagenome for an mOTU that passed instrain read filtering

#a. produced from instrain with --detailed_mapping_info flag
#line,read_pair,scaffold,mm,insert_dist,mapq,length,reads,start,stop

#b. grep tax_abbr_name from file - gets all reads mapped to it
grep "$mOTU" $detailed_mapping_file/$detailed_mapping > filtered_$detailed_mapping


#c. filter for just paired reads
awk -F, '$8==2 {print}' filtered_$detailed_mapping > temp_file && mv -f temp_file filtered_$detailed_mapping


#d. filter for reads with less than 95% ANI
# read length * 95% = 7.5 different bps = < 8 mm
# remove reads with less than 8mm
awk -F, '$4<8 {print}' filtered_$detailed_mapping > temp_file && mv -f temp_file filtered_$detailed_mapping


#e. filter for min/max insert size, min = => 50, max =< 3* mean insert size
mean=`awk -F, '{sum += $5} END {if (NR > 0) print sum / NR}' filtered_$detailed_mapping`
max=$(echo "$mean * 3" | bc)
awk -v var="$max" -F, '$5<var {print}' filtered_$detailed_mapping > temp_file && mv -f temp_file filtered_$detailed_mapping
awk -F, '$5>49 {print}' filtered_$detailed_mapping > temp_file && mv -f temp_file filtered_$detailed_mapping


#f. print out list of reads that belong in metag to screen in bam/fastq files

cut -d "," -f2 filtered_$detailed_mapping > ${mOTU}_${sampleID}_filtered.reads
awk '{print "@" $1}' ${mOTU}_${sampleID}_filtered.reads > temp_file && mv -f temp_file ${mOTU}_${sampleID}_filtered.reads

#g. final output is a list of reads that passed instrain filter for an mOTU for each metagenome


### Generating filtered bam files from filtered reads file
SAMTOOLS=/home/glbrc.org/kjkibler/miniconda3/bin/samtools
#version 1.20

mkdir bt2/

# filter bam file based on generated ^^^
$SAMTOOLS view -b $maps -N filtered_$detailed_mapping > bt2/${mOTU}_${sampleID}_filtered.sorted.bam
#    -0 /dev/null -s /dev/null -n
