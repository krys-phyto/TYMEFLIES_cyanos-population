#! /bin/bash

#################################################################
###                  Krys Kibler 2025-07-24                   ###
### Purpose: microdiversity in the cyano_diamond_50           ###
### https://github.com/kussell-lab/mcorr                      ###
#################################################################





# start interactive condor shell

# activate mcorr
apptainer shell /mnt/bigdata/processed_data/computational_biology/containers/mcorr/mcorr.sif

for f in bt2/*.sorted.bam
do

output=${f%_filtered.sorted.bam}
# correlation profile
/usr/local/bin/mcorr-bam APHAN_134.gff bt2/APHAN_134_rd0665_filtered.bam APHAN_134_rd0665_mcorr-bam

# fit the correlation profile
/usr/local/bin/mcorr-fit APHAN_134_rd0665_mcorr-bam.csv APHAN_134_rd0665_mcorr-fit
