#! /bin/bash


#################################################################
###                  Krys Kibler 2024-07-22                   ###
###                    Condor BASH SCRIPT                     ###
### Purpose: run anvio on the diamond Mags                    ###
### https://merenlab.org/2016/06/22/anvio-tutorial-v2/        ###
#################################################################

source /home/GLBRCORG/kjkibler/miniconda3/etc/profile.d/conda.sh
export PATH=//home/GLBRCORG/kjkibler/miniconda3/bin:$PATH
unset PYTHONPATH

PYTHONPATH=""
PERL5LIB=""


# activate anvio env
conda activate anvio-7.1


# anvio code

# Make/Move to working directory #
cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/anvio


# Generate Contig Database #
#$ anvi-gen-contigs-database -f contigs.fa -o contigs.db -n 'An example contigs database'
anvi-gen-contigs-database -f /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/bowtie2/cyano_diamond_76.fna \
  -o cyano_diamond_76.db -n 'Diamond-76 Cyano Mags' --split-length 0

# Run Hmms
#$ anvi-run-hmms -c contigs.db

anvi-run-hmms -c cyano_diamond_76.db --num-threads 8

# NCBI cogs/functions
#$ anvi-run-ncbi-cogs -c contigs.db

anvi-run-ncbi-cogs -c cyano_diamond_76.db --num-threads 8

# get those annotation
#anvi-export-functions -c contigs-db

anvi-export-functions -c cyano_diamond_76.db -o gene_functions.txt

# get gene calls too
#anvi-export-gene-calls -c contigs-db --gene-caller Prodigal --skip-sequence-reporting -o gene-calls-txt

anvi-export-gene-calls -c cyano_diamond_76.db \
                       --gene-caller prodigal \
                       --skip-sequence-reporting \
                       -o gene-calls-txt


# im getting sick of ncbi cogs... pivot to the keggs
cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/anvio
anvi-run-kegg-kofams -c cyano_diamond_76.db -T 8

anvi-export-functions -c cyano_diamond_76.db -o gene_functions-kegg.txt






#####################################################
### Generating stats
for f in *.db
do
anvi-profile -c $f  \
             -o ${f%.db}_PROFILE \
             -S blank \
            --blank-profile
done

for f in *.db
do
anvi-script-add-default-collection -c $f \
                                   -p ${f%.db}_PROFILE/PROFILE.db
done

anvi-show-collections-and-bins -p VULCA_96_PROFILE

for f in *.db
do
anvi-summarize -c $f \
               -p ${f%.db}_PROFILE/PROFILE.db \
               -C DEFAULT \
               -o $f.SUMMARY
done

for f in *SUMMARY
do
  cat $f/bins_summary.txt >> SUMMARY.txt
