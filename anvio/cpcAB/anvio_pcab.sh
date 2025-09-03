#! /bin/bash


#################################################################
###                  Krys Kibler 2024-09-12                   ###
###                    Condor BASH SCRIPT                     ###
### Purpose: run anvio to get phycocyan gene                  ###
### https://anvio.org/help/7.1/#anvio-artifacts               ###
#################################################################

cd
condor_submit --interactive anvio-interactive.submit

source /home/GLBRCORG/kjkibler/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
PERL5LIB=""

conda activate anvio-7.1



#Phycocyanin alpha chain, CpcA/CpeB (CpcA) (PDB:1JBO)
#Phycocyanin/phycoerythrin beta chain, CpcB/CpeB (CpcB) (PDB:1JBO)

cd cyanoSTRONG/krys/analysis/anvio/

### Create a separate .db for each cyano mag
for f in ../STRONG/EVAL/Genome/*.tmp
do
  mag=`basename $f`
  anvi-gen-contigs-database -f $f -o mags/$mag.db  --split-length 0
done

cd mags/
for f in *.db
do
anvi-run-hmms -c $f --num-threads 8
anvi-run-ncbi-cogs -c $f --num-threads 8
done

# Anvio flank target locus and print out sequence
for f in mags/*.db
do
mag=`basename $f`
anvi-export-locus -c $f \
                   --flank-mode \
                   -o locus_output \
                   -O cpcA_to_cpcB-$mag \
                   --search-term "Phycocyanin alpha chain","Phycocyanin/phycoerythrin beta chain"
done
