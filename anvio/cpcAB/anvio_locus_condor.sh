#! /bin/bash


#################################################################
###                  Krys Kibler 2024-09-16                   ###
###                    Condor BASH SCRIPT                     ###
### Purpose: run anvio to get phycocyan gene                  ###
### https://anvio.org/help/7.1/#anvio-artifacts               ###
#################################################################

# arguments
genome=$1
cpcA=$2
cpcB=$3


# conda
source /home/GLBRCORG/kjkibler/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
PERL5LIB=""

conda activate anvio-7.1

cd /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/anvio/mags

#code
mag=`basename $genome | sed 's/-/_/'`

anvi-export-locus -c $genome \
                   --flank-mode \
                   -o ../cpcAB/locus_output/ \
                   -O ${cpcA}_to_${cpcB}_${mag%_0seq.tmp.db} \
                   --gene-caller-ids $cpcA,$cpcB


anvi-get-sequences-for-gene-calls -c $genome \
                                  -o ../cpcAB/cpcA/${cpcA}_${mag%_0seq.tmp.db} \
                                  --gene-caller-ids ${cpcA}

anvi-get-sequences-for-gene-calls -c $genome \
                                  -o ../cpcAB/cpcB/${cpcB}_${mag%_0seq.tmp.db} \
                                  --gene-caller-ids ${cpcB} 
