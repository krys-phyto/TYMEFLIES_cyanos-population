################################
### Krys Kibler
### Purpose: Use anvio to run tree
### Resource: https://gitpub.wei.wisc.edu/htcondor/strong-strain-resolution-docker


################################
cd
condor_submit --interactive anvio-interactive.submit
conda activate anvio-7.1

source /home/GLBRCORG/kjkibler/miniconda3/etc/profile.d/conda.sh
export PATH=//home/GLBRCORG/kjkibler/miniconda3/bin:$PATH
unset PYTHONPATH

PYTHONPATH=""
PERL5LIB=""


cd cyanoSTRONG/krys/analysis/anvio
mkdir tree

################################
echo name$'\t'contigs_db_path > external_genomes.txt
for f in *.db
do
  name=${f%_0seq.tmp.db}
  name="${name//-/_}"
  echo $name$'\t'$f >> external_genomes.txt
done


################################
anvi-get-sequences-for-hmm-hits --external-genomes mags/external_genomes.txt \
                                --hmm-source Bacteria_71 \
                                --list-available-gene-names

#Karthiks list of 16 ribosomal proteins to make a tree
#L2, L3, L4, L5, L6, L14, L16, L18, L22, L24, S3, S8, S10, S17 and S19
#Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6,Ribosomal_L14,Ribosomal_L16,Ribosomal_L18p,Ribosomal_L22,ribosomal_L24,Ribosomal_S3_C,Ribosomal_S8,Ribosomal_S10,Ribosomal_S17,Ribosomal_S19

anvi-get-sequences-for-hmm-hits --external-genomes mags/external_genomes.txt \
                                --hmm-source Bacteria_71 \
                                --gene-names Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6,Ribosomal_L14,Ribosomal_L16,Ribosomal_L18p,Ribosomal_L22,ribosomal_L24,Ribosomal_S3_C,Ribosomal_S8,Ribosomal_S10,Ribosomal_S17,Ribosomal_S19 \
                                --get-aa-sequences \
                                --concatenate-genes \
                                --return-best-hit \
                                -o tree/16ribosomal-fasta.fna


################################
anvi-gen-phylogenomic-tree -f tree/16ribosomal-fasta.fna \
                           -o tree/16ribosomal-fasta_tree.txt
