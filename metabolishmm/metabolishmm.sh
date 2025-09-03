
source /home/GLBRCORG/kjkibler/miniconda3/etc/profile.d/conda.sh
export PATH=//home/GLBRCORG/kjkibler/miniconda3/bin:$PATH
unset PYTHONPATH

PYTHONPATH=""
PERL5LIB=""


# creating metadata sheet
for f in *.fna
do
name=${f%.fna}
echo $name,$name >> /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/metabolishmm/metadata.csv
done

conda activate metabolishmm
summarize-metabolism --input /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/final_50_tree/cyano_mOTUs/ \
      --output /home/glbrc.org/kjkibler/cyanoSTRONG/krys/analysis/metabolishmm/output/ \
      --metadata metadata.csv \
      --summary metabolishmm-summary.csv
