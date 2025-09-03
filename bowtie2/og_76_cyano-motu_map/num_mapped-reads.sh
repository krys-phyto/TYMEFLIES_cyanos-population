
source /home/GLBRCORG/kjkibler/miniconda3/etc/profile.d/conda.sh
export PATH=//home/GLBRCORG/kjkibler/miniconda3/bin:$PATH
unset PYTHONPATH

PYTHONPATH=""
PERL5LIB=""

echo "sampleID,num_mapped" > num_mapped-reads.csv
for f in *[0-9].bam
do
echo $f
num=`samtools view -c -F 260 $f`
echo $f,$num >> num_mapped-reads.csv
done
