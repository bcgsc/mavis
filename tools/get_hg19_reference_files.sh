mkdir reference_inputs
cd reference_inputs

ENV_FILE=hg19_env.sh

echo "downloading the reference genome file"
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar -xvzf chromFa.tar.gz

CWD=$( pwd )

# concatenate the chromosome fa files into a single file
for fname in chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}.fa
do
	cat $fname >> hg19.fa
done

# Clean up the non concatenated and alt chromosome files
rm -f chr*.fa
echo export MAVIS_REFERENCE_GENOME="${CWD}/hg19.fa" >> $ENV_FILE

echo "downloading the gene annotations file"
wget http://www.bcgsc.ca/downloads/mavis/ensembl69_hg19_annotations.json
echo export MAVIS_ANNOTATIONS="${CWD}/ensembl69_hg19_annotations.json" >> $ENV_FILE

echo "downloading the masking file"
wget http://www.bcgsc.ca/downloads/mavis/hg19_masking.tab
echo export MAVIS_MASKING="${CWD}/hg19_masking.tab" >> $ENV_FILE

echo "downloading the dgv annotation file"
wget http://www.bcgsc.ca/downloads/mavis/dgv_hg19_variants.tab
echo export MAVIS_DGV_ANNOTATION="${CWD}/dgv_hg19_variants.tab" >> $ENV_FILE

echo "downloading the aligner reference file"
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
echo export MAVIS_ALIGNER_REFERENCE="${CWD}/hg19.2bit" >> $ENV_FILE

echo "downloading the template metadata file"
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
gunzip cytoBand.txt.gz
echo export MAVIS_TEMPLATE_METADATA="${CWD}/cytoBand.txt" >> $ENV_FILE

echo "Source $CWD/$ENV_FILE prior to running MAVIS to set MAVIS default arguments"
