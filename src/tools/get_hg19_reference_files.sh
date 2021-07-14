set -euo pipefail

echo "downloading the reference genome file"
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar -xvzf chromFa.tar.gz

# concatenate the chromosome fa files into a single file
for fname in chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}.fa
do
	cat $fname >> hg19.fa
done

# Clean up the non concatenated and alt chromosome files
rm -f chr*.fa
rm -f chromeFa.tar.gz

echo "downloading the gene annotations file"
wget http://www.bcgsc.ca/downloads/mavis/ensembl69_hg19_annotations.json

echo "downloading the masking file"
wget http://www.bcgsc.ca/downloads/mavis/hg19_masking.tab

echo "downloading the dgv annotation file"
wget http://www.bcgsc.ca/downloads/mavis/dgv_hg19_variants.tab

echo "downloading the aligner reference file"
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

echo "downloading the template metadata file"
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
gunzip cytoBand.txt.gz
