set -euo pipefail

echo "downloading the reference genome (no alt) file"
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

echo "downloading the gene annotations file"
wget http://www.bcgsc.ca/downloads/mavis/v3/ensembl79_hg38_annotations.v3.json.gz
gunzip ensembl79_hg38_annotations.v3.json.gz

echo "downloading the masking file"
wget http://www.bcgsc.ca/downloads/mavis/GRCh38_masking.tab

echo "downloading the dgv annotation file"
wget http://www.bcgsc.ca/downloads/mavis/dgv_hg38_variants.tab

echo "downloading the aligner reference file"
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit

echo "downloading the template metadata file"
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
gunzip cytoBand.txt.gz
