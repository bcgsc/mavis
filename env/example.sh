export MAVIS_TEMPLATE_METADATA='/home/creisle/git/mavis/tests/integration/data/cytoBand.txt'
export MAVIS_REFERENCE_GENOME='/projects/seqref/genomes/Homo_sapiens/GRCh37/1000genomes/bwa_ind/genome/GRCh37-lite.fa'
export MAVIS_ANNOTATIONS='/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv'
export MAVIS_MASKING='/projects/tumour_char/analysis_scripts/SVIA/delly/reference_data/GRCh37/human_nspan.hg19.excl.with_header.tsv'
export MAVIS_ALIGNER_REFERENCE='/home/pubseq/genomes/Homo_sapiens/GRCh37/blat/hg19.2bit'
export MAVIS_DGV_ANNOTATION='/home/cchoo/git/mavis/references/GRCh37_hg19_variants_2016-05-15.input.txt'
export MAVIS_MAX_FILES=100
export MAVIS_MIN_CLUSTERS_PER_FILE=30
export PYTHONUNBUFFERED='True'

#Add paths for samtools, blat and git
export PATH=/projects/trans_scratch/transabyss/trans-ABySS/v1.4.10/bin/:/gsc/software/linux-x86_64-centos6/git-2.12.0/bin/:$PATH
