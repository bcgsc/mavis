#!/bin/bash

# add the ensembl api modules to the path
PATH=$(pwd):$PATH
PERL5LIB=${PERL5LIB}:/home/creisle/applications/ensembl_79/bioperl-live
PERL5LIB=${PERL5LIB}:/home/creisle/applications/ensembl_79/ensembl/modules
PERL5LIB=${PERL5LIB}:/home/creisle/applications/ensembl_79/ensembl-compara/modules
PERL5LIB=${PERL5LIB}:/home/creisle/applications/ensembl_79/ensembl-variation/modules
PERL5LIB=${PERL5LIB}:/home/creisle/applications/ensembl_79/ensembl-funcgen/modules
export PERL5LIB

# default perl
PATH=/projects/trans_scratch/software/perl/perl-5.20.3/bin:$PATH

# required data files
export HUGO_ENSEMBL_MAPPING=/projects/tumour_char/analysis_scripts/databases/processed_files/drug_target_tables/current_gene_drug_pathway.hg38.tsv
export BEST_TRANSCRIPTS=/home/creisle/svn/ensembl_flatfiles/ens69_best_transcript.txt

# connection information for the ensembl local server
export ENSEMBL_HOST='ensembl02'
export ENSEMBL_PASS='ensembl'
export ENSEMBL_USER='ensembl'
export ENSEMBL_PORT=3306
