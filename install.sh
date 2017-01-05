svn export https://svn.bcgsc.ca/svn/SVIA/vocab/tags/v0.2.0/vocab.py vocab.py
svn export https://svn.bcgsc.ca/svn/SVIA/TSV/tags/v2.0.1/TSV.py TSV.py
echo "generating the documentation"
export PATH=/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin:$PATH

cd docs
bash create.sh
