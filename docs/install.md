# Install Instructions

Once the install steps are complete [MAVIS](http://mavis.bcgsc.ca) is ready to be run.
See the MAVIS [tutorial](https://mavis.readthedocs.io/en/latest/tutorials/mini) to learn about running MAVIS.

For either install option you will want to install the main Snakefile. It is best to use a tag to
specify the version of interest but you can download the latest version from the master branch

```bash
wget https://raw.githubusercontent.com/bcgsc/mavis/master/Snakefile -O Snakefile
```

## Install for Docker/Singularity

The simplest way to use MAVIS is via Singularity. The MAVIS docker container used
by singularity will take care of installing the aligner as well.

```bash
pip install -U setuptools pip
pip install mavis_config  # also installs snakemake
```

Now you will run mavis via Snakemake as follows

```bash
snakemake \
    -j <MAX JOBS> \
    --configfile <YOUR CONFIG> \
    --use-singularity \
    -s Snakefile
```

## Install (Python Only)

MAVIS can also be run with just python. However you will need to install the aligner(s) required
by MAVIS separately and ensure they are availble on the default PATH variable when MAVIS is run

### 1. Install Aligner

In addition to the python package dependencies, [MAVIS](http://mavis.bcgsc.ca) also requires an aligner to be installed.
Currently the only aligners supported are [blat](https://mavis.readthedocs.io/en/latest/glossary/#blat) and [bwa mem](https://mavis.readthedocs.io/en/latest/glossary/#bwa).
For MAVIS to run successfully the aligner must be installed and accessible on the path.
If you have a non-standard install you may find it useful to edit the PATH environment variable. For example

``` bash
export PATH=/path/to/directory/containing/blat/binary:$PATH
```

[blat](http://mavis.bcgsc.ca/docs/latest/glossary.html#term-blat) is the default aligner. To configure MAVIS to use [bwa mem](http://mavis.bcgsc.ca/docs/latest/glossary.html#term-bwa) it must be specified
in the [config](https://mavis.readthedocs.io/en/latest/configuration/settings/) JSON file.

After this has been installed MAVIS itself can be installed through [pip](https://pypi.org/project/mavis/)

### 2. Install MAVIS

#### Install using pip

The easiest way to install [MAVIS](http://mavis.bcgsc.ca) is through the python package manager, pip. If you do not have python3 installed it can be found [here](https://www.python.org/downloads)

Ensuring you have a recent version of pip and setuptools will improve the install experience. Older versions of pip and setuptools may have issues with obtaining some of the mavis python dependencies

``` bash
pip install --upgrade pip setuptools
```

or (for Anaconda users)

``` bash
conda update pip setuptools
```

If this is not a clean/new python install it may be useful to set up mavis in a [virtual python environment](https://docs.python.org/3/tutorial/venv.html)

Then install mavis itself

``` bash
pip install mavis
```

This will install mavis and its python dependencies.

#### Install using Buildout

Alternatively you can use the [bootstrap/buildout](http://www.buildout.org/en/latest/) to install mavis into bin/mavis

``` bash
git clone https://github.com/bcgsc/mavis.git
cd mavis
pip install zc.buildout
python bootstrap.py
bin/buildout
```

This will install mavis and its python dependencies into eggs inside the cloned mavis directory which can be used by simply running bin/mavis

Finally you will need to Build/Download the necessary reference files

## Build or Download Reference Files

After [MAVIS](http://mavis.bcgsc.ca) is installed the [reference files](https://mavis.readthedocs.io/en/latest/inputs/reference) must be generated (or downloaded) before it can be run. A simple bash script to download the hg19 reference files is provided under mavis/tools for convenience.

### Download Hg19 Files

``` bash
cd /path/to/where/you/want/to/put/the/files
wget https://raw.githubusercontent.com/bcgsc/mavis/master/src/tools/get_hg19_reference_files.sh
bash get_hg19_reference_files.sh
```

You should now see the reference files in the current directory

```text
.
|-- cytoBand.txt
|-- dgv_hg19_variants.tab
|-- ensembl69_hg19_annotations.json
|-- get_hg19_reference_files.sh
|-- hg19.2bit
|-- hg19.fa
`-- hg19_masking.tab
```

### Download Hg38 Files

``` bash
cd /path/to/where/you/want/to/put/the/files
wget https://raw.githubusercontent.com/bcgsc/mavis/master/src/tools/get_hg38_reference_files.sh
bash get_hg19_reference_files.sh
```

You should now see the reference files in the current directory

```text
.
|-- cytoBand.txt
|-- dgv_hg38_variants.tab
|-- ensembl79_hg38_annotations.json
|-- get_hg38_reference_files.sh
|-- GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
|-- GRCh38_masking.tab
`-- hg38.2bit
```
