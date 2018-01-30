.. _reference-input:

Reference Input Files
======================

There are several reference files that are required for full functionality of the MAVIS pipeline. If the same
reference file will be reused often then the user may find it helpful to set reasonable defaults. Default values
for any of the reference file arguments can be :ref:`configured through environment variables <config-environment>`.

To improve the install experience for the users, different configurations of the MAVIS annotations file have been made available. These files can be downloaded below, or if the required configuration is not available, :ref:`instructions on generating the annotations file <generate-reference-annotations>` can be found below.

.. list-table::
    :header-rows: 1
    
    *   - File Name (Type/Format)
        - Environment Variable
        - Download
    *   - :ref:`reference genome <reference-files-reference-genome>` (:term:`fasta`)
        - ``MAVIS_REFERENCE_GENOME``
        - .. raw:: html
            
            <a class='download-button btn btn-neutral' href='http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz' download>
                <img src='_static/Ic_cloud_download_48px.svg'>GRCh37/Hg19
            </a><br>
            <a class='download-button btn btn-neutral' href='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' download>
                <img src='_static/Ic_cloud_download_48px.svg'>GRCh38
            </a>
    *   - :ref:`annotations <reference-files-annotations>` (:term:`JSON`)
        - ``MAVIS_ANNOTATIONS``
        - .. raw:: html
    
            <a class='download-button btn btn-neutral' href='http://www.bcgsc.ca/downloads/mavis/ensembl69_hg19_annotations.json' download>
                <img src='_static/Ic_cloud_download_48px.svg'>GRCh37/Hg19 + Ensembl69
            </a><br>
            <a class='download-button btn btn-neutral' href='http://www.bcgsc.ca/downloads/mavis/ensembl79_hg38_annotations.json' download>
                <img src='_static/Ic_cloud_download_48px.svg'>GRCh38 + Ensembl79
            </a>
    *   - :ref:`masking <reference-files-masking>` (text/tabbed)
        - ``MAVIS_MASKING``
        - .. raw:: html
    
            <a class='download-button btn btn-neutral' href='http://www.bcgsc.ca/downloads/mavis/hg19_masking.tab' download>
                <img src='_static/Ic_cloud_download_48px.svg'>GRCh37/Hg19
            </a><br>
            <a class='download-button btn btn-neutral' href='http://www.bcgsc.ca/downloads/mavis/GRCh38_masking.tab' download>
                <img src='_static/Ic_cloud_download_48px.svg'>GRCh38
            </a>
    *   - :ref:`template metadata <reference-files-template-metadata>` (text/tabbed)
        - ``MAVIS_TEMPLATE_METADATA``
        - .. raw:: html
    
            <a class='download-button btn btn-neutral' href='http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz' download>
                <img src='_static/Ic_cloud_download_48px.svg'>GRCh37/Hg19
            </a><br>
            <a class='download-button btn btn-neutral' href='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz' download>
                <img src='_static/Ic_cloud_download_48px.svg'>GRCh38
            </a>
    *   - :ref:`DGV annotations <reference-files-dgv-annotations>` (text/tabbed)
        - ``MAVIS_DGV_ANNOTATION``
        - .. raw:: html
    
            <a class='download-button btn btn-neutral' href='http://www.bcgsc.ca/downloads/mavis/dgv_hg19_variants.tab' download>
                <img src='_static/Ic_cloud_download_48px.svg'>GRCh37/Hg19
            </a><br>
            <a class='download-button btn btn-neutral' href='http://www.bcgsc.ca/downloads/mavis/dgv_hg38_variants.tab' download>
                <img src='_static/Ic_cloud_download_48px.svg'>GRCh38
            </a>
    *   - :ref:`aligner reference <reference-files-aligner-reference>`
        - ``MAVIS_ALIGNER_REFERENCE``
        - .. raw:: html
            
            <a class='download-button btn btn-neutral' href='http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit' download>
                <img src='_static/Ic_cloud_download_48px.svg'>GRCh37/Hg19 2bit (blat)
            </a><br>
            <a class='download-button btn btn-neutral' href='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit' download>
                <img src='_static/Ic_cloud_download_48px.svg'>GRCh38 2bit (blat)
            </a>



If the environment variables above are set they will be used as the default values when any step of the pipeline
script is called (including generating the template config file)


.. _reference-files-reference-genome:

Reference Genome
,,,,,,,,,,,,,,,,,,,,,,,

These are the sequence files in fasta format that are used in aligning and generating the fusion sequences.


.. _reference-files-template-metadata:

Template Metadata
,,,,,,,,,,,,,,,,,,,,,,,,

This is the file which contains the band information for the chromosomes. This is only used during visualization.


The structure of the file should look something like this

.. code-block:: text

    chr1    0       2300000 p36.33  gneg
    chr1    2300000 5400000 p36.32  gpos25
    chr1    5400000 7200000 p36.31  gneg
    chr1    7200000 9200000 p36.23  gpos25
    chr1    9200000 12700000        p36.22  gneg

.. _reference-files-masking:

Masking File
,,,,,,,,,,,,,,,,,,,,,,,

The masking file is a tab delimited file which contains regions that we should ignore calls in. 
This can be used to filter out regions with known false positives, bad mapping, centromeres, telomeres etc. 
An example of the expected format is shown below. The file should have four columns: chr, start, end and name.

.. code-block:: text

    #chr    start   end     name
    chr1    0       2300000 centromere
    chr1    9200000 12700000        telomere

The pre-built masking files in the downloads table above are telomere regions, centromere regions (based on the cytoband file), 
and nspan regions (computed with tools/find_repeats.py).

Masking is not required (can provide a header-only file), but is recommended as it will improve performance and specificity.

.. _reference-files-annotations:

Annotations
,,,,,,,,,,,,,,,,,,,,,,,

This is a custom file format. It is a :term:`JSON` file which contains the gene, transcript, exon,
translation and protein domain positional information

Pre-built annotation files can be downloaded above. The 'best transcript' flag is based on an in-house model, as are the
aliases (custom one-to-one hugo gene name mapping)

.. warning::

    the :func:`~mavis.annotate.file_io.load_reference_genes` will
    only load valid translations. If the cds sequence in the annotation is not
    a multiple of :attr:`~mavis.constants.CODON_SIZE` or if a
    reference genome (sequences) is given and the cds start and end are not
    M and * amino acids as expected the translation is not loaded

Example of the :term:`JSON` file structure can be seen below

.. code-block:: javascript

    [
        {
            "name": string,
            "start": int,
            "end": int
            "aliases": [string, string, ...],
            "transcripts": [
                {
                    "name": string,
                    "start": int,
                    "end": int,
                    "exons": [
                        {"start": int, "end": int, "name": string},
                        ...
                    ],
                    "cdna_coding_start": int,
                    "cdna_coding_end": int,
                    "domains": [
                        {
                            "name": string,
                            "regions": [
                                {"start" aa_start, "end": aa_end}
                            ],
                            "desc": string
                        },
                        ...
                    ]
                },
                ...
            ]
        },
        ...
    }

The provided files were generated with :ref:`Ensembl <Yates-2016>`, however it can be generated from any database with the 
necessary information so long as the above :term:`JSON` structure is respected.

.. _generate-reference-annotations:

Generating the Annotations from :ref:`Ensembl <Yates-2016>`
...............................................................

There is a helper script included with mavis to facilitate generating the custom annotations
file from an instance of the :ref:`Ensembl <Yates-2016>` database. This uses the :ref:`Ensembl <Yates-2016>` perl api to connect and
pull information from the database. This has been tested with both Ensembl69 and Ensembl79.

Instructions for downloading and installing the perl api can be found on the `ensembl site <http://www.ensembl.org/info/docs/api/api_installation.html>`_

1. **Make sure the ensembl perl api modules are added to the PERL5LIB environment variable**

.. code-block:: bash

   PERL5LIB=${PERL5LIB}:$HOME/ensembl_79/bioperl-live
   PERL5LIB=${PERL5LIB}:$HOME/ensembl_79/ensembl/modules
   PERL5LIB=${PERL5LIB}:$HOME/ensembl_79/ensembl-compara/modules
   PERL5LIB=${PERL5LIB}:$HOME/ensembl_79/ensembl-variation/modules
   PERL5LIB=${PERL5LIB}:$HOME/ensembl_79/ensembl-funcgen/modules
   export PERL5LIB

2. **Configure the environment variables to set defaults for the perl script**

.. code-block:: bash

   # required data files
   export HUGO_ENSEMBL_MAPPING=/path/to/mapping/file
   export BEST_TRANSCRIPTS=/path/to/transcripts/file

   # connection information for the ensembl local (or external) server
   export ENSEMBL_HOST=HOSTNAME
   export ENSEMBL_PASS=PASSWORD
   export ENSEMBL_USER=USERNAME
   export ENSEMBL_PORT=PORT_NUMBER

3. **Run the perl script**

you can view the help menu by running

.. code-block:: bash

    perl generate_ensembl_json.pl

you can override the default input file parameters (configured in the above step) by providing arguments
to the script itself

.. code-block:: bash

    perl generate_ensembl_json.pl --best_transcript_file /path/to/best/transcripts/file --output /path/to/output/json/file.json

or if you have configured the environment variables as given in step 2, then simply provide the output path

.. code-block:: bash

    perl generate_ensembl_json.pl --output /path/to/output/json/file.json


.. _reference-files-dgv-annotations:

:ref:`DGV (Database of Genomic Variants) <Macdonald-2014>` Annotations
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

The DGV annotations file contains regions corresponding to what is found in the database of genomic variants. This is
used to annotate events that are found in healthy control samples and therefore may not be of interest
if looking for somatic events. 

The above (downloads table) files were generated from from `DGV <http://dgv.tcag.ca/dgv/app/download>`_
and reformatted to have 4 columns after download. We used awk to convert the raw file

.. code-block:: bash

    awk '{print $2"\t"$3"\t"$4"\t"$1} GRCh37_hg19_variants_2016-05-15.txt > dgv_hg19_variants.tab

Note in hg19 the column is called "name" and in hg38 the column is called "variantaccession".
An example is shown below

.. code-block:: text

    #chr     start   end     name
    1       1       2300000 nsv482937
    1       10001   22118   dgv1n82
    1       10001   127330  nsv7879

.. _reference-files-aligner-reference:

Aligner Reference
,,,,,,,,,,,,,,,,,,,,,,,

The aligner reference file is the reference genome file used by the aligner during the validate stage. For example,
if :term:`blat` is the aligner then this will be a :term:`2bit` file.

