# Reference Input Files

There are several reference files that are required for full
functionality of the MAVIS pipeline. If the same reference file will be
reused often then the user may find it helpful to set reasonable
defaults. Default values for any of the reference file arguments can be
[configured through environment variables](../../configuration/general/#environment-variables)

To improve the install experience for the users, different
configurations of the MAVIS annotations file have been made available.
These files can be downloaded below, or if the required configuration is
not available,
(instructions on generating the annotations file)[/inputs/reference/#generating-the-annotations-from-ensembl] can be found below.

| File Name (Type/Format)                                                                       | Environment Variable      | Download                                                                                                                                                                                                                                                      |
| --------------------------------------------------------------------------------------------- | ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [reference genome](../../inputs/reference/#reference-genome) ([fasta](../../glossary/#fasta)) | `MAVIS_REFERENCE_GENOME`  | [![](../images/get_app-24px.svg) GRCh37/Hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz) <br> [![](../images/get_app-24px.svg) GRCh38](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.tar.gz)                        |
| [annotations](../../inputs/reference/#annotations) ([JSON](../../glossary/#json))             | `MAVIS_ANNOTATIONS`       | [![](../images/get_app-24px.svg) GRCh37/Hg19 + Ensembl69](http://www.bcgsc.ca/downloads/mavis/ensembl69_hg19_annotations.json) <br> [![](../images/get_app-24px.svg) GRCh38 + Ensembl79](http://www.bcgsc.ca/downloads/mavis/ensembl79_hg38_annotations.json) |
| [masking](../../inputs/reference/#masking-file) (text/tabbed)                                 | `MAVIS_MASKING`           | [![](../images/get_app-24px.svg) GRCh37/Hg19](http://www.bcgsc.ca/downloads/mavis/hg19_masking.tab)<br>[![](../images/get_app-24px.svg) GRCh38](http://www.bcgsc.ca/downloads/mavis/GRCh38_masking.tab)                                                       |
| [template metadata](../../inputs/reference/#template-metadata) (text/tabbed)                  | `MAVIS_TEMPLATE_METADATA` | [![](../images/get_app-24px.svg) GRCh37/Hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz)<br>[![](../images/get_app-24px.svg) GRCh38](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz)                   |
| [DGV annotations](../../inputs/reference/#dgv-database-of-genomic-variants) (text/tabbed)     | `MAVIS_DGV_ANNOTATION`    | [![](../images/get_app-24px.svg) GRCh37/Hg19](http://www.bcgsc.ca/downloads/mavis/dgv_hg19_variants.tab)<br>[![](../images/get_app-24px.svg) GRCh38](http://www.bcgsc.ca/downloads/mavis/dgv_hg38_variants.tab)                                               |
| [aligner reference](../../inputs/reference/#aligner-reference)                                | `MAVIS_ALIGNER_REFERENCE` | [![](../images/get_app-24px.svg) GRCh37/Hg19 2bit (blat)](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit)<br>[![](../images/get_app-24px.svg) GRCh38 2bit (blat)](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit)         |


If the environment variables above are set they will be used as the
default values when any step of the pipeline script is called (including
generating the template config file)

## Reference Genome**

These are the sequence files in fasta format that are used in aligning
and generating the fusion sequences.

## Template Metadata

This is the file which contains the band information for the
chromosomes. This is only used during visualization.

The structure of the file should look something like this

    chr1    0       2300000 p36.33  gneg
    chr1    2300000 5400000 p36.32  gpos25
    chr1    5400000 7200000 p36.31  gneg
    chr1    7200000 9200000 p36.23  gpos25
    chr1    9200000 12700000        p36.22  gneg

## Masking File

The masking file is a tab delimited file which contains regions that we
should ignore calls in. This can be used to filter out regions with
known false positives, bad mapping, centromeres, telomeres etc. An
example of the expected format is shown below. The file should have four
columns: chr, start, end and name.

    #chr    start   end     name
    chr1    0       2300000 centromere
    chr1    9200000 12700000        telomere

The pre-built masking files in the downloads table above are telomere
regions, centromere regions (based on the cytoband file), and nspan
regions (computed with tools/find\_repeats.py).

Masking is not required (can provide a header-only file), but is
recommended as it will improve performance and specificity.

## Annotations

This is a custom file format. It is a [JSON](../../glossary/#json) file which contains the gene, transcript, exon, translation
and protein domain positional information

Pre-built annotation files can be downloaded above. The 'best
transcript' flag is based on an in-house model. We have also pre-built
the ensembl annotations file including non-coding transcripts below.

!!! warning
    It is worth noting that using the reference annotation file including
    the non-coding genes will require an increase in the default amount of
    memory for the annotation step due to the increased size of the
    annotations file. On our standard COLO829 we increased the default
    memory for the annotation step from 12G to 18G.

[![](../images/get_app-24px.svg) GRCh37/Hg19 + Ensembl69 (includes non-coding genes)](http://www.bcgsc.ca/downloads/mavis/ensembl69_hg19_annotations_with_ncrna.json)


!!! warning
    the `mavis.annotate.file_io.load_reference_genes`{.interpreted-text
    role="func"} will only load valid translations. If the cds sequence in
    the annotation is not a multiple of
    `mavis.constants.CODON_SIZE`{.interpreted-text role="attr"} or if a
    reference genome (sequences) is given and the cds start and end are not
    M and \* amino acids as expected the translation is not loaded

Example of the [JSON](../../glossary/#json) file structure can
be seen below

```text
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
```

The provided files were generated with
[Ensembl](../../background/citations/#yates-2016), however it can be
generated from any database with the necessary information so long as
the above [JSON](../../glossary/#json) structure is respected.

### Generating the Annotations from [Ensembl](../../background/citations/#yates-2016)

There is a helper script included with mavis to facilitate generating
the custom annotations file from an instance of the
[Ensembl](../../background/citations/#yates-2016) database. This uses
the [Ensembl](../../background/citations/#yates-2016) perl api to
connect and pull information from the database. This has been tested
with both Ensembl69 and Ensembl79.

Clone and install the mavis tools dependencies

```bash
git clone https://github.com/bcgsc/mavis.git
cd mavis
python3 -m venv venv
source venv/bin/activate
pip install .[tools]  # could also install just the dependencies for the tools
```

Run the Python script choosing your species and Ensembl version

```bash
python tools/generate_ensembl_json.py -s human -r 75 -o ensembl_human_v75.json
```

This will produce the JSON file required as input by MAVIS


## DGV (Database of Genomic Variants)

The DGV annotations file contains regions corresponding to what is found
in the database of genomic variants. This is used to annotate events
that are found in healthy control samples and therefore may not be of
interest if looking for somatic events.

The above (downloads table) files were generated from from
[DGV](http://dgv.tcag.ca/dgv/app/download) and reformatted to have 4
columns after download. We used awk to convert the raw file

```bash
awk '{print $2"\t"$3"\t"$4"\t"$1} GRCh37_hg19_variants_2016-05-15.txt > dgv_hg19_variants.tab
```

Note in hg19 the column is called "name" and in hg38 the column is
called "variantaccession". An example is shown below

    #chr     start   end     name
    1       1       2300000 nsv482937
    1       10001   22118   dgv1n82
    1       10001   127330  nsv7879

## Aligner Reference

The aligner reference file is the reference genome file used by the
aligner during the validate stage. For example, if
[blat](../../glossary#blat) is the aligner then this will be a
[2bit](../../glossary#2bit) file.
