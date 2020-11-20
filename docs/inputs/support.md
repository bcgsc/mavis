# Supported Dependencies

MAVIS integrates with
[SV callers](../../inputs/support/#sv-callers),
[job schedulers](../../inputs/support/#job-schedulers), and
[aligners](../../inputs/support/#aligners). While some of
these dependencies are optional, all currently supported options are
detailed below. The versions column in the tables below list all the
versions which were tested for each tool. Each version listed is known
to be compatible with MAVIS.

## Job Schedulers

MAVIS can be run locally without a job scheduler
(`MAVIS_SCHEDULER=LOCAL`) however, due to the computational resources
generally required, it is recommended that you use one of the supported
schedulers listed below.

| Name                             | Version(s)  | Environment Setting      |
| -------------------------------- | ----------- | ------------------------ |
| [TORQUE](../../glossary/#torque) | `6.1.2`     | `MAVIS_SCHEDULER=TORQUE` |
| [SGE](../../glossary/#sge)       | `8.1.8`     | `MAVIS_SCHEDULER=SGE`    |
| [SLURM](../../glossary/#slurm)   | `17.02.1-2` | `MAVIS_SCHEDULER=SLURM`  |

Users requiring support for other schedulers may make a request by
[submitting an issue to our github
page](https://github.com/bcgsc/mavis/issues). Additionally, developers
looking to extend the functionality may submit a pull request (Please
see the
[guidelines for contributors](../../development/)

MAVIS running locally uses the python
`concurrent.futures` library to manage
jobs.
##

## Aligners

Two aligners are supported [bwa](../../glossary/#bwa) and
[blat](../../glossary/#blat) (default).

| Name                                           | Version(s)              | Environment Setting       |
| ---------------------------------------------- | ----------------------- | ------------------------- |
| [blat](../../glossary/#blat)                   | `36x2` `36`             | `MAVIS_ALIGNER=blat`      |
| [bwa mem <bwa>](../../glossary/#bwa mem <bwa>) | `0.7.15-r1140` `0.7.12` | `MAVIS_ALIGNER='bwa mem'` |

!!! note
    When setting the aligner you will also need to set the
    [aligner_reference](../../configuration/settings/#aligner_reference) to match

## SV Callers

MAVIS supports output from a wide-variety of [SV](../../glossary/#sv) callers. Assumptions are made for each tool based on
interpretation of the output and the publications for each tool. The
tools and versions currently supported are given below. Versions listed
indicate the version of the tool for which output files have been tested
as input into MAVIS

MAVIS also supports a [general VCF input](../../inputs/support/#general-vcf-inputs).
It should be noted however that the tool tracked will only be listed as
'vcf' then.

| Name                                       | Version(s)       | MAVIS input                                   | Publication                                                 |
| ------------------------------------------ | ---------------- | --------------------------------------------- | ----------------------------------------------------------- |
| [BreakDancer](../../glossary/#breakdancer) | `1.4.5`          | `Tools main output file(s)`                   | [Chen-2009](../../background/citations#chen-2009)           |
| [BreakSeq](../../glossary/#breakseq)       | `2.2`            | `work/breakseq.vcf.gz`                        | [Abyzov-2015](../../background/citations#abyzov-2015)       |
| [Chimerascan](../../glossary/#chimerascan) | `0.4.5`          | `*.bedpe`                                     | [Iyer-2011](../../background/citations#Iyer-2011)           |
| [CNVnator](../../glossary/#cnvnator)       | `0.3.3`          | `Tools main output file(s)`                   | [Abyzov-2011](../../background/citations#abyzov-2011)       |
| [DeFuse](../../glossary/#defuse)           | `0.6.2`          | `results/results.classify.tsv`                | [McPherson-2011](../../background/citations#mcpherson-2011) |
| [DELLY](../../glossary/#delly)             | `0.6.1` `0.7.3`  | `combined.vcf` (converted from bcf)           | [Rausch-2012](../../background/citations#rausch-2012)       |
| [Manta](../../glossary/#manta)             | `1.0.0`          | `{diploidSV,somaticSV}.vcf`                   | [Chen-2016](../../background/citations#chen-2016)           |
| [Pindel](../../glossary/#pindel)           | `0.2.5b9`        | `Tools main output file(s)`                   | [Ye-2009](../../background/citations#ye-2009)               |
| [Trans-ABySS](../../glossary/#trans-abyss) | `1.4.8 (custom)` | `{indels/events_novel_exons,fusions/*}.tsv`   | [Robertson-2010](../../background/citations#robertson-2010) |
| [Strelka](../../glossary/#strelka)         | `1.0.6`          | `passed.somatic.indels.vcf`                   | [Saunders-2012](../../background/citations#saunders-2012)   |
| [STAR-Fusion](../../glossary/#star-fusion) | `1.4.0`          | `star-fusion.fusion_predictions.abridged.tsv` | [Haas-2017](../../background/citations#haas-2017)           |

!!! note
    [Trans-ABySS](../../glossary/#trans-abyss): The trans-abyss version
    used was an in-house dev version. However the output columns are
    compatible with 1.4.8 as that was the version branched from.
    Additionally, although indels can be used from both genome and
    transcriptome outputs of Trans-ABySS, it is reccommended to only use the
    genome indel calls as the transcriptome indels calls (for versions
    tested) introduce a very high number of false positives. This will slow
    down validation. It is much faster to simply use the genome indels for
    both genome and transcriptome.

### [DELLY](../../glossary/#delly) Post-processing

Some post-processing on the delly output files is generally done prior
to input. The output BCF files are converted to a VCF file

```bash
bcftools concat -f /path/to/file/with/vcf/list --allow-overlaps --output-type v --output combined.vcf
```

### Writing A Custom Conversion Script

#### Logic Example - [Chimerascan](../../glossary/#chimerascan)

The following is a description of how the conversion script for
[Chimerascan](../../background/citations/#iyer-2011) was generated.
While this is a built-in conversion command now, the logic could also
have been put in an external script. As mentioned above, there are a
number of assumptions that had to be made about the tools output to
convert it to the
[standard mavis format](../../inputs/standard/). Assumptions were then verified by reviewing at a series of
called events in [IGV](../../glossary/#igv). In the current
example, [Chimerascan](../../background/citations/#iyer-2011) output
has six columns of interest that were used in the conversion

- start3p
- end3p
- strand3p
- start5p
- end5p
- strand5p

The above columns describe two segments which are joined. MAVIS requires
the position of the join. It was assumed that the segments are always
joined as a [sense fusion](../../glossary/#sense-fusion). Using this
assumption there are four logical cases to determine the position of the
breakpoints.

i.e. the first case would be: If both strands are positive, then the end
of the five-prime segment (end5p) is the first breakpoint and the start
of the three-prime segment is the second breakpoint

### Calling A Custom Conversion Script

Custom conversion scripts can be specified during
[automatic config generation](../../configuration/general/#pipeline-configuration-file)
using the `--external_conversion` option.

!!! note
    Any external conversion scripts must take a `-o` option which requires a
    single outputfile argument. This outputfile must be the converted file
    output by the script. Additionally, the conversion script must be
    specified by its full path name and have executeable permissions.

In the following example the user has created a custom conversion script
`my_convert_script.py` which they are passing an input file named
`my_input1.txt`.

```bash
mavis config --external_conversion my_converted_input1 "my_convert_script.py my_input1.txt ... "
```

This will then be called during the pipeline step as

```bash
my_convert_script.py my_input1.txt ... -o /path/to/output/dir/converted_inputs/my_converted_input1.tab
```

You can also re-use the same conversion script if you have multiple
inputs to convert simply by specifying a different alias

```bash
mavis config \
    --external_conversion my_converted_input1 "my_convert_script.py my_input1.txt" \
    --external_conversion my_converted_input2 "my_convert_script.py my_input2.txt"
```

## General VCF inputs

Assuming that the tool outputting the VCF file follows standard
conventions, then it is possible to use a
[general VCF conversion](../../package/mavis/tools/vcf)
that is not tool-specific. Given the wide variety in content for VCF files,
MAVIS makes a number of assumptions and the VCF conversion may not work
for all VCFs. In general MAVIS follows the [VCF 4.2
specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf). If the
input tool you are using differs, it would be better to use a
[custom conversion script](#calling-a-custom-conversion-script).

**Assumptions on non-standard INFO fields**

- `PRECISE` if given, Confidence intervals are ignored if given in favour of exact breakpoint calls using pos and END as the breakpoint positions
- `CT` values if given are representative of the breakpoint orientations.
- `CHR2` is given for all interchromosomal events

### Translating BND type Alt fields

There are four possible configurations for the alt field of a BND type structural variant
based on the VCF specification. These correspond 1-1 to the orientation types for MAVIS
translocation structural variants.

```text
r = reference base/seq
u = untemplated sequence/alternate sequence
p = chromosome:position
```

| alt format | orients |
| ---------- | ------- |
| `ru[p[`    | LR      |
| `[p[ur`    | RR      |
| `]p]ur`    | RL      |
| `ru]p]`    | LL      |
