# SV Callers

MAVIS supports output from a wide-variety of SV callers. Assumptions are made for each tool based on interpretation of the output and the publications for each tool.

## Configuring Conversions

Adding a conversion step to your MAVIS run is as simple as adding that section to the input JSON config.

The general structure of this section is as follows

```jsonc
{
    "convert": {
        "<ALIAS>": {
            "file_type": "<TOOL OUTPUT TYPE>",
            "name": "<TOOL NAME>",  // optional field for supported tools
            "inputs": [
                "/path/to/tool/output/file"
            ]
        }
    }
}
```

A full version of the input configuration file specification can be found in the [configuration](../configuration/general.md) section.

## Supported Tools

The tools and versions currently supported are given below. Versions listed indicate the version of the tool for which output files have been tested as input into MAVIS. MAVIS also supports a [general VCF input](#general-vcf-inputs).

| SV Caller                                                                   | Version(s) Tested | Files used as MAVIS input                     |
| --------------------------------------------------------------------------- | ----------------- | --------------------------------------------- |
| [BreakDancer (Chen, 2009)](../../background/citations#chen-2009)            | `1.4.5`           | `Tools main output file(s)`                   |
| [BreakSeq (Abyzov, 2015)](../../background/citations#abyzov-2015)           | `2.2`             | `work/breakseq.vcf.gz`                        |
| [Chimerascan (Iyer, 2011)](../../background/citations#iyer-2011)            | `0.4.5`           | `*.bedpe`                                     |
| [CNVnator (Abyzov, 2011)](../../background/citations#abyzov-2011)           | `0.3.3`           | `Tools main output file(s)`                   |
| [CuteSV (Jiang, 2020)](../../background/citations#jiang-2020)               | `1.0.10`          | `*.vcf`                                       |
| [DeFuse (McPherson. 2011)](../../background/citations#mcpherson-2011)       | `0.6.2`           | `results/results.classify.tsv`                |
| [DELLY (Rausch, 2012)](../../background/citations#rausch-2012)              | `0.6.1` `0.7.3`   | `combined.vcf` (converted from bcf)           |
| [Manta (Chen, 2016)](../../background/citations#chen-2016)                  | `1.0.0`           | `{diploidSV,somaticSV}.vcf`                   |
| [Pindel (Ye, 2009)](../../background/citations#ye-2009)                     | `0.2.5b9`         | `Tools main output file(s)`                   |
| [Sniffles (Sedlazeck, 2018)](../../background/citations#sedlazeck-2018)     | `1.0.12b`         | `*.vcf`                                       |
| [STAR-Fusion (Haas, 2017)](../../background/citations#haas-2017)            | `1.4.0`           | `star-fusion.fusion_predictions.abridged.tsv` |
| [Straglr (Chiu, 2021)](../../background/citations#chiu-2021)                |                   |
| [Strelka (Saunders, 2012)](../../background/citations#saunders-2012)        | `1.0.6`           | `passed.somatic.indels.vcf`                   |
| [Trans-ABySS (Robertson, 2010)](../../background/citations/#robertson-2010) | `1.4.8 (custom)`  | `{indels/events_novel_exons,fusions/*}.tsv`   | `<output_prefix>.bed` |

!!! note
    [Trans-ABySS](../../glossary/#trans-abyss): The trans-abyss version
    used was an in-house dev version. However the output columns are
    compatible with 1.4.8 as that was the version branched from.
    Additionally, although indels can be used from both genome and
    transcriptome outputs of Trans-ABySS, it is recommended to only use the
    genome indel calls as the transcriptome indels calls (for versions
    tested) introduce a very high number of false positives. This will slow
    down validation. It is much faster to simply use the genome indels for
    both genome and transcriptome.

## [DELLY](../../glossary/#delly) Post-processing

Some post-processing on the delly output files is generally done prior
to input. The output BCF files are converted to a VCF file

```bash
bcftools concat -f /path/to/file/with/vcf/list --allow-overlaps --output-type v --output combined.vcf
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
[custom conversion script](#custom-conversions).

Using the general VCF tool with a non-standard tool can be done as follows

```json
{
    "convert": {
        "my_tool_alias": {
            "file_type": "vcf",
            "name": "my_tool",
            "inputs": ["/path/to/my_tool/output.vcf"]
        }
    }
}
```

### Assumptions on non-standard INFO fields

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

## Custom Conversions

If there is a tool that is not yet supported by MAVIS and you would like it to be, you can either add a [feature request](https://github.com/bcgsc/mavis/issues) to our GitHub page or tackle writing the conversion script yourself. Either way there are a few things you will need

- A sample output from the tool in question
- Tool metadata for the citation, version, etc

### Logic Example - [Chimerascan](../../glossary/#chimerascan)

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

### Calling a Custom Conversion Script

Since MAVIS v3+ is run using [snakemake](https://snakemake.readthedocs.io/en/stable/) the simplest way to incorporate your custom conversion scripts is to modify the Snakefile and add them as rules.
