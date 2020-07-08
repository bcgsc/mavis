# Glossary

## flanking read pair

A pair of reads where one read maps to one side of a set of
breakpoints and its mate maps to the other.

## split read

A read which aligns next to a breakpoint and is softclipped at one
or more sides.

## spanning read

Applies primarily to small structural variants. Reads which span
both breakpoints.

## half-mapped read

A read whose mate is unaligned. Generally this refers to reads in
the evidence stage that are mapped next to a breakpoint.

## breakpoint

A breakpoint is a genomic position (interval) on some
reference/template/chromosome which has a strand and orientation.
The orientation describes the portion of the reference that is
retained.

## event

Used interchangeably with [structural variant](#structural-variant).

## event type

Classification for a structural variant. see
[event_type](#event_type).

## structural variant

A genomic alteration that can be described by a pair of breakpoints
and an [event type](#event type). The two
breakpoints represent regions in the genome that are broken apart
and reattached together.

## breakpoint pair

Basic definition of a [structural variant](#structural variant). Does not automatically imply a classification/type.

## bed

File format specification. See
[https://genome.ucsc.edu/FAQ/FAQformat#format1](https://genome.ucsc.edu/FAQ/FAQformat#format1).

## IGV batch file

This is a file format type defined by [IGV](#IGV) see [running IGV with a batch
file](https://software.broadinstitute.org/software/igv/batch).

## BAM

File format specification. See
[https://genome.ucsc.edu/FAQ/FAQformat#format5.1](https://genome.ucsc.edu/FAQ/FAQformat#format5.1).

## 2bit

File format specification. See
[https://genome.ucsc.edu/FAQ/FAQformat#format7](https://genome.ucsc.edu/FAQ/FAQformat#format7).

## fasta

File format specification. See
[https://genome.ucsc.edu/FAQ/FAQformat#format18](https://genome.ucsc.edu/FAQ/FAQformat#format18).

## psl

File format specification. See
[https://genome.ucsc.edu/FAQ/FAQformat#format2](https://genome.ucsc.edu/FAQ/FAQformat#format2).

## pslx

Extended format of a [psl](#psl).

## SVG

SVG (Scalable vector graph) is an image format. see
[https://www.w3schools.com/graphics/svg_intro.asp](https://www.w3schools.com/graphics/svg_intro.asp).

## JSON

JSON (JavaScript Object Notation) is a data file format. see
[https://www.w3schools.com/js/js_json_intro.asp](https://www.w3schools.com/js/js_json_intro.asp).

## blat

Alignment tool. see [https://genome.ucsc.edu/FAQ/FAQblat.html#blat3](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3)
for instructions on download and install.

## IGV

Integrative Genomics Viewer is a visualization tool. see
[http://software.broadinstitute.org/software/igv](http://software.broadinstitute.org/software/igv).

## SGE

Sun Grid Engine (SGE) is a job scheduling system for cluster
management see
[http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html](http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html).

## SLURM

SLURM is a job scheduling system for cluster management see
[https://slurm.schedmd.com/archive/slurm-17.02.1](https://slurm.schedmd.com/archive/slurm-17.02.1).

## TORQUE

TORQUE is a job scheduling system for cluster management see
[http://www.adaptivecomputing.com/products/open-source/torque/](http://www.adaptivecomputing.com/products/open-source/torque/).

## BWA

BWA is an alignement tool. See [https://github.com/lh3/bwa](https://github.com/lh3/bwa) for
install instructions.

## HGVS

Community based standard of reccommendations for variant notation.
See [http://varnomen.hgvs.org/](http://varnomen.hgvs.org/)

## BreakDancer

BreakDancer is an SV caller. Source for BreakDancer can be found
[here](https://github.com/genome/breakdancer) [Chen-2009](../background/citations#chen-2009)

## BreakSeq

BreakSeq is an SV caller. Source for BreakSeq can be found
[here](https://github.com/bioinform/breakseq2) [Abyzov-2015](../background/citations#abyzov-2015)

## Chimerascan

Chimerascan is an SV caller. Source for Chimerascan can be found
[here](https://code.google.com/archive/p/chimerascan)
[Iyer-2011](../background/citations#iyer-2011)

## CNVnator

CNVnator is an SV caller. Source for CNVnator can be found
[here](https://github.com/abyzovlab/CNVnator) [Abyzov-2011](../background/citations#abyzov-2011)

## DeFuse

DeFuse is an SV caller. Source for DeFuse can be found
[here](https://bitbucket.org/dranew/defuse) [McPherson-2011](../background/citations#mcpherson-2011)

## DELLY

DELLY is an SV caller. Source for DELLY can be found
[here](https://github.com/dellytools/delly) [Rausch-2012](../background/citations#rausch-2012)

## Manta

Manta is an SV caller. Source for Manta can be found
[here](https://github.com/Illumina/manta) [Chen-2016](../background/citations#chen-2016)

## STAR-Fusion

STAR-Fusion is an SV caller. Source for STAR-Fusion can be found
[here](https://github.com/STAR-Fusion/STAR-Fusion) [Haas-2017](../background/citations#haas-2017)

## Strelka

Strelka is an SNV and small indel caller. Only the small indels can
be processed, since SNVs are not currently suported. Source for
Strelka can be found [here](https://github.com/Illumina/strelka)
[Saunders-2012](../background/citations#saunders-2012)

## Pindel

Pindel is an SV caller. Source for Pindel can be found
[here](https://github.com/genome/pindel) [Ye-2009](../background/citations#ye-2009)

## Trans-ABySS

Trans-ABySS is an SV caller. Source for Trans-ABySS can be found
[here](https://github.com/bcgsc/transabyss) [Robertson-2010](../background/citations#robertson-2010)

## SV

Structural Variant
