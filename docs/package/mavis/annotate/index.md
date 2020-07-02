# Sub-package Documentation


## Types of Output Files

| expected name/suffix           | file type/format | content                                  |
| ------------------------------ | ---------------- | ---------------------------------------- |
| ``annotations.tab``            | text/tabbed      | annotated events                         |
| ``annotations.fusion-cdna.fa`` | :term:`fasta`    | putative fusion unspliced cDNA sequences |
| ``drawings/*.svg``             | :term:`SVG`      | diagrams                                 |
| ``drawings/*.legend.json``     | :term:`JSON`     | diagram legend/metadata                  |


## Algorithm Overview

see :ref:`theory - annotating events <theory-annotating-events>`

- read in breakpoint pairs
- generate strand-specific annotations (one annotation per strand, multiple if multiple genes/transcripts in the region)
- try building fusion transcripts for bp-specific calls
- generate :term:`SVG` diagrams

## Levels of Annotations

![](../../../images/feature_levels.svg)


## Overview of Class Relationships

![](../../../images/annotation_model.svg)

The Annotation sub-package has objects for genetic annotations and related calculations. The basic layout of the
package is shown above. IS-A relationships are given by the blue arrows. HAS-A relationships are shown in black.
And reference_object/parent
type relationships are shown in red. mavis.annotate.genomic.Gene is a gene. Start and end are
genomic positions wrt to the template/chr. mavis.annotate.genomic.PreTranscript is the
unspliced transcript. Start and end are genomic positions wrt to the template/chr.
mavis.annotate.genomic.Transcript: is the spliced transcript. Start and end coordinates are
1 to the length of the spliced product in base pairs.
mavis.annotate.protein.Translation: is the translation of the spliced transcript. Start and
end are cdna positions wrt the 5' end of the spliced transcript. The start and end here describe the start and end
of the coding sequence
