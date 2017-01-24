"""

Overview
-----------

.. figure:: _static/svmerge_annotation_model.svg
    :width: 100%

    The Annotation subpackage has objects for genetic annotations and related calculations. The basic layout of the 
    package is shown above. IS-A relationships are given by the blue arrows. HAS-A relationships are shown in black. 
    And reference_object/parent
    type relationships are shown in red. :class:`~structural_variant.annotate.genomic.Gene` is a gene. Start and end are 
    genomic positions wrt to the template/chr. :class:`~structural_variant.annotate.genomic.usTranscript` is the 
    unspliced transcript. Start and end are genomic positions wrt to the template/chr. 
    :class:`~structural_variant.annotate.genomic.Transcript`: is the spliced transcript. Start and end coordinates are 
    1 to the length of the spliced product in base pairs. 
    :class:`~structural_variant.annotate.protein.Translation`: is the translation of the spliced transcript. Start and 
    end are cdna positions wrt the 5' end of the spliced transcript. The start and end here describe the start and end 
    of the coding sequence
"""

from .file_io import *
