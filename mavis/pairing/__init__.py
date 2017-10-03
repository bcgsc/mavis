"""
Sub-package Documentation
============================

This is the package responsible for pairing/grouping calls between libraries. In many cases
this will be where somatic vs germline is determined or genomic only vs expressed.

Output Files
--------------

+----------------------------+------------------+------------------------------------------------------------+
| expected name/suffix       | file type/format | content                                                    |
+============================+==================+============================================================+
| ``mavis_paired_*.tab``     | text/tabbed      |  call information and pairing information using product id |
+----------------------------+------------------+------------------------------------------------------------+


Algorithm Overview
---------------------

- pairwise comparison of breakpoint pairs between libraries

    - fail if orientations do not match
    - fail if template/chromosomes do not match
    - if the protocols are mixed

        - pass if the fusion products match at the sequence level
        - pass if the breakpoint predicted from the genome matches the transcriptome breakpoint

    - if the protocols are the same

        - pass if the breakpoints are co-located

- filter matches based on annotations

    - if both breakpoints have the same gene annotation, they must also both have the same transcript annotation

"""
