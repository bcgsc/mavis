"""
Sub-package Documentation
============================

This is the package responsible for summarizing the calls between libraries. In many cases
this will be where somatic vs germline is determined or genomic only vs expressed.

Output Files
--------------

+----------------------------+------------------+------------------------------------------------------------+
| expected name/suffix       | file type/format | content                                                    |
+============================+==================+============================================================+
| ``mavis_summary_*.tab``    | text/tabbed      |  ?                                                         |
+----------------------------+------------------+------------------------------------------------------------+


Algorithm Overview
---------------------
TODO

"""
from .summary import annotate_dgv, filter_by_annotations, filter_by_call_method, filter_by_evidence, group_events
