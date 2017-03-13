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


"""
from .pairing import equivalent_events
from .main import main
