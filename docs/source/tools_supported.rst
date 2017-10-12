.. _supported-sv-callers:


Supported SV Callers
================================

MAVIS supports output from a wide-variety of SV callers. Assumptions are made for each tool based on interpretation of
the output and the publications for each tool. The tools and versions currently supported are given below. Versions listed
indicate the version of the tool for which output files have been tested as input into MAVIS


+-------------+-------------------+---------------------------------------------------------------------+--------------------+
| Name        | Publication       | Source                                                              | Versions Supported |
+=============+===================+=====================================================================+====================+
| Chimerascan | [Iyer-2011]_      | `code.google.com <https://code.google.com/archive/p/chimerascan>`_  | 0.4.5              |
+-------------+-------------------+---------------------------------------------------------------------+--------------------+
| DeFuse      | [McPherson-2011]_ | `bitbucket <https://bitbucket.org/dranew/defuse>`_                  | 0.6.2              |
+-------------+-------------------+---------------------------------------------------------------------+--------------------+
| DELLY       | [Rausch-2012]_    | `github <https://github.com/dellytools/delly>`_                     | 0.7.3              |
+-------------+-------------------+---------------------------------------------------------------------+--------------------+
| Manta       | [Chen-2016]_      | `github <https://github.com/Illumina/manta>`_                       | 1.0.0              |
+-------------+-------------------+---------------------------------------------------------------------+--------------------+
| Pindel      | [Ye-2009]_        | `github <https://github.com/genome/pindel>`_                        | 0.2.5b9            |
+-------------+-------------------+---------------------------------------------------------------------+--------------------+
| Trans-ABySS | [Robertson-2010]_ | `github <https://github.com/bcgsc/transabyss>`_                     | 1.4.8 (custom)     |
+-------------+-------------------+---------------------------------------------------------------------+--------------------+

.. note:: 

    The trans-abyss version used was an in-house dev version. However the output columns are compatible with 1.4.8 as that
    was the version branched from
