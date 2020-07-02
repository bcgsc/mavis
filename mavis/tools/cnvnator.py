"""
from cnvnator: https://github.com/abyzovlab/CNVnator

    CNV_type coordinates CNV_size normalized_RD e-val1 e-val2 e-val3 e-val4 q0

    normalized_RD -- normalized to 1.
    e-val1        -- is calculated using t-test statistics.
    e-val2        -- is from the probability of RD values within the region to be in
    the tails of a gaussian distribution describing frequencies of RD values in bins.
    e-val3        -- same as e-val1 but for the middle of CNV
    e-val4        -- same as e-val2 but for the middle of CNV
    q0            -- fraction of reads mapped with q0 quality
"""
import re


def convert_row(row):
    """

    Args:
        row (Dict[str]): dict representing the row output from cnvnator

    Returns:
        dict: transformed row using mavis starndard column names
    """
    result = {k: v for k, v in row.items() if k != 'coordinates'}
    chrom, start, end = re.split(r'[-:]', row['coordinates'])
    result['break1_chromosome'] = result['break2_chromosome'] = chrom
    result['break1_position_start'] = result['break1_position_end'] = start
    result['break2_position_start'] = result['break2_position_end'] = end
    return result
