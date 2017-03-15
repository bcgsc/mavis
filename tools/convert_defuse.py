#! /projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v4.3.0/envs/python3.4/bin/python
import argparse
import sys, os
import pandas as pd
import time

versions = ['0.6.1']

def make_tsv(pogid, tsv, libraryid,version=None, outputDir=""):
    '''
    given   a pogid, libraryid and its deFuse results (results.filtered.tsv)
    return  another tsv with subset columns & orientations
    '''
    df = pd.read_table(tsv, sep='\t')
    nrows = len(df)
    outfile = outputDir + pogid + '.defuse_svmerge.tsv'

    # columns to subset from deFuse results tsv
    info = ['span_count', 'splitr_count', 'cluster_id', 'probability']
    constants = ['tool_version', 'protocol', 'library', 'stranded']

    # determine the orientation from the strands of the two genes
    df_o = {'orientation1':[], 'orientation2':[],'opposing_strand':[]}
    for i in range(nrows):
        o1 = 'L' if df.iloc[i]['genomic_strand1'] == '+' else 'R'
        o2 = 'L' if df.iloc[i]['genomic_strand2'] == '+' else 'R'

        df_o['orientation1'].append(o1)
        df_o['orientation2'].append(o2)

        #use the opposite here since it is based on genes (which act similar to reads) rather than contigs
        if o1 == o2:
            df_o['opposing_strand'].append(True)
        else:
            df_o['opposing_strand'].append(False)
    df_o = pd.DataFrame(df_o)

    # add tool and tool_version columns
    df_end = pd.DataFrame({'tool_version':['deFuse_v{}'.format(version)]*nrows,
                           'protocol':['transcriptome']*nrows,
                           'library':[libraryid]*nrows,
                           'stranded':[False]*nrows})

    # combine columns and write out tsv
    frames = [df['gene_chromosome1'], df['genomic_break_pos1'], df['genomic_break_pos1'], df_o['orientation1'],
              df['gene_chromosome2'], df['genomic_break_pos2'], df['genomic_break_pos2'], df_o['orientation2'],
              df_o['opposing_strand'],
              df_end[constants],
              df[info]]
    result = pd.concat(frames, axis=1)
    result.columns = ['break1_chromosome', 'break1_position_start', 'break1_position_end', 'start_orientation',
                      'break2_chromosome', 'break2_position_start', 'break2_position_end','end_orientation',
                      'opposing_strands',
                      'tools', 'protocol', 'library', 'stranded',
                      'spanning_read_count', 'split_read_count', 'cluster_id','probability']

    with open(outfile, 'w') as f:
        f.write('## inputs: {}\n'.format(" ".join(sys.argv)))
        f.write('## file generated on {}\n'.format(time.strftime('%B %d, %Y')))
        result.to_csv(f, sep='\t', index=False)

    print("Wrote {} gene fusion events to {}".format(len(result.index), outfile))

def __main__():
    parser = argparse.ArgumentParser(description="Pulls gene breakpoint coordinates, strands, and orientations from deFuse result tsv and generates the sv merge input file", add_help=False,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('path_to_result', type=str,
                        help='absolute path to results.filtered.tsv from deFuse')
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--id', help='The id used in the deFuse run.', type=str, required=True)
    required.add_argument('-l', '--library_id', help='The library id for the input bam file.', required=True)

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-h','--help', action='help', help='Show this help message and exit')
    optional.add_argument('-v', '--version', help='the version of defuse that was used in the analysis', default='0.6.1')
    optional.add_argument('-o', '--outdir', help='the directory where the output will be placed', default='')

    args = parser.parse_args()

    result = args.path_to_result
    pogid = args.id
    libraryid = args.library_id
    outputDir = args.outdir
    version = args.version
    if version not in versions:
        print("ERROR: The defuse version given is not valid")
        sys.exit()

    if os.path.isfile(result):
        make_tsv(pogid, result, libraryid, version, outputDir)
    else:
        print("ERROR: Cannot find file: " + result)

if __name__ == "__main__":
    __main__()
