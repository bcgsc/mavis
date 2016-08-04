import re
import TSV
from constants import ORIENT
from constants import STRAND
from sv import Breakpoint, BreakpointPair

def load_input_file(filename):
    header, rows = TSV.read_file(filename,
            retain=[
                'start_chromosome', 
                'end_chromosome', 
                'start_orientation', 
                'end_orientation', 
                'start_strand', 
                'end_strand',
                'protocol',
                'tool_version'],
            split={
                'start_position': ('^(\d+)-(\d+)$', ['start_pos1', 'start_pos2']),
                'end_position': ('^(\d+)-(\d+)$', ['end_pos1', 'end_pos2']),
                },
            cast={'start_pos1': 'int', 'start_pos2': 'int', 'end_pos1': 'int', 'end_pos2': 'int'},
            validate={
                'start_orientation': '^{0}$'.format('|'.join([ re.escape(x) for x in ORIENT.values()])),
                'end_orientation': '^{0}$'.format('|'.join([ re.escape(x) for x in ORIENT.values()])),
                'start_strand': '^{0}$'.format('|'.join([ re.escape(x) for x in STRAND.values()])),
                'end_strand': '^{0}$'.format('|'.join([ re.escape(x) for x in STRAND.values()])),
                'tool_version': '^.+_v\d+\.\d+\.\d+$',
                'protocol': '^(genome|transcriptome)$',
                'libraries': '^[\w-]+(;[\w-]+)*$'
                },
            row_index = 'row_index'
            )
    breakpoints = []
    for row in rows:
        label = row['tool_version'] + '_' + str(row['row_index'])
        b1 = Breakpoint(
                row['start_chromosome'], 
                row['start_pos1'], 
                row['start_pos2'], 
                row['start_orientation'], 
                row['start_strand'],
                label = label)
        b2 = Breakpoint(
                row['end_chromosome'], 
                row['end_pos1'], 
                row['end_pos2'], 
                row['end_orientation'], 
                row['end_strand'],
                label = label)
        breakpoints.append(BreakpointPair(b1, b2))
    return breakpoints

