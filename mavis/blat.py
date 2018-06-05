"""


::

    In general the coordinates in psl files are “zero based half open.” The first base in a sequence is numbered
    zero rather than one. When representing a range the end coordinate is not included in the range. Thus the first
    100 bases of a sequence are represented as 0-100, and the second 100 bases are represented as 100-200. There is
    another little unusual feature in the .psl format. It has to do with how coordinates are handled on the
    negative strand. In the qStart/qEnd fields the coordinates are where it matches from the point of view of the forward
    strand (even when the match is on the reverse strand). However on the qStarts[] list, the coordinates are reversed.

-- http://wiki.bits.vib.be/index.php/Blat

"""
import logging
import math
import re

import tab

from .align import query_coverage_interval
from .bam import cigar as _cigar
from .bam.cigar import QUERY_ALIGNED_STATES
from .bam.read import SamRead
from .constants import CIGAR, DNA_ALPHABET, NA_MAPPING_QUALITY, PYSAM_READ_FLAGS, reverse_complement, STRAND
from .util import LOG
from .interval import Interval


class Blat:
    """
    """
    @staticmethod
    def millibad(row, is_protein=False, is_mrna=True):
        """
        this function is used in calculating percent identity
        direct translation of the perl code
        # https://genome.ucsc.edu/FAQ/FAQblat.html#blat4
        """
        size_mul = 1 if not is_protein else 3
        if is_protein and is_mrna:
            raise AttributeError('cannot be a protein AND and mRNA')
        q_ali_size = size_mul * (row['qend'] - row['qstart'])
        t_ali_size = row['tend'] - row['tstart']
        ali_size = min(q_ali_size, t_ali_size)
        if ali_size <= 0:
            return 0

        size_dif = q_ali_size - t_ali_size
        if size_dif < 0:
            if is_mrna:
                size_dif = 0
            else:
                size_dif = abs(size_dif)

        insert_factor = row['qgap_count']
        if not is_mrna:
            insert_factor += row['tgap_count']

        total = size_mul * (row['match'] + row['repmatch'] + row['mismatch'])
        if total != 0:
            millibad = 0
            round_away_from_zero = 3 * math.log(1 + size_dif)
            if round_away_from_zero < 0:
                round_away_from_zero = int(round_away_from_zero - 0.5)
            else:
                round_away_from_zero = int(round_away_from_zero + 0.5)
            millibad = (1000 * (row['mismatch'] * size_mul + insert_factor + round_away_from_zero)) / total
            return millibad
        else:
            return 0

    @staticmethod
    def score(row, is_protein=False):
        """
        direct translation from ucsc guidelines on replicating the web blat score
        https://genome.ucsc.edu/FAQ/FAQblat.html#blat4

        below are lines from the perl code i've re-written in python

        ::

            my $sizeMul = pslIsProtein($blockCount, $strand, $tStart, $tEnd, $tSize, $tStarts, $blockSizes);
            sizmul = 1 for DNA
            my $pslScore = $sizeMul * ($matches + ($repMatches >> 1) ) - $sizeMul * $misMatches - $qNumInsert - $tNumIns
                ert)
        """

        size_mul = 1 if not is_protein else 3
        score = size_mul * (row['match'] + (row['repmatch'] >> 1)) \
            - size_mul * row['mismatch'] \
            - row['qgap_count'] \
            - row['tgap_count']
        return score

    @staticmethod
    def percent_identity(row, is_protein=False, is_mrna=True):
        return 100 - int(Blat.millibad(row, is_protein, is_mrna)) * 0.1

    @staticmethod
    def read_pslx(filename, seqid_to_sequence_mapping, is_protein=False, verbose=True):
        pslx_header = [
            'match', 'mismatch', 'repmatch', 'ncount',
            'qgap_count', 'qgap_bases',
            'tgap_count', 'tgap_bases',
            'strand',
            'qname', 'qsize', 'qstart', 'qend',
            'tname', 'tsize', 'tstart', 'tend',
            'block_count', 'block_sizes',
            'qstarts', 'tstarts',
            'qseqs', 'tseqs'
        ]

        def split_csv_trailing_seq(x):
            return [s.upper() for s in re.sub(',$', '', x).split(',')]

        def split_csv_trailing_ints(x):
            return [int(s) for s in re.sub(',$', '', x).split(',')]

        header, rows = tab.read_file(
            filename,
            header=pslx_header,
            cast={
                'match': int,
                'mismatch': int,
                'repmatch': int,
                'ncount': int,
                'qgap_count': int,
                'qgap_bases': int,
                'tgap_count': int,
                'tgap_bases': int,
                'qsize': int,
                'qstart': int,
                'qend': int,
                'tsize': int,
                'tstart': int,
                'tend': int,
                'block_count': int,
                'tname': lambda x: re.sub('^chr', '', x),
                'block_sizes': split_csv_trailing_ints,
                'qstarts': split_csv_trailing_ints,
                'tstarts': split_csv_trailing_ints,
                'qseqs': split_csv_trailing_seq,
                'tseqs': split_csv_trailing_seq
            },
            validate={
                'strand': r'^[\+-]$'
            }
        )

        final_rows = []
        for row in rows:
            try:
                row['score'] = Blat.score(row, is_protein=is_protein)
                row['percent_ident'] = Blat.percent_identity(row, is_protein=is_protein)
                qseq = seqid_to_sequence_mapping[row['qname']]
                row['qseq_full'] = qseq

                for x in [
                    'qgap_count', 'qgap_bases', 'tgap_count',
                    'tgap_bases', 'qsize', 'tsize', 'ncount',
                    'match', 'mismatch', 'repmatch'
                ]:
                    if row[x] < 0 and verbose:
                        raise AssertionError(
                            'Blat error: blat returned a negative number, which are not allowed: {}={}'.format(
                                x, row[x]))
                final_rows.append(row)
            except AssertionError as err:
                LOG(type(err), ':', str(err), level=logging.DEBUG)
        return header, final_rows

    @staticmethod
    def pslx_row_to_pysam(row, bam_cache, reference_genome):
        """
        given a 'row' from reading a pslx file. converts the row to a BlatAlignedSegment object

        Args:
            row (dict of str): a row object from the 'read_pslx' method
            bam_cache (BamCache): the bam file/cache to use as a template for creating reference_id from chr name
            reference_genome (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`):
              dict of reference sequence by template/chr name

        """
        chrom = bam_cache.reference_id(row['tname'])
        query_sequence = row['qseq_full']
        if row['strand'] == STRAND.NEG:
            query_sequence = reverse_complement(query_sequence)
            temp = [q + b for q, b in zip(row['qstarts'], row['block_sizes'])]
            temp = [len(query_sequence) - q for q in temp][::-1]

        # note: converting to inclusive range [] vs end-exclusive [)
        reference_sequence = reference_genome[row['tname']].seq if reference_genome else None
        query_ranges = [Interval(x, x + y - 1) for x, y in zip(row['qstarts'], row['block_sizes'])]
        ref_ranges = [Interval(x, x + y - 1) for x, y in zip(row['tstarts'], row['block_sizes'])]

        for i in range(1, len(query_ranges)):
            # first check for blat errors
            if query_ranges[i].start <= query_ranges[i - 1].end or ref_ranges[i].start <= ref_ranges[i - 1].end:
                raise AssertionError('block ranges overlap in row', row)

        # try extending by consuming from the next aligned portion
        if reference_sequence:
            i = 0
            new_query_ranges = []
            new_ref_ranges = []

            while i < len(query_ranges):
                qpos = query_ranges[i][1] + 1
                rpos = ref_ranges[i][1] + 1
                shift = 0
                while qpos + shift < len(query_sequence) and rpos + shift < len(reference_sequence):
                    if DNA_ALPHABET.match(query_sequence[qpos + shift], reference_sequence[rpos + shift]):
                        shift += 1
                    else:
                        break

                if shift == 0:
                    new_query_ranges.append(query_ranges[i])
                    new_ref_ranges.append(ref_ranges[i])
                    i += 1
                    continue

                new_query_ranges.append((query_ranges[i][0], query_ranges[i][1] + shift))
                new_ref_ranges.append((ref_ranges[i][0], ref_ranges[i][1] + shift))

                next_index = i + 1
                while shift > 0 and next_index < len(query_ranges):
                    size = query_ranges[next_index][1] - query_ranges[next_index][0] + 1
                    if size > shift:
                        new_query_ranges.append((query_ranges[next_index][0] + shift, query_ranges[next_index][1]))
                        new_ref_ranges.append((ref_ranges[next_index][0] + shift, ref_ranges[next_index][1]))
                        shift = 0
                    else:
                        shift -= size
                    next_index += 1
                i = next_index
            query_ranges = new_query_ranges
            ref_ranges = new_ref_ranges
        seq = ''
        cigar = []
        for i in range(0, len(query_ranges)):
            rcurr = ref_ranges[i]
            qcurr = query_ranges[i]
            size = qcurr[1] - qcurr[0] + 1
            if i > 0:  # will be an ins/del depending on the distance from the last block
                # append based on the prev range
                rprev = ref_ranges[i - 1]
                qprev = query_ranges[i - 1]
                qjump = qcurr[0] - qprev[1]
                rjump = rcurr[0] - rprev[1]
                if rjump == 1:  # reference is consecutive
                    if qjump > 1:  # query range skipped. insertion to the reference sequence
                        cigar.append((CIGAR.I, qjump - 1))
                        # adds the inserted seq for the pysam read
                        seq += query_sequence[qprev[1] + 1:qcurr[0]]
                elif qjump == 1:  # query is consecutive
                    if rjump > 1:  # reference range skipped. deletion of the reference sequence
                        cigar.append((CIGAR.D, rjump - 1))
                else:  # indel
                    seq += query_sequence[qprev[1] + 1:qcurr[0]]
                    cigar.append((CIGAR.I, qjump - 1))
                    cigar.append((CIGAR.D, rjump - 1))
            # compute the match/mismatch for the current block
            if not reference_sequence:
                cigar.append((CIGAR.M, size))
            else:
                for ref_seq, query_seq in zip(reference_sequence[rcurr[0]:rcurr[1] + 1], query_sequence[qcurr[0]:qcurr[1] + 1]):
                    if DNA_ALPHABET.match(ref_seq, query_seq):
                        cigar.append((CIGAR.EQ, 1))
                    else:
                        cigar.append((CIGAR.X, 1))
            seq += query_sequence[qcurr[0]:qcurr[1] + 1]
        # add initial soft-clipping
        if query_ranges[0][0] > 0:  # first block starts after the query start
            temp = query_sequence[0:query_ranges[0][0]]
            seq = temp + seq
            cigar.insert(0, (CIGAR.S, len(temp)))
        if query_ranges[-1][1] < len(query_sequence) - 1:
            temp = query_sequence[query_ranges[-1][1] + 1:]
            seq += temp
            cigar.append((CIGAR.S, len(temp)))
        read = SamRead(reference_name=row['tname'], alignment_score=row['score'])
        read.query_sequence = seq
        read.reference_start = row['tstarts'][0]
        read.reference_id = chrom
        read.cigar = _cigar.join(cigar)
        read.query_name = row['qname']
        read.mapping_quality = NA_MAPPING_QUALITY
        if row['strand'] == STRAND.NEG:
            read.flag = read.flag | PYSAM_READ_FLAGS.REVERSE
            # read.cigar = read.cigar[::-1] # DON't REVERSE b/c blat reports on the positive strand already
        if read.query_sequence != row['qseq_full'] and read.query_sequence != reverse_complement(row['qseq_full']):
            raise AssertionError(
                'read sequence should reproduce input sequence',
                read.cigar,
                read.query_sequence,
                row['qseq_full'],
                reverse_complement(row['qseq_full'])
            )
        qcons = sum([v for c, v in read.cigar if c in QUERY_ALIGNED_STATES])
        assert len(read.query_sequence) == qcons
        try:
            query_coverage_interval(read)
        except (AttributeError, ValueError) as err:
            raise err
        return read


def process_blat_output(
        input_bam_cache,
        query_id_mapping,
        reference_genome,
        aligner_output_file='aligner_out.temp',
        blat_min_percent_of_max_score=0.8,
        blat_min_identity=0.7,
        blat_limit_top_aln=25,
        is_protein=False):
    """
    converts the blat output pslx (unheadered file) to bam reads
    """
    if is_protein:
        raise NotImplementedError('currently does not support aligning protein sequences')

    try:
        _, rows = Blat.read_pslx(aligner_output_file, query_id_mapping, is_protein=is_protein)
    except tab.tab.EmptyFileError:
        rows = []

    # split the rows by query id
    rows_by_query = {}
    for row in rows:
        rows_by_query.setdefault(row['qname'], []).append(row)

    reads_by_query = {}
    sequences = set(query_id_mapping.values())
    for seq in sequences:
        reads_by_query[seq] = []
    for query_id, rows in rows_by_query.items():
        query_seq = query_id_mapping[query_id]

        reads = []
        for row in rows:
            try:
                read = Blat.pslx_row_to_pysam(row, input_bam_cache, reference_genome)
            except KeyError as err:
                LOG('warning: reference template name not recognized', str(err), level=logging.DEBUG)
            except AssertionError as err:
                LOG('warning: invalid blat alignment', repr(err), level=logging.DEBUG)
            else:
                reads.append((row, read))

        filtered_rows = [(row, read) for row, read in reads if round(row['percent_ident'], 0) >= blat_min_identity]
        # filter on score
        filtered_rows.sort(key=lambda x: x[0]['score'], reverse=True)
        # filter on percent id
        score_ranks = {}
        for count, score in enumerate(sorted([row['score'] for row, read in reads], reverse=True)):
            score_ranks[score] = count
        min_rank = min(list(score_ranks.values()) + [0])

        filtered_reads = []
        for count, (row, read) in enumerate(sorted(reads, key=lambda x: x[0]['score'], reverse=True)):
            if count >= blat_limit_top_aln:
                break
            row['rank'] = score_ranks[row['score']]
            if row['rank'] > min_rank:
                read.mapping_quality = 0
            read.alignment_rank = row['rank']
            read.set_tag(PYSAM_READ_FLAGS.BLAT_SCORE, row['score'], value_type='i')
            read.set_tag(PYSAM_READ_FLAGS.BLAT_ALIGNMENTS, len(filtered_rows), value_type='i')
            read.set_tag(PYSAM_READ_FLAGS.BLAT_PMS, blat_min_percent_of_max_score, value_type='f')
            read.set_tag(PYSAM_READ_FLAGS.BLAT_RANK, row['rank'], value_type='i')
            read.set_tag(PYSAM_READ_FLAGS.BLAT_PERCENT_IDENTITY, row['percent_ident'], value_type='f')
            filtered_reads.append(read)
        reads_by_query[query_seq] = filtered_reads
    return reads_by_query
