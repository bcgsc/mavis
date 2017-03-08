"""

"In general the coordinates in psl files are “zero based half open.” The first base in a sequence is numbered
zero rather than one. When representing a range the end coordinate is not included in the range. Thus the first
100 bases of a sequence are represented as 0-100, and the second 100 bases are represented as 100-200. There is
another little unusual feature in the .psl format. It has to do with how coordinates are handled on the
negative strand. In the qStart/qEnd fields the coordinates are where it matches from the point of view of the forward
strand (even when the match is on the reverse strand). However on the qStarts[] list, the coordinates are reversed."
--- http://wiki.bits.vib.be/index.php/Blat


"""
import pysam
import itertools
import math
import subprocess
import warnings
import re
import os
import TSV
from .constants import *
from .bam import cigar as cigar_tools
from .bam import read as read_tools
import tempfile
from .interval import Interval


class BlatAlignedSegment(pysam.AlignedSegment):
    """
    """
    def __init__(self, reference_name=None, blat_score=None):
        """
        Args:
            row (:class:`dict` by :class:`str`): a row dictionary from the Blat.read_pslx method
        """
        pysam.AlignedSegment.__init__(self)
        if reference_name is None:
            self._reference_name = pysam.AlignedSegment.reference_name(self)
        else:
            self._reference_name = reference_name
        self.blat_score = blat_score

    def query_coverage_interval(self):
        """
        Returns:
            :class:`~structural_variant.interval.Interval`: The portion of the original query sequence that is aligned by this read
        """
        seq = self.query_sequence
        s = 0
        t = len(seq) - 1
        if self.cigar[0][0] == CIGAR.S:
            s += self.cigar[0][1]
        if self.cigar[-1][0] == CIGAR.S:
            t -= self.cigar[-1][1]
        return Interval(s, t)

    @property
    def reference_name(self):
        return self._reference_name


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
    def read_pslx(filename, seqid_to_sequence_mapping, is_protein=False):
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

        split_csv_trailing_seq = lambda x: [s.upper() for s in re.sub(',$', '', x).split(',')]
        split_csv_trailing_ints = lambda x: [int(s) for s in re.sub(',$', '', x).split(',')]

        header, rows = TSV.read_file(
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
                'strand': '^[\+-]$'
            }
        )

        for row in rows:
            row['score'] = Blat.score(row, is_protein=is_protein)
            row['percent_ident'] = Blat.percent_identity(row, is_protein=is_protein)
            qseq = seqid_to_sequence_mapping[row['qname']]
            row['qseq_full'] = qseq
        return header, rows

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
        seq = ''
        cigar = []
        #print(query_ranges, ref_ranges)
        #print([query_sequence[q.start:q.end + 1] for q in query_ranges])
        #print('ref', [str(reference_sequence[r.start:r.end + 1]) for r in ref_ranges])
        #print('ref 0-20', reference_sequence[0:20])

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

                n = i + 1
                while shift > 0 and n < len(query_ranges):
                    size = query_ranges[n][1] - query_ranges[n][0] + 1
                    if size > shift:
                        new_query_ranges.append((query_ranges[n][0] + shift, query_ranges[n][1]))
                        new_ref_ranges.append((ref_ranges[n][0] + shift, ref_ranges[n][1]))
                        shift = 0
                    else:
                        shift -= size
                    n += 1
                i = n
            query_ranges = new_query_ranges
            ref_ranges = new_ref_ranges
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
                for r, q in zip(reference_sequence[rcurr[0]:rcurr[1] + 1], query_sequence[qcurr[0]:qcurr[1] + 1]):
                    if DNA_ALPHABET.match(r, q):
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
        read = BlatAlignedSegment(row['tname'], row['score'])
        read.query_sequence = seq
        read.reference_start = row['tstarts'][0]
        read.reference_id = chrom
        read.cigar = cigar_tools.join(cigar)
        read.query_name = row['qname']
        read.mapping_quality = NA_MAPPING_QUALITY
        if row['strand'] == STRAND.NEG:
            read.flag = read.flag | PYSAM_READ_FLAGS.REVERSE
            #read.cigar = read.cigar[::-1] # DON't REVERSE b/c blat reports on the positive strand already
        if read.query_sequence != row['qseq_full'] and read.query_sequence != reverse_complement(row['qseq_full']):
            raise AssertionError(
                'read sequence should reproduce input sequence',
                read.cigar,
                read.query_sequence,
                row['qseq_full'],
                reverse_complement(row['qseq_full'])
            )
        return read


def paired_alignment_score(read1, read2=None):
    score = read_tools.calculate_alignment_score(read1) * read1.query_coverage_interval().length()
    if read2 is not None:
        qci1 = read1.query_coverage_interval()
        qci2 = read2.query_coverage_interval()
        total_len = qci1.length() + qci2.length()
        s1 = read_tools.calculate_alignment_score(read1)
        s2 = read_tools.calculate_alignment_score(read2)
        avg = s1 * qci1.length() / total_len + s2 * qci2.length() / total_len
        score = avg * len(qci1 | qci2)
    return score


def blat_contigs(
        evidence,
        INPUT_BAM_CACHE,
        reference_genome,
        blat_pslx_output_file='blat_out.pslx',
        blat_fa_input_file='blat_in.fa',
        blat_2bit_reference='/home/pubseq/genomes/Homo_sapiens/GRCh37/blat/hg19.2bit',
        blat_min_percent_of_max_score=0.8,
        blat_min_identity=0.7,
        blat_min_query_consumption=0.5,
        is_protein=False,
        min_extend_overlap=10,
        pair_scoring_function=paired_alignment_score,
        clean_files=True,
        log=lambda *pos, **kwargs: None,
        **kwargs):
    """
    given a set of contigs, call blat from the command line and adds the results to the contigs
    associated with each Evidence object

    Args:
        evidence (list of Evidence): the iterable container of of Evidence object which has associated contigs
        INPUT_BAM_CACHE (BamCache): the bam to use as a template in generating bam-like reads
        reference_genome (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence by template/chr name
        blat_2bit_reference (str): path to the 2bit file for blat
        blat_min_percent_of_max_score (float): ignores all alignments with a score less
        blat_min_identity (float): minimum percent identity
        is_protein (bool): is the sequence an amino acid sequence (used in the blat calculations)
        min_extend_overlap (int): minimum amount of non-shared coverage of the template sequence required to pair alignments
        =blat_options (list of str): optional, can specify alternate blat parameters to use

    .. todo::
        add support for blatting protein sequences
        allow one alignment to be lower score if its partner is higher score? or recompute a combined score?
    """
    if is_protein:
        raise NotImplementedError('currently does not support blatting protein sequences')
    blat_min_identity *= 100
    blat_options = kwargs.pop('blat_options',
                              ["-stepSize=5", "-repMatch=2253", "-minScore=0", "-minIdentity={0}".format(blat_min_identity)])

    try:
        # write the input sequences to a fasta file
        query_id_mapping = {}
        count = 1
        sequences = set()
        ev_by_seq = {}
        for e in evidence:
            for c in e.contigs:
                sequences.add(c.seq)
                ev_by_seq.setdefault(c.seq, []).append(e.data.get(COLUMNS.cluster_id, None))
        with open(blat_fa_input_file, 'w') as fh:
            for seq in sequences:
                n = 'seq{}_{}'.format(count, '_'.join(sorted([x for x in ev_by_seq[seq] if x is not None])))
                query_id_mapping[n] = seq
                fh.write('>' + n + '\n' + seq + '\n')
                count += 1
        if len(sequences) == 0:
            return
        log('will blat', len(sequences), 'unique sequences', time_stamp=False)
        # call the blat subprocess
        # will raise subprocess.CalledProcessError if non-zero exit status
        # parameters from https://genome.ucsc.edu/FAQ/FAQblat.html#blat4
        # print(["blat", blat_2bit_reference, fasta_name, psl.name, '-out=pslx', '-noHead'] + blat_options)
        subprocess.check_output(["blat", blat_2bit_reference,
            blat_fa_input_file, blat_pslx_output_file, '-out=pslx', '-noHead'] + blat_options)

        header, rows = Blat.read_pslx(blat_pslx_output_file, query_id_mapping, is_protein=is_protein)

        # split the rows by query id
        rows_by_query = {}
        for row in rows:
            if row['qname'] not in rows_by_query:
                rows_by_query[row['qname']] = []
            rows_by_query[row['qname']].append(row)

        reads_by_query = {}
        for s in sequences:
            reads_by_query[s] = []
        for query_id, rows in rows_by_query.items():
            query_seq = query_id_mapping[query_id]
            # filter on percent id
            filtered_rows = [row for row in rows if round(row['percent_ident'], 0) >= blat_min_identity]

            # filter on score
            scores = sorted([r['score'] for r in rows])
            max_score = scores[-1]
            filtered_rows.sort(key=lambda x: x['score'], reverse=True)
            reads = []
            for rank, row in enumerate(filtered_rows):
                try:
                    read = Blat.pslx_row_to_pysam(row, INPUT_BAM_CACHE, reference_genome)
                    read.set_tag(PYSAM_READ_FLAGS.BLAT_SCORE, row['score'], value_type='i')
                    read.set_tag(PYSAM_READ_FLAGS.BLAT_ALIGNMENTS, len(filtered_rows), value_type='i')
                    read.set_tag(PYSAM_READ_FLAGS.BLAT_PMS, blat_min_percent_of_max_score, value_type='f')
                    read.set_tag(PYSAM_READ_FLAGS.BLAT_RANK, rank, value_type='i')
                    read.set_tag(PYSAM_READ_FLAGS.BLAT_PERCENT_IDENTITY, row['percent_ident'], value_type='f')
                    reads.append(read)
                except KeyError as e:
                    warnings.warn(
                        'warning: reference template name not recognized {0}'.format(e))
            reads_by_query[query_seq] = reads

        # now for each evidence assign an alignment to each contig
        for e in evidence:
            for contig in e.contigs:
                aln = reads_by_query.get(contig.seq, [])
                putative_alignments = []
                combo_prohibited = set()
                putative_event_types = e.putative_event_types()

                if e.break1.chr == e.break2.chr and not e.opposing_strands:
                    for read in aln:
                        # if it covers both breakpoints add to putative alignments
                        temp = Interval(read.reference_start, read.reference_end - 1)
                        if INPUT_BAM_CACHE.chr(read) == e.break1.chr \
                                and Interval.overlaps(e.outer_window1, temp) \
                                and Interval.overlaps(e.outer_window2, temp):
                            # split the continuous alignment, assume ins/dup or indel
                            ins = sum([v for c, v in read.cigar if c == CIGAR.I] + [0])
                            dln = sum([v for c, v in read.cigar if c in [CIGAR.D, CIGAR.N]] + [0])
                            consume = len(read.query_coverage_interval()) / len(read.query_sequence)
                            if consume < blat_min_query_consumption:
                                continue
                            for event_type in e.putative_event_types():
                                if event_type in [SVTYPE.INS, SVTYPE.DUP] and ins > 0 and ins > dln:
                                    putative_alignments.append((read, None))
                                elif event_type == SVTYPE.DEL and dln > 0 and dln > ins:
                                    putative_alignments.append((read, None))

                combo_prohibited = [x for x, y in putative_alignments]
                for a1, a2 in itertools.combinations([x for x in aln if x not in combo_prohibited], 2):
                    # do they overlap both breakpoints
                    if a1.reference_id > a2.reference_id or \
                            (a1.reference_id == a2.reference_id and a1.reference_start > a2.reference_start):
                        a1, a2 = (a2, a1)

                    q1 = a1.query_coverage_interval()
                    q2 = a2.query_coverage_interval()

                    if a1.is_reverse != a1.is_reverse:
                        l = len(a1.query_sequence) - 1
                        q2 = Interval(l - q2.end, l - q2.start)

                    if not e.interchromosomal and a1.reference_end > a2.reference_end:
                        continue

                    if a2.query_sequence != a1.query_sequence:  # alignments aligned to opposite strands
                        if not e.opposing_strands:
                            continue
                        q2 = Interval(len(a2.query_sequence) - q2.end, len(a2.query_sequence) - q2.start)
                    elif e.opposing_strands:
                        continue
                    union = q1 | q2
                    consume = len(union)
                    if not Interval.overlaps(q1, q2):
                        consume = len(q1) + len(q2)
                    if consume / len(a1.query_sequence) < blat_min_query_consumption:
                        continue

                    if len(union) - len(q1) < min_extend_overlap or len(union) - len(q2) < min_extend_overlap:
                        continue
                    if INPUT_BAM_CACHE.chr(a1) == e.break1.chr \
                            and Interval.overlaps(e.outer_window1, (a1.reference_start, a1.reference_end - 1)) \
                            and INPUT_BAM_CACHE.chr(a2) == e.break2.chr \
                            and Interval.overlaps(e.outer_window2, (a2.reference_start, a2.reference_end - 1)):
                        putative_alignments.append((a1, a2))
                    elif INPUT_BAM_CACHE.chr(a2) == e.break1.chr \
                            and Interval.overlaps(e.outer_window1, (a2.reference_start, a2.reference_end - 1)) \
                            and INPUT_BAM_CACHE.chr(a1) == e.break2.chr \
                            and Interval.overlaps(e.outer_window2, (a1.reference_start, a1.reference_end - 1)):
                        putative_alignments.append((a2, a1))
                if len(putative_alignments) == 0:
                    continue
                score_by_alignments = {}
                for read1, read2 in putative_alignments:
                    score = pair_scoring_function(read1, read2)
                    score_by_alignments[(read1, read2)] = score
                max_score = max(list(score_by_alignments.values()))
                for pair in putative_alignments:
                    score = score_by_alignments[pair]
                    if score / max_score >= blat_min_percent_of_max_score:
                        contig.alignments.append(pair)
    finally:
        # clean up the temporary files
        if clean_files:
            for f in [blat_pslx_output_file, blat_fa_input_file]:
                if os.path.exists(f):
                    try:
                        os.remove(f)
                    except OSError as e:
                        warnings.warn(repr(e))
