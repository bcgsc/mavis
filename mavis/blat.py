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
from .bam.cigar import QUERY_ALIGNED_STATES
from .bam import read as read_tools
from .interval import Interval
from .error import InvalidRearrangement
from .breakpoint import BreakpointPair
from .util import devnull


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
    
    def __repr__(self):
        return '{}({}:{}, {}, {})'.format(
            self.__class__.__name__, self.reference_name, self.reference_start,
            self.cigar, self.query_sequence
        )

    def __copy__(self):
        cp = BlatAlignedSegment(self.reference_name, self.blat_score)
        cp.query_sequence = self.query_sequence
        cp.reference_start = self.reference_start
        cp.reference_id = self.reference_id
        cp.cigar = self.cigar
        cp.query_name = self.query_name
        cp.mapping_quality = self.mapping_quality
        cp.set_tags(self.get_tags())
        cp.flag = self.flag
        cp.next_reference_id = self.next_reference_id
        cp.next_reference_start = self.next_reference_start
        return cp
    
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
                if verbose:
                    warnings.warn(repr(err))
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
        assert(len(read.query_sequence) == qcons)
        try:
            read.query_coverage_interval()
        except (AttributeError, ValueError) as err:
            print(row)
            raise err
        return read


def get_blat_version():
    proc = subprocess.getoutput(['blat'])
    for line in proc.split('\n'):
        m = re.search('blat - Standalone BLAT v. (\d+(x\d+)?)', line)
        if m:
            return m.group(1)
    raise ValueError("unable to parse blat version number from:'{}'".format(proc))


def paired_alignment_score(read1, read2=None):
    score = read_tools.calculate_alignment_score(read1) * BlatAlignedSegment.query_coverage_interval(read1).length()
    if read2 is not None:
        qci1 = BlatAlignedSegment.query_coverage_interval(read1)
        qci2 = BlatAlignedSegment.query_coverage_interval(read2)
        total_len = qci1.length() + qci2.length()
        s1 = read_tools.calculate_alignment_score(read1)
        s2 = read_tools.calculate_alignment_score(read2)
        avg = s1 * qci1.length() / total_len + s2 * qci2.length() / total_len
        score = avg * len(qci1 | qci2)
    return score


def select_paired_alignments(
    bpp, aligned_contigs, min_query_consumption, min_extend_overlap, max_event_size, min_anchor_size
):
    """
    give a breakpoint pair and a set of alignments for contigs associated with the given pair,
    alignments are paired (some events cannot be represented with a single bamfile alignment)
    and the most appropriate alignments supporting the breakpoint pair are selected and returned
    """
    # now for each bpp assign an alignment to each contig
    putative_alignments = []
    putative_event_types = set(bpp.putative_event_types())
    if {SVTYPE.INS, SVTYPE.DUP} & putative_event_types:
        putative_event_types = putative_event_types | {SVTYPE.INS, SVTYPE.DUP}
    
    # for events on the same template and strand we expect to find a single contig alignment
    if not bpp.interchromosomal and not bpp.opposing_strands:
        for read in aligned_contigs:
            # if it covers both breakpoints add to putative alignments
            read_cover = Interval(read.reference_start, read.reference_end - 1)
            if all([
                read.reference_name == bpp.break1.chr,
                Interval.overlaps(bpp.outer_window1, read_cover),
                Interval.overlaps(bpp.outer_window2, read_cover)
            ]):
                # split the continuous alignment, assume ins/dup or indel
                ins = sum([v for c, v in read.cigar if c == CIGAR.I] + [0])
                dln = sum([v for c, v in read.cigar if c in [CIGAR.D, CIGAR.N]] + [0])
                consume = len(BlatAlignedSegment.query_coverage_interval(read)) / len(read.query_sequence)
                if consume < min_query_consumption:
                    continue
                
                for event_type in putative_event_types:
                    if event_type in {SVTYPE.INS, SVTYPE.DUP} and ins > 0 and ins > dln:
                        putative_alignments.append((read, None))
                    elif event_type in {SVTYPE.DEL, SVTYPE.INV} and dln > 0 and dln > ins:
                        putative_alignments.append((read, None))
    
    # don't use reads in combined alignments if they have already been assigned in a single alignment
    combo_prohibited = [x for x, y in putative_alignments]
    
    for read1, read2 in itertools.combinations([x for x in aligned_contigs if x not in combo_prohibited], 2):
        # do they overlap both breakpoints
        if any([
            read1.reference_name > read2.reference_name,
            read1.reference_name == read2.reference_name and read1.reference_start > read2.reference_start
        ]):
            read1, read2 = read2, read1

        if read1.reference_name != bpp.break1.chr or read2.reference_name != bpp.break2.chr:
            continue
        read1 = read_tools.convert_events_to_softclipping(
            read1, bpp.break1.orient, max_event_size=max_event_size, min_anchor_size=min_anchor_size)
        read2 = read_tools.convert_events_to_softclipping(
            read2, bpp.break2.orient, max_event_size=max_event_size, min_anchor_size=min_anchor_size)
        # check that the coverage intervals overlap
        if any([
            not Interval.overlaps((read1.reference_start + 1, read1.reference_end), bpp.outer_window1),
            not Interval.overlaps((read2.reference_start + 1, read2.reference_end), bpp.outer_window2)
        ]):
            continue
        # reads should have unique reference overlap
        if not bpp.interchromosomal and read1.reference_end > read2.reference_end:
            continue
        # check that the combination extends the amount of the initial query sequence we consume
        query_cover1 = BlatAlignedSegment.query_coverage_interval(read1)
        query_cover2 = BlatAlignedSegment.query_coverage_interval(read2)

        if read2.query_sequence != read1.query_sequence:  # alignments aligned to opposite strands
            if not bpp.opposing_strands:
                continue
            if read2.query_sequence != reverse_complement(read1.query_sequence):
                continue
            l = len(read1.query_sequence) - 1
            query_cover2 = Interval(l - query_cover2.end, l - query_cover2.start)
        elif bpp.opposing_strands:
            continue
        
        consume = len(query_cover1 | query_cover2)
        if not Interval.overlaps(query_cover1, query_cover2):
            consume = len(query_cover1) + len(query_cover2)
        
        if consume / len(read1.query_sequence) < min_query_consumption:
            continue
        if consume - len(query_cover1) < min_extend_overlap or consume - len(query_cover2) < min_extend_overlap:
            continue
        
        try:
            call = BreakpointPair.call_breakpoint_pair(read1, read2)
            if not set(BreakpointPair.classify(call)) & putative_event_types:
                continue
        except (InvalidRearrangement, AssertionError) as err:
            continue
        putative_alignments.append((read1, read2))
    return putative_alignments


def blat_contigs(
        evidence,
        INPUT_BAM_CACHE,
        reference_genome,
        blat_pslx_output_file='blat_out.pslx',
        blat_fa_input_file='blat_in.fa',
        blat_2bit_reference='/home/pubseq/genomes/Homo_sapiens/GRCh37/blat/hg19.2bit',
        blat_min_percent_of_max_score=0.8,
        blat_min_identity=0.7,
        contig_aln_min_query_consumption=0.5,
        contig_aln_max_event_size=50,
        contig_aln_min_anchor_size=50,
        is_protein=False,
        min_extend_overlap=10,
        pair_scoring_function=paired_alignment_score,
        clean_files=True,
        log=devnull,
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
    blat_options = kwargs.pop(
        'blat_options', ["-stepSize=5", "-repMatch=2253", "-minScore=0", "-minIdentity={0}".format(blat_min_identity)])

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
                n = 'seq{}'.format(count)
                log(n, [x for x in ev_by_seq[seq] if x is not None])
                query_id_mapping[n] = seq
                fh.write('>' + n + '\n' + seq + '\n')
                count += 1
        if len(sequences) == 0:
            return
        log('will blat', len(sequences), 'unique sequences', time_stamp=False)
        # call the blat subprocess
        # will raise subprocess.CalledProcessError if non-zero exit status
        # parameters from https://genome.ucsc.edu/FAQ/FAQblat.html#blat4
        log(['blat', blat_2bit_reference,
            blat_fa_input_file, blat_pslx_output_file, '-out=pslx', '-noHead'] + blat_options)
        subprocess.check_output([
            'blat', blat_2bit_reference,
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
                putative_alignments = select_paired_alignments(
                    e, aln,
                    min_extend_overlap=min_extend_overlap,
                    min_query_consumption=contig_aln_min_query_consumption,
                    min_anchor_size=contig_aln_min_anchor_size,
                    max_event_size=contig_aln_max_event_size
                )
                if len(putative_alignments) == 0:
                    continue
                score_by_alignments = {}
                for read1, read2 in sorted(putative_alignments, key=lambda x: pair_scoring_function(x[0], x[1])):
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
