import pysam
import itertools
import math
import subprocess
import warnings
import re
import os
import TSV
from structural_variant.constants import *
from structural_variant.read_tools import CigarTools
from Bio.Seq import Seq
import tempfile
from structural_variant.interval import Interval


class BlatAlignedSegment(pysam.AlignedSegment):
    """
    """
    def __init__(self, reference_name=None, blat_score=None):
        """
        Args:
            row (Dict[str,]): a row dictionary from the Blat.read_pslx method
        """
        pysam.AlignedSegment.__init__(self)
        if reference_name is None:
            self._reference_name = pysam.AlignedSegment.reference_name(self)
        else:
            self._reference_name = reference_name
        self.blat_score = blat_score

    def query_coverage_interval(self):
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
            row (dict[str,]): a row object from the 'read_pslx' method
            bam_cache (BamCache): the bam file/cache to use as a template for creating reference_id from chr name

        """
        chrom = bam_cache.reference_id(row['tname'])
        query_sequence = row['qseq_full']
        if row['strand'] == STRAND.NEG:
            temp = Seq(query_sequence, DNA_ALPHABET)
            temp = temp.reverse_complement()
            query_sequence = str(temp)

        # note: converting to inclusive range [] vs end-exclusive [)
        reference_sequence = reference_genome[row['tname']].seq if reference_genome else None
        query_ranges = [(x, x + y - 1) for x, y in zip(row['qstarts'], row['block_sizes'])]
        ref_ranges = [(x, x + y - 1) for x, y in zip(row['tstarts'], row['block_sizes'])]
        seq = ''
        cigar = []

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
        read.cigar = CigarTools.join(cigar)
        read.query_name = row['qname']
        read.mapping_quality = NA_MAPPING_QUALITY

        if row['strand'] == STRAND.NEG:
            read.flag = read.flag | PYSAM_READ_FLAGS.REVERSE
        if read.query_sequence != row['qseq_full'] and read.query_sequence != reverse_complement(row['qseq_full']):
            raise AssertionError(
                'read sequence should reproduce input sequence',
                read.cigar,
                read.query_sequence,
                row['qseq_full'],
                reverse_complement(row['qseq_full'])
            )
        return read


def blat_contigs(
        evidence,
        INPUT_BAM_CACHE,
        reference_genome,
        ref_2bit='/home/pubseq/genomes/Homo_sapiens/GRCh37/blat/hg19.2bit',
        min_percent_of_max_score=0.8,
        min_identity=0.95,
        is_protein=False,
        MIN_EXTEND_OVERLAP=10,
        **kwargs):
    """
    given a set of contigs, call blat from the commandline and adds the results to the contigs
    associated with each Evidence object

    Args:
        evidence (List[Evidence]): the iterable container of of evidence object which has associated contigs
        bam_cache (BamCache): the bam to use as a template in generating bam-like reads
        ref (str): path to the reference genome 2bit file for blat
        min_percent_of_max_score (float): filters all results with a score of a lower fraction of the best score
        min_identity (float): minimum percent identity
        is_protein (boolean): used in blat calculations
    """
    min_identity *= 100
    blat_options = kwargs.pop('blat_options',
                              ["-stepSize=5", "-repMatch=2253", "-minScore=0", "-minIdentity={0}".format(min_identity)])

    tempfiles = []
    try:
        # write the input sequences to a fasta file

        fasta = tempfile.NamedTemporaryFile(mode='w', delete=False)
        tempfiles.append(fasta.name)
        fasta_name = fasta.name
        query_id_mapping = {}
        count = 1
        sequences = set()
        for e in evidence:
            for c in e.contigs:
                sequences.add(c.seq)
        for seq in sequences:
            n = 'seq{0}'.format(count)
            query_id_mapping[n] = seq
            fasta.write('>' + n + '\n' + seq + '\n')
            count += 1
        fasta.close()

        # call the blat subprocess
        psl = tempfile.NamedTemporaryFile(delete=False)
        tempfiles.append(psl.name)
        # will raise subprocess.CalledProcessError if non-zero exit status
        # parameters from https://genome.ucsc.edu/FAQ/FAQblat.html#blat4
        # print(["blat", ref_2bit, fasta_name, psl.name, '-out=pslx', '-noHead'] + blat_options)
        subprocess.check_output(["blat", ref_2bit, fasta_name, psl.name, '-out=pslx', '-noHead'] + blat_options)
        psl.close()

        header, rows = Blat.read_pslx(psl.name, query_id_mapping, is_protein=is_protein)

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
            filtered_rows = [row for row in rows if round(row['percent_ident'], 0) >= min_identity]

            # filter on score
            scores = sorted([r['score'] for r in rows])
            max_score = scores[-1]
            second_score = scores[-2]
            min_score = max_score * min_percent_of_max_score
            filtered_rows = [row for row in filtered_rows if row['score'] >= min_score or row['score'] == second_score]

            filtered_rows.sort(key=lambda x: x['score'], reverse=True)
            reads = []
            for rank, row in enumerate(filtered_rows):
                try:
                    read = Blat.pslx_row_to_pysam(row, INPUT_BAM_CACHE, reference_genome)
                    read.set_tag('bs', row['score'], value_type='i')
                    read.set_tag('ba', len(filtered_rows), value_type='i')
                    read.set_tag('bp', min_percent_of_max_score, value_type='f')
                    read.set_tag('br', rank, value_type='i')
                    read.set_tag('bi', row['percent_ident'], value_type='f')
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

                if e.break1.chr == e.break2.chr and not e.opposing_strands:
                    for read in aln:
                        # if it covers both breakpoints add to putative alignments
                        temp = Interval(read.reference_start, read.reference_end - 1)
                        if INPUT_BAM_CACHE.chr(read) == e.break1.chr \
                                and Interval.overlaps(e.window1, temp) \
                                and Interval.overlaps(e.window2, temp):
                            # split the continuous alignment, assume ins/dup or indel
                            putative_alignments.append((read, None))

                for a1, a2 in itertools.combinations([x for x in aln if (x, None) not in putative_alignments], 2):
                    # do they overlap both breakpoints
                    if a1.reference_id > a2.reference_id or \
                            (a1.reference_id == a2.reference_id and a1.reference_start > a2.reference_start):
                        a1, a2 = (a2, a1)

                    if (a1.is_reverse != a2.is_reverse) != e.opposing_strands:
                        continue

                    union = Interval.union(a1.query_coverage_interval(),
                                           a2.query_coverage_interval())
                    if len(union) - len(a1.query_coverage_interval()) < MIN_EXTEND_OVERLAP \
                            or len(union) - len(a2.query_coverage_interval()) < MIN_EXTEND_OVERLAP:
                        continue

                    if INPUT_BAM_CACHE.chr(a1) == e.break1.chr \
                            and Interval.overlaps(e.window1, (a1.reference_start, a1.reference_end - 1)) \
                            and INPUT_BAM_CACHE.chr(a2) == e.break2.chr \
                            and Interval.overlaps(e.window2, (a2.reference_start, a2.reference_end - 1)):
                        putative_alignments.append((a1, a2))
                    elif INPUT_BAM_CACHE.chr(a2) == e.break1.chr \
                            and Interval.overlaps(e.window1, (a2.reference_start, a2.reference_end - 1)) \
                            and INPUT_BAM_CACHE.chr(a1) == e.break2.chr \
                            and Interval.overlaps(e.window2, (a1.reference_start, a1.reference_end - 1)):
                        putative_alignments.append((a2, a1))
                contig.alignments = putative_alignments
    finally:
        # clean up the temporary files
        for f in tempfiles:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError as e:
                    warnings.warn(repr(e))
