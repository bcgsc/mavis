import pysam
import math
import subprocess
import warnings
import re
import os
import TSV
from structural_variant.constants import *
from structural_variant.align import reverse_complement
from Bio.Seq import Seq
import tempfile
from structural_variant.interval import Interval
from structural_variant.align import CigarTools


class BlatAlignedSegment(pysam.AlignedSegment):

    def __init__(self, row):
        pysam.AlignedSegment.__init__(self)
        self.blat = row

    def blat_score(self):  # convenience
        return self.blat['score']

    def query_coverage_interval(self):
        query_ranges = [(x, x + y - 1) for x, y in zip(self.blat['qstarts'], self.blat['block_sizes'])]
        u = Interval.union(*query_ranges)
        if not self.is_reverse:
            return u
        else:
            l = len(self.blat['qseq_full']) - 1
            return Interval(l - u[1], l - u[0])


class Blat:

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
        
        Perl:
            my $sizeMul = pslIsProtein($blockCount, $strand, $tStart, $tEnd, $tSize, $tStarts, $blockSizes);
            sizmul = 1 for DNA
            my $pslScore = $sizeMul * ($matches + ($repMatches >> 1) ) - $sizeMul * $misMatches - $qNumInsert - $tNumInsert)
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
                'match': 'int',
                'mismatch': 'int',
                'repmatch': 'int',
                'ncount': 'int',
                'qgap_count': 'int',
                'qgap_bases': 'int',
                'tgap_count': 'int',
                'tgap_bases': 'int',
                'qsize': 'int',
                'qstart': 'int',
                'qend': 'int',
                'tsize': 'int',
                'tstart': 'int',
                'tend': 'int',
                'block_count': 'int'
            },
            validate={
                'strand': '^[\+-]$'
            },
            transform={
                'tname': lambda x: re.sub('^chr', '', x),
                'block_sizes': split_csv_trailing_ints,
                'qstarts': split_csv_trailing_ints,
                'tstarts': split_csv_trailing_ints,
                'qseqs': split_csv_trailing_seq,
                'tseqs': split_csv_trailing_seq
            }
        )

        for row in rows:
            row['score'] = Blat.score(row, is_protein=is_protein)
            row['percent_ident'] = Blat.percent_identity(row, is_protein=is_protein)
            qseq = seqid_to_sequence_mapping[row['qname']]
            row['qseq_full'] = qseq
        return header, rows

    @staticmethod
    def pslx_row_to_pysam(row, bam_cache):
        """
        given a 'row' from reading a pslx file. converts the row to a BlatAlignedSegment object
        
        Args:
            row (dict[str,]): a row object from the 'read_pslx' method
            bam_cache (BamCache): the bam file/cache to use as a template for creating reference_id from chr name

        """
        chrom = bam_cache.reference_id(row['tname'])
        qseq = row['qseq_full']
        if row['strand'] == STRAND.NEG:
            temp = Seq(qseq, DNA_ALPHABET)
            temp = temp.reverse_complement()
            qseq = str(temp)

        # note: converting to inclusive range [] vs end-exclusive [)
        query_ranges = [(x, x + y - 1) for x, y in zip(row['qstarts'], row['block_sizes'])]
        ref_ranges = [(x, x + y - 1) for x, y in zip(row['tstarts'], row['block_sizes'])]

        cigar = []
        seq = ''

        # add initial soft-clipping
        if query_ranges[0][0] > 0:  # first block starts after the query start
            temp = qseq[0:query_ranges[0][0]]
            seq += temp
            cigar.append((CIGAR.S, len(temp)))
        for i in range(0, len(query_ranges)):
            rcurr = ref_ranges[i]
            qcurr = query_ranges[i]
            if i > 0:
                # append based on the prev range
                rprev = ref_ranges[i - 1]
                qprev = query_ranges[i - 1]
                qjump = qcurr[0] - qprev[1]
                rjump = rcurr[0] - rprev[1]

                if rjump == 1:  # reference is consecutive
                    if qjump > 1:  # query range skipped. insertion to the reference sequence
                        cigar.append((CIGAR.I, qjump - 1))
                        # adds the inserted seq for the pysam read
                        seq += qseq[qprev[1] + 1:qcurr[0]]
                elif qjump == 1:  # query is consecutive
                    if rjump > 1:  # reference range skipped. deletion of the reference sequence
                        cigar.append((CIGAR.D, rjump - 1))
                else:  # indel
                    seq += qseq[qprev[1] + 1:qcurr[0]]
                    cigar.append((CIGAR.I, qjump - 1))
                    cigar.append((CIGAR.D, rjump - 1))
            # add the current range of matches
            temp, offset = CigarTools.compute(row['tseqs'][i], row['qseqs'][i])
            if temp[0][0] == CIGAR.S:
                temp[0] = (CIGAR.X, temp[0][1])
            if temp[-1][0] == CIGAR.S:
                temp[-1] = (CIGAR.X, temp[-1][1])
            cigar = CigarTools.join(cigar, temp)
            seq += qseq[qcurr[0]:qcurr[1] + 1]

        if query_ranges[-1][1] < len(qseq) - 1:
            temp = qseq[query_ranges[-1][1] + 1:]
            seq += temp
            cigar.append((CIGAR.S, len(temp)))
        read = BlatAlignedSegment(row)
        read.query_sequence = seq
        read.reference_start = row['tstart']
        read.reference_id = chrom
        read.cigar = cigar
        read.query_name = row['qname']
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
        sequences,
        bam_cache,
        ref='/home/pubseq/genomes/Homo_sapiens/GRCh37/blat/hg19.2bit',
        min_percent_of_max_score=0.8,
        min_identity=0.95,
        is_protein=False,
        **kwargs):
    """
    given a set of contigs, call blat from the commandline and return the results
    """
    min_identity *= 100
    blat_options = kwargs.pop('blat_options',
                              ["-stepSize=5", "-repMatch=2253", "-minScore=0", "-minIdentity={0}".format(min_identity)])
    
    tempfiles = []
    try:
        # write the input sequences to a fasta file
        fasta = tempfile.NamedTemporaryFile(mode='w', delete=False)
        fasta_name = fasta.name
        query_id_mapping = {}
        count = 1
        for seq in sequences:
            n = 'seq{0}'.format(count)
            query_id_mapping[n] = seq
            fasta.write('>' + n + '\n' + seq + '\n')
            count += 1
        fasta.close()

        # call the blat subprocess
        psl = tempfile.NamedTemporaryFile(delete=False)
        # will raise subprocess.CalledProcessError if non-zero exit status
        # parameters from https://genome.ucsc.edu/FAQ/FAQblat.html#blat4
        print(["blat", ref, fasta_name, psl.name, '-out=pslx', '-noHead'] + blat_options)
        subprocess.check_output(["blat", ref, fasta_name, psl.name, '-out=pslx', '-noHead'] + blat_options)
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
            filtered_rows = []
            max_score = max([r['score'] for r in rows])

            # filter by alignment quality
            for row in rows:
                if row['score'] >= max_score * min_percent_of_max_score \
                        and round(row['percent_ident'], 0) >= min_identity:
                    filtered_rows.append(row)

            filtered_rows.sort(key=lambda x: x['score'], reverse=True)
            reads = []
            for rank, row in enumerate(filtered_rows):
                try:
                    read = Blat.pslx_row_to_pysam(row, bam_cache)
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
        return reads_by_query
    finally:
        for f in tempfiles:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError as e:
                    warnings.warn(repr(e))
