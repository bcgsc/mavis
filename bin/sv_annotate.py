
"""
About
------

This is the third step in the svmerge pipeline. It is responsible for annotating breakpoint pairs with reference
formation and drawing visualizations. Outputs are written to the annotation subfolder in the following pattern

::

    <output_dir_name>/
    |-- clustering/
    |-- validation/
    |-- annotation/
    |   `--<library>_<protocol>/
    |       |-- annotations.tab
    |       |-- annotations.fusions-cdna.fa
    |       `-- drawings/
    |           |-- <annotation_id>.legend.json
    |           `-- <annotation_id>.svg
    |-- pairing/
    `-- summary/

General Process
----------------

- Breakpoint pairs are first annotated by what transcripts are at each breakpoint (or lack thereof). All combinations
  are kept going forward.
- The related gene annotations are collected.
- The putative protein products are predicted.
- A Fusion transcript is built with different splicing possibilities according to the splicing model.
- For each splicing model a splice transcript is built.
- ORFs are computed for the spliced transcript and translated to create the putative AA sequence
- From the original transcript(s). The amino acid sequences of the domains is gathered and aligned to the new AA
  sequence
- Each new 'protein' product is drawn and those without products are drawn without a fusion track

.. todo::

    allow multiple duplicates between input files to be filtered before annotating

"""
import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from structural_variant.breakpoint import read_bpp_from_input_file, BreakpointPair
from structural_variant.annotate import load_reference_genes, load_reference_genome, load_templates
from structural_variant.annotate.variant import gather_annotations, FusionTranscript, determine_prime
from structural_variant.error import DrawingFitError, NotSpecifiedError
from structural_variant import __version__
from structural_variant.illustrate.settings import DiagramSettings
from structural_variant.illustrate.draw import draw_sv_summary_diagram
import TSV
from structural_variant.constants import PROTOCOL, SVTYPE, COLUMNS, sort_columns, PRIME, build_batch_id, log
import re
import json
from datetime import datetime
import warnings
import glob


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='Outputs the version number'
    )
    parser.add_argument(
        '-f', '--force_overwrite', default=False, type=bool,
        help='set flag to overwrite existing reviewed files')
    parser.add_argument(
        '--no_draw', default=True, action='store_false',
        help='set flag to suppress svg drawings of putative annotations')
    parser.add_argument(
        '-o', '--output',
        help='path to the output directory', required=True
    )
    parser.add_argument(
        '-n', '--input',
        help='path to the input file(s)', required=True, nargs='*'
    )
    parser.add_argument(
        '--low_memory', default=False, type=bool,
        help='when working on a machine with less memory this is sacrifice time for memory where possible')
    g = parser.add_argument_group('reference files')
    g.add_argument(
        '-a', '--annotations',
        default='/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv',
        help='path to the reference annotations of genes, transcript, exons, domains, etc.'
    )
    g.add_argument(
        '-r', '--reference_genome',
        default='/home/pubseq/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa',
        help='path to the human reference genome in fa format'
    )
    g.add_argument(
        '--template_metadata', default=os.path.join(os.path.dirname(__file__), 'cytoBand.txt'),
        help='file containing the cytoband template information')
    parser.add_argument(
        '-p', '--max_proximity', default=5000, type=int,
        help='The maximum distance away from breakpoints to look for proximal genes')
    parser.add_argument(
        '--min_orf_size', default=120, type=int, help='minimum size for putative ORFs')
    parser.add_argument(
        '--max_orf_cap', default=3, type=int, help='keep the n longest orfs')
    parser.add_argument(
        '--min_domain_mapping_match', default=0.8, type=float,
        help='minimum percent match for the domain to be considered aligned')
    g = parser.add_argument_group('visualization options')
    g.add_argument(
        '-d', '--domain_regex_filter', default='^PF\d+$',
        help='only show domains which names (external identifiers) match the given pattern')
    args = parser.parse_args()
    return args

class Index:
    def __init__(self, start=1, prefix=None):
        self.curr = start
        self.prefix = prefix

    def __call__(self):
        temp = self.curr
        self.curr += 1
        if self.prefix:
            return '{}{}'.format(self.prefix, temp)
        return temp


def main():
    # load the file
    args = parse_arguments()

    DRAWINGS_DIRECTORY = os.path.join(args.output, 'drawings')
    TABBED_OUTPUT_FILE = os.path.join(args.output, 'annotations.tab')
    FA_OUTPUT_FILE = os.path.join(args.output, 'annotations.fusion-cdna.fa')

    glob_checks = [
        os.path.join(args.output, '*.fa'),
        os.path.join(args.output, '*.tab'),
        os.path.join(DRAWINGS_DIRECTORY, '*.svg'),
        os.path.join(DRAWINGS_DIRECTORY, '*.json')
    ]
    for g in glob_checks:
        if len(glob.glob(g)) > 0:
            warnings.warn('existing files will be overwritten and directories will not be cleaned')
            break

    if not os.path.exists(DRAWINGS_DIRECTORY):
        os.mkdir(DRAWINGS_DIRECTORY)

    log('input arguments listed below')
    for arg, val in sorted(args.__dict__.items()):
        log(arg, '=', val, time_stamp=False)

    # test that the sequence makes sense for a random transcript
    bpps = []
    for f in args.input:
        log('loading:', f)
        bpps.extend(
            read_bpp_from_input_file(
                f,
                cast={
                    COLUMNS.stranded: TSV.tsv_boolean
                },
                _in={
                    COLUMNS.protocol: PROTOCOL,
                    COLUMNS.event_type: SVTYPE
                },
                add={
                    COLUMNS.cluster_id: Index(1),
                    COLUMNS.validation_id: Index(1)
                },
                simplify=False
            ))
    log('read {} breakpoint pairs'.format(len(bpps)))

    log('loading:' if not args.low_memory else 'indexing:', args.reference_genome)
    REFERENCE_GENOME = load_reference_genome(args.reference_genome, low_mem=args.low_memory)

    log('loading:', args.template_metadata)
    TEMPLATES = load_templates(args.template_metadata)

    log('loading:', args.annotations)
    REFERENCE_ANNOTATIONS = load_reference_genes(args.annotations, REFERENCE_GENOME=REFERENCE_GENOME)

    annotations = []
    for bpp in bpps:
        log('gathering annotations for', bpp)
        try:
            ann = gather_annotations(
                REFERENCE_ANNOTATIONS,
                bpp,
                event_type=bpp.data[COLUMNS.event_type],
                proximity=args.max_proximity
            )
            annotations.extend(ann)
        except KeyError:
            pass
        log('generated', len(ann), 'annotations', time_stamp=False)

    for bpp in annotations:
        if bpp.data[COLUMNS.event_type] not in BreakpointPair.classify(bpp):
            raise AssertionError(
                'input type does not fit with breakpoint pair description:',
                bpp.data[COLUMNS.event_type], BreakpointPair.classify(bpp))

    id_prefix = 'annotation_{}-'.format(build_batch_id())
    rows = []  # hold the row information for the final tsv file
    fa_sequences = {}
    for i, ann in enumerate(annotations):
        annotation_id = id_prefix + str(i + 1)
        ann.data[COLUMNS.annotation_id] = annotation_id
        row = ann.flatten()
        row[COLUMNS.break1_strand] = ann.transcript1.get_strand()
        row[COLUMNS.break2_strand] = ann.transcript2.get_strand()
        row[COLUMNS.fusion_sequence_fasta_file] = FA_OUTPUT_FILE

        log('current annotation', annotation_id, ann.transcript1, ann.transcript2, ann.event_type)

        # try building the fusion product
        ann_rows = []
        ft = None
        try:
            ft = FusionTranscript.build(
                ann, REFERENCE_GENOME,
                min_orf_size=args.min_orf_size,
                max_orf_cap=args.max_orf_cap,
                min_domain_mapping_match=args.min_domain_mapping_match
            )
            # add fusion information to the current row
            for t in ft.transcripts:
                fusion_fa_id = '{}_{}'.format(annotation_id, t.splicing_pattern.splice_type)
                if fusion_fa_id in fa_sequences:
                    raise AssertionError('should not be duplicate fa sequence ids', fusion_fa_id)
                fa_sequences[fusion_fa_id] = ft.get_cdna_seq(t.splicing_pattern)

                # duplicate the row for each translation
                for tl in t.translations:
                    nrow = dict()
                    nrow.update(row)
                    nrow[COLUMNS.fusion_splicing_pattern] = tl.transcript.splicing_pattern.splice_type
                    nrow[COLUMNS.fusion_cdna_coding_start] = tl.start
                    nrow[COLUMNS.fusion_cdna_coding_end] = tl.end
                    nrow[COLUMNS.fusion_sequence_fasta_id] = fusion_fa_id

                    domains = []
                    for dom in tl.domains:
                        m, t = dom.score_region_mapping()
                        temp = {
                            "name": dom.name,
                            "sequences": dom.get_seqs(),
                            "regions": [
                                {"start": dr.start, "end": dr.end} for dr in sorted(dom.regions, key=lambda x: x.start)
                            ],
                            "mapping_quality": round(m * 100 / t, 0),
                            "matches": m
                        }
                        domains.append(temp)
                    nrow[COLUMNS.fusion_mapped_domains] = json.dumps(domains)
                    ann_rows.append(nrow)
        except NotSpecifiedError as err:
            print(repr(err))
        except AttributeError as err:
            print(repr(err))
        except NotImplementedError as err:
            print(repr(err))

        # now try generating the svg
        DS = DiagramSettings()
        DS.DOMAIN_NAME_REGEX_FILTER = args.domain_regex_filter
        drawing = None
        retry_count = 0
        while drawing is None:  # continue if drawing error and increase width
            try:
                canvas, legend = draw_sv_summary_diagram(
                    DS, ann, ft, REFERENCE_GENOME=REFERENCE_GENOME, templates=TEMPLATES)

                gene_aliases1 = 'NA'
                gene_aliases2 = 'NA'
                try:
                    if len(ann.transcript1.gene.aliases) > 0:
                        gene_aliases1 = '-'.join(ann.transcript1.gene.aliases)
                    if ann.transcript1.is_best_transcript:
                        gene_aliases1 = 'b-' + gene_aliases1
                except AttributeError:
                    pass
                try:
                    if len(ann.transcript2.gene.aliases) > 0:
                        gene_aliases2 = '-'.join(ann.transcript2.gene.aliases)
                    if ann.transcript2.is_best_transcript:
                        gene_aliases2 = 'b-' + gene_aliases2
                except AttributeError:
                    pass
                try:
                    if determine_prime(ann.transcript1, ann.break1) == PRIME.THREE:
                        gene_aliases1, gene_aliases2 = gene_aliases2, gene_aliases1
                except NotSpecifiedError:
                    pass

                name = 'chr{}_chr{}.{}_{}.{}'.format(
                    ann.break1.chr, ann.break2.chr,
                    gene_aliases1, gene_aliases2, ann.data[COLUMNS.annotation_id])

                drawing = os.path.join(DRAWINGS_DIRECTORY, name + '.svg')
                l = os.path.join(DRAWINGS_DIRECTORY, name + '.legend.json')
                for r in ann_rows + [row]:
                    r[COLUMNS.annotation_figure] = drawing
                    r[COLUMNS.annotation_figure_legend] = l
                log('generating svg:', drawing, time_stamp=False)
                canvas.saveas(drawing)

                log('generating legend:', l, time_stamp=False)
                with open(l, 'w') as fh:
                    json.dump(legend, fh)
                break
            except DrawingFitError as err:
                log('extending width:', DS.WIDTH, DS.WIDTH + 500, time_stamp=False)
                DS.WIDTH += 500
                retry_count += 1
                if retry_count > 10:
                    raise err
        if len(ann_rows) == 0:
            rows.append(row)
        else:
            rows.extend(ann_rows)

    with open(TABBED_OUTPUT_FILE, 'w') as fh:
        log('writing:', TABBED_OUTPUT_FILE)
        header = set()

        for row in rows:
            header.update(row.keys())

        header = sort_columns([str(c) for c in header if not str(c).startswith('_')])
        fh.write('\t'.join([str(c) for c in header]) + '\n')

        for i, row in enumerate(rows):
            fh.write('\t'.join([str(row.get(c, None)) for c in header]) + '\n')

        log('generated {} annotations'.format(len(annotations)))

    with open(FA_OUTPUT_FILE, 'w') as fh:
        log('writing:', FA_OUTPUT_FILE)
        for name, seq in sorted(fa_sequences.items()):
            fh.write('> {}\n'.format(name))
            fh.write('{}\n'.format(seq))


if __name__ == '__main__':
    main()
