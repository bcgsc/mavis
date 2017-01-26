
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
    |       |-- annotation-#.failed
    |       |-- annotation-#.passed
    |       `-- annotation-#_drawings/
    |           |-- legend.svg
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

"""
import argparse
from structural_variant.breakpoint import read_bpp_from_input_file
from structural_variant.annotate import load_reference_genes, load_reference_genome, load_templates
from structural_variant.annotate.variant import gather_annotations, FusionTranscript
from structural_variant.error import DiscontinuousMappingError, DrawingFitError, NotSpecifiedError
from structural_variant import __version__
from structural_variant.draw import Diagram
import TSV
from structural_variant.constants import PROTOCOL, SVTYPE, COLUMNS, sort_columns
import re
import json
import os
from datetime import datetime


def log(*pos, time_stamp=True):
    if time_stamp:
        print('[{}]'.format(datetime.now()), *pos)
    else:
        print(' ' * 28, *pos)


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='Outputs the version number'
    )
    parser.add_argument(
        '-f', '--overwrite',
        action='store_true', default=False,
        help='set flag to overwrite existing reviewed files'
    )
    parser.add_argument(
        '--no_draw', default=True, action='store_false',
        help='set flag to suppress svg drawings of putative annotations')
    parser.add_argument(
        '-o', '--output',
        help='path to the output directory', required=True
    )
    parser.add_argument(
        '-n', '--input',
        help='path to the input file', required=True
    )
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
    parser.add_argument(
        '-p', '--max_proximity', default=5000,
        help='The maximum distance away from breakpoints to look for proximal genes')
    parser.add_argument(
        '--min_orf_size', default=120, type=int, help='minimum size for putative ORFs')
    parser.add_argument(
        '--max_orf_cap', default=3, type=int, help='keep the n longest orfs')
    parser.add_argument(
        '--min_domain_mapping_match', default=0.8, type=float, 
        help='minimum percent match for the domain to be considered aligned')
    parser.add_argument(
        '--template_metadata', default='/home/creisle/svn/svmerge/trunk/cytoBand.txt',
        help='file containing the cytoband template information')
    args = parser.parse_args()
    return args


def main():
    # load the file
    args = parse_arguments()
    log('input arguments listed below')
    for arg, val in sorted(args.__dict__.items()):
        log(arg, '=', val, time_stamp=False)
    FILENAME_PREFIX = re.sub('\.(tab|tsv|txt)$', '', os.path.basename(args.input))
    
    log('loading:', args.template_metadata)
    TEMPLATES = load_templates(args.template_metadata)
    log('loading:', args.annotations)
    REFERENCE_ANNOTATIONS = load_reference_genes(args.annotations)
    
    
    log('loading:', args.reference_genome)
    REFERENCE_GENOME = load_reference_genome(args.reference_genome)
    # test that the sequence makes sense for a random transcript
    log('loading:', args.input)
    bpps = read_bpp_from_input_file(
        args.input,
        require=[COLUMNS.cluster_id, COLUMNS.validation_id],
        cast={
            COLUMNS.stranded.name: TSV.tsv_boolean
        },
        _in={
            COLUMNS.protocol: PROTOCOL,
            COLUMNS.event_type: SVTYPE
        },
        simplify=False)
    log('read {} breakpoint pairs'.format(len(bpps)))

    annotations = []
    for bpp in bpps:
        log('gathering annotations for', bpp)
        ann = gather_annotations(
            REFERENCE_ANNOTATIONS,
            bpp,
            event_type=bpp.data[COLUMNS.event_type],
            proximity=args.max_proximity
        )
        annotations.extend(ann)
        log('generated', len(ann), 'annotations', time_stamp=False)

    id_prefix = 'annotation_{}-'.format(re.sub(':', '-', re.sub(' ', '_', str(datetime.now()))))
    rows = []  # hold the row information for the final tsv file
    for i, ann in enumerate(annotations):
        ann.data[COLUMNS.annotation_id] = id_prefix + str(i + 1)
        row = ann.flatten()
        # try building the fusion product
        d = Diagram()
        ann_rows = []
        drawing = None
        retry_count = 0
        while drawing is None:  # continue if drawing error and increase width
            try:
                ft = None
                try:
                    ft = FusionTranscript.build(
                        ann, REFERENCE_GENOME,
                        min_orf_size=args.min_orf_size,
                        max_orf_cap=args.max_orf_cap,
                        min_domain_mapping_match=args.min_domain_mapping_match
                    )
                    # add fusion information to the current row
                    # duplicate the row for each translation
                    for tl in ft.translations:
                        nrow = dict()
                        nrow.update(row)
                        nrow[COLUMNS.fusion_splicing_pattern] = tl.transcript.splicing_pattern.splice_type
                        nrow[COLUMNS.fusion_cdna_coding_start] = tl.start
                        nrow[COLUMNS.fusion_cdna_coding_end] = tl.end
                        nrow[COLUMNS.fusion_cdna_sequence] = ft.get_cdna_sequence(tl.transcript.splicing_pattern)
                        domains = []
                        for dom in tl.domains:
                            m, t = dom.score_region_mapping()
                            temp = {
                                "name": dom.name,
                                "sequences": dom.get_sequences(),
                                "regions": [{"start": dr.start, "end": dr.end} for dr in sorted(dom.regions, key=lambda x: x.start)],
                                "mapping_quality": round(m * 100 / t, 0),
                                "matches": m
                            }
                            domains.append(temp)
                        nrow[COLUMNS.fusion_mapped_domains] = json.dumps(domains)
                        ann_rows.append(nrow)
                except NotSpecifiedError:
                    pass
                except AttributeError as err:
                    pass
                canvas = d.draw(ann, ft, REFERENCE_GENOME=REFERENCE_GENOME, draw_template=True, templates=TEMPLATES)
                drawing = os.path.join(args.output, FILENAME_PREFIX + '.' + ann.data[COLUMNS.annotation_id] + '.svg')
                for r in ann_rows:
                    r[COLUMNS.annotation_figure] = drawing
                log('generating svg:', drawing)
                canvas.saveas(drawing)
                break
            except (NotImplementedError, AttributeError, DiscontinuousMappingError) as err:
                print(repr(err))
                raise err
            except DrawingFitError as err:
                print(repr(err))
                print('extending width:', d.WIDTH, d.WIDTH + 500)
                d.WIDTH += 500
                retry_count += 1
                if retry_count > 3:
                    raise err
        if len(ann_rows) == 0:
            rows.append(row)
        else:
            rows.extend(ann_rows)


    of = os.path.join(args.output, FILENAME_PREFIX + '.annotation.tab')
    with open(of, 'w') as fh:
        log('writing:', of)
        header = set()

        for row in rows:
            header.update(row.keys())

        header = sort_columns([str(c) for c in header if not str(c).startswith('_')])
        fh.write('\t'.join([str(c) for c in header]) + '\n')

        for i, row in enumerate(rows):
            fh.write('\t'.join([str(row.get(c, None)) for c in header]) + '\n')

        log('generated {} annotations'.format(len(annotations)))


if __name__ == '__main__':
    main()
