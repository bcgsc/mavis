import os
from typing import Dict, List, Tuple

from . import annotate as _annotate
from . import util as _util
from .annotate.file_io import ReferenceFile
from .error import DrawingFitError
from .illustrate.constants import DiagramSettings
from .illustrate.diagram import draw_multi_transcript_overlay
from .illustrate.scatter import bam_to_scatter


def check_overlay_args(args, parser):
    """
    parse the overlay options and check the formatting
    """
    # check complex options
    for marker in args.markers:
        if len(marker) < 3:
            marker.append(marker[-1])
        try:
            marker[1] = int(marker[1])
            marker[2] = int(marker[2])
        except ValueError:
            parser.error('argument --marker: start and end must be integers: {}'.format(marker))

    defaults = [None, None, 0.5, None, True]
    bam_file, density, ymax, stranded = range(1, 5)

    for plot in args.read_depth_plots:
        for i, d in enumerate(defaults):
            if i >= len(plot):
                plot.append(d)
        if not os.path.exists(plot[bam_file]):
            parser.error(
                'argument --read_depth_plots: the bam file given does not exist: {}'.format(
                    plot[bam_file]
                )
            )
        try:
            plot[density] = float(plot[density])
            if plot[density] < 0 or plot[density] > 1:
                raise ValueError()
        except ValueError:
            parser.error(
                'argument --read_depth_plots: density must be an float between 0 and 1: {}'.format(
                    plot[density]
                )
            )
        try:
            if str(plot[ymax]).lower() in ['null', 'none']:
                plot[ymax] = None
            else:
                plot[ymax] = int(plot[ymax])
        except ValueError:
            parser.error(
                'argument --read_depth_plots: ymax must be an integer: {}'.format(plot[ymax])
            )
        try:
            plot[stranded] = _util.cast_boolean(plot[stranded])
        except TypeError:
            parser.error(
                'argument --read_depth_plots: stranded must be an boolean: {}'.format(
                    plot[stranded]
                )
            )
    return args


def main(
    gene_name: str,
    output: str,
    config: Dict,
    buffer_length: int,
    read_depth_plots,
    markers: List[Tuple[str, int, int]],
    ymax_color='#FF0000',
    **kwargs,
):
    """
    generates an overlay diagram
    """
    annotations = ReferenceFile.load_from_config(config, 'annotations')
    annotations.load()
    drawing_width_iter_increase = config['illustrate.drawing_width_iter_increase']
    max_drawing_retries = config['illustrate.max_drawing_retries']
    min_mapping_quality = config['validate.min_mapping_quality']
    # check options formatting
    gene_to_draw = None

    for chrom in annotations.content:
        for gene in annotations.content[chrom]:
            if gene_name in gene.aliases or gene_name == gene.name:
                gene_to_draw = gene
                _util.logger.info(
                    f'Found target gene: {gene.name}(aka. {gene.aliases}) {gene.chr}:{gene.start}-{gene.end}'
                )
                break
    if gene_to_draw is None:
        raise KeyError('Could not find gene alias or id in annotations file', gene_name)

    settings = DiagramSettings(**kwargs)

    genomic_min = max(gene_to_draw.start - buffer_length, 1)
    genomic_max = gene_to_draw.end + buffer_length

    plots = []
    for axis_name, bam_file, density, ymax, stranded in read_depth_plots:
        # one plot per bam
        plots.append(
            bam_to_scatter(
                bam_file,
                gene_to_draw.chr,
                genomic_min,
                genomic_max,
                strand=gene_to_draw.get_strand() if stranded else None,
                ymax=ymax,
                density=density,
                axis_name=axis_name,
                min_mapping_quality=min_mapping_quality,
                ymax_color=ymax_color,
            )
        )

    vmarkers = []

    for i, (marker_name, marker_start, marker_end) in enumerate(markers):
        vmarkers.append(
            _annotate.base.BioInterval(gene_to_draw.chr, marker_start, marker_end, name=marker_name)
        )

    canvas = None
    attempts = 1
    while True:
        try:
            canvas = draw_multi_transcript_overlay(
                settings,
                gene_to_draw,
                vmarkers=vmarkers,
                plots=plots,
                window_buffer=buffer_length,
            )
            break
        except DrawingFitError as err:
            if attempts > max_drawing_retries:
                raise err
            _util.logger.info(f'Drawing fit: extending window {drawing_width_iter_increase}')
            settings.width += drawing_width_iter_increase
            attempts += 1

    svg_output_file = os.path.join(output, '{}_{}_overlay.svg'.format(gene_to_draw.name, gene_name))
    _util.logger.info(f'writing: {svg_output_file}')

    canvas.saveas(svg_output_file)
