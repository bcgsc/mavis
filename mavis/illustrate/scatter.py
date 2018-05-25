import os

from ..bam.read import sequenced_strand, pileup
from ..util import LOG, DEVNULL
from ..interval import Interval
from ..validate.constants import DEFAULTS as VALIDATION_DEFAULTS


def bam_to_scatter(bam_file, chrom, start, end, density, strand=None, axis_name=None, ymax=None, min_mapping_quality=0, ymax_color='#FF0000'):
    """
    pull data from a bam file to set up a scatter plot of the pileup

    Args:
        bam_file (str): path to the bam file
        chrom (str): chromosome name
        start (int): genomic start position for the plot
        end (int): genomic end position for the plot
        bin_size (int): number of genomic positions to group together and average to reduce data
        strand (STRAND): expected strand
        axis_name (str): axis name
        ymax (int): maximum value to plot the y axis
        min_mapping_quality (int): minimum mapping quality for reads to be considered in the plot

    Returns:
        ScatterPlot: the scatter plot representing the bam pileup
    """
    import pysam
    if not axis_name:
        axis_name = os.path.basename(bam_file)
    # one plot per bam
    LOG('reading:', bam_file)
    plot = None
    samfile = pysam.AlignmentFile(bam_file, 'rb')

    def read_filter(read):
        if read.mapping_quality < min_mapping_quality:
            return True
        if strand is None:
            return False
        try:
            return sequenced_strand(read, VALIDATION_DEFAULTS.strand_determining_read) != strand
        except ValueError:
            return True

    try:
        points = []
        try:
            for refpos, count in pileup(samfile.fetch(chrom, start, end), filter_func=read_filter):
                if refpos <= end and refpos >= start:
                    points.append((refpos, count))
        except ValueError:  # chrom not in bam
            pass

        LOG('scatter plot {} has {} points'.format(axis_name, len(points)))
        plot = ScatterPlot(
            points, axis_name,
            ymin=0,
            ymax=max([y for x, y in points] + [100]) if ymax is None else ymax,
            density=density,
            ymax_color=ymax_color
        )
    finally:
        samfile.close()
    return plot


class ScatterPlot:
    """
    holds settings that will go into matplotlib after conversion using the mapping system
    """

    def __init__(
        self, points, y_axis_label,
        ymax=None, ymin=None, xmin=None, xmax=None, hmarkers=None, height=100, point_radius=2,
        title='', yticks=None, colors=None, density=1, ymax_color='#FF0000'
    ):
        self.hmarkers = hmarkers if hmarkers is not None else []
        self.yticks = yticks if yticks is not None else []
        self.colors = colors if colors else {}
        self.ymin = ymin
        self.ymax = ymax
        self.points = points
        if self.ymin is None and (yticks or points):
            self.ymin = min([y for x, y in points] + yticks)
        if self.ymax is None and (yticks or points):
            self.ymax = max([y for x, y in points] + yticks)
        self.xmin = xmin
        self.xmax = xmax
        if self.xmin is None and points:
            self.xmin = min([x for x, y in points])
        if self.xmax is None and points:
            self.xmax = max([x for x, y in points])
        self.y_axis_label = y_axis_label
        self.height = height
        self.point_radius = point_radius
        self.title = title
        self.ymax_color = ymax_color
        self.density = density


def draw_scatter(ds, canvas, plot, xmapping, log=DEVNULL):
    """
    given a xmapping, draw the scatter plot svg group

    Args:
        ds (DiagramSettings): the settings/constants to use for building the svg
        canvas (svgwrite.canvas): the svgwrite object used to create new svg elements
        plot (ScatterPlot): the plot to be drawn
        xmapping (:class:`dict` of :class:`Interval` by :class:`Interval`):
            dict used for conversion of coordinates in the xaxis to pixel positions
    """
    from shapely.geometry import Point as sPoint
    # generate the y coordinate mapping
    plot_group = canvas.g(class_='scatter_plot')

    yratio = plot.height / (abs(plot.ymax - plot.ymin))
    px_points = []
    circles = []
    for x_pos, y_pos in plot.points:
        try:
            x_px = Interval.convert_ratioed_pos(xmapping, x_pos, forward_to_reverse=False)
            y_px = Interval(plot.height - abs(min(y_pos, plot.ymax) - plot.ymin) * yratio)
            current_circle = sPoint(x_px.center, y_px.center).buffer(ds.scatter_marker_radius)
            if circles:
                ratio = circles[-1].intersection(current_circle).area / current_circle.area
                if ratio > plot.density:
                    continue
            circles.append(current_circle)
            px_points.append((x_px, y_px, plot.ymax_color if y_pos > plot.ymax else plot.colors.get((x_pos, y_pos), '#000000')))
        except IndexError:
            pass
    log('drew {} of {} points (density={})'.format(len(circles), len(plot.points), plot.density), time_stamp=False)

    for x_px, y_px, color in px_points:
        if x_px.length() > ds.scatter_marker_radius:
            plot_group.add(canvas.line(
                (x_px.start, y_px.center),
                (x_px.end, y_px.center),
                stroke='#000000',
                stroke_width=ds.scatter_error_bar_stroke_width
            ))
        if y_px.length() > ds.scatter_marker_radius:
            plot_group.add(canvas.line(
                (x_px.center, y_px.start),
                (x_px.center, y_px.end),
                stroke='#000000',
                stroke_width=ds.scatter_error_bar_stroke_width
            ))
        plot_group.add(canvas.circle(
            center=(x_px.center, y_px.center),
            fill=color,
            r=ds.scatter_marker_radius
        ))

    xmax = Interval.convert_ratioed_pos(xmapping, plot.xmax, forward_to_reverse=False).end
    for py in plot.hmarkers:
        py = plot.height - abs(py - plot.ymin) * yratio
        plot_group.add(
            canvas.line(
                start=(0, py),
                end=(xmax, py),
                stroke='blue'
            )
        )
    # draw left y axis
    plot_group.add(canvas.line(
        start=(0, 0), end=(0, plot.height), stroke='#000000'
    ))
    ytick_labels = [0]
    # draw start and end markers on the y axis
    for y in plot.yticks:
        ytick_labels.append(len(str(y)))
        py = plot.height - abs(y - plot.ymin) * yratio
        plot_group.add(
            canvas.line(
                start=(0 - ds.scatter_yaxis_tick_size, py),
                end=(0, py),
                stroke='#000000'
            ))
        plot_group.add(
            canvas.text(
                str(y),
                insert=(
                    0 - ds.scatter_yaxis_tick_size - ds.padding,
                    py + ds.scatter_ytick_font_size * ds.font_central_shift_ratio),
                fill=ds.label_color,
                style=ds.font_style.format(font_size=ds.scatter_ytick_font_size, text_anchor='end')
            ))

    shift = max(ytick_labels)
    x = 0 - ds.padding * 2 - ds.scatter_axis_font_size - ds.scatter_yaxis_tick_size - \
        ds.scatter_ytick_font_size * ds.font_width_height_ratio * shift
    y = plot.height / 2
    yaxis = canvas.text(
        plot.y_axis_label,
        insert=(x, y),
        fill=ds.label_color,
        style=ds.font_style.format(font_size=ds.scatter_axis_font_size, text_anchor='start'),
        class_='y_axis_label'
    )
    plot_group.add(yaxis)
    cx = len(plot.y_axis_label) * ds.font_width_height_ratio * ds.scatter_axis_font_size / 2
    yaxis.rotate(270, (x + cx, y))
    yaxis.translate(0, 0)

    y = plot.height
    setattr(plot_group, 'height', y)
    return plot_group
