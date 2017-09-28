from ..interval import Interval


class ScatterPlot:
    """
    holds settings that will go into matplotlib after conversion using the mapping system
    """

    def __init__(
        self, points, y_axis_label,
        ymax=None, ymin=None, xmin=None, xmax=None, hmarkers=None, height=100, point_radius=2,
        title='', yticks=None, colors=None
    ):
        self.hmarkers = hmarkers if hmarkers is not None else []
        self.yticks = yticks if yticks is not None else []
        self.colors = colors if colors else {}
        self.ymin = ymin
        self.ymax = ymax
        self.points = points
        if self.ymin is None:
            self.ymin = min([y.start for x, y in points] + yticks)
        if self.ymax is None:
            self.ymax = max([y.end for x, y in points] + yticks)
        self.xmin = xmin
        self.xmax = xmax
        if self.xmin is None:
            self.xmin = min([x.start for x, y in points])
        if self.xmax is None:
            self.xmax = max([x.end for x, y in points])
        self.y_axis_label = y_axis_label
        self.height = 100
        self.point_radius = 2
        self.title = title


def draw_scatter(ds, canvas, plot, xmapping):
    """
    given a xmapping, draw the scatter plot svg group

    Args:
        ds (DiagramSettings): the settings/constants to use for building the svg
        canvas (svgwrite.canvas): the svgwrite object used to create new svg elements
        plot (ScatterPlot): the plot to be drawn
        xmapping (:class:`dict` of :class:`Interval` by :class:`Interval`):
            dict used for conversion of coordinates in the xaxis to pixel positions
    """
    # generate the y coordinate mapping
    plot_group = canvas.g(class_='scatter_plot')

    yratio = plot.height / (abs(plot.ymax - plot.ymin))
    ypx = []
    xpx = []
    for xpo, ypo in plot.points:
        try:
            temp = Interval.convert_ratioed_pos(xmapping, xpo.start)
            xp = Interval.convert_ratioed_pos(xmapping, xpo.end)
            xp = xp | temp
            xpx.append((xp, xpo))
            temp = plot.height - abs(ypo.start - plot.ymin) * yratio
            yp = Interval(plot.height - abs(ypo.end - plot.ymin) * yratio, temp)
            ypx.append((yp, ypo))
        except IndexError:
            pass

    for x, y in zip(xpx, ypx):
        xp, xpo = x
        yp, ypo = y
        if xp.length() > ds.scatter_marker_radius:
            plot_group.add(canvas.line(
                (xp.start, yp.center),
                (xp.end, yp.center),
                stroke='#000000',
                stroke_width=ds.scatter_error_bar_stroke_width
            ))
        if yp.length() > ds.scatter_marker_radius:
            plot_group.add(canvas.line(
                (xp.center, yp.start),
                (xp.center, yp.end),
                stroke='#000000',
                stroke_width=ds.scatter_error_bar_stroke_width
            ))
        plot_group.add(canvas.circle(
            center=(xp.center, yp.center),
            fill=plot.colors.get((xpo, ypo), '#000000'),
            r=ds.scatter_marker_radius
        ))

    xmax = Interval.convert_ratioed_pos(xmapping, plot.xmax).end
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
