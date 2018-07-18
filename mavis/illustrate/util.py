from colour import Color
import svgwrite

from ..error import DrawingFitError
from ..interval import Interval, IntervalMapping


MIN_PIXEL_ACCURACY = 1


def dynamic_label_color(color):
    """
    calculates the luminance of a color and determines if a black or white label will be more contrasting
    """
    color = Color(color)
    if color.get_luminance() < 0.5:
        return '#FFFFFF'
    return '#000000'


class LabelMapping:

    def __init__(self, **kwargs):
        self._mapping = dict()
        self._reverse_mapping = dict()
        for attr, val in kwargs.items():
            self[attr] = val

    def __setitem__(self, key, value):
        if key in self._mapping:
            raise KeyError('duplicate key: keys must be unique', key)
        if value in self._reverse_mapping:
            raise KeyError('duplicate value: values must be unique', value)
        self._mapping[key] = value
        self._reverse_mapping[value] = key

    def items(self):
        return self._mapping.items()

    def __getitem__(self, key):
        return self._mapping[key]

    def __len__(self):
        return len(self._mapping.keys())

    def get_key(self, value):
        return self._reverse_mapping[value]

    def set_key(self, key, value):
        if key in self._mapping:
            current_value = self._mapping[key]
            if value == current_value:
                return
            elif value in self._reverse_mapping:
                raise KeyError('duplicate value: values must be unique', value)
            del self._mapping[key]
            del self._reverse_mapping[current_value]
        elif value in self._reverse_mapping:
            raise KeyError('duplicate value: values must be unique', value)
        self[key] = value

    def add(self, value, prefix=''):
        if value in self._reverse_mapping:
            return self._reverse_mapping[value]
        i = 1
        while True:
            key = '{}{}'.format(prefix, i)
            if key not in self._mapping:
                self[key] = value
                break
            i += 1
        return self._reverse_mapping[value]


def split_intervals_into_tracks(intervals):
    tracks = [[]]
    for itvl in sorted(intervals, key=lambda x: x[0]):
        added = False
        for track in tracks:
            overlaps = False
            for track_itvl in track:
                if Interval.overlaps(itvl, track_itvl):
                    overlaps = True
                    break
            if not overlaps:
                added = True
                track.append(itvl)
                break
        if not added:
            tracks.append([itvl])
    return tracks


def generate_interval_mapping(
        input_intervals, target_width, ratio, min_width,
        buffer_length=None, start=None, end=None, min_inter_width=None,
        min_pixel_accuracy=MIN_PIXEL_ACCURACY):
    min_inter_width = min_width if min_inter_width is None else min_inter_width
    if all([x is not None for x in [start, end, buffer_length]]):
        raise AttributeError('buffer_length is a mutually exclusive argument with start/end')
    if not input_intervals and (start is None or end is None):
        raise AttributeError('must specify an interval or start/end at minimum to generate an interval mapping')
    intervals = []
    for i in Interval.min_nonoverlapping(*input_intervals):
        if not intervals or abs(Interval.dist(intervals[-1], i)) > 1:
            intervals.append(i)
        else:
            intervals[-1] = intervals[-1] | i
    # break up the intervals by any intervals of length 1
    for itvl_in in input_intervals:
        if len(itvl_in) > 1:
            continue
        # try splitting all current interval
        temp = []
        for itvl in intervals:
            split = itvl - itvl_in
            if split is not None:
                temp.extend(split)
        intervals = temp
    for itvl_in in input_intervals:
        if len(itvl_in) == 1:
            intervals.append(itvl_in)
    # now split any intervals by start/end
    breaks = {}
    for i in intervals:
        # split by input intervals
        breaks[i] = set([i.start, i.end])
        for ii in input_intervals:
            if ii.start >= i.start and ii.start <= i.end:
                breaks[i].add(ii.start)
            if ii.end >= i.start and ii.end <= i.end:
                breaks[i].add(ii.end)
    temp = []
    for itvl, breakpoints in breaks.items():
        breakpoints.add(itvl.start)
        breakpoints.add(itvl.end)
        pos = sorted(breakpoints)
        if len(pos) == 1:
            temp.append(Interval(pos[0]))
        else:
            # remove all the single intervals to start?
            pos[0] -= 1
            for i in range(1, len(pos)):
                temp.append(Interval(pos[i - 1] + 1, pos[i]))
    intervals = sorted(temp, key=lambda x: x.start)

    if buffer_length is None:
        buffer_length = 0

    if start is None:
        start = max(intervals[0].start - buffer_length, 1)
    elif start <= 0:
        raise AttributeError('start must be a natural number', start)

    if end is None:
        end = intervals[-1].end + buffer_length
    elif end <= 0:
        raise AttributeError('end must be a natural number', end)

    if not intervals:  # if no input intervals are given, then use the start/end of the entire range as the focus
        intervals = [Interval(start, end)]

    total_length = end - start + 1
    genic_length = sum([len(i) for i in intervals])
    intergenic_length = total_length - genic_length
    intermediate_intervals = 0
    if start < intervals[0].start:
        intermediate_intervals += 1
    if end > intervals[-1].end:
        intermediate_intervals += 1

    for i in range(1, len(intervals)):
        if intervals[i].start > intervals[i - 1].end + 1:
            intermediate_intervals += 1
    width = target_width - intermediate_intervals * min_inter_width - len(intervals) * min_width  # reserved width

    if width < 0:
        raise DrawingFitError('width cannot accommodate the number of expected objects')

    intergenic_width = width // (ratio + 1) if intergenic_length > 0 else 0
    genic_width = width - intergenic_width

    def intergenic_unit(x):
        return x * intergenic_width / intergenic_length

    def genic_unit(x):
        return x * genic_width / genic_length
    mapping = []

    pos = 0
    # do the intergenic region prior to the first genic region
    if start < intervals[0].start:
        ifrom = Interval(start, intervals[0].start - 1)
        s = max(intergenic_unit(len(ifrom)), 0)
        ito = Interval(pos, pos + min_inter_width + s)
        mapping.append((ifrom, ito))
        pos = ito.end

    for i, curr in enumerate(intervals):
        if i > 0 and intervals[i - 1].end + 1 < curr.start:  # add between the intervals

            prev = intervals[i - 1]
            ifrom = Interval(prev.end + 1, curr.start - 1)
            s = max(intergenic_unit(len(ifrom)), 0)
            ito = Interval(pos, pos + min_inter_width + s)
            mapping.append((ifrom, ito))
            pos = ito.end

        s = max(genic_unit(len(curr)), 0)
        assert s >= 0
        ito = Interval(pos, pos + min_width + s)
        mapping.append((curr, ito))
        pos = ito.end

    # now the last intergenic region will make up for the rounding error
    if end > intervals[-1].end:
        ifrom = Interval(intervals[-1].end + 1, end)
        s = max(intergenic_unit(len(ifrom)), 0)
        ito = Interval(pos, pos + min_inter_width + s)
        mapping.append((ifrom, ito))
        pos = ito.end
    # mapping[-1][1].end = target_width  # min(int(target_width), mapping[-1][1].end)
    if abs(mapping[-1][1].end - target_width) > min_pixel_accuracy:
        raise AssertionError(
            'end is off by more than the expected pixel allowable error',
            mapping[-1][1].end, target_width, min_pixel_accuracy)
    mapping = {k: v for k, v in mapping}
    # assert that that mapping is correct
    for ifrom, ito in mapping.items():
        if ifrom not in input_intervals:
            continue
        if ito.length() < min_width and abs(ito.length() - min_width) > min_pixel_accuracy:  # precision error allowable
            raise AssertionError(
                'interval mapping should not map any intervals to less than the minimum required width. Interval {}'
                ' was mapped to a pixel interval of length {} but the minimum width is {}'.format(
                    ifrom, ito.length(), min_width), mapping,
                input_intervals, target_width, ratio, min_width, buffer_length, start, end, min_inter_width)
    return IntervalMapping(mapping)


class Tag(svgwrite.base.BaseElement):

    def __init__(self, elementname, content='', **kwargs):
        self.elementname = elementname
        super(Tag, self).__init__(**kwargs)
        self.content = content

    def get_xml(self):
        xml = super(Tag, self).get_xml()
        xml.text = self.content
        return xml
