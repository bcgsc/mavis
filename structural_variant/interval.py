import numpy as np

from structural_variant.constants import *
from structural_variant.error import *


class Interval:
    """
    Intervals are inclusive
    """

    def __init__(self, start, end=None, freq=1):
        self.start = int(start)
        self.end = int(end) if end is not None else self.start
        if self.start > self.end:
            raise AttributeError('interval start > end is not allowed')
        self.freq = int(freq)
        if self.freq <= 0:
            raise AttributeError('Interval frequency must be a natural number')

    def __getitem__(self, index):
        try:
            index = int(index)
        except ValueError:
            raise IndexError('index input accessor must be an integer', index)
        if index == 0:
            return self.start
        elif index == 1:
            return self.end
        elif index == 2:
            return self.freq
        raise IndexError(
            'index input accessor is out of bounds: 1 or 2 only', index)

    def overlaps(self, other):
        """
        >>> x, y, z = ( Interval(1, 4), Interval(-1, 0), Interval(0, 3) )
        >>> x.overlap(y)
        False
        >>> x.overlap(z)
        True
        >>> y.overlap(x)
        False
        >>> y.overlap(z)
        True
        """
        if self - other == 0:
            return True
        return False

    def __len__(self):
        return self[1] - self[0] + 1

    def __repr__(self):
        temp = str(self[0])
        if self[1] != self[0]:
            temp += '-' + str(self[1])
        if self.freq != 1:
            temp += 'x' + str(self.freq)
        return self.__class__.__name__ + '(' + temp + ')'

    @property
    def center(self):
        """
        >>> x, y, z = ( Interval(1, 5), Interval(-1, 0), Interval(1, 2) )
        >>> x.center
        3.0
        >>> y.center
        -0.5
        >>> z.center
        1.5
        """
        return self[0] + (len(self) - 1) / 2

    @classmethod
    def weighted_mean(cls, intervals):
        """
        returns the weighted mean for a set of intervals
        the weight is the inverse of the size of the interval
        so that broader intervals are weighted less than
        more specific/tighter intervals
        """
        if len(intervals) == 0:
            raise AttributeError('input list cannot be empty')
        first = next(iter(intervals))
        centers = []
        weights = []

        for i in intervals:
            for temp in range(0, i.freq):
                centers.append(i.center)
                weights.append(1 / len(i))

        return np.average(centers, weights=weights)

    def combine(self, other):
        """
        adding two intervals returns the minimum interval that covers both input intervals
        >>> x, y, z = ( Interval(1, 4), Interval(-1, 0), Interval(1, 2) )
        >>> x.combine(y)
        Interval(-1-4)
        >>> x.combine(z)
        Interval(1-4)
        >>> y.combine(z)
        Interval(-1-2)
        """
        return Interval.union(self, other)

    def __eq__(self, other):
        if not hasattr(other, 'start') \
                or not hasattr(other, 'end') \
                or not hasattr(other, 'freq') \
                or self[0] != other[0] \
                or self[1] != other[1] \
                or self[2] != other[2]:
            return False
        return True

    def __lt__(self, other):
        if self[0] < other[0]:
            return True
        elif self[0] == other[0]:
            if self[1] < other[1]:
                return True
            elif self[1] == other[1]:
                if self.freq < other.freq:
                    return True
        return False

    def __gt__(self, other):
        if self[0] > other[0]:
            return True
        elif self[0] == other[0]:
            if self[1] > other[1]:
                return True
            elif self[1] == other[1]:
                if self.freq > other.freq:
                    return True
        return False

    def __sub__(self, other):
        """returns the minimum distance between intervals

        >>> x, y, z = ( Interval(1, 4), Interval(-1, 0), Interval(0, 3) )
        >>> x - y
        1
        >>> y - x
        -1
        >>> x - z
        0
        >>> z - x
        0
        """
        if self[1] < other[0]:
            return self[1] - other[0]
        elif self[0] > other[1]:
            return self[0] - other[1]
        else:
            return 0

    def __hash__(self):
        return hash((self[0], self[1], self.freq))

    @classmethod
    def paired_weighted_means(cls, intervals):
        int_a = Interval.weighted_mean([x[0] for x in intervals])
        int_b = Interval.weighted_mean([x[1] for x in intervals])
        return int_a, int_b

    @classmethod
    def position_in_range(cls, segments, pos):
        """
        >>> b = (12, 12)
        >>> Interval.position_in_range([(1, 2), (3, 6), (7, 15)], b)
        (2, False)
        >>> Interval.position_in_range([(1, 2), (3, 6), (7, 10), (14, 16)], b)
        (3, True)
        >>> Interval.position_in_range([(1, 2), (3, 6), (7, 10)], b)
        (3, False)
        >>> Interval.position_in_range([(15, 16), (17, 19)], b)
        (0, True)
        """
        if len(segments) == 0:
            raise AttributeError('cannot compute on an empty list')
        num = 0
        found_inbetween_segment = False

        segments = sorted(segments)

        while num < len(segments):
            current = segments[num]

            if pos[0] >= current[0] \
                    and pos[1] <= current[1]:
                # pos range is fully contained in the current segment
                break
            elif num == 0:  # first segment
                if pos[1] < current[0]:
                    # before the first segment
                    found_inbetween_segment = True
                    break
            else:
                # check the previous segment
                previous = segments[num - 1]
                if pos[0] > previous[1] and pos[1] < current[0]:
                    found_inbetween_segment = True
                    break
            num += 1
        return num, found_inbetween_segment

    @classmethod
    def convert_pos(cls, mapping, pos):
        """ convert any given position given a mapping of intervals to another range

        @param mapping \a required (type: \b Dict<Interval, Interval>) a mapping of a set of continuous intervals
        @param pos \a required (type: \b int) a position in the first coordinate system
        @param front_to_end \a optional (type: \b boolean) start-start & end-end or start-end & end-start

        @return (type: \b int) the position in the alternate coordinate system given the input mapping

        @except AttributeError if the input position is outside the set of input segments
        @except DiscontiuousMappingError if the input position cannot be converted to the output system

        >>> mapping = {(1, 10): (101, 110), (21, 30): (201, 210), (41, 50): (301, 310)}
        >>> Interval.convert_pos(mapping, 5)
        105
        >>> Interval.convert_pos(mapping, 15)
        Traceback (most recent call last):
        ...
        structural_variant.error.DiscontiuousMappingError: DiscontiuousMappingError<between=(110, 201), outside mapped range>
        >>> Interval.convert_pos(mapping, 0)
        Traceback (most recent call last):
        ...
        structural_variant.error.DiscontiuousMappingError: DiscontiuousMappingError<before=101, outside mapped range>
        >>> Interval.convert_pos(mapping, 80)
        Traceback (most recent call last):
        ...
        structural_variant.error.DiscontiuousMappingError: DiscontiuousMappingError<after=310, outside mapped range>
        >>> mapping = {(41, 50): (101, 110), (21, 30): (201, 210), (1, 10): (301, 310)}
        >>> Interval.convert_pos(mapping, 5)
        306
        """
        if len(mapping.keys()) < 2:
            raise AttributeError(
                'mapping is insufficient to determine orientation')

        # order the input intervals
        input_intervals = sorted(mapping.keys())

        # input checking
        for curr in input_intervals:
            if curr[1] - curr[0] != mapping[curr][1] - mapping[curr][0]:
                raise AttributeError(
                    'input mappings must have segments of equal length')

        front_to_end = None
        for i, curr in enumerate(input_intervals):
            if i == 0:
                continue
            curr = mapping[curr]
            prev = mapping[input_intervals[i - 1]]
            if front_to_end is None:
                if prev[1] < curr[0]:
                    front_to_end = False
                elif curr[1] < prev[0]:
                    front_to_end = True
                else:
                    raise AttributeError('mapping must be two non-overlapping sequences mapped one-to-one forward or '
                                         'reversed', prev, curr)
            elif not front_to_end:
                if prev[1] >= curr[0]:
                    raise AttributeError('mapping must be two non-overlapping sequences mapped one-to-one forward or '
                                         'reversed prev[1] >= curr[0]', prev, curr)
            else:
                if curr[1] >= prev[0]:
                    raise AttributeError('mapping must be two non-overlapping sequences mapped one-to-one forward or '
                                         'reversed curr[1] >= prev[0]', prev, curr)

        i, previous_flag = Interval.position_in_range(
            input_intervals, (pos, pos))  # get the input position
        if i == len(input_intervals):
            curr = input_intervals[i - 1]
            if not front_to_end:
                raise DiscontiuousMappingError('outside mapped range', after=mapping[curr][1])
            else:
                raise DiscontiuousMappingError('outside mapped range', before=mapping[curr][0])
        elif previous_flag:
            curr = input_intervals[i]
            if i == 0:
                if not front_to_end:
                    raise DiscontiuousMappingError('outside mapped range', before=mapping[curr][0])
                else:
                    raise DiscontiuousMappingError('outside mapped range', after=mapping[curr][1])
            else:  # between two segments
                prev = input_intervals[i - 1]
                if not front_to_end:
                    raise DiscontiuousMappingError('outside mapped range', between=(mapping[prev][1], mapping[curr][0]))
                else:
                    raise DiscontiuousMappingError('outside mapped range', between=(mapping[curr][1], mapping[prev][0]))
        else:
            # fell into a mapped region
            curr = input_intervals[i]
            if not front_to_end:
                shift = pos - curr[0]
                return mapping[curr][0] + shift
            else:
                shift = curr[1] - pos
                return mapping[curr][0] + shift

    @classmethod
    def paired_set_distance(cls, intervals, other_intervals):
        """
        for two sets of interval pairs (as tuples) computes the weighted mean
        of each interval set (a total of four) and then returns the 
        distance between the sets of pairs as the cumulative distance
        of the weighted means of their pairs
        """
        int_a, int_b = cls.paired_weighted_means(intervals)
        oint_a, oint_b = cls.paired_weighted_means(other_intervals)
        return abs(int_a - oint_a) + abs(int_b - oint_b)

    @classmethod
    def redundant_ordered_hierarchical_clustering(cls, clusters, r=None):
        """
        for an input set of of clusters, do hierarchical clustering
        redundant b/c we allow clusters to be grouped more than once
        into either of their immediate neighbours
        """
        r = int(r)
        if r < 0:
            raise AttributeError('r must be a positive integer')
        # order the clusters by weighted mean
        complete = []
        queue = sorted(clusters, key=lambda x: cls.paired_weighted_means(x))

        while len(queue) > 0:
            temp_queue = []
            for i in range(0, len(queue)):
                curr = queue[i]
                joined = False

                if i > 0:  # try joining your previous neighbor
                    dist = cls.paired_set_distance(curr, queue[i - 1])
                    if dist <= r:
                        temp_queue.append(curr.union(queue[i - 1]))
                        joined = True
                if i < len(queue) - 1:
                    dist = cls.paired_set_distance(curr, queue[i + 1])
                    if dist <= r:
                        temp_queue.append(curr.union(queue[i + 1]))
                        joined = True
                if not joined:
                    complete.append(curr)
            queue = temp_queue
        for c in clusters:
            for n in c:
                found = False
                for cc in complete:
                    for cn in cc:
                        if n == cn:
                            found = True
                            break
                if not found:
                    raise AssertionError('error node is no longer assigned', clusters, complete, n)
        return complete

    @classmethod
    def union(cls, *intervals):
        """
        returns the union of the set of input intervals

        >>> l = [Interval(1, 10), Interval(5, 7), Interval(7)]
        >>> Interval.union(l)
        Interval(1-10)
        >>> l.append(Interval(11))
        >>> Interval.union(l)
        Interval(1-11)
        """
        if len(intervals) < 1:
            raise AttributeError('cannot compute the union of an empty set of intervals')
        return Interval(min([i[0] for i in intervals]), max([i[1] for i in intervals]))

    @classmethod
    def intersection(cls, intervals):
        """
        returns None if there is no intersection

        >>> l = [Interval(1, 10), Interval(5, 7), Interval(7)]
        >>> Interval.intersection(l)
        Interval(7)
        >>> l.append(Interval(11))
        >>> Interval.intersection(l)
        """
        if len(intervals) < 1:
            raise AttributeError('cannot compute the intersection of an empty set of intervals')
        curr = next(iter(intervals))
        low = curr[0]
        high = curr[1]

        for i in intervals:
            if high < i[0] or i[1] < low:
                return None
            low = max(low, i[0])
            high = min(high, i[1])
        return Interval(low, high)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
