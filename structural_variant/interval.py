import numpy as np
import itertools
import networkx as nx
from structural_variant.constants import *
from structural_variant.error import *


class Interval:
    """
    """
    def __init__(self, start, end=None, freq=1):
        """
        Args:
            start (int): the start of the interval (inclusive)
            end (int, default=start): the end of the interval (inclusive)
            freq (int, default=1): the frequency or weight of the interval
        """
        self.start = int(start)
        self.end = int(end) if end is not None else self.start
        if self.start > self.end:
            raise AttributeError('interval start > end is not allowed')
        self.freq = int(freq)
        if self.freq <= 0:
            raise AttributeError('Interval frequency must be a natural number')

    def __sub__(self, other):  # difference
        """the difference of two intervals

        Example:
            >>> Interval(1, 10) - Interval(5, 50)
            [Interval(1, 4)]
            >>> Interval(1, 2) - Interval(10, 11)
            [Interval(1, 2)]
            >>> Interval(1, 2) - Interval(-1, 10)
            []
        """
        if Interval.overlaps(self, other):
            if other[0] <= self[0]:
                if other[1] >= self[1]:
                    return []
                else:
                    return [Interval(other[1] + 1, self[1])]
            elif other[1] >= self[1]:
                return [Interval(self[0], other[0] - 1)]
            else:
                return [Interval(self[0], other[0] - 1), Interval(other[1] + 1, self[1])]
        else:
            return [Interval(self[0], self[1])]

    def __and__(self, other):  # intersection
        """the intersection of two intervals

        Example:
            >>> Interval(1, 10) & Interval(5, 50)
            Interval(5, 10)
            >>> Interval(1, 2) & Interval(10, 11)
            None
        """
        return Interval.intersection(self, other)

    def __or__(self, other):  # union
        """the union of two intervals

        Example:
            >>> Interval(1, 10) | Interval(5, 50)
            Interval(1, 50)
            >>> Interval(1, 2) | Interval(10, 11)
            Interval(1, 11)
        """
        return Interval.union(self, other)

    def __xor__(self, other):
        return (self - other) + (other - self)

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

    @classmethod
    def overlaps(cls, self, other):
        """
        checks if two intervals have any portion of their given ranges in common

        Example:
            >>> Interval.overlaps(Interval(1, 4), Interval(5, 7))
            False
            >>> Interval.overlaps(Interval(1, 10), Interval(10, 11))
            True
            >>> Interval.overlaps((1, 10), (10, 11))
            True
        """
        if self[1] < other[0]:
            return False
        elif self[0] > other[1]:
            return False
        else:
            return True

    def __len__(self):
        """
        the length of the interval

        Example:
            >>> len(Interval(1, 11))
            12
        """
        return self[1] - self[0] + 1

    def __lt__(self, other):
        if self[0] < other[0]:
            return True
        return False

    def __gt__(self, other):
        if self[1] > other[1]:
            return True
        return False

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
        the middle of the interval

        Example:
            >>> Interval(1, 10).center
            5.5
            >>> Interval(1, 11).center
            6
        """
        return self[0] + (len(self) - 1) / 2

    def __eq__(self, other):
        if self[0] != other[0] or self[1] != other[1] or self[2] != other[2]:
            return False
        return True

    def __contains__(self, other):
        if other[0] >= self[0] and other[1] <= self[1]:
            return True
        return False

    @classmethod
    def dist(cls, self, other):
        """returns the minimum distance between intervals
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
    def weighted_mean(cls, intervals):
        """
        Args:
            intervals (List[Interval]): a list of intervals

        Returns:
            Interval: the weighted mean interval of the input intervals

        Raises:
            AttributeError: if the input list is empty

        Example:
            >>> Interval.weighted_mean([(1, 2), (1, 9), (2, 10)])
            Interval(1, 4)
            >>> Interval.weighted_mean([(1, 1), (10, 10)])
            Interval(6)
        """
        centers = []
        weights = []
        lengths = []
        if len(intervals) == 0:
            raise AttributeError('cannot compute the weighted mean interval of an empty set of intervals')
        for i in intervals:
            if not isinstance(i, Interval):
                i = Interval(i[0], i[1])
            for temp in range(0, i.freq):
                centers.append(i.center)
                weights.append(1 / len(i))
                lengths.append(len(i))

        center = np.average(centers, weights=weights)
        size = np.average(lengths, weights=weights) - 1
        return Interval(round(center - size / 2, 0), round(center + size / 2, 0))

    @classmethod
    def position_in_range(cls, segments, pos):
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

        Args:
            mapping (Dict[Interval, Interval]): a mapping of a set of continuous intervals
            pos (int): a position in the first coordinate system

        Returns:
            int: the position in the alternate coordinate system given the input mapping

        Raises:
            AttributeError: if the input position is outside the set of input segments
            DiscontiuousMappingError: if the input position cannot be converted to the output system

        Example:
            >>> mapping = {(1, 10): (101, 110), (11, 20): (555, 564)}
            >>> Interval.convert_pos(mapping, 5)
            5
            >>> Interval.convert_pos(mapping, 15)
            559
        """
        if len(mapping.keys()) < 2:
            raise AttributeError(
                'mapping is insufficient to determine orientation')

        # order the input intervals
        input_intervals = sorted(mapping.keys())
        mapped_to_intervals = [mapping[i] for i in input_intervals]
        forward_to_reverse = None
        
        # input check interval ranges are non-overlapping and increasing/decreasing
        for i in range(0, len(input_intervals)):
            if input_intervals[i][1] - input_intervals[i][0] != mapped_to_intervals[i][1] - mapped_to_intervals[i][0]:
                raise AttributeError('input mappings must have segments of equal length')
            if i < 1:
                continue
            if Interval.overlaps(input_intervals[i - 1], input_intervals[i]):
                raise AttributeError(
                    'input intervals cannot be overlapping',
                    input_intervals[i], input_intervals[i - 1]
                )
            if Interval.overlaps(mapped_to_intervals[i - 1], mapped_to_intervals[i]):
                raise AttributeError(
                    'mapped_to intervals cannot be overlapping',
                    mapped_to_intervals[i], mapped_to_intervals[i - 1]
                )
            if mapped_to_intervals[i][0] > mapped_to_intervals[i - 1][1]:
                if forward_to_reverse is None:
                    forward_to_reverse = False
                elif forward_to_reverse:
                    raise AttributeError('direction of mapped intervals is not consistent')
            elif mapped_to_intervals[i][1] < mapped_to_intervals[i - 1][0]:
                if forward_to_reverse is None:
                    forward_to_reverse = True
                elif not forward_to_reverse:
                    raise AttributeError('direction of mapped intervals is not consistent')

        i, previous_flag = Interval.position_in_range(
            input_intervals, (pos, pos))  # get the input position
        if i == len(input_intervals):
            curr = input_intervals[i - 1]
            if not forward_to_reverse:
                raise DiscontiuousMappingError('outside mapped range', after=mapping[curr][1])
            else:
                raise DiscontiuousMappingError('outside mapped range', before=mapping[curr][0])
        elif previous_flag:
            curr = input_intervals[i]
            if i == 0:
                if not forward_to_reverse:
                    raise DiscontiuousMappingError('outside mapped range', before=mapping[curr][0])
                else:
                    raise DiscontiuousMappingError('outside mapped range', after=mapping[curr][1])
            else:  # between two segments
                prev = input_intervals[i - 1]
                if not forward_to_reverse:
                    raise DiscontiuousMappingError(
                        'outside mapped range',
                        before=mapping[prev][1],
                        after=mapping[curr][0]
                    )
                else:
                    raise DiscontiuousMappingError(
                        'outside mapped range',
                        before=mapping[curr][1],
                        after=mapping[prev][0]
                    )
        else:
            # fell into a mapped region
            curr = input_intervals[i]
            if not forward_to_reverse:
                shift = pos - curr[0]
                return mapping[curr][0] + shift
            else:
                shift = curr[1] - pos
                return mapping[curr][0] + shift

    @classmethod
    def union(cls, *intervals):
        """
        returns the union of the set of input intervals
        """
        if len(intervals) < 1:
            raise AttributeError('cannot compute the union of an empty set of intervals')
        return Interval(min([i[0] for i in intervals]), max([i[1] for i in intervals]))

    @classmethod
    def intersection(cls, *intervals):
        """
        returns None if there is no intersection
        """
        if len(intervals) < 1:
            raise AttributeError('cannot compute the intersection of an empty set of intervals')
        low = max([i[0] for i in intervals])
        high = min([i[1] for i in intervals])
        if low > high:
            return None
        return Interval(low, high)

    @classmethod
    def min_nonoverlapping(cls, *intervals):
        """
        for a list of intervals, orders them and merges any overlap to return a list of non-overlapping intervals
        O(n^2)
        """
        g = nx.Graph()
        for i in intervals:
            g.add_node(i)

        for x, y in itertools.combinations(intervals, 2):
            if Interval.overlaps(x, y):
                g.add_edge(x, y)

        merges = []
        for c in nx.connected_components(g):
            if len(c) == 1:
                merges.append(c[0])
            else:
                merges.append(Interval.union(*c))
        return merges
