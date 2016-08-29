from structural_variant.constants import *

class Interval:
    """
    Intervals are inclusive
    """
    def __init__(self, start, end = None, **kwargs):
        self.start = int( start )
        self.end = int( end ) if end is not None else self.start
        if self.start > self.end:
            raise AttributeError('interval start > end is not allowed')
        self.freq = int(kwargs.pop('freq', 1))
        if self.freq <= 0:
            raise AttributeError('Interval frequency must be a natural number')
    
    def overlap(self, other):
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
        return self.end - self.start + 1

    def __repr__(self):
        temp = str(self.start)
        if self.end != self.start:
            temp += '-' + str(self.end)
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
        return self.start + (len(self) - 1)/2
    
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

        return np.average(centers,weights=weights)

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
        return Interval(min(self.start, other.start), max(self.end, other.end))
    
    def __eq__(self, other):
        if not hasattr(other, 'start') \
                or not hasattr(other, 'end') \
                or not hasattr(other, 'weight') \
                or self.start != other.start \
                or self.end != other.end \
                or self.freq != other.freq:
            return False
        return True
    
    def __lt__(self, other):
        if self.start < other.start:
            return True
        elif self.start == other.start:
            if self.end < other.end:
                return True
            elif self.end == other.end:
                if self.freq < other.freq:
                    return True
        return False

    def __gt__(self, other):
        if self.start > other.start:
            return True
        elif self.start == other.start:
            if self.end > other.end:
                return True
            elif self.end == other.end:
                if self.freq > other.freq:
                    return True
        return False

    def __sub__(self, other):
        """ 
        returns the minimum distance between intervals
        
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
        if self.end < other.start:
            return self.end - other.start
        elif self.start > other.end:
            return self.start - other.end
        else:
            return 0
    
    def __hash__(self):
        return hash((self.start, self.end, self.freq))
    
    @classmethod
    def paired_weighted_means(cls, intervals):
        int_a = Interval.weighted_mean( [ x[0] for x in intervals ] )
        int_b = Interval.weighted_mean( [ x[1] for x in intervals ] )
        return int_a, int_b
    
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
    def redundant_ordered_hierarchical_clustering(cls, clusters, **kwargs):
        """
        for an input set of of clusters, do hierarchical clustering
        redundant b/c we allow clusters to be grouped more than once
        into either of their immediate neighbours
        """
        r = int(kwargs.pop('r'))
        if kwargs:
            raise AttributeError('invalid parameter', kwargs)
        if r < 0:
            raise AttributeError('r must be a positive integer')
        # order the clusters by weighted mean
        complete = []
        queue = sorted(clusters, key=lambda x: cls.paired_weighted_means(x) )
    
        
        while len(queue) > 0:
            temp_queue = []
            
            for i in range(0, len(queue)):
                curr = queue[i]
                joined = False
                
                if i > 0:
                    dist = cls.paired_set_distance(curr, clusters[i - 1])
                    if dist <= r:
                        joined = True
                if i < len(queue) - 1:
                    dist = cls.paired_set_distance(curr, clusters[i + 1])
                    if dist <= r:
                        temp_queue.append(curr.union(clusters[i + 1]))
                        joined = True 
                if not joined:
                    complete.append(curr)
            queue = temp_queue
        return complete

    @classmethod
    def union(cls, intervals):
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
        curr = next(iter(intervals))
        low = curr.start
        high = curr.end
        
        for i in intervals:
            if i.start < low:
                low = i.start
            if i.end > high:
                high = i.end
        return Interval(low, high)

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
        low = curr.start
        high = curr.end
        
        for i in intervals:
            if high < i.start or i.end < low:
                return None
            low = max(low, i.start)
            high = min(high, i.end)
        return Interval(low, high)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
