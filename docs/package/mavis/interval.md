# mavis.interval

## class Interval







### Interval.overlaps()

checks if two intervals have any portion of their given ranges in common

```python
@classmethod
def overlaps(cls, first, other):
```

**Args**

- first (`Interval`): an interval to be compared
- other (`Interval`): an interval to be compared

**Examples**

```python
Interval.overlaps(Interval(1, 4), Interval(5, 7))
False
Interval.overlaps(Interval(1, 10), Interval(10, 11))
True
Interval.overlaps((1, 10), (10, 11))
True
```







### Interval.center()

the middle of the interval

```python
@property
def center(self):
```

**Args**

- self

**Examples**

```python
Interval(1, 10).center
5.5
Interval(1, 11).center
6
```




### Interval.dist()

returns the minimum distance between intervals

```python
@classmethod
def dist(cls, first, other):
```

**Args**

- first
- other

**Examples**

```python
Interval.dist((1, 4), (5, 7))
-1
Interval.dist((5, 7), (1, 4))
1
Interval.dist((5, 8), (7, 9))
0
```





### Interval.convert\_ratioed\_pos()

convert any given position given a mapping of intervals to another range

```python
@classmethod
def convert_ratioed_pos(cls, mapping, pos, forward_to_reverse=None):
```

**Args**

- mapping (`:class:`dict` of :class:`Interval` and :class:`Interval``): a mapping of a set of continuous intervals
- pos (`int`): a position in the first coordinate system
- forward_to_reverse

**Returns**

- `Interval`: the position in the alternate coordinate system given the input mapping

**Raises**

- `AttributeError`: if the input position is outside the set of input segments
- `IndexError`: if the input position cannot be converted to the output system

**Examples**

```python
mapping = {(1, 10): (101, 110), (11, 20): (555, 564)}
Interval.convert_pos(mapping, 5)
5
Interval.convert_pos(mapping, 15)
559
```


### Interval.union()

returns the union of the set of input intervals

```python
@classmethod
def union(cls, *intervals):
```

**Examples**

```python
Interval.union((1, 2), (4, 6), (4, 9), (20, 21))
Interval(1, 21)
```


### Interval.intersection()

returns None if there is no intersection

```python
@classmethod
def intersection(cls, *intervals):
```

**Examples**

```python
Interval.intersection((1, 10), (2, 8), (7, 15))
Interval(7, 8)
Interval.intersection((1, 2), (5, 9))
None
```


### Interval.min\_nonoverlapping()

for a list of intervals, orders them and merges any overlap to return a list of non-overlapping intervals
O(nlogn)

```python
@classmethod
def min_nonoverlapping(cls, *intervals):
```

**Examples**

```python
Interval.min_nonoverlapping((1, 10), (7, 8), (6, 14), (17, 20))
[Interval(1, 14), Interval(17, 20)]
```



### Interval.split\_overlap()

for a given set of intervals splits any overlap so that the result is a new list of
intervals with a single bp overlap

Warning:
this method is meant for integer intervals only

```python
@classmethod
def split_overlap(cls, *intervals, weight_mapping={}):
```


## class IntervalMapping

mapping between coordinate systems using intervals.
source intervals cannot overlap but no such assertion is enforced on the target intervals

### IntervalMapping.\_\_init\_\_()

```python
def __init__(self, mapping=None, opposing=None):
```

**Args**

- mapping
- opposing




### IntervalMapping.convert\_ratioed\_pos()

convert any given position given a mapping of intervals to another range

```python
def convert_ratioed_pos(self, pos):
```

**Args**

- pos (`Interval`): a position in the first coordinate system

**Returns**

: the position in the alternate coordinate system given the input mapping - int: if simplify is True - Interval: if simplify is False

**Raises**

- `IndexError`: if the input position is not in any of the mapped intervals

**Examples**

```python
mapping = IntervalMapping(mapping={(1, 10): (101, 110), (11, 20): (555, 564)})
mapping.convert_pos(5)
5
mapping.convert_pos(15)
559
```


### IntervalMapping.convert\_pos()

convert any given position given a mapping of intervals to another range

```python
def convert_pos(self, pos):
```

**Args**

- pos (`int`): a position in the first coordinate system

**Returns**

: the position in the alternate coordinate system given the input mapping - int: if simplify is True - Interval: if simplify is False

**Raises**

- `IndexError`: if the input position is not in any of the mapped intervals

**Examples**

```python
mapping = IntervalMapping(mapping={(1, 10): (101, 110), (11, 20): (555, 564)})
mapping.convert_pos(5)
5
mapping.convert_pos(15)
559
```
