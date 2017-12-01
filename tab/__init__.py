"""
# About

The tab module is a python module used for easy transformation of tab files
It is fairly basic and does not support quoting text or escaping delimiters

## Order of Operations

There are a number of different transformations which can be applied to the rows. The order in which they are applied
is as follows

1. add
2. add_default
3. require
4. validate
5. rename
6. split
7. combine
8. cast
9. in_
10. drop
11. simplify

### add

this adds a new column with a default value. If the column name already exists in the input header. The existing column
is overwritten and replaced with the default value

### add_default

this adds a new column with a default value. If the column name already exists in the input header the existing value
is retained instead

### require

checks a list of column names to ensure they exist

### validate

checks a given column name against a regular expression

### rename

renames an input column to one or more new column names

### split

based on named capture groups in a regular expression, splits an existing column into one or more new columns

### combine

based on python template strings. Combines the values of 1 or more columns into a new column

### cast

applies any cast function to the column value

### in_

check the column value for membership of a specified object

### drop

deletes columns with a given name

### simplify

drops any input (not new) column names not specified in the require option


# Use-Cases

1. reading a tab file with no transformations

```
>>> header, rows = tab.read_file(filename, suppress_index=True)
```

2. reading a tab file and getting the line numbers

```
>>> header, rows = tab.read_file(filename)
>>> for row in rows:
>>>    print('row number:', row['_index'])
'row number:' 1
```

3. split a column with an expected pattern into multiple columns

```
>>> header, rows = tab.read_file(filename, split={'colname': r'^(?P<chr>\w+):(?P<pos>\d+)$'})
>>> print(header)
['colname', 'chr', 'pos']
```

4. drop specific unwanted columns

```
>>> header, rows = tab.read_file(filename, drop=['colname'])
```

5. drop all but specific columns

```
>>> header, rows = tab.read_file(filename, require=['colname'], simplify=True)
```

6. add a column with a default value

```
>>> header, rows = tab.read_file(filename, add={'colname': 'default value'})
```

7. combine columns into a new column

```
>>> header, rows = tab.read_file(filename, combine={'new_colname': '{colname1}_{colname2}'})
```

8. cast a column to a specific type

```
>>> header, rows = tab.read_file(filename, cast={'colname': int})
>>> header, rows = tab.read_file(filename, cast={'colname': tab.cast_boolean})
```
"""
from .tab import FileTransform, cast_boolean, cast_null, read_file, VERBOSE
