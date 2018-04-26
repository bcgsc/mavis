#!/usr/bin/env python3
"""
Order of transform operations

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
"""

from __future__ import division

import re
import string
import warnings


VERBOSE = False  # Output extra logging information (useful in debugging)


def cast_boolean(input_value):
    value = str(input_value).lower()
    if value in ['t', 'true', '1', 'y', 'yes', '+']:
        return True
    elif value in ['f', 'false', '0', 'n', 'no', '-']:
        return False
    raise TypeError('casting to boolean failed', input_value)


def cast_null(input_value):
    value = str(input_value).lower()
    if value in ['none', 'null']:
        return None
    raise TypeError('casting to null/None failed', input_value)


def null(input_value):
    warnings.warn('null is deprecated in favor of cast_null', DeprecationWarning, stacklevel=2)
    return cast_null(input_value)


class EmptyFileError(Exception):
    pass


class FileTransform:
    """
    Holds a set of operations which define the transform_line function.
    Generally a single FileTransform object is required per file as lines are expected to have the same format
    """
    def __init__(self, header, **kwargs):
        """
        Args:
            header (list of str): the header from the file as a list of column names (in-order)
            require (list of str): list of columns that must be in the input header
            rename (dict of str and list of str): mapping of old to new column(s)
            drop (list of str): list of columns in the old input header to drop
            add_default (dict of str): mapping of new column names to default values (if the column does not exist already)
            cast (dict of str and func): mapping of new/final columns to the type to cast them to
            split (dict of str and str):
                a dictionary mapping original column names to regex groups to create as the new column names
            combine (dict of str and str):
                a dictionary of the final column name to the format string. The field names in the format
                string must correspond to existing column names
            simplify (bool): drop all columns not created or retained
            validate (dict of str and str): mapping of old columns to regex they must satisfy

        Returns:
            FileTransform: an object with the validated rules for transforming lines in an input file
        """
        self.input = header[:]
        self.require = kwargs.pop('require', [])
        self.rename = kwargs.pop('rename', {})
        self.drop = kwargs.pop('drop', [])
        self.add = kwargs.pop('add', {})
        self.add_default = kwargs.pop('add_default', {})
        self.split = kwargs.pop('split', {})
        self.combine = kwargs.pop('combine', {})
        self.validate = kwargs.pop('validate', {})
        self.cast = kwargs.pop('cast', {})
        self.simplify = kwargs.pop('simplify', False)
        self.in_ = kwargs.pop('in_', {})
        self.header = []                               # holds the new header after the transform

        if kwargs:
            raise TypeError('invalid argument(s)', list(kwargs.keys()))

        current_columns = set(header)
        cant_simplify = set()  # columns that are restricted against being dropped in simplify

        if VERBOSE:
            print('input header:', header)

        # check that the header columns are unique
        if len(set(header)) != len(header):
            raise KeyError('duplicate input col: column names in input header must be unique', header)

        for col in self.add:
            current_columns.add(col)
            cant_simplify.add(col)
        # add_default: add_default new columns with default values if not already present
        for col in self.add_default:
            current_columns.add(col)
            cant_simplify.add(col)

        # 1. require: check that the required columns exist in the input header
        for col in self.require:
            if col not in current_columns:
                raise KeyError('cannot require: column not found in the input header', col, current_columns)
            cant_simplify.add(col)

        # 2. validate: check that the input column matches the expected pattern
        for col, regex in self.validate.items():
            if col not in current_columns:
                raise KeyError('cannot validate: column not found in the input header', col, current_columns)
            cant_simplify.add(col)

        # 4. rename: rename a column to one or more new column names
        for col, new_names in self.rename.items():
            if col not in current_columns:
                raise KeyError('cannot rename column. column not found in header', col, current_columns)
            for new_name in new_names:
                if new_name in current_columns:
                    raise KeyError('duplicate column name', new_name, current_columns)
                current_columns.add(new_name)
                cant_simplify.add(new_name)

        # 5. split: split a column into a set of new columns
        for col, regex in self.split.items():
            robj = re.compile(regex)
            new_columns = robj.groupindex.keys()
            if col not in current_columns:
                raise KeyError('cannot split column. column not found in header', col, current_columns)
            for new_col in new_columns:
                if new_col in current_columns:
                    raise KeyError('duplicate column name', new_col, current_columns)
                current_columns.add(new_col)
                cant_simplify.add(new_col)

        # 6. combine:
        for ncol, format_string in self.combine.items():
            old_column_names = [t[1] for t in list(string.Formatter().parse(format_string))]
            if ncol in current_columns:
                raise KeyError('duplicate column name', ncol, current_columns)
            current_columns.add(ncol)
            cant_simplify.add(ncol)
            for col in old_column_names:
                if col not in current_columns:
                    raise KeyError('cannot combine column. column not found in header', col, current_columns)

        # 7. cast: apply some callable
        for col, func in self.cast.items():
            if col not in current_columns:
                raise KeyError('cannot cast column. column not found in header', col, current_columns)
            if not callable(func):
                raise TypeError('function applied to column must be callable', col, func)
            cant_simplify.add(col)

        # 8. in_: check for satisfying some controlled vocab
        for col, item in self.in_.items():
            if col not in current_columns:
                raise KeyError('cannot check membership column. column not found in header', col, current_columns)
            if None in item:
                pass
            cant_simplify.add(col)

        # 9. drop: drop any columns from the original input IF EXIST
        for col in self.drop:
            if col in self.require:
                raise AssertionError('cannot both drop and retain a column', col)
            current_columns.discard(col)        # 8. simplify: drop any columns that are not new, added, or retained

        if self.simplify:
            for col in list(current_columns):
                if col not in cant_simplify:
                    current_columns.discard(col)

        # retain the original input order except for new columns
        order = {}
        for col in current_columns:
            order[col] = len(header)
        for i, col in enumerate(header):
            if col in current_columns:
                order[col] = i

        self.header = [m for m, n in sorted(order.items(), key=lambda x: (x[1], x[0]))]

        if VERBOSE:
            print('output header:', self.header)

    def transform_line(self, line, allow_short=False):
        """
        transforms the input line into a hash of the new/final column names with the transform rules applied

        Args:
            line (list of str): list of values for a row with the same input header as the transform
        Raises:
            exception exceptions occur if validation, split or combine fails

        Returns:
            dict of str: the hash representation of the new row
        """
        if any([
            not allow_short and len(self.input) != len(line),
            allow_short and len(self.input) < len(line)
        ]):
            raise AssertionError('length of input list {0} does not match length of the expected header {1}: '.format(
                len(line), len(self.input)) + re.sub('\n', '\\n', '\\t'.join(line)), self.input)

        row = {}
        cant_simplify = set()

        for i in range(0, len(self.input)):
            row[self.input[i]] = line[i] if i < len(line) else None

        for col, default in self.add.items():
            row[col] = default
            cant_simplify.add(col)

        # add_default: add new columns with default values if not already present
        for col, default in self.add_default.items():
            row.setdefault(col, default)
            cant_simplify.add(col)

        # 1. require: check that the required columns exist in the input header
        cant_simplify.update(self.require)

        # 2. validate: check that the input column matches the expected pattern
        for col, regex in self.validate.items():
            cant_simplify.add(col)
            if not re.match(regex, row[col]):
                raise UserWarning('validation failed', col, regex, row[col])

        # 4. rename: rename a column to one or more new column names
        for col, new_names in self.rename.items():
            for new_name in new_names:
                row[new_name] = row[col]
                cant_simplify.add(new_name)

        # 5. split: split a column into a set of new columns
        for col, regex in self.split.items():
            robj = re.compile(regex)
            new_columns = robj.groupindex.keys()
            match = robj.match(row[col])
            if not match:
                raise UserWarning('split of column failed', col, regex, row[col])
            for new_col in new_columns:
                row[new_col] = match.group(new_col)
                cant_simplify.add(new_col)

        # 6. combine:
        for ncol, format_string in self.combine.items():
            old_column_names = [t[1] for t in list(string.Formatter().parse(format_string))]
            cant_simplify.add(ncol)
            substitutions = {}
            for col in old_column_names:
                substitutions[col] = row[col]
            row[ncol] = format_string.format(**substitutions)

        # 7. cast: apply some callable
        for col, func in self.cast.items():
            try:
                row[col] = func(row[col])
            except Exception as err:
                raise type(err)('error in casting column: {}. {}'.format(col, str(err)))
            cant_simplify.add(col)

        # 8. in_: check for satisfying some controlled vocab
        for col, item in self.in_.items():
            if row[col] not in item:
                raise KeyError('failed in_ check', col, row[col], item)
            cant_simplify.add(col)

        # 9. drop: drop any columns from the original input IF EXIST
        for col in self.drop:
            row.pop(col, None)

        # 10. simplify: drop any columns that are not new, added, or retained
        if self.simplify:
            for col in list(row):
                if col not in cant_simplify:
                    row.pop(col, None)

        return row


def read_file(inputfile, delimiter='\t', header=None, strict=True, suppress_index=False, allow_short=False, **kwargs):
    """
    Args:
        inputfile (str): the path to the inputfile
        header (list of str): for non-headered files
        delimiter (str): the delimiter (what to split on)
        strict (bool): if false will ignore lines that fail transform
        suppress_index (bool): do not create an index
    Returns:
        list of str and dict of str: header and the row dictionaries
    """
    if VERBOSE:
        print("read_file(", inputfile, ", ", kwargs, ")")

    new_header = None
    is_file_handle = True if hasattr(inputfile, 'readlines') else False
    index = '_index'
    objects = []
    line_count = 0

    fh = inputfile if is_file_handle else open(inputfile, 'r')

    # first grab the header and skip comments
    lines = fh.readlines()
    if not lines:
        raise EmptyFileError('empty file has no lines to read')
    current_line_index = 0
    line = re.sub(r'[\r\n]*$', '', lines[current_line_index])
    while current_line_index < len(lines):
        if not re.match(r'^\s*##', lines[current_line_index]):  # skip comment lines
            break
        current_line_index += 1

    # first line is the header unless a header was input
    if not header:
        if current_line_index >= len(lines):
            raise EmptyFileError('no lines beyond comments to read as header')
        line = re.sub(r'(^#)|([\r\n\s]*$)', '', lines[current_line_index])  # clean the header
        current_line_index += 1
        header = line.split(delimiter) if line else []
    if not header:
        raise EmptyFileError('header is empty', inputfile)
    # create the file transform object
    transform = FileTransform(header, **kwargs)
    new_header = transform.header

    if not suppress_index and index in new_header:
        raise AttributeError('column name {0} is reserved and cannot be used as an input'.format(repr(index)))

    # now go through the lines in the file
    while current_line_index < len(lines):
        line_count += 1
        line = re.sub(r'[\r\n]*$', '', lines[current_line_index])  # clean the line
        try:
            row = line.split(delimiter)
            row = transform.transform_line(row, allow_short=allow_short)
            if not suppress_index:
                row[index] = current_line_index
            objects.append(row)
        except Exception as error:  # General b/c will be re-raised unless strict mode is off
            if strict:
                print('error at line', current_line_index)
                raise type(error)('{0} happens at line {1}'.format(error, current_line_index))
            elif VERBOSE:
                print('[ERROR]', str(error))
        current_line_index += 1

    if not is_file_handle:
        fh.close()
    return (new_header, objects)
