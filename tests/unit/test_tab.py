import unittest
from tab import FileTransform, cast_boolean, cast_null


class MockFileTransform:
    def __init__(self, h, **kwargs):
        self.input = h
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

    def transform_line(self, *pos, **kwargs):
        return FileTransform.transform_line(self, *pos, **kwargs)


class TestCast(unittest.TestCase):

    def test_cast_boolean_true(self):
        self.assertEqual(True, cast_boolean('+'))
        self.assertEqual(True, cast_boolean('T'))
        self.assertEqual(True, cast_boolean('true'))
        self.assertEqual(True, cast_boolean('y'))
        self.assertEqual(True, cast_boolean(1))

    def test_cast_boolean_false(self):
        self.assertEqual(False, cast_boolean('-'))
        self.assertEqual(False, cast_boolean('f'))
        self.assertEqual(False, cast_boolean('false'))
        self.assertEqual(False, cast_boolean('n'))
        self.assertEqual(False, cast_boolean(0))

    def test_cast_boolean_error(self):
        with self.assertRaises(TypeError):
            cast_boolean(2)

    def test_cast_null_ok(self):
        self.assertEqual(None, cast_null('none'))
        self.assertEqual(None, cast_null(None))

    def test_cast_null_error(self):
        with self.assertRaises(TypeError):
            cast_null('f')


class TestFileTransform(unittest.TestCase):

    def test_simplify(self):
        h = ['a', 'b', 'c']
        ft = FileTransform(header=h)
        self.assertEqual(h, ft.input)
        self.assertEqual(h, ft.header)
        ft = FileTransform(h, simplify=True)
        self.assertEqual(h, ft.input)
        self.assertEqual([], ft.header)

    def test_require_simplify(self):
        h = ['a', 'b', 'c']
        ft = FileTransform(header=h, require=['a'], simplify=True)
        self.assertEqual(h, ft.input)
        self.assertEqual(['a'], ft.header)

    def test_require_error(self):
        h = ['a', 'b', 'c']
        with self.assertRaises(KeyError):
            FileTransform(header=h, require=['k'])

    def test_rename(self):
        h = ['a', 'b', 'c']
        ft = FileTransform(h, rename={'a': ['k', 'm']})
        self.assertEqual(h, ft.input)
        self.assertEqual(['a', 'b', 'c', 'k', 'm'], ft.header)

    def test_rename_error(self):
        h = ['a', 'b', 'c']
        with self.assertRaises(KeyError):
            FileTransform(header=h, rename={'k': ['t']})

    def test_cast_error(self):
        h = ['b', 'c']
        with self.assertRaises(KeyError):
            FileTransform(h, cast={'a': int})

    def test_add(self):
        h = ['a', 'b', 'c']
        ft = FileTransform(h, add_default={'k': 1})
        self.assertEqual(h, ft.input)
        self.assertEqual(['a', 'b', 'c', 'k'], ft.header)

    def test_require__in(self):
        h = ['a', 'b', 'c']
        ft = FileTransform(h, require=['c'], in_={'a': []}, simplify=True)
        self.assertEqual(h, ft.input)
        self.assertEqual(['a', 'c'], ft.header)

    def test_combine(self):
        h = ['a', 'b', 'c']
        ft = FileTransform(h, combine={'k': '{a}{b}{c}'})
        self.assertEqual(ft.input, h)
        self.assertEqual(ft.header, h + ['k'])

    def test_combine_error_name_conflict(self):
        h = ['a', 'b', 'c']
        with self.assertRaises(KeyError):
            FileTransform(h, combine={'b': '{a}{b}{c}'})

    def test_combine_keyerror(self):
        h = ['a', 'b', 'c']
        with self.assertRaises(KeyError):
            FileTransform(h, combine={'k': '{m}{b}{c}'})

    def test_duplicate_input_column(self):
        with self.assertRaises(KeyError):
            FileTransform(['a', 'a'])

    def test_validate_missing_column(self):
        with self.assertRaises(KeyError):
            FileTransform(['a', 'b'], validate={'c': ''})

    def test_drop_and_require_error(self):
        with self.assertRaises(AssertionError):
            FileTransform(['a'], require=['a'], drop=['a'])

    def test_membership_of_missing_column_error(self):
        with self.assertRaises(KeyError):
            FileTransform(['a'], in_={'x': []})

    def test_membership_bad_object(self):
        with self.assertRaises(TypeError):
            FileTransform(['a'], in_={'a': 1})

    def test_cast_noncallable_error(self):
        FileTransform(['a'], cast={'a': int})
        with self.assertRaises(TypeError):
            FileTransform(['a'], cast={'a': 1})

    def test_split_missing_column_error(self):
        FileTransform(['a'], split={'a': r'^(?P<thing>\w+)'})
        with self.assertRaises(KeyError):
            FileTransform(['a'], split={'x': r'^(?P<thing>\w+)'})

    def test_split_duplicate_column_error(self):
        FileTransform(['a', 'b'], split={'a': r'^(?P<thing>\w+)'})
        with self.assertRaises(KeyError):
            FileTransform(['a', 'b'], require=['b'], split={'a': r'^(?P<b>\w+)'})

    def test_add(self):
        ft = FileTransform(['a', 'b'], add={'c': 1})
        self.assertEqual(['a', 'b', 'c'], ft.header)

    def test_add_default(self):
        ft = FileTransform(['a', 'b'], add_default={'c': 1})
        self.assertEqual(['a', 'b', 'c'], ft.header)
        ft = FileTransform(['a', 'b'], add_default={'b': 1})
        self.assertEqual(['a', 'b'], ft.header)

    def test_require(self):
        ft = FileTransform(['a', 'b'], require=['a'])
        self.assertEqual(['a', 'b'], ft.header)

        ft = FileTransform(['a', 'b'], require=['a'], simplify=True)
        self.assertEqual(['a'], ft.header)

    def test_invalid_option(self):
        with self.assertRaises(TypeError):
            FileTransform(['a', 'b'], require=['a'], blargh=1)


class TestTransformLine(unittest.TestCase):

    def test_add(self):
        h = ['a', 'b', 'c']
        ft = MockFileTransform(h, add={'a': 'blargh'})
        row = ft.transform_line(['1', '2', '3'])
        self.assertEqual('blargh', row['a'])

        h = ['a', 'b', 'c']
        ft = MockFileTransform(h, add={'x': 'blargh'})
        row = ft.transform_line(['1', '2', '3'])
        self.assertEqual('blargh', row['x'])

    def test_combine(self):
        h = ['a', 'b', 'c']
        ft = MockFileTransform(h, combine={'k': '{a}{b}{c}'})
        row = ft.transform_line(['1', '2', '3'])
        self.assertEqual('123', row['k'])

    def test_combine_then_cast(self):
        h = ['a', 'b', 'c']
        ft = MockFileTransform(h, combine={'k': '{a}{b}{c}'}, cast={'k': int})
        row = ft.transform_line(['1', '2', '3'])
        self.assertEqual(123, row['k'])

    def test_cast_to_cast_boolean(self):
        h = ['a', 'b', 'c']
        ft = MockFileTransform(h, cast={'a': cast_boolean, 'b': cast_boolean})
        row = ft.transform_line(['1', '0', '3'])
        self.assertEqual(True, row['a'])
        self.assertEqual(False, row['b'])

    def test_split_combine_cast(self):
        h = ['a', 'b', 'c']
        ft = MockFileTransform(h, split={'a': r'^(?P<a1>\d+)_(?P<a2>\d+)$'}, combine={'k': '{a1}{b}{c}'}, cast={'k': int})
        row = ft.transform_line(['1_10', '2', '3'])
        self.assertEqual(123, row['k'])

    def test_add_default(self):
        h = ['a', 'b', 'c']
        ft = MockFileTransform(h, add_default={'k': 1})
        line = ['1', '2', '3']
        row = ft.transform_line(line)
        self.assertEqual(1, row['k'])
        self.assertEqual('1', row['a'])
        self.assertEqual('2', row['b'])
        self.assertEqual('3', row['c'])

    def test_add_default_override_default(self):
        h = ['a', 'b', 'c']
        ft = MockFileTransform(h, add_default={'a': 8}, in_={'a': ['1']})
        line = ['1', '2', '3']
        row = ft.transform_line(line)
        self.assertEqual('1', row['a'])
        self.assertEqual('2', row['b'])
        self.assertEqual('3', row['c'])

    def test_validate(self):
        h = ['a']
        ft = MockFileTransform(h, validate={'a': r'^[t]+$'})
        line = ['ttttt']
        row = ft.transform_line(line)
        self.assertEqual('ttttt', row['a'])

    def test_rename(self):
        h = ['a']
        ft = MockFileTransform(h, rename={'a': ['b', 'c']})
        line = ['ttttt']
        row = ft.transform_line(line)
        self.assertEqual('ttttt', row['a'])
        self.assertEqual('ttttt', row['b'])
        self.assertEqual('ttttt', row['c'])

    def test_length_mismatch_error(self):
        h = ['a', 'b']
        ft = MockFileTransform(h)
        line = ['ttttt']
        with self.assertRaises(AssertionError):
            ft.transform_line(line)

    def test_rename_drop_original(self):
        h = ['a']
        ft = MockFileTransform(h, rename={'a': ['b', 'c']}, drop=['a'])
        line = ['ttttt']
        row = ft.transform_line(line)
        self.assertTrue('a' not in row)
        self.assertEqual('ttttt', row['b'])
        self.assertEqual('ttttt', row['c'])

        ft = MockFileTransform(h, rename={'a': ['b', 'c']}, simplify=True)
        row = ft.transform_line(line)
        self.assertTrue('a' not in row)
        self.assertEqual('ttttt', row['b'])
        self.assertEqual('ttttt', row['c'])

    def test_split(self):
        h = ['a']
        ft = MockFileTransform(
            h,
            split={'a': r'^(?P<a1>\d+)[_]+(?P<a2>\d+)$'}
        )
        row = ft.transform_line(['1_2'])
        self.assertEqual('1', row['a1'])
        self.assertEqual('2', row['a2'])
        row = ft.transform_line(['1__2'])
        self.assertEqual('1', row['a1'])
        self.assertEqual('2', row['a2'])
        with self.assertRaises(UserWarning):
            ft.transform_line(['_1__4'])


if __name__ == '__main__':
    unittest.main()
