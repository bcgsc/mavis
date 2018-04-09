import os
import pycodestyle
import unittest


BASEPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


class TestCodeFormat(unittest.TestCase):

    def remove_expected_error(self, result):
        if result.messages.get('E741', '') == 'ambiguous variable name \'I\'':
            del result.messages['E741']
            result.counters['E741'] -= 1
            result.file_errors -= 1
            result.total_errors -= 1

    def test_conformance(self):
        """Test that we conform to PEP-8."""
        style = pycodestyle.StyleGuide(config_file=os.path.join(BASEPATH, 'setup.cfg'))
        result = style.check_files([os.path.join(BASEPATH, 'mavis'), os.path.join(BASEPATH, 'tests')])
        self.remove_expected_error(result)
        self.assertEqual(result.total_errors, 0, "Found code style errors (and warnings).")
