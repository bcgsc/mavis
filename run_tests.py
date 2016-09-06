import sys
import os
import unittest
import pkgutil

basedir = os.path.abspath( os.path.realpath( os.path.dirname( __file__ )))

import tests

pkgpath = os.path.dirname(tests.__file__)
testmodules = [ tests.__name__ + '.' + name for _, name, _ in pkgutil.iter_modules([pkgpath])]

suite = unittest.TestSuite()

for t in testmodules:
    new_tests = unittest.defaultTestLoader.loadTestsFromName(t)
    suite.addTest(new_tests)

unittest.TextTestRunner(verbosity=1).run(suite)
