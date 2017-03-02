import os
import sys
import subprocess
import glob
import re

d = os.path.dirname(os.path.abspath(__file__))

sys.path.insert(0, os.path.join(d, '..'))

# auto build the other documentation
subprocess.check_call('sphinx-apidoc -f -P -o {} {} --separate'.format(
    os.path.join(d, 'source/auto'),
    os.path.join(d, './../mavis')), shell=True)

subprocess.check_call('sphinx-apidoc -f -P -o {} {}'.format(
    os.path.join(d, 'source/auto'),
    os.path.join(d, './../bin')), shell=True)

# now we need to add showing only select special members
for f in glob.glob(os.path.join(d, 'source/auto/*.rst')):
    # now open and read the file
    lines = []
    with open(f, 'r') as fh:
        lines = fh.readlines()
    with open(f, 'w') as fh:
        saw_automodule = False
        for line in lines:
            if re.match('^\.\.\s+automodule::\s+.*$', line):
                fh.write(line)
                fh.write('    :special-members: __and__, __or__, __xor__, __len__, __sub__, __add__\n')
            elif re.match('(\S+)\.(\S+)\s+(module|package)', line):
                m = re.match('(\S+)\.(\S+)\s+(module|package)', line)
                if m.group(1) == 'package':
                    line = re.sub('(\S+)\.(\S+)\s+(package)', '\g<2> package', line)
                else:
                    line = re.sub('(\S+)\.(\S+)\s+(module)', '\g<2> module', line)
                fh.write(line)
            else:
                fh.write(line)


# copy the README file to the source directory
# subprocess.check_call('cp {} {}'.format(os.path.join(d, './../README.rst'), os.path.join(d, 'source')), shell=True)
