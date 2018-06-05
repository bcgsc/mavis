"""
module for pre-processing code into rst files which are then built into html/latex/etc.
"""
import os
import sys
import subprocess
import glob
import re
from mavis.config import REFERENCE_DEFAULTS
from mavis.schedule.constants import OPTIONS as SUBMIT_OPTIONS
from mavis.summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from mavis.pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from mavis.validate.constants import DEFAULTS as VALIDATION_DEFAULTS
from mavis.annotate.constants import DEFAULTS as ANNOTATION_DEFAULTS
from mavis.cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from mavis.illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS
from mavis.util import ENV_VAR_PREFIX

dirname = os.path.dirname(os.path.abspath(__file__))

sys.path.insert(0, os.path.join(dirname, '..'))

# auto build the other documentation
subprocess.check_call('sphinx-apidoc -f -P -o {} {} --separate'.format(os.path.join(dirname, 'source', 'auto'), os.path.join(dirname, '..', 'mavis')), shell=True)

# now we need to add showing only select special members
for filename in glob.glob(os.path.join(dirname, 'source', 'auto', '*.rst')):
    # now open and read the file
    lines = []
    with open(filename, 'r') as fh:
        lines = fh.readlines()
    with open(filename, 'w') as fh:
        saw_automodule = False
        for line in lines:
            if re.match(r'^\.\.\s+automodule::\s+.*$', line):
                fh.write(line)
                fh.write('    :special-members: __and__, __or__, __xor__, __len__, __sub__, __add__\n')
            elif re.match(r'(\S+)\.(\S+)\s+(module|package)', line):
                m = re.match(r'(\S+)\.(\S+)\s+(module|package)', line)
                if m.group(1) == 'package':
                    line = re.sub(r'(\S+)\.(\S+)\s+(package)', r'\g<2> package', line)
                else:
                    line = re.sub(r'(\S+)\.(\S+)\s+(module)', r'\g<2> module', line)
                fh.write(line)
            else:
                fh.write(line)

fname = os.path.join(dirname, 'source', 'config_settings_glossary.rst')
print('writing:', fname)
with open(fname, 'w') as fh:
    fh.write('Configurable Settings\n')
    fh.write('+' * 50)
    tab = ' ' * 4
    fh.write('\n\n.. glossary::\n{}:sorted:\n\n'.format(tab))
    glossary = {}
    CUSTOM_TYPES = {
        'cast_boolean': 'bool',
        'float_fraction': '~mavis.constants.float_fraction',
        'ChrListString': '~mavis.util.ChrListString'
    }
    for namespace in [
        SUBMIT_OPTIONS,
        REFERENCE_DEFAULTS,
        SUMMARY_DEFAULTS,
        PAIRING_DEFAULTS,
        ANNOTATION_DEFAULTS,
        VALIDATION_DEFAULTS,
        CLUSTER_DEFAULTS,
        ILLUSTRATION_DEFAULTS
    ]:
        for term, value in namespace.items():
            typ = namespace.type(term).__name__
            typ = CUSTOM_TYPES.get(typ, typ)
            defn = ':class:`{}` - {}. The corresponding environment variable is ``{}{}`` and the default value is ``{}``'.format(
                typ, re.sub(r'\.?$', '', namespace.define(term, '')).capitalize(), ENV_VAR_PREFIX, term.upper(), repr(value))
            try:
                defn += '. Accepted values include: {}'.format(', '.join(['``{}``'.format(repr(v)) for v in namespace.type(term).values()]))
            except AttributeError:
                pass
            glossary[term] = defn
    for term, defn in sorted(glossary.items()):
        fh.write('\n{}{}\n'.format(tab, term))
        fh.write('{}{}{}\n'.format(tab, tab, defn))
