import sys
import os
import subprocess
import glob
import re

import mavis
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
fname = os.path.join(dirname, 'glossary.md')
print('writing:', fname)
with open(fname, 'a') as fh:
    fh.write('\n\n## Configurable Settings\n')
    glossary = {}
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
            # typ = CUSTOM_TYPES.get(typ, typ)
            defn = '`{}` - {}. The corresponding environment variable is `{}{}` and the default value is `{}`'.format(
                typ, re.sub(r'\.?$', '', namespace.define(term, '')).capitalize(), ENV_VAR_PREFIX, term.upper(), repr(value))
            try:
                defn += '. Accepted values include: {}'.format(', '.join(['`{}`'.format(repr(v)) for v in namespace.type(term).values()]))
            except AttributeError:
                pass
            glossary[term] = defn
    for term, defn in sorted(glossary.items()):
        fh.write(f'### {term}\n\n')
        fh.write(f'{defn}\n\n')
