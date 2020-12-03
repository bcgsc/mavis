import os
import re

from markdown_refdocs.main import extract_to_markdown
from mavis.annotate.constants import DEFAULTS as ANNOTATION_DEFAULTS
from mavis.cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from mavis.config import REFERENCE_DEFAULTS
from mavis.illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS
from mavis.pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from mavis.schedule.constants import OPTIONS as SUBMIT_OPTIONS
from mavis.summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from mavis.util import ENV_VAR_PREFIX
from mavis.validate.constants import DEFAULTS as VALIDATION_DEFAULTS


def generate_settings_doc():
    dirname = os.path.dirname(os.path.abspath(__file__))

    for (filepath, title, namespaces) in [
        (
            'configuration/settings.md',
            'Configurable Settings',
            [
                SUBMIT_OPTIONS,
                REFERENCE_DEFAULTS,
                SUMMARY_DEFAULTS,
                PAIRING_DEFAULTS,
                ANNOTATION_DEFAULTS,
                VALIDATION_DEFAULTS,
                CLUSTER_DEFAULTS,
                ILLUSTRATION_DEFAULTS,
            ],
        ),
    ]:
        fname = os.path.join(dirname, filepath)
        print('writing:', fname)
        with open(fname, 'w') as fh:
            fh.write(f'\n\n# {title}\n')
            glossary = {}
            for namespace in namespaces:
                for term, value in namespace.items():
                    typ = namespace.type(term).__name__
                    # typ = CUSTOM_TYPES.get(typ, typ)
                    desc = re.sub(r"\.?$", "", namespace.define(term, "")).capitalize()
                    accepted = ''
                    try:
                        accepted = '\n\n**accepted values**: {}\n'.format(
                            ', '.join(['`{}`'.format(repr(v)) for v in namespace.type(term).values()])
                        )
                    except AttributeError:
                        pass
                    defn = f'''## {term}

**type**: `#!python {typ}`

**environment variable**: `{ENV_VAR_PREFIX}{term.upper()}`

**default**: `#!python {repr(value)}`{accepted}

{desc}
        '''
                    glossary[term] = defn
            for term, defn in sorted(glossary.items()):
                fh.write(f'{defn}\n\n')


def build_package_docs(config):
    generate_settings_doc()
    package_dir = os.path.join(os.path.dirname(__file__), '../mavis')
    output_dir = os.path.join(os.path.dirname(__file__), 'package')
    extract_to_markdown(
        [package_dir],
        output_dir,
        link=True,
        hide_private=True,
        hide_undoc=True,
        hide_undoc_args=True,
        namespace_headers=True,
    )
