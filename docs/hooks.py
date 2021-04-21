import json
import os
import re
from textwrap import dedent

from markdown_refdocs.main import extract_to_markdown
from mavis.schemas import DEFAULTS
from mavis.util import ENV_VAR_PREFIX


def json_to_pytype(record):
    input_type = record
    try:
        input_type = record['type']
    except TypeError:
        pass
    types = {'string': 'str', 'integer': 'int', 'float': 'float'}

    if input_type == 'array':
        try:
            sub_type = json_to_pytype(record['items']['type'])
            return f'List[{sub_type}]'
        except TypeError:
            return 'List'

    if isinstance(input_type, list):
        # Union
        types = ', '.join([json_to_pytype(t) for t in input_type])
        return f'Union[{types}]'
    return types.get(input_type, input_type)


def generate_settings_doc(schema_file):
    with open(schema_file, 'r') as fh:
        schema = json.load(fh)
    dirname = os.path.dirname(os.path.abspath(__file__))
    filepath = 'configuration/settings.md'
    title = 'Configurable Settings'

    fname = os.path.join(dirname, filepath)
    print('writing:', fname)
    with open(fname, 'w') as fh:
        fh.write(f'\n\n# {title}\n')
        glossary = {}
        for term, defn in schema['properties'].items():
            if term in ['libraries', 'convert']:
                continue
            typ = json_to_pytype(defn)
            desc = defn.get('description', '')
            default_value = defn.get('default')
            schema_defn = json.dumps(
                {k: v for k, v in defn.items() if k not in ['description', 'default']},
                sort_keys=True,
                indent='    ',
            )
            schema_defn = f'**schema definition**:\n```json\n{schema_defn}\n```\n'

            lines = [
                f'## {term}',
                f'**type**: `#!python {typ}`',
                f'**default**: `#!python {repr(default_value)}`',
                desc,
                schema_defn,
            ]
            glossary[term] = '\n\n'.join(lines)
        for term, defn in sorted(glossary.items()):
            fh.write(f'{defn}\n\n')


def build_package_docs(config):
    schema_file = os.path.join(os.path.dirname(__file__), '../mavis/schemas/config.json')
    generate_settings_doc(schema_file)
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
