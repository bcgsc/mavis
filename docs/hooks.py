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
    types = {
        'string': 'str',
        'integer': 'int',
        'float': 'float',
        'boolean': 'bool',
        'number': 'float',
    }

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


def list_properties(schema, skip_terms=tuple()):
    glossary = {}
    for term, defn in schema['properties'].items():
        if term in skip_terms:
            continue
        typ = json_to_pytype(defn)
        desc = defn.get('description', '')
        default_value = defn.get('default')
        schema_fields = {k: v for k, v in defn.items() if k not in ['description', 'default']}

        if len(schema_fields) > 1:
            schema_defn = json.dumps(
                schema_fields,
                sort_keys=True,
                indent='    ',
            )
            schema_defn = f'**schema definition**:\n```json\n{schema_defn}\n```\n'
        else:
            schema_defn = ''

        lines = [
            f'### {term}',
            f'**type**: `#!python {typ}`',
            f'**default**: `#!python {repr(default_value)}`' if default_value is not None else '',
            desc,
            schema_defn,
        ]
        glossary[term] = '\n\n'.join(lines)
    return [v for k, v in sorted(glossary.items())]


def generate_settings_doc(schema_file):
    with open(schema_file, 'r') as fh:
        schema = json.load(fh)
    dirname = os.path.dirname(os.path.abspath(__file__))
    filepath = 'configuration/settings.md'
    title = 'Configurable Settings'

    fname = os.path.join(dirname, filepath)

    result = [f'\n\n# {title}\n']
    result.append(
        dedent(
            '''\
            ## Defining Samples/Libraries

            The `libraries` property of the mavis config is required to run the snakemake
            workflow. This is the section that defines what inputs to use, and what types of
            samples are available.

            ```json
            {
                "libraries": {
                    "<LIBRARY_NAME>": { }  // mapping of library name to library settings
                }
            }
            ```

            The library specific settings are listed below
            '''
        )
    )
    result.extend(list_properties(schema['properties']['libraries']['additionalProperties']))
    result.append(
        dedent(
            '''\
            ## Defining Conversions

            If the input to MAVIS is raw tool output and has not been pre-converted to the
            standard tab delimited format expected by MAVIS then you will need to add
            a section to the config to tell mavis how to perform the required conversions

            ```json
            {
                "convert": {
                    "<ALIAS>": { }  // mapping of alias to conversion settings
                }
            }
            ```

            The conversion specific settings are listed below
            '''
        )
    )
    result.extend(list_properties(schema['properties']['convert']['additionalProperties']))
    result.append('\n## General Settings\n')
    result.extend(list_properties(schema, ('libraries', 'convert')))

    print('writing:', fname)
    with open(fname, 'w') as fh:
        fh.write('\n\n'.join(result) + '\n')


def build_package_docs(config):
    schema_file = os.path.join(os.path.dirname(__file__), '../src/mavis/schemas/config.json')
    generate_settings_doc(schema_file)
    package_dir = os.path.join(os.path.dirname(__file__), '../src/mavis')
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
