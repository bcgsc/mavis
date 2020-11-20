import os

from markdown_refdocs.main import extract_to_markdown


def build_package_docs(config):
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
