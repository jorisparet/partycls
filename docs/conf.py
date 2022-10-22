# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'partycls'
copyright = '2022, Joris Paret'
author = 'Joris Paret, Daniele Coslovich'

# The full version, including alpha/beta/rc tags
with open('../partycls/core/_version.py') as f:
    exec(f.read())
release = __version__
version = release


# -- General configuration ---------------------------------------------------

def setup(app):
    app.add_css_file('custom.css')

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
	'sphinx.ext.autodoc',
	'sphinx.ext.coverage',
	'sphinx.ext.napoleon',
	'sphinx.ext.viewcode',
	'sphinx.ext.intersphinx',
	'sphinx.ext.mathjax',
	'sphinx_copybutton',
	'sphinx_rtd_theme',
	'sphinxcontrib.bibtex',
	'nbsphinx'
]

# include __init__ methods
autoclass_content = 'both'

# bibtex files
bibtex_bibfiles = [
	'joss.bib',
	'references.bib'
]

# latex directives
latex_elements = {
    'inputenc': '',
    'utf8extra': ''
}


# Add any paths that contain templates here, relative to this directory.
#templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
	'logo_only': True,
	'display_version': True,
	'vcs_pageview_mode': 'edit'
}

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#
html_logo = "./_static/img/logo.svg"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
	'python': ('https://docs.python.org/3', None),
	'numpy': ('https://numpy.org/doc/stable/', None),
	'matplotlib' : ('https://matplotlib.org/stable/', None),
	'scikit-learn': ('https://scikit-learn.org/stable/', None)
}

# nbsphinx options (for notebooks)
# This is processed by Jinja2 and inserted before each notebook
nbsphinx_prolog = r"""
{% set docname = env.doc2path(env.docname, base='doc') %}

.. tip::
	
	Run this notebook online: |binder|

	.. |binder| image:: https://mybinder.org/badge_logo.svg
			:target: https://mybinder.org/v2/gh/jorisparet/partycls/{{ env.config.release }}?filepath={{ docname }}
			:alt: Binder badge
"""
