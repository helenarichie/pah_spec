# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import pah_spec

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "pah_spec"
copyright = "2026, Helena Richie"
author = "Helena Richie"
release = pah_spec.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx_inline_tabs",
]
source_suffix = [".rst", ".md"]

napoleon_numpy_docstring = True
html_theme = "sphinx_rtd_theme"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
automodule_members = True

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "module-doc-separator": "",
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_static_path = ["_static"]

# -- Options for intersphinx -------------------------------------------------
intersphinx_mapping = {
    "CPython": ("https://docs.python.org/3", None),
    "astropy": ("http://docs.astropy.org/en/stable", None),
    "pooch": ("https://www.fatiando.org/pooch/latest/", None),
}

# Recommended by readthedocsadding the following config value.
# Sphinx defaults to automatically resolve *unresolved* labels using Intersphinx
# mappings. This behavior has unintended side-effects, namely that documentations local
# references can suddenly resolve to an external location.
# See also:
# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html#confval-intersphinx_disabled_reftypes
intersphinx_disabled_reftypes = ["*"]

# -- Options for MyST --------------------------------------------------------

myst_enable_extensions = ["colon_fence", "fieldlist"]
