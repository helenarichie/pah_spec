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

extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon", "sphinx.ext.viewcode", "myst_parser"]
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

#html_static_path = ["_static"]
