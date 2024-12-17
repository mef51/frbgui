# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'FRBGui'
copyright = '2023, Mohammed A. Chamma'
author = 'Mohammed A. Chamma'
version = '0.10.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx_copybutton',
    'sphinxnotes.strike',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
suppress_warnings = ['epub.unknown_project_files']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_favicon = '../frbgui.ico'
html_logo = 'imgs/frbgui.png'
html_context = {
    "display_github": True, # Integrate GitHub
    "github_user": "mef51", # Username
    "github_repo": "frbgui", # Repo name
    "github_version": "main", # Version
    "conf_py_path": "/docs/", # Path in the checkout to the docs root
}

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'sticky_navigation': True
}
html_static_path = ['_static']
html_extra_path = ['tutorials/files']
html_css_files = ['style.css']

pygments_style = 'sphinx'
