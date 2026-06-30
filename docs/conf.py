# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2017 Mikkel Schubert <mikkelsch@gmail.com>
#
# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os

# General information about the project.
project = "AdapterRemoval"
copyright = "%Y, Mikkel Schubert; Stinus Lindgreen"
author = "Mikkel Schubert; Stinus Lindgreen"

# The short X.Y version (`|version|`)
version = "3.0"
# The full version, including alpha/beta/rc tags (`|release|`)
release = "3.0.1"


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

source_suffix = ".rst"
master_doc = "index"
language = "en"

# Prevent conversion of -- to emdashes
smartquotes = False


# -- Themes ----------------------------------------------

if os.environ.get("READTHEDOCS") == "True":
    # https://sphinx-rtd-theme.readthedocs.io/
    html_theme = "sphinx_rtd_theme"
else:
    # http://alabaster.readthedocs.io/
    html_theme = "alabaster"

    # This is required for the alabaster theme
    html_sidebars = {
        "**": [
            "about.html",
            "searchfield.html",
            "navigation.html",
            "relations.html",
            "donate.html",
        ]
    }
