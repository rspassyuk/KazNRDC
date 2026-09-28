"""Documentation builds without importing the GUI or nuclear databases."""
project = "KazNRDC"
author = "KazNRDC contributors"
extensions = ["myst_parser", "sphinx.ext.mathjax"]
myst_enable_extensions = ["dollarmath", "colon_fence"]
myst_heading_anchors = 3
html_theme = "sphinx_rtd_theme"
html_title = "KazNRDC documentation"
exclude_patterns = ["_build"]
source_suffix = {".md": "markdown"}
root_doc = "index"
