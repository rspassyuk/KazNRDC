"""Documentation builds without importing the GUI or nuclear databases."""
project = "NuMatRx"
author = "NuMatRx contributors"
extensions = ["myst_parser", "sphinx.ext.mathjax"]
myst_enable_extensions = ["dollarmath", "colon_fence"]
myst_heading_anchors = 3
html_theme = "sphinx_rtd_theme"
html_title = "NuMatRx documentation"
exclude_patterns = ["_build"]
source_suffix = {".md": "markdown"}
root_doc = "index"

html_logo = "../model/resources/branding/numatrx_logo_light.png"
