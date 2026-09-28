# Publishing source and documentation

## Repository contents

Keep source, tests, application assets, documentation and small intentional examples in Git. `.gitignore` excludes external `xsdir/`, Python caches, virtual environments, `.kaz` sessions and generated documentation. Ignoring a path does not untrack a previously committed file or remove it from history.

Keep a separate, versioned distribution of nuclear data with source attribution, access instructions and checksums. The application currently has no automatic database downloader. Do not advertise a self-contained installation without supplying those instructions.

## Validate locally

Run from the repository root:

```bash
python -B -m unittest discover -s tests -p "test_*.py"
python -m pip install -r docs/requirements.txt
python -m sphinx -W --keep-going -b html docs docs/_build/html
```

Open `docs/_build/html/index.html`. Review the installation steps against a fresh environment and the data you plan to distribute. Unit tests do not establish agreement with observed stellar abundances.

## GitHub

Review Source Control in VS Code before committing. Exclude datasets, credentials and personal sessions. Commit source and documentation changes, then push to the intended repository. The documentation workflow builds HTML and uploads it as an Actions artifact; it does not publish a website by itself.

Before a tagged public release, the authors need to supply a project license, a version identifier and a data-access location. The current documentation does not invent those release decisions.

## Read the Docs

1. Sign in to Read the Docs and connect the GitHub account that owns this repository.
2. Add/import the repository and select the intended default branch.
3. Keep `.readthedocs.yaml` in the repository root.
4. Trigger a build and inspect its log; enable the versions you intend to publish.
5. After a successful build, add the assigned documentation URL to the repository description and README.

The configuration uses Ubuntu 24.04, Python 3.12, Sphinx and the pinned packages in `docs/requirements.txt`. Documentation builds do not need PyQt5, ENDF or `xsdir`. The configuration follows the [Read the Docs configuration reference](https://docs.readthedocs.com/platform/stable/config-file/v2.html).

Do not add a badge claiming a successful hosted build until that project exists and builds successfully.

## Editing policy

User-facing UI text and source-code comments are English. The maintained English guide is in `docs/user-guide.md`. Update method descriptions when solver behavior changes. Keep development history out of the user workflow.
