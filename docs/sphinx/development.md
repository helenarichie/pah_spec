# Development

## Testing

We use the [pytest](https://docs.pytest.org/en/stable/) framework to write tests.

To install the appropriate version of this tool, you can invoke the following command from the root of the repository.

```sh
pip install --group test
```

To actually execute the tests (``pah_spec`` must be installed), you should invoke ``pytest`` (or ``python -m pytest``) from the root of the repository.


## Style Formatting

The ``pah_spec`` repository is configured with tools to enforce checks on the various files in the repository.
The primary checks apply the Ruff formatter and Ruff linter to the python files.

Because we realize that these tools may seem overwhelming or complicated, we have configured the repository to try to simplify the experience to the greatest extent possible.
All of these checks are managed by the [pre-commit](https://pre-commit.com/) software.
We discuss how to automatically invoke these checks as part of a PR (no local installation is required) or on your local machine down below.

:::{note}
It's worth clarifying there are essentially 3 distinct entities named "pre-commit":

1. the [pre-commit](https://pre-commit.com/) software.
   Contributors **only** need to know about this if they want to apply the check enforcement tools locally.
2. the [pre-commit.ci](https://pre-commit.ci/) continuous integration service.
   This is named because the service simply executes the pre-commit software.
   (All contributors will encounter this).
3. the ["pre-commit" git hook](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks).
   This is one of multiple different "hooks" offered by git.
   The pre-commit software is named after this hook because it was originally designed to be used with this hook.
   We do **NOT** currently recommend that contributors use the pre-commit software with the pre-commit hook (unless they fully understand the choice that they are making).
:::

### Automatic Checks (No installation required)

At a basic level, you don't have to worry about manually invoking any of these tools on your machine.
In fact, you are free to entirely ignore the existence of these tools until it comes time to submit a Pull Request.
When you submit a Pull Request (and whenever you update it), various forms of continuous integration are triggered.

For the present discussion, the [pre-commit.ci](https://pre-commit.ci/) continuous integration tool is of primary relevance.
This tool executes all of the formatting tools and if your submission doesn't satisfy all of the requirements, it will fail and report each problem.

If you don't want to install anything locally, there are 2 approaches for addressing problems:

1. You can leave a comment on the Pull Request that simply states

   > pre-commit.ci autofix

   and pre-commit.ci will contribute push a commit to your branch that fixes as many issues as possible.
   **This is the recommended way to fix all code-formatting issues.**

2. You can also manually fix the issues locally (and then push your changes).
   **You will need to do this to address most linting errors.**


### Running the Checks Locally

The easiest way to run the checks locally is to install the pre-commit software and use pre-commit to invoke the checks. The pre-commit software is written in python and can be installed with ``pip`` (the [installation instructions](https://pre-commit.com/#installation) provides further details and alternative approaches).

Once you have installed ``pre-commit``, you can enforce the checks by invoking the following command from the root of your ``pah_spec`` repository:

```sh
pre-commit run --all-files
```

The above command does 2 things:

1. First, it ensure that local copies of the correct versions of the required enforcement tools are installed.
   These local copies are only accessed by pre-commit and won't affect other parts of your system.
   These copies are also cached (so that the tools don't need to be reinstalled on every invocation).

2. Then the command applies the enforcement tools on the files in your repository (tool-specific exclusions, like files listed by ``tool.ruff.format.exclude`` in **pyproject.toml**, are obviously respected).

:::{caution}
The above command will modify the files in your repository (after all, that's the whole point of the command).
The pre-commit software does not provide a way to reverse this change.
:::

### Summary of Code Checks

#### Ruff Formatter (Python code formatting)

Python code is formatted by the [Ruff Formatter](https://docs.astral.sh/ruff/formatter/) tool.
This is provided as part of the popular [Ruff](https://github.com/astral-sh/ruff) tool (i.e. it's invoked with ``ruff format``).
It performs a similar role to the [Black code formatter](https://black.readthedocs.io/en/stable/).

#### Ruff Linter (Python code linting)

Python code is linted by the [Ruff Linter](https://docs.astral.sh/ruff/linter/) tool.
This is also provided as part of the popular [Ruff](https://github.com/astral-sh/ruff) tool (i.e. it's invoked via ``ruff check``).
It performs a similar role to the [Flake8 linter](https://pypi.org/project/flake8/).


#### Miscellaneous checks
Some miscellaneous checks are also implemented by a set of miscellaneous enforcement scripts provided by the authors of pre-commit.

## Writing and Building Docs

The ``pah_spec`` documentation is written in a mixture of ReStructuredText and [MyST Markdown](https://myst-parser.readthedocs.io/en/latest/index.html).
It is converted into a website with the [Sphinx](https://www.sphinx-doc.org/) documentation generator.

### What is MyST Markdown?

There are *MANY* flavors of Markdown with various extensions.
[MyST](https://myst-parser.readthedocs.io/en/latest/index.html) ("Markedly Structured Text") extends the CommonMark Markdown specification to address the deficiencies in the core Markdown language that make the core language a poor choice for writing extensive and modern documentation.

The [myst-parser](https://myst-parser.readthedocs.io/en/latest/index.html) is intended to integrate with the [Sphinx](https://www.sphinx-doc.org/) documentation generator.
Sphinx was originally designed to generate documentation from the ReStructuredText markup language, which implements extensions in terms of roles and directives.
Consequently, the MyST Markdown focuses on implements syntax extensions for implementing roles and directives so that it achieves parity with RestructuredText.

### Building a local copy of the documentation

To locally build the documentation, you should invoke the following commands from the root of the ``pah_spec`` repository (we suggest that you do this in a virtual python environment).
We're assuming that you already installed ``pah_spec``.

```shell-session
# ensure that pip is up to date (we need version 25.1+)
$ python -m pip install --upgrade pip

# install the python-dependencies
$ python -m pip install --group dev

# actually build the documentation
$ python -m sphinx -M html docs/sphinx "_build" -W
```

At this point, you can render the documentation.
The precise command depends on your system and browser.

::::{tab} Linux
:::{tab} chrome
```shell-session
$ google-chrome _build/html/index.html
```
:::
:::{tab} firefox
```shell-session
$ firefox _build/html/index.html
```
:::
::::

::::{tab} macOS
:::{tab} chrome
```shell-session
$ open -a "Google Chrome" _build/html/index.html
```
:::
:::{tab} firefox
```shell-session
$ open -a firefox _build/html/index.html
```
:::
::::

