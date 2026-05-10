# Development

## Testing

To run the tests, you should install ``pah_spec`` and then execute the following command from the root of the repository:

```sh
python -m unittest discover -s tests -p "pah_spec_test.py"
```

## Writing and Building Docs

The ``pah_spec`` documentation is written in a mixture of ReStructuredText and [MyST Markdown](https://myst-parser.readthedocs.io/en/latest/index.html).
It is converted into a website with the [Sphinx](https://www.sphinx-doc.org/) documentation generator.

## What is MyST Markdown?

There are *MANY* flavors of Markdown with various extensions.
[MyST](https://myst-parser.readthedocs.io/en/latest/index.html) ("Markedly Structured Text") extends the CommonMark Markdown specification to address the deficiencies in the core Markdown language that make the core language a poor choice for writing extensive and modern documentation.

The [myst-parser](https://myst-parser.readthedocs.io/en/latest/index.html) is intended to integrate with the [Sphinx](https://www.sphinx-doc.org/) documentation generator.
Sphinx was originally designed to generate documentation from the ReStructuredText markup language, which implements extensions in terms of roles and directives.
Consequently, the MyST Markdown focuses on implements syntax extensions for implementing roles and directives so that it achieves parity with RestructuredText.

## Building a local copy of the documentation

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

