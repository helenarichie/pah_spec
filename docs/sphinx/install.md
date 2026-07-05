# Installation

## Installing ``pah_spec``

Currently, you **MUST** install this package from source.
We describe two approaches down below.

Each approach will automatically install ``pah_spec``'s dependencies (please open an issue if you encounter a version incompatability).
For context, most development has used python 3.13.7, h5py 3.14.0, astropy 6.1.2, numpy 2.3.2, and scipy 1.16.1.

### Easy Approach (for normal users)

For normal users, the easiest way to do this is to invoke:

```sh
python -m pip install --user git+https://github.com/helenarichie/pah_spec
```

Behind the scenes, pip will download the full repository, install the python package, and delete the repository.

### Editable Mode (if you want to modify the code)

If you are interested in developing ``pah_spec``, you probably want to install the package in editable mode.
To do this, you need to clone the git repository and invoke the following command from the root of the repository (i.e. from the same directory as the pyproject.toml file)

```sh
python -m pip install --user -e .
```

In this mode, pip effectively installs links to the source files.
Consequently, if you mutate one of the source files, the installed package will also see those updates.

## Installing the data files

``pah_spec`` makes extensive use of external data files.
For proper packaging, this data can't be shipped within the package and must be installed separately.
For user convenience, ``pah_spec`` provides functionality to easily download this data.

By default, this package assumes that files are stored in a data-cache.
Internally, we make use of the [pooch](https://www.fatiando.org/pooch) package to implement this logic.
- The default location of this cache depends on OS conventions.
- Care is taken so that when multiple versions of ``pah_spec`` are installed on a single machine different versions of cached data files won't have conflicting paths.

### Kinds of datafiles

There are 2 kinds of datafiles

:::{object} internal_data

   This consists of assorted data totalling about 2 MB in size and that are downloaded from ``pah_spec``'s git repository.
   The version number associated with these files is the same as ``pah_spec``'s version number.

   You can overwrite the path to the cache directory holding these files by setting the :envvar:`!PAHSPEC_INTERNALDATA_DIR` environment variable.
:::

:::{object} sample_basis

   This consists of sample basis spectra totalling about 2 GB in size.
   They are downloaded from a [data repository on the Harvard Dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/LUUXEJ).
   The version number associated with these files is separate from ``pah_spec``'s version number.
   Thus, different versions of ``pah_spec`` may use a single version of these files (care is taken to avoid unnecessarily duplicating these files in the cache directory).

   You can overwrite the path to the cache directory holding these files by setting the :envvar:`!PAHSPEC_SAMPLEBASIS_DIR` environment variable.
:::

### Approaches for installing data

We highlight the approaches to install this data below:

:::{tab} command-line (simple)
If you just want to install all files to the data cache, you can simply invoke the following from the command line:

```shell-session
$ python -m "pah_spec" download --to-cache all
```
:::
:::{tab} command-line (custom)
If you want to install a particular kind of data to the cache, you can invoke

```shell-session
$ python -m "pah_spec" download --to-cache KIND
```

If you want to install a particular kind of data to an arbitrary directory (denoted as ``<PATH/TO/DIR>``, you can invoke:

```shell-session
$ python -m "pah_spec" download --to=<PATH/TO/DIR> <KIND>
```

In both cases ``<KIND>`` should be either ``sample_basis`` or ``internal_data``.

:::
:::{tab} programmatically

You can use the {py:func}`pah_spec.retrieve_internal_data` and {py:func}`pah_spec.retrieve_sample_basis` functions
By default, each function will try to install the data to the appropriate data cache.
You can overwrite the destination with the ``path`` kwarg.
:::

No matter what approach you take,

- a progressbar will be shown if the [tqdm](https://tqdm.github.io/) python package is installed.
- the installation logic will not initiate a download of a file if a file already exists at the output location (and errors will be raised an existing file has the wrong checksum).

:::{warning}
You should generally try to avoid using the capability to specify an arbitrary path to manually specify the path to the default cache directory (in other words, if you want to use the default cache directory, don't specify an arbitrary path).
:::

:::{note}
While we allow ``internal_data`` to an arbitrary location, that is mostly for parity with ``sample_basis``.
You should generally prefer to install ``internal_data`` to the cache.
:::

## Uninstalling ``pah_spec``

In order to fully uninstall ``pah_spec``, you need to manually delete the data file cache directories.

You can find these directories by invoking the following command from the commandline

```shell-session
$ python -m "pah_spec" show-cache
```

Of course, to uninstall ``pah_spec`` you can invoke

```shell-session
$ python -m pip uninstall pah_spec
```
