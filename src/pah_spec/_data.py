"""
Logic for retrieving the data files
"""

from ._version import __version__

import enum
import functools
import os
import pathlib
from types import MappingProxyType
from typing import Any, Optional, TYPE_CHECKING


import pooch

if TYPE_CHECKING:
    import argparse.Namespace as argparse_Namespace
else:
    argparse_Namespace = Any

SubParserBuilder = Any  # <- type object returned by ArgumentParser.add_subparsers


class DataKind(enum.Enum):
    """Identifies the kinds of data sets known to ``pah_spec``."""

    # this contains all data that doesn't belong to another kind
    # - this includes absorption coefficient data as well as default choices
    #   for the grain size-distributions & interstellar radiation field
    # - if this just contained text-files under 100 kB, it would make sense to
    #   directly embed the files in the packages uploaded to pypi (sdist and
    #   wheels) using importlib.resources
    # - However, the 1.8 MB binary FITs file should clearly needs to be
    #   separately downloaded (rather than embedded). If we are downloading one
    #   file, then we may as well download all of them.
    # - since the total size of this dataset is ~2 MB, it makes sense to
    #   download everything at once
    INTERNAL_DATA = "internal_data"

    # this is the standard pre-computed set of basis spectra
    # - this is about 2 GB of data
    # - users may want to track this data separately (and they may want more
    #   control over where the data is stored)
    # - the remote location of this data is actually outside of pah_spec's
    #   GitHub repository and actually has its own version number
    SAMPLE_BASIS = "sample_basis"


# this fetches the appropriate (OS-specific) location for storing cache locations
_CACHE_BASE_PATH: pathlib.Path = pooch.os_cache("pah_spec")

# we primarily use pooch.Pooch as a collection of data
# - we access its path and registry attributes and use its get_url method
# - For context all the other methods are very focused on lazily downloading the data
#   only as it becomes necessary.
#   - While that strategy is sensible for providing sample data, it doesn't provide us
#     with any benefits (we generally need all of none of the files in a dataset).
#   - This behavior can be surprising to users
#   - This behavior may be undesirable when performing parallelized calculations,
#     especially on a cluster with a shared file-system
# - Nevertheless, its still convenient to use pooch.Pooch because we are basically
#   following all of its documented conventions for actually specifying what data is
#   downloaded, where it comes from, and where it is saved
_POOCH_MAP: MappingProxyType[DataKind, pooch.Pooch] = MappingProxyType(
    {
        DataKind.INTERNAL_DATA: pooch.create(
            # default storage location for caching files
            path=_CACHE_BASE_PATH / "internal_data",
            # we don't currently take advanatage of the ability to change urls based on
            # the version kwarg
            base_url="https://raw.githubusercontent.com/helenarichie/pah_spec/99044bf52e0c48de4819017818ec16b06dc28e65/data/",
            # for simplicity, we tie these files' version to pah_spec's version
            version=__version__,
            # users can use this environment variable to specify a custom location that
            # overwrides the path kwarg
            env="PAHSPEC_INTERNALDATA_DIR",
            # files (and checksums) managed by the constructed instance
            registry={
                "c_abs_data/graphite_cabs.fits": "sha256:390fa0963b5687bebe591c1d99af83069b9411976f934e0ff51f2ca316a25fb8",
                "c_abs_data/draine21_Table4.dat": "sha256:a184608c735479a1d403fcd88e187d23ed2214dedc1bf39f0db230659068568a",
                "defaults/isrf_mmpisrf_0.00": "sha256:50002023a13c4343e9c3e3a450d9b9985c34344a453ea7cb05a87679702206b1",
                "defaults/pahspec_dnda.out_st_std": "sha256:063f27a24300e704016ea5eb6fb60bc6af11799acfba80b1689d6d8a617c6b59",
            },
            # we're being a little conservative here (I suspect this choice is better
            # for parallelized calculations on shared file systems), but I don't think
            # we actually call methods where this is relevant
            allow_updates=False,
        ),
        DataKind.SAMPLE_BASIS: pooch.create(
            path=_CACHE_BASE_PATH / "sample_basis",
            # we explicitly require that a custom url is set for each file (since the
            # urls don't end in the name of the file)
            #
            # While pooch accepts a DOI for the url (in our case, it would be
            # `doi:10.7910/DVN/LUUXEJ`) that would be a HUGE mistake since
            # -> dataverse doesn't change the doi when you upload a newer version of
            #    files and pooch always tries to fetch the newest version of the files
            # -> thus, uploading an updated version of any of the data files (even if
            #    we were just adding some extra metadata) would change the checksum and
            #    older versions of pah_spec (using the old checksum) would break
            base_url="",
            # if we need to update the doi/hashes, update the version number so that it
            # remains synchronized with the version-number on Harvard Dataverse
            # - if you forget to update the version number, then cache paths may be
            #   accidently reused for different versions of a file (if you try to
            #   download data using different versions of pah_spec)
            version="4.0",
            # users can use this environment variable to specify a custom location that
            # overwrides the path kwarg
            env="PAHSPEC_SAMPLEBASIS_DIR",
            # files (and checksums) managed by the constructed instance
            registry={
                "basis_ion.h5": "sha256:8ac9e685b56afa15bb80185b5fb46ec9688056a32d63ac10e929924e36303adb",
                "basis_neu.h5": "sha256:3f1a02cd9dae135d6cd9c0b49a93b6b84c39b888bdaa9a89a69b9183d01fc3b3",
            },
            # here we specify custom URLs for each file
            urls={
                "basis_ion.h5": "https://dataverse.harvard.edu/api/access/datafile/13967678",
                "basis_neu.h5": "https://dataverse.harvard.edu/api/access/datafile/13967679",
            },
            # we're being a little conservative here (I suspect this choice is better
            # for parallelized calculations on shared file systems), but I don't think
            # we actually call methods where this is relevant
            allow_updates=False,
        ),
    }
)


def _retrieve_data(
    *, path: Optional[os.PathLike], kind: DataKind, progressbar: Any = None
):
    """
    Internal function that tries to retrieve all files in the specified data collection

    Parameters
    ----------
    path: ``os.PathLike``, optional
        When specified, all files are downloaded to this directory. When None,
        the default, files are downloaded to the cache directory.
    kind: ``DataKind``
        Specifies the kind of data to download.
    progressbar: bool or arbitrary progressbar, optional
        Forwarded onto :py:func:`pooch.retrieve`. When ``None`` (the default),
        the value is inferred based on whether the
        `tqdm <https://tqdm.github.io/>`__ package is installed
    """
    if progressbar is None:
        try:
            import tqdm  # noqa: F401

            progressbar = True
        except ImportError:
            progressbar = False

    pooch_obj = _POOCH_MAP[kind]
    for fname, cksum in pooch_obj.registry.items():
        pooch.retrieve(
            url=pooch_obj.get_url(fname),
            known_hash=cksum,
            fname=fname,
            path=pooch_obj.path if path is None else path,
            progressbar=progressbar,
        )


_DOCSTRING_TEMPLATE = """
Ensure that all {descr} has been downloaded.

If the data has already been downloaded and has the appropriate checksum
it won't be redownloaded.

Parameters
----------
path: ``os.PathLike``, optional
    When specified, all files are downloaded to this directory. When None,
    the default, files are downloaded to the cache directory.
progressbar: bool or arbitrary progressbar, optional
    Forwarded onto :external:py:func:`pooch.retrieve`. When ``None`` (the default),
    the value is inferred based on whether the
    `tqdm <https://tqdm.github.io/>`__ package is installed
"""

retrieve_internal_data = functools.partial(
    _retrieve_data, path=None, kind=DataKind.INTERNAL_DATA
)
retrieve_internal_data.__doc__ = _DOCSTRING_TEMPLATE.format(descr="internal data")

retrieve_sample_basis = functools.partial(
    _retrieve_data, path=None, kind=DataKind.SAMPLE_BASIS
)
retrieve_sample_basis.__doc__ = _DOCSTRING_TEMPLATE.format(descr="sample basis spectra")


def build_data_file_path(
    kind: DataKind, fname: str, override_path: Optional[os.PathLike] = None
) -> pathlib.Path:
    """
    Internal function that constructs the path to a file based

    Parameters
    ----------
    kind
        Specifies the kind of dataset that the file is part of.
    fname: str
        The name of the file to try to download
    override_path: path-like, optional
        When not None, this encodes a user-specified path

    Returns
    -------
    pathlib.Path
        Full concatenated path to ``fname``
    """
    pooch_obj = _POOCH_MAP[kind]
    if fname not in pooch_obj.registry:
        raise ValueError(f"{fname} is not a known file of type {kind._value_}")

    if override_path is None:
        out_path = pathlib.Path(pooch_obj.path) / fname
        if not os.path.exists(out_path):
            raise FileNotFoundError(
                f"{out_path} is missing from the data cache. Did you forget to "
                f"prefetch all files of type {kind._value_}?"
            )
    else:
        out_path = pathlib.Path(override_path) / fname
        if not os.path.exists(out_path):
            raise FileNotFoundError(f"{override_path!s} does not contain {fname}")
    return out_path


def add_data_cli_subcommands(subparsers: Any):
    """
    Add command line subcommands pertaining to data management

    Parameters
    ----------
    subparsers
        Object returned by :py:func:`argparse.ArgumentParser.add_subparsers`
    """

    def show_cache_cmd(args: argparse_Namespace) -> int:
        for member in DataKind:
            name = member._value_
            pooch_obj = _POOCH_MAP[member]
            print(f"{name}: {pooch_obj.path}")
        return 0

    show_cache_parser = subparsers.add_parser(
        "show-cache",
        help=(
            "Show paths to cache directories where pah_spec will store each kind of "
            "data. This is primarily used to help clean up cached data when "
            "uninstalling the package."
        ),
    )
    show_cache_parser.set_defaults(fn=show_cache_cmd)

    def download_cmd(args: argparse_Namespace) -> int:
        path = args.download_dst_path
        if args.kind == "all":
            if path is not None:
                raise ValueError("can only download `all` data if saving to cache")
            kinds = list(DataKind)
        else:
            kinds = [DataKind(args.kind)]
        for kind in kinds:
            print(
                f"\nchecking for {kind._value_} data "
                "(it will be downloaded, if not present)"
            )
            _retrieve_data(path=path, kind=kind)
        return 0

    download_parser = subparsers.add_parser(
        "download", help="Download data used by pah_spec"
    )
    dst_grp = download_parser.add_mutually_exclusive_group(required=True)
    dst_grp.add_argument(
        "--to-cache",
        action="store_const",
        const=None,
        dest="download_dst_path",
        help="data is downloaded to the cache location",
    )
    dst_grp.add_argument(
        "--to",
        dest="download_dst_path",
        help="path to the directory where data is saved",
    )
    download_parser.add_argument(
        "kind",
        choices=["all"] + [member._value_ for member in DataKind],
        help="specify the kind of data to download",
    )
    download_parser.set_defaults(fn=download_cmd)
