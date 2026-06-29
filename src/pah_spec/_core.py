"""
Implements core logic
"""

import os
from astropy.constants import c, h, k_B
from astropy.io import fits
import astropy.units as u
import datetime
import h5py
import numpy as np
from numpy import exp
import pandas as pd
from scipy.integrate import trapezoid
from scipy import interpolate

from ._data import DataKind, build_data_file_path
from ._version import version as __version__

__all__ = ["PahSpec", "calc_pah_energy", "calc_pah_mode_energies"]


_DELTA_LAMBDA = 0.01

_THETA_OP_CC = 863  # C-C out-of-plane bending mode Debye temperature in Kelvin
_THETA_IP_CC = 2500  # C-C in-plane bending mode Debye temperature in Kelvin
_THETA_OP_CH = 1275  # C-H out-of-plane beinding mode Debye temperature in Kelvin
_THETA_IP_CH = 1670  # C-H in-plane bending mode Debye temperature in Kelvin
_THETA_STR_CH = 4360  # C-H stretching mode Debye temperature in Kelvin

_C_CGS = 2.99792e10  # speed of light
_H_CGS = 6.6261e-27  # Placnk constant
_KB_CGS = 1.380649e-16  # Boltzmann constant


class PahSpec:
    """Class for generating PAH emission spectra with the single photon approximation

    Parameters
    ----------
    basis_dir: path-like, optional
        Specifies the path to the directory containing "basis_ion.h5" and
        "basis_neu.h5" which store the desired basis spectra for the ionized
        and neutral PAHs, respectively. When None (the default), basis spectra
        are read from the appropriate data cache directory.
    internal_data_dir: path-like, optional
        Specifies the path to the directory holding internal data files. When
        None (the default), these files are read from the appropriate data cache
        directory. This is mostly provided for testing purposes (and should NOT
        be considered part of the public API)
    c_abs_func : callable, optional
        Function to calculate absorption cross-sections. Must have the same signature as
        PahSpec.calc_c_abs(wavelength_arr, radius_arr). Defaults to PahSpec.calc_c_abs.
    pah_energy_func : callable, optional
        Function to calculate PAH vibrational energy as a function of temperature. Must
        have the same signature as calc_pah_energy(grain_radius, temp_arr). Defaults to
        calc_pah_energy.
    """

    def __init__(
        self,
        *,
        basis_dir=None,
        internal_data_dir=None,
        c_abs_func=None,
        pah_energy_func=None,
    ):
        # Load the basis spectra into memory
        _data = {}
        _dset_unit_pairs = [
            ("lambda_em", "um"),
            ("lambda_abs", "um"),
            ("basis_spectra", "erg/(s*cm)"),
            ("grain_sizes", "Angstrom"),
        ]
        _header_attrs = {}
        _header_keywords = ["pah_energy_func", "c_abs_func"]
        for fname in ["basis_ion.h5", "basis_neu.h5"]:
            path = build_data_file_path(
                kind=DataKind.SAMPLE_BASIS, fname=fname, override_path=basis_dir
            )
            _data[fname] = {"path": path}
            _header_attrs[fname] = {}
            with h5py.File(path, "r") as f:
                for dset, unit in _dset_unit_pairs:
                    _data[fname][dset] = u.Quantity(f[dset][...], unit=unit)
                for kw in _header_keywords:
                    _header_attrs[fname][kw] = f.attrs[kw]

        for dset in ("lambda_em", "lambda_abs", "grain_sizes"):
            if not np.all(_data["basis_ion.h5"][dset] == _data["basis_neu.h5"][dset]):
                raise ValueError(
                    f"the {dset} datasets in {_data['basis_ion.h5']['path']!s} and "
                    f"{_data['basis_neu.h5']['path']!s} are not identical"
                )

        for kw in _header_keywords:
            if not np.all(
                _header_attrs["basis_ion.h5"][kw] == _header_attrs["basis_neu.h5"][kw]
            ):
                raise ValueError(
                    f"the {kw} header attributes in {_data['basis_ion.h5']['path']!s} and "
                    f"{_data['basis_neu.h5']['path']!s} are not identical"
                )

        self.emission_wavelengths = _data["basis_ion.h5"]["lambda_em"]
        self.photon_wavelengths = _data["basis_ion.h5"]["lambda_abs"]
        self.grain_sizes = _data["basis_ion.h5"]["grain_sizes"]
        self.basis_spectra_ion = _data["basis_ion.h5"]["basis_spectra"]
        self.basis_spectra_neu = _data["basis_neu.h5"]["basis_spectra"]

        energy_funcs = [
            _header_attrs[f]["pah_energy_func"]
            for f in ["basis_ion.h5", "basis_neu.h5"]
        ]
        c_abs_funcs = [
            _header_attrs[f]["c_abs_func"] for f in ["basis_ion.h5", "basis_neu.h5"]
        ]

        # If basis spectra were generated using a custom energy/c_abs function, ensure PahSpec is initialized with custom energy/c_abs function (and vice versa)
        if "custom" in energy_funcs and pah_energy_func is None:
            raise ValueError(
                "Basis spectra were computed with a custom PAH energy function, but pah_energy_func was not specified."
            )
        if "custom" in c_abs_funcs and c_abs_func is None:
            raise ValueError(
                "Basis spectra were computed with a custom cross section function, but c_abs_func was not specified."
            )
        if "default" in energy_funcs and pah_energy_func is not None:
            raise ValueError(
                "To use custom pah_energy_func, basis spectra must be computed with the same custom pah_energy_func."
            )
        if "default" in c_abs_funcs and c_abs_func is not None:
            raise ValueError(
                "To use custom c_abs_func, basis spectra must be computed with the same custom c_abs_func."
            )

        self.pah_energy_func_str = energy_funcs[0]
        self.c_abs_func_str = c_abs_funcs[0]

        # Set PAH energy and cross section functions
        if pah_energy_func is None:
            self.pah_energy_func = calc_pah_energy
        else:
            self.pah_energy_func = pah_energy_func
            print(
                "PahSpec object successfully initialized with custom PAH energy function"
            )

        if c_abs_func is None:
            self.c_abs_func = self.calc_c_abs
        else:
            self.c_abs_func = c_abs_func
            print(
                "PahSpec object successfully initialized with custom cross section function"
            )

        # Load the cross section data into memory
        # _lamj_tab units are um
        draine_table_path = build_data_file_path(
            kind=DataKind.INTERNAL_DATA,
            fname="c_abs_data/draine21_Table4.dat",
            override_path=internal_data_dir,
        )
        (
            self._lamj_tab,
            self._gamj_tab,
            self._sigj_neu_tab,
            self._sigj_ion_tab,
            self._hc_tab,
        ) = np.genfromtxt(fname=draine_table_path, unpack=True)
        self._sigj_neu_tab *= 1.0e-20  # units of cm
        self._sigj_ion_tab *= 1.0e-20  # units of cm

        fits_path = build_data_file_path(
            kind=DataKind.INTERNAL_DATA,
            fname="c_abs_data/graphite_cabs.fits",
            override_path=internal_data_dir,
        )
        with fits.open(fits_path, mode="readonly") as hdul:
            # (1) rad_graphite (um), (2) wav_graphite (um), (3) cabs_graphite (cm**2)
            self._cabs_graphite_spl = interpolate.RectBivariateSpline(
                hdul[1].data, hdul[2].data, hdul[3].data
            )

        # Load the default size distribution and ionization function into memory (std. dn/da, st. f_ion;
        # Draine et al. 2021)
        _, self.size_dist_neu, self.size_dist_ion = _read_size_dist(internal_data_dir)

        # Load the default radiation field into memory (U=1 mMMP ISRF; Draine 2011)
        self.wavelength_u_arr, self.u_lambda_arr = _read_radiation_field(
            internal_data_dir
        )

    def generate_spectrum(
        self,
        wavelength_arr=None,
        u_lambda_arr=None,
        size_dist_neu=None,
        size_dist_ion=None,
    ):
        """Scale the basis spectra for ionized and neutral PAHs to an input radiation field and integrate
        over the size distribution.

        Parameters
        ----------
        wavelength_arr : astropy.units.Quantity (array_like), optional
            Wavelength array for the radiation field u_lambda (in u.um)
        u_lambda_arr : astropy.units.Quantity (array_like), optional
            Array of length len(wavelength_arr) with the radiation field (in u.erg / u.cm ** 4)
        size_dist_neu : array_like, optional
            PAH0 size distribution and ionization function for each grain size in self.grain_sizes
        size_dist_ion : array_like, optional
            PAH+ size distribution and ionization function for each grain size in self.grain_sizes

        Returns
        -------
        spectrum_neu : astropy.units.Quantity (array_like)
            PAH0 size- and ionization-integrated spectrum for grains heated by the input u_lambda
            (in u.erg / (u.cm * u.s))
        spectrum_ion : astropy.units.Quantity (array_like)
            PAH+ size- and ionization-integrated spectrum for grains heated by the input u_lambda
            (in u.erg / (u.cm * u.s))
        """
        # If the radiation field, size distribution, or ionization function is not specified, then use the defaults
        if wavelength_arr is None:
            wavelength_arr = self.wavelength_u_arr
        if u_lambda_arr is None:
            u_lambda_arr = self.u_lambda_arr
        if size_dist_neu is None:
            size_dist_neu = self.size_dist_neu
        if size_dist_ion is None:
            size_dist_ion = self.size_dist_ion

        _check_param(wavelength_arr, u.um)
        _check_param(u_lambda_arr, u.erg / u.cm**4)

        spectra_neu_a, spectra_ion_a = self.generate_spectra_per_grain(
            wavelength_arr, u_lambda_arr
        )

        spectrum_neu = _size_integrate(self.grain_sizes, spectra_neu_a, size_dist_neu)
        spectrum_ion = _size_integrate(self.grain_sizes, spectra_ion_a, size_dist_ion)

        return spectrum_neu, spectrum_ion

    def generate_spectra_per_grain(self, wavelength_arr, u_lambda_arr):
        """Scale the basis spectra for ionized and neutral PAHs to an input radiation field

        Parameters
        ----------
        wavelength_arr : astropy.units.Quantity (array_like), optional
            Wavelength array for the radiation field u_lambda (in u.um)
        u_lambda_arr : astropy.units.Quantity (array_like), optional
            Array of length len(wavelength_arr) with the radiation field (in u.erg / u.cm ** 4)

        Returns
        -------
        spectra_neu_a : astropy.units.Quantity (array-like)
            Array of spectra for neutral PAHs in response to radiation field, with dimensions
            (len(grain_sizes), len(emission_wavelengths)) in u.erg / (u.cm * u.s)
        spectra_ion_a : astropy.units.Quantity (array-like)
            Array of spectra for ionized PAHs in response to radiation field, with dimensions
            (len(grain_sizes), len(emission_wavelengths)) in u.erg / (u.cm * u.s)
        """
        _check_param(wavelength_arr, u.um)
        _check_param(u_lambda_arr, u.erg / u.cm**4)

        basis_spectrum_unit = u.erg / (u.cm * u.s)

        spectra_neu_a = (
            np.zeros((len(self.grain_sizes), len(self.emission_wavelengths)))
            * basis_spectrum_unit
        )
        spectra_ion_a = (
            np.zeros((len(self.grain_sizes), len(self.emission_wavelengths)))
            * basis_spectrum_unit
        )

        # Check that that inputs for _scale_basis_spectra() have the expected units
        _check_param(self.photon_wavelengths, u.um, iterable=True)
        _check_param(self.emission_wavelengths, u.um, iterable=True)
        _check_param(wavelength_arr, u.um, iterable=True)
        _check_param(u_lambda_arr, u.erg / (u.cm**4), iterable=True)
        _check_param(self.basis_spectra_neu, basis_spectrum_unit, iterable=True)
        _check_param(self.basis_spectra_ion, basis_spectrum_unit, iterable=True)

        for i, grain_radius in enumerate(self.grain_sizes):
            # Remove units from arguments to call _scale_basis_spectra() (for performance reasons)
            spectra_neu_a[i] = (
                np.sum(
                    self._scale_basis_spectra(
                        self.photon_wavelengths.value,
                        _DELTA_LAMBDA,
                        self.emission_wavelengths.value,
                        grain_radius.to(u.AA).value,
                        wavelength_arr.value,
                        u_lambda_arr.value,
                        self.basis_spectra_neu[i].value,
                        ion=False,
                    ),
                    axis=0,
                )
            ) * basis_spectrum_unit  # manually set units of the function output
            spectra_ion_a[i] = (
                np.sum(
                    self._scale_basis_spectra(
                        self.photon_wavelengths.value,
                        _DELTA_LAMBDA,
                        self.emission_wavelengths.value,
                        grain_radius.to(u.AA).value,
                        wavelength_arr.value,
                        u_lambda_arr.value,
                        self.basis_spectra_ion[i].value,
                        ion=True,
                    ),
                    axis=0,
                )
            ) * basis_spectrum_unit  # manually set units of the function output

        return spectra_neu_a, spectra_ion_a

    def generate_basis_spectra(
        self,
        grain_sizes,
        emission_wavelengths,
        output_directory="./",
        ion=False,
        lambda_min=0.0912 * u.um,
        lambda_max=10 * u.um,
        pah_mode_energy_func=None,
    ):
        """Generates basis spectra file for input grain sizes for an ionized or neutral PAHs.

        If basis spectra file already exists in output_directory, appends to exisitng file.

        Parameters
        ----------
        grain_sizes : astropy.units.Quantity (float or array_like)
            Array of dust grain radii to calculate basis spectra for
        emission_wavelengths : astropy.units.Quantity (array_like)
            Array of emission wavelengths to define basis spectra over
        output_directory : str, optional
            Directory to output basis spectra to, default is ./
        ion : Bool, optional
            PAH ionization to run basis spectra for, default is False
        lambda_min : astropy.units.Quantity, optional
            Lowest lambda_abs wavelength, recommended default is 912 A
        lambda_max : astropy.units.Quantity, optional
            Highest lambda_abs wavelength, recommended default is 10 um
        pah_mode_energy_func : callable, optional
            Function to calculate PAH vibrational mode energies as a function of grain size.
            Must accept the grain radius as an argument and return an ordered list of PAH vibrational
            mode energies (in ascending order) in ergs. Defaults to calc_pah_mode_energies.


        Returns
        -------
        None
        """
        grain_sizes = _check_param(grain_sizes, u.AA, force_iterable=True)
        _check_param(lambda_min, u.um)
        _check_param(lambda_max, u.um)
        _check_param(emission_wavelengths, u.um, iterable=True)

        if self.pah_energy_func_str == "custom" and pah_mode_energy_func is None:
            raise AttributeError(
                "pah_mode_energy_func must be specified when using a custom pah_energy_func."
            )

        if not os.path.exists(output_directory):
            print(f"Path {output_directory} does not exist.")
            return None

        # generate array of lambda_abs absorbed photon wavelengths with constant spacing of _DELTA_LAMBDA
        # (dlambda / lambda = _DELTA_LAMBDA)
        photon_wavelengths = _generate_photon_wavelengths(
            _DELTA_LAMBDA, lambda_min, lambda_max
        )

        ionization = None
        if ion:
            ionization = "ion"
        else:
            ionization = "neu"

        filename = os.path.join(output_directory, f"basis_{ionization}.h5")

        # Create HDF5 data file if it doesn't already exist
        if not os.path.isfile(filename):
            with h5py.File(filename, "w") as f:
                f.attrs["description"] = (
                    f"PAH basis spectra ({ionization}) generated by pah_spec"
                )
                f.attrs["ionization_state"] = ionization
                f.attrs["c_abs_func"] = self.c_abs_func_str
                f.attrs["pah_energy_func"] = self.pah_energy_func_str
                f.attrs["pah_spec_version"] = __version__
                f.attrs["created"] = datetime.datetime.now().isoformat()
                f.attrs["reference"] = "Richie & Hensley (in prep.); arXiv:2510.16861v2"

                f.create_dataset("lambda_abs", data=photon_wavelengths.value)
                f.create_dataset("lambda_em", data=emission_wavelengths.value)
                f.create_dataset("grain_sizes", data=grain_sizes.to(u.AA).value)
                f.create_dataset(
                    "basis_spectra",
                    data=np.zeros(
                        (
                            len(grain_sizes),
                            len(photon_wavelengths),
                            len(emission_wavelengths),
                        )
                    ),
                )

                f["basis_spectra"].attrs["units"] = "erg / (s cm)"
                f["basis_spectra"].attrs[
                    "dimensions"
                ] = "(n_grains, n_lambda_abs, n_lambda_em)"

                f["lambda_abs"].attrs["units"] = "microns"
                f["lambda_abs"].attrs["description"] = "Absorbed photon wavelengths"

                f["lambda_em"].attrs["units"] = "microns"
                f["lambda_em"].attrs["description"] = "Emission wavelengths"

                f["grain_sizes"].attrs["units"] = "Angstrom"
                f["grain_sizes"].attrs["description"] = "PAH grain radii"
        else:
            print(f"Basis spectra file already exists at: {filename}")
            print("Appending to existing file.")

        for a, grain_radius in enumerate(grain_sizes):
            if not ion:
                c_abs = self.c_abs_func(
                    emission_wavelengths, [grain_radius.value] * grain_radius.unit
                )[1][0]
                print(f"C_abs computed for PAH0 of size {grain_radius:.2f}")
            else:
                c_abs = self.c_abs_func(
                    emission_wavelengths, [grain_radius.value] * grain_radius.unit
                )[0][0]
                print(f"C_abs computed for PAH+ of size {grain_radius:.2f}")

            temp_arr = np.linspace(0, 5e3, 10000) * u.K
            energy_arr = self.pah_energy_func(grain_radius, temp_arr)

            if pah_mode_energy_func is None:
                pah_mode_energy_func = calc_pah_mode_energies

            # Get the PAH fundmanetal mode energies (in order of increasing energy)
            pah_energy_modes = pah_mode_energy_func(grain_radius) * u.erg

            basis_spectra_a = (
                np.zeros((len(photon_wavelengths), len(emission_wavelengths)))
                * u.erg
                / (u.cm * u.s)
            )

            for i, lambda_abs in enumerate(photon_wavelengths):
                basis_spectra_a[i] = _compute_basis_spectrum(
                    lambda_abs,
                    grain_radius.to(u.AA),
                    emission_wavelengths,
                    c_abs,
                    temp_arr,
                    energy_arr,
                    self.pah_energy_func,
                    pah_energy_modes,
                )

            # write 2D basis spectra array for grain size to data file
            with h5py.File(
                os.path.join(output_directory, f"basis_{ionization}.h5"), "r+"
            ) as f:
                f["basis_spectra"][a] = basis_spectra_a.value

    def calc_c_abs(self, wavelength_arr, radius_arr):
        """Calculate the absorption cross-section, C_abs, for input grain sizes and wavelengths based on method from
        Draine et al. (2021).

        Parameters
        ----------
        wavelength_arr : astropy.units.Quantity (array_like)
            Array of wavelengths to calculate C_abs for
        radius_arr : astropy.units.Quantity (float or array_like)
            Array of dust grain radii to calculate C_abs for

        Returns
        -------
        c_abs_ion_out : astropy.units.Quantity (array_like)
            Array with C_abs values for ionized grains (in u.cm ** 2)
        c_abs_neu_out : astropy.units.Quantity (array_like)
            Array with C_abs values for neutral grains (in u.cm ** 2)
        """
        # TODO: check that grain does not exceed maximum allowed size

        wavelength_unit, radius_unit = u.um, u.AA

        # ensure that wavelength and grain radius variables are array-like with correct units
        wavelength_arr = _check_param(
            wavelength_arr, wavelength_unit, force_iterable=True
        )
        radius_arr = _check_param(radius_arr, radius_unit, force_iterable=True)

        c_abs_ion_out, c_abs_neu_out = self._calc_c_abs(
            wavelength_arr.value, radius_arr.value
        )

        return c_abs_ion_out * u.cm**2, c_abs_neu_out * u.cm**2

    def _calc_c_abs(self, wavelength_arr, radius_arr):
        """Internals for calc_c_abs(). Note that Astropy Units have been stripped for performance reasons.

        Parameters
        ----------
        wavelength_arr : numpy.ndarray
            Array of wavelengths to calculate C_abs for (in um)
        radius_arr : numpy.ndarray
            Array of dust grain radii to calculate C_abs for (in Angstroms)

        Returns
        -------
        c_abs_ion_out : numpy.ndarray
            Array with C_abs values for ionized grains (in cm**2)
        c_abs_neu_out : numpy.ndarray
            Array with C_abs values for neutral grains (in cm**2)
        """

        nc = _calc_nc(radius_arr * u.AA)
        hc = 0.5 * np.ones_like(nc)  # Note that this is from Equation 4 of DL07
        hc[(nc >= 25) & (nc <= 100)] = 0.5 * np.sqrt(25 / nc[(nc >= 25) & (nc <= 100)])
        hc[nc > 100] = 0.25
        xi_gra = np.zeros_like(radius_arr)
        xi_gra[radius_arr < 50] = 0.01
        xi_gra[radius_arr >= 50] = 0.01 + 0.99 * (
            1.0 - (50 / radius_arr[radius_arr >= 50]) ** 3
        )

        def s_func(lam, lamj, gamma, sigma):
            num = 2 * gamma * lamj * sigma
            denom = np.pi * (((lam / lamj) - (lamj / lam)) ** 2 + gamma**2)
            return num / denom

        def c_func(y):
            return np.arctan(1.0e3 * (y - 1.0) ** 3 / y) / np.pi + 0.5

        # Cutoff function parameters
        m_ring = 0.3 * nc
        m_ring[nc > 40] = 0.4 * nc[nc > 40]
        lamc_neu = 0.951 / (1.0 + 3.616 / np.sqrt(m_ring))  # um
        lamc_ion = 1.125 / (1.0 + 2.567 / np.sqrt(m_ring))  # um

        x = 1.0 / wavelength_arr  # um
        c_abs_ion_out = np.zeros((len(radius_arr), len(wavelength_arr)))  # cm**2
        c_abs_neu_out = np.zeros((len(radius_arr), len(wavelength_arr)))  # cm**2

        for i in range(len(radius_arr)):
            # cabs_graphite_spl expects units of um for both radius and wavelength parameters
            graph_cont = self._cabs_graphite_spl.ev(
                radius_arr[i] * 0.0001, wavelength_arr
            )  # cm**2
            c_abs_ion_out[i, :] += xi_gra[i] * graph_cont  # cm**2
            c_abs_neu_out[i, :] += xi_gra[i] * graph_cont  # cm**2

            # PAH contribution
            c_fac_neu = c_func(
                (wavelength_arr / lamc_neu[i]) ** -1
            )  # Note error in Equation A7, should be C(lambda_c/lambda)
            c_fac_ion = c_func(
                (wavelength_arr / lamc_ion[i]) ** -1
            )  # Note error in Equation A7, should be C(lambda_c/lambda)
            s_mat_neu = np.zeros((len(self._lamj_tab), len(wavelength_arr)))  # cm**2
            s_mat_ion = np.zeros((len(self._lamj_tab), len(wavelength_arr)))  # cm**2
            for j in range(len(self._lamj_tab)):
                if self._hc_tab[j] == 0:
                    # s_func param units: um -> cm, um -> cm, dimensionless, cm; units of s_func are cm**2
                    s_mat_neu[j, :] = s_func(
                        wavelength_arr * 1e-4,
                        self._lamj_tab[j] * 1e-4,
                        self._gamj_tab[j],
                        self._sigj_neu_tab[j],
                    )
                    s_mat_ion[j, :] = s_func(
                        wavelength_arr * 1e-4,
                        self._lamj_tab[j] * 1e-4,
                        self._gamj_tab[j],
                        self._sigj_ion_tab[j],
                    )
                else:
                    s_mat_neu[j, :] = s_func(
                        wavelength_arr * 1e-4,
                        self._lamj_tab[j] * 1e-4,
                        self._gamj_tab[j],
                        self._sigj_neu_tab[j] * hc[i],
                    )
                    s_mat_ion[j, :] = s_func(
                        wavelength_arr * 1e-4,
                        self._lamj_tab[j] * 1e-4,
                        self._gamj_tab[j],
                        self._sigj_ion_tab[j] * hc[i],
                    )

            idx = np.where((10 < x) & (x < 15))
            c_abs_ion_out[i, idx] += (
                (s_mat_ion[0, idx] + (1.35 * x[idx] - 3.0) * 1.0e-18)
                * nc[i]
                * (1.0 - xi_gra[i])
            )
            c_abs_neu_out[i, idx] += (
                (s_mat_neu[0, idx] + (1.35 * x[idx] - 3.0) * 1.0e-18)
                * nc[i]
                * (1.0 - xi_gra[i])
            )

            idx = np.where((7.7 < x) & (x <= 10))
            c_abs_ion_out[i, idx] += (
                (66.302 - 24.367 * x[idx] + 2.950 * x[idx] ** 2 - 0.1057 * x[idx] ** 3)
                * 1.0e-18
                * nc[i]
                * (1.0 - xi_gra[i])
            )
            c_abs_neu_out[i, idx] += (
                (66.302 - 24.367 * x[idx] + 2.950 * x[idx] ** 2 - 0.1057 * x[idx] ** 3)
                * 1.0e-18
                * nc[i]
                * (1.0 - xi_gra[i])
            )

            idx = np.where((5.9 < x) & (x <= 7.7))
            c0 = 1.8687e-18  # cm**2
            c1 = 1.905e-19  # cm**2
            c2 = 4.175e-19  # cm**2
            c3 = 4.37e-20  # cm**2, Note Equations A11 and A12 are mislabeled
            c_abs_ion_out[i, idx] += (
                (
                    s_mat_ion[1, idx]
                    + c0
                    + c1 * x[idx]
                    + c2 * (x[idx] - 5.9) ** 2
                    + c3 * (x[idx] - 5.9) ** 3
                )
                * nc[i]
                * (1.0 - xi_gra[i])
            )
            c_abs_neu_out[i, idx] += (
                (
                    s_mat_neu[1, idx]
                    + c0
                    + c1 * x[idx]
                    + c2 * (x[idx] - 5.9) ** 2
                    + c3 * (x[idx] - 5.9) ** 3
                )
                * nc[i]
                * (1.0 - xi_gra[i])
            )

            idx = np.where((3.3 < x) & (x <= 5.9))
            c_abs_ion_out[i, idx] += (
                (s_mat_ion[1, idx] + c0 + c1 * x[idx]) * nc[i] * (1.0 - xi_gra[i])
            )
            c_abs_neu_out[i, idx] += (
                (s_mat_neu[1, idx] + c0 + c1 * x[idx]) * nc[i] * (1.0 - xi_gra[i])
            )

            idx = np.where(x <= 3.3)
            c_abs_ion_out[i, idx] += (
                34.58e-18
                * 10 ** (-3.431 / x[idx])
                * c_fac_ion[idx]
                * nc[i]
                * (1.0 - xi_gra[i])
            )
            c_abs_neu_out[i, idx] += (
                34.58e-18
                * 10 ** (-3.431 / x[idx])
                * c_fac_neu[idx]
                * nc[i]
                * (1.0 - xi_gra[i])
            )

            for j in range(len(self._lamj_tab)):
                if j > 1:
                    c_abs_ion_out[i, idx] += (
                        s_mat_ion[j, idx] * nc[i] * (1.0 - xi_gra[i])
                    )
                    c_abs_neu_out[i, idx] += (
                        s_mat_neu[j, idx] * nc[i] * (1.0 - xi_gra[i])
                    )

            c_abs_ion_out[i, idx] += (
                3.5e-19
                * 10 ** (-1.45 / x[idx])
                * np.exp(-((0.1 * x[idx]) ** 2))
                * nc[i]
                * (1.0 - xi_gra[i])
            )  # Note: does not appear in D21+

        return c_abs_ion_out, c_abs_neu_out

    def _scale_basis_spectra(
        self,
        photon_wavelength_arr,
        dlambda,
        wavelength_arr,
        grain_radius,
        wavelength_arr_u,
        u_lambda_arr,
        p_lambda_arr,
        ion,
    ):
        """Calculate the energy conservation normalization to scale basis vectors to the input radiation field.
        Note that Astropy Units attributes should be stripped from input parameters (for performance reasons).

        Parameters
        ----------
        photon_wavelength_arr : numpy.ndarray
            Array of absorbed photon wavelengths (in microns)
        dlambda : float
            Width of the "monochromatic" radiation field, defined as a percentage of lambda_abs
        wavelength_arr : numpy.ndarray
            Array of emission wavelengths (in microns)
        grain_radius : float
            Dust grain radius (in Angstroms)
        wavelength_arr_u : numpy.ndarray
            Wavelength array for the radiation field u_lambda (in microns)
        u_lambda_arr : numpy.ndarray
            Array of length len(wavelengths_u) with the radiation field (in erg / cm**4)
        p_lambda_arr : numpy.ndarray
            Basis vector array of size (len(photon_wavelength_arr), len(wavelength_arr)) for a given grain
            (in erg / (cm * s))
        ion : bool
            Specifies whether or not the PAH is ionized, used to calculate the cross-section

        Returns
        -------
        scaled_spectrum : array_like
            Integrated spectrum for the input grain and radiation field (in erg / (cm * s)) with length
            len(wavelength_arr))
        """

        # call calc_c_abs once for all wavelength values needed by calc_normalization
        lambda1_max = photon_wavelength_arr[-1] + photon_wavelength_arr[-1] * dlambda
        if ion:
            c_abs_arr_u = self.c_abs_func(
                wavelength_arr_u * u.um, np.array([grain_radius]) * u.AA
            )[0][0].value
            c_abs_phot_arr = np.concatenate(
                [
                    self.c_abs_func(
                        photon_wavelength_arr * u.um, np.array([grain_radius]) * u.AA
                    )[0][0].value,
                    self.c_abs_func(
                        np.array([lambda1_max]) * u.um, np.array([grain_radius]) * u.AA
                    )[0][0].value,
                ],
            )
        else:
            c_abs_arr_u = self.c_abs_func(
                wavelength_arr_u * u.um, np.array([grain_radius]) * u.AA
            )[1][0].value
            c_abs_phot_arr = np.concatenate(
                [
                    self.c_abs_func(
                        photon_wavelength_arr * u.um, np.array([grain_radius]) * u.AA
                    )[1][0].value,
                    self.c_abs_func(
                        np.array([lambda1_max]) * u.um, np.array([grain_radius]) * u.AA
                    )[1][0].value,
                ],
            )

        normalizations = np.zeros(len(photon_wavelength_arr))

        for i, lambda_abs in enumerate(photon_wavelength_arr):
            # define wavelengths of monochromatic radiation field (MRF)
            lambda0, lambda1 = lambda_abs, lambda_abs + lambda_abs * dlambda
            wh_mrf = np.where(
                np.logical_and(wavelength_arr_u > lambda0, wavelength_arr_u < lambda1)
            )
            wavelength_arr_mrf = np.concatenate(
                [[lambda0], wavelength_arr_u[wh_mrf], [lambda1]]
            )

            # get cross-section values of MRF
            wh0 = np.argmin(np.abs(lambda0 - photon_wavelength_arr))
            wh1 = np.argmin(np.abs(lambda1 - photon_wavelength_arr))
            c_abs_arr_mrf = np.concatenate(
                [[c_abs_phot_arr[wh0]], c_abs_arr_u[wh_mrf], [c_abs_phot_arr[wh1]]]
            )

            # interpolate radiation field to the MRF wavelength range
            u_lambda_arr_mrf = np.interp(
                wavelength_arr_mrf, wavelength_arr_u, u_lambda_arr
            )

            # integrate to determine the power of the radiation field in this wavelength range
            numerator = trapezoid(
                u_lambda_arr_mrf * c_abs_arr_mrf * _C_CGS,
                x=wavelength_arr_mrf
                * 1e-4,  # convert wavelength_arr_mrf from u.um to u.cm
            )

            # integrate to determine the power radiated by the grain over all wavelengths
            denominator = trapezoid(
                p_lambda_arr[i],
                x=wavelength_arr * 1e-4,  # convert wavelength_arr from u.um to u.cm
            )

            normalizations[i] = numerator / denominator

        return p_lambda_arr * normalizations[:, np.newaxis]


############################## general pah_spec functions ##############################


def calc_pah_energy(grain_radius, temp_arr):
    """Calculate PAH vibrational energy as a function of temperature.

    Parameters
    ----------
    grain_radius : astropy.units.Quantity (float)
        The PAH effective radius
    temp_arr : astropy.units.Quantity (array-like)
        Array of temperature to calculate energies for

    Returns
    -------
    energy_arr : astropy.units.Quantity (array-like)
        Resulting PAH energy array (in u.erg)
    """
    radius_unit, temp_unit = u.AA, u.K
    _check_param(grain_radius, radius_unit)
    _check_param(temp_arr, temp_unit, iterable=True)

    # after ensuring units are correct, strip astropy units and call private version of energy function, which
    # is much faster
    energy_arr = _calc_pah_energy(grain_radius.value, temp_arr.value) * u.erg

    return energy_arr


def _calc_pah_energy(grain_radius, temp_arr):
    """Calculate PAH vibrational energy as a function of temperature.

    Parameters
    ----------
    grain_radius : float
        The PAH effective radius in Angstrom
    temp_arr : numpy.ndarray
        Array of temperature to calculate energies for in Kelvin

    Returns
    -------
    energy_arr : numpy.ndarray
        Resulting PAH energy array in ergs
    """
    nc, nh, nm_cc_op, nm_cc_ip = _calc_pah_structure(grain_radius * u.AA)

    # determines if grain energy will be calculated using Debye spectrum method or discrete mode method
    nc_cutoff = 7360

    energy_arr = None

    # For grains with more carbon atoms than nc_cutoff, calculate C-C energies using the Debye spectrum approximation
    # Eq. 33 of Draine & Li (2001)
    if nc > nc_cutoff:
        energy_arr = _calc_pah_energy_debye(temp_arr, nh, nm_cc_ip, nm_cc_op)

    # For grains with fewer carbon atoms than nc_cutoff, calculate energies by summing contributions from individual
    # modes
    # Eq. 2 of Draine & Li (2001)
    if nc <= nc_cutoff:
        energy_arr = _calc_pah_energy_modes(
            temp_arr, nc, nh, nm_cc_ip, nm_cc_op, grain_radius
        )

    return energy_arr


def _calc_pah_energy_debye(temp_arr, nh, nm_cc_ip, nm_cc_op):
    """Calculate PAH energy using eq. 33 of Draine & Li (2001).

    Parameters
    ----------
    temp_arr : numpy.ndarray
        Array of temperature to calculate energies for in Kelvin
    nm_cc_ip : int
        Number of C-C in-plane stretching modes
    nm_cc_op : int
        Number of C-C out-of-plane stretching modes
    nh : int
        Number of hydrogen atoms

    Returns
    -------
    energy_arr : numpy.ndarray
        Resulting PAH energy array in ergs
    """
    # contributions from individual C-H modes, as implemented in lines 157-181 of pah_spec_heat.f
    energy_ch = np.zeros(len(temp_arr))  # erg

    x = _THETA_OP_CH / temp_arr
    y = exp(x)
    tmin = 32  # Kelvin
    energy_ch[temp_arr > tmin] += (
        nh
        * (_KB_CGS * temp_arr[temp_arr > tmin])
        * (x[temp_arr > tmin] / (y[temp_arr > tmin] - 1))
    )

    x = _THETA_IP_CH / temp_arr
    y = exp(x)
    tmin = 42  # Kelvin
    energy_ch[temp_arr > tmin] += (
        nh
        * (_KB_CGS * temp_arr[temp_arr > tmin])
        * (x[temp_arr > tmin] / (y[temp_arr > tmin] - 1))
    )

    x = _THETA_STR_CH / temp_arr
    y = exp(x)
    tmin = 109  # Kelvin
    energy_ch[temp_arr > tmin] += (
        nh
        * (_KB_CGS * temp_arr[temp_arr > tmin])
        * (x[temp_arr > tmin] / (y[temp_arr > tmin] - 1))
    )

    # contributions from C-C modes (approximated as Debye spectrum)
    energy_cc_op = np.zeros(len(temp_arr))  # erg
    energy_cc_op += _debye_2(_THETA_OP_CC / temp_arr) * _KB_CGS * temp_arr

    energy_cc_ip = np.zeros(len(temp_arr))  # erg
    energy_cc_ip += _debye_2(_THETA_IP_CC / temp_arr) * _KB_CGS * temp_arr

    # total energy
    energy_arr = energy_ch + nm_cc_op * energy_cc_op + nm_cc_ip * energy_cc_ip

    return energy_arr


def calc_pah_mode_energies(grain_radius):
    """Calculate vibrational mode energies for a PAH

    Parameters
    ----------
    grain_radius : astropy.units.Quantity (float)
        PAH radius in u.AA

    Returns
    -------
    emode_arr : astropy.units.Quantity (array-like)
        Array of PAH vibrational energies in u.erg, sorted in order of increasing energy
    """
    _check_param(grain_radius, u.AA)

    emode_arr = _calc_pah_mode_energies(grain_radius.value) * u.erg

    return emode_arr


def _calc_pah_energy_modes(temp_arr, nc, nh, nm_cc_ip, nm_cc_op, grain_radius):
    """Calculate PAH energy using eq. 2 of Draine & Li (2001).

    Parameters
    ----------
    temp_arr : numpy.ndarray
        Array of temperature to calculate energies for in Kelvin
    nc : int
        Number of carbon atoms
    nh : int
        Number of hydrogen atoms
    nm_cc_ip : int
        Number of C-C in-plane stretching modes
    nm_cc_op : int
        Number of C-C out-of-plane stretching modes
    grain_radius : float
        PAH radius in Angstroms

    Returns
    -------
    energy_arr : astropy.units.Quantity (array_like)
        Resulting PAH energy array (in u.erg)
    emode_arr : astropy.units.Quantity (array_like)
        Energies of all PAH vibrational modes in order of increasing energy (in u.erg)
    """

    emode_arr = _calc_pah_mode_energies(grain_radius)

    nmodes = len(emode_arr)

    energy_arr = np.zeros(len(temp_arr))  # ergs

    exp_cutoff = 100  # ignore contributions from exp(x) when x is large

    # emodes array contains mode energies from C-C in-plane and out-of-plane bending and C-H stretching and in-plane
    # in-plane and out-of-plane bending
    for j in range(0, nmodes):
        x = emode_arr[j] / (_KB_CGS * temp_arr)
        y = exp(x)
        energy_arr[x < exp_cutoff] += (
            (x[x < exp_cutoff] / (y[x < exp_cutoff] - 1))
            * _KB_CGS
            * temp_arr[x < exp_cutoff]
        )

    return energy_arr


def _calc_pah_cooling_discrete(
    lambda_abs,
    grain_radius,
    wavelength_arr,
    c_abs_arr,
    temp_arr,
    energy_arr,
    pah_energy_modes,
):
    """Calculate the temperature evolution of a PAH following a single-photon absorption.

    Parameters
    ----------
    lambda_abs : astropy.units.Quantity (float)
        Wavelength of absorbed photon
    grain_radius : astropy.units.Quantity (float)
        Dust grain radius
    wavelength_arr : astropy.units.Quantity (array_like)
        Emission wavelengths to integrate over
    c_abs_arr : astropy.units.Quantity (array_like)
        Array of length len(wavelength_arr) with C_abs values for a grain of size grain_radius
    temp_arr : astropy.units.Quantity (array_like)
        Array of temperatures corresponding to PAH vibrational energies for a grain of size grain_radius
    energy_arr : astropy.units.Quantity (array_like)
        Array of PAH vibrational energies as a function of temperature for a grain of size grain_radius
    pah_energy_modes : astropy.units.Quantity (array_like)
        Array of PAH fundamental mode energies for a grain of size grain_radius

    Returns
    -------
    dt_arr_out : astropy.units.Quantity (array_like)
        Array of time-steps used to solve for T(t) (in u.s)
    time_arr_out : astropy.units.Quantity (array_like)
        Array with time values for T(t) (in u.s)
    temp_arr_out : astropy.units.Quantity (array_like)
        Array with temperature values for T(t) (at the energy bin edges, in u.K)
    """
    wavelength_unit, radius_unit, c_abs_unit, energy_unit, temp_unit = (
        u.um,
        u.AA,
        u.cm**2,
        u.erg,
        u.K,
    )
    _check_param(lambda_abs, wavelength_unit)
    _check_param(grain_radius, radius_unit)
    _check_param(wavelength_arr, wavelength_unit, iterable=True)
    _check_param(c_abs_arr, c_abs_unit, iterable=True)
    _check_param(energy_arr, energy_unit, iterable=True)
    _check_param(temp_arr, temp_unit, iterable=True)

    nu_arr = c.cgs / wavelength_arr.to(u.cm)

    energy_abs = c.cgs * h.cgs / lambda_abs.to(u.cm)
    temp_abs = np.interp(energy_abs.value, energy_arr.value, temp_arr.value) * u.K

    # print(f"Photon wavelength: {lambda_abs:.2f}, initial temperature: {temp_abs:.2f}")

    time_i, energy_i, temp_i = 0 * u.s, energy_abs, temp_abs
    time_arr_out, temp_arr_out, dt_arr_out = [0], [temp_i.value], []
    dt_unit, time_unit, temp_unit = None, None, None
    # The change in energy per timestep, results are converged for dE < 3%
    dE_max = 0.03

    # Define first 10 energy bins following method in DL01 Appendix B
    nbins = 10

    energy_j_min = np.zeros(nbins) * u.erg
    energy_j_max = np.zeros(nbins) * u.erg

    # Eq. B1
    energy_j_min[1] = 3 / 2 * pah_energy_modes[1] - 1 / 2 * pah_energy_modes[2]
    energy_j_max[0] = energy_j_min[1]

    # Eq. B2
    for j in range(0, 2):
        bin_edge = 1 / 2 * (pah_energy_modes[j] + pah_energy_modes[j + 1])
        energy_j_max[j] = bin_edge
        energy_j_min[j + 1] = bin_edge

    # Eq. B3
    for j in range(2, nbins):
        bin_edge = 1 / 2 * (pah_energy_modes[2 * j - 1] + pah_energy_modes[2 * j])
        energy_j_max[j] = bin_edge
        if j != nbins - 1:
            energy_j_min[j + 1] = bin_edge

    # exclude emission from frequencies exceeding the frequency of the absorbed photon
    weights_sawteeth = np.ones(len(nu_arr))
    weights_sawteeth[nu_arr > c.cgs / lambda_abs.to(u.cm)] = 0

    # Cache variables outside of loop to modify the final temperature bin before manually implementing temperature bins
    dE, dt, dE_dt = None, None, None
    while energy_i > energy_j_max[-1]:
        dE_dt = -trapezoid(
            4
            * np.pi
            * _planck_function_nu(nu_arr, temp_i)
            * c_abs_arr
            * weights_sawteeth,
            x=nu_arr,
        )

        dt = dE_max * energy_i / dE_dt
        dE = dE_dt * dt

        energy_i -= dE
        time_i += dt

        temp_i = np.interp(energy_i.value, energy_arr.value, temp_arr.value) * u.K

        dt_arr_out.append(dt.value)
        time_arr_out.append(time_i.value)
        temp_arr_out.append(temp_i.value)
        dt_unit, time_unit, temp_unit = dt.unit, time_i.unit, temp_i.unit

    # Ensure that the final temperature bin lines up with the pre-defined energy bins
    # Undo the final timestep
    energy_i += dE
    time_i -= dt

    # Calculate the intended timestep
    dE_correction = energy_i - energy_j_max[-1]
    dt_correction = dE_correction / dE_dt

    # Apply the correction
    energy_i -= dE_correction  # should line up with energy_j_max[-1]
    time_i += dt_correction
    temp_i = np.interp(energy_i, energy_arr, temp_arr)

    # Overwrite the final values
    dt_arr_out[-1] = dt_correction.value
    time_arr_out[-1] = time_i.value
    temp_arr_out[-1] = temp_i.value

    for j in range(nbins - 1, -1, -1):
        dE_dt = -trapezoid(
            4
            * np.pi
            * _planck_function_nu(nu_arr, temp_i)
            * c_abs_arr
            * weights_sawteeth,
            x=nu_arr,
        )

        dE = energy_j_max[j] - energy_j_min[j]
        dt = dE / dE_dt

        energy_i -= dE
        time_i += dt
        temp_i = np.interp(energy_i.value, energy_arr.value, temp_arr.value) * u.K

        dt_arr_out.append(dt.value)
        time_arr_out.append(time_i.value)
        temp_arr_out.append(temp_i.value)

        dt_unit, time_unit, temp_unit = dt.unit, time_i.unit, temp_i.unit

    # ensure that zero energy corresponds to zero temperature (np.interp sometimes returns nonzero/NaN values)
    temp_arr_out[-1] = 0

    dt_arr_out = np.array(dt_arr_out) * dt_unit
    time_arr_out = np.array(time_arr_out) * time_unit
    temp_arr_out = np.array(temp_arr_out) * temp_unit

    # print(f"Final temperature of {temp_arr_out[-1]:.2f} at time {time_arr_out[-1]:.2e}")

    return dt_arr_out, time_arr_out, temp_arr_out


def _calc_basis_vector(
    grain_radius,
    lambda_abs,
    wavelength_arr,
    c_abs_arr,
    weighting_arr,
    temp_arr,
    pah_energy_func,
):
    """Calculate the basis vector for a single-photon absorption.

    Parameters
    ----------
    grain_radius : float
        Dust grain radius in Angstroms
    lambda_abs : float
        Wavelength of absorbed photon in um
    wavelength_arr : numpy.ndarray
        Array of emission wavelengths in um
    c_abs_arr : numpy.ndarray
        Array of length len(wavelength_arr) with C_abs values for a given grain in cm**2
    weighting_arr : numpy.ndarray
        Dimensionless array of length len(temp_arr) with values to weight temperatures by
    temp_arr : numpy.ndarray
        Array with temperature values for T(t) (upper energy bin edges only, in Kelvin)
    pah_energy_func : callable
        Function for calculating the PAH energy as a function of grain size and temperature

    Returns
    -------
    basis_vector : numpy.ndarray
        Array of floats of length len(wavelength_arr) with basis vector for a given grain (in erg / (cm * s))
    """
    weighting_energy_arr = pah_energy_func(
        grain_radius=grain_radius * u.AA, temp_arr=temp_arr * u.K
    ).value  # erg

    basis_vector = np.zeros(len(wavelength_arr))
    for i, lambda_i in enumerate(wavelength_arr):
        weights_sawteeth = np.ones(len(weighting_energy_arr))
        # Condition in Eq. 56 of DL01 excluding emission from lambda_em such that E(lambda_em) > E(lambda_abs)
        weights_sawteeth[_H_CGS * _C_CGS / lambda_i > weighting_energy_arr] = 0
        p_lambda_i = np.sum(
            4
            * np.pi
            * _planck_function_lambd(lambda_i * 1e-4, temp_arr)
            * weights_sawteeth
            * weighting_arr
            * c_abs_arr[i]
        )
        basis_vector[i] = p_lambda_i

    return basis_vector


def _compute_basis_spectrum(
    lambda_abs,
    grain_radius,
    emission_wavelengths,
    c_abs_arr,
    temp_arr,
    energy_arr,
    pah_energy_func,
    pah_energy_modes,
):
    dt_arr, _, temp_arr_t = _calc_pah_cooling_discrete(
        lambda_abs,
        grain_radius.to(u.AA),
        emission_wavelengths,
        c_abs_arr,
        temp_arr,
        energy_arr,
        pah_energy_modes,
    )
    temp_arr_t = temp_arr_t[0:-1]  # temperatures of the top edges of each energy bin
    temp_weights = dt_arr / np.sum(dt_arr)

    _check_param(grain_radius, u.AA)
    _check_param(lambda_abs, u.um)
    _check_param(emission_wavelengths, u.um, iterable=True)
    _check_param(temp_arr, u.K, iterable=True)
    _check_param(c_abs_arr, u.cm**2, iterable=True)
    _check_param(temp_weights, u.dimensionless_unscaled)

    return (
        _calc_basis_vector(
            grain_radius.value,
            lambda_abs.value,
            emission_wavelengths.value,
            c_abs_arr.value,
            temp_weights.value,
            temp_arr_t.value,
            pah_energy_func,
        )
        * u.erg
        / (u.cm * u.s)
    )


def _calc_nh(nc):
    """Eq. 8 of Draine & Li (2001)"""
    if nc <= 25:
        return round(0.5 * nc + 0.5)
    if nc > 25 and nc <= 100:
        return round(2.5 * np.sqrt(nc) + 0.5)
    if nc > 100:
        return round(0.25 * nc + 0.5)

    return None


def _calc_nc(a):
    """Eq. 8 of Draine & Li (2021)"""
    _check_param(a, u.AA)
    nc = 418 * (a / (10 * u.AA)) ** 3
    return nc.value.round(0).astype(int)


def _calc_pah_structure(grain_radius):
    _check_param(grain_radius, u.AA)
    # 5 types of vibration: out-of-plane C-C modes, in-plane C-C modes, out-of-plane C-H bending, in-plane C-H bending,
    # and C-H stretching.
    nc = _calc_nc(grain_radius)  # number of carbon atoms
    nh = _calc_nh(nc)  # number of hydrogen atoms
    # TODO: consider turning these into functions
    nm_cc_op = nc - 2  # total number of C-C out-of-plane modes
    nm_cc_ip = 2 * (nc - 2)  # total number of C-C in-plane modes

    return nc, nh, nm_cc_op, nm_cc_ip


def _calc_delta_j(j):
    """Eq. 6 and 5 of Draine & Li (2001). Adjusts location of first 3 modes to bring agreement with mode spectrum of
    coronene C_{24}H_{12}.
    """
    if j == 2 or j == 3:
        return 1
    return 1 / 2


def _calc_beta(nc, nm):
    """Eq. 7 of Draine & Li (2001)"""
    if nc <= 54:
        return 0
    if nc > 54 and nc <= 102:
        return (1 / (2 * nm - 1)) * ((nc - 54) / 52)
    if nc > 102:
        return (1 / (2 * nm - 1)) * (((nc - 2) / 52) * (102 / nc) ** (2 / 3) - 1)

    return None


def _calc_pah_mode_energies(grain_radius):
    nc, nh, nm_cc_op, nm_cc_ip = _calc_pah_structure(grain_radius * u.AA)

    mode_arr_cc_op = _calc_cc_mode_energies(nc, nm_cc_op, _THETA_OP_CC)
    mode_arr_cc_ip = _calc_cc_mode_energies(nc, nm_cc_ip, _THETA_IP_CC)
    mode_arr_ch_op = _calc_ch_mode_energies(nh, _THETA_OP_CH)
    mode_arr_ch_ip = _calc_ch_mode_energies(nh, _THETA_IP_CH)
    mode_arr_ch_str = _calc_ch_mode_energies(nh, _THETA_STR_CH)

    emode_arr = sorted(
        np.concatenate(
            (
                mode_arr_cc_op,
                mode_arr_cc_ip,
                mode_arr_ch_op,
                mode_arr_ch_ip,
                mode_arr_ch_str,
            )
        )
    )  # ergs

    return emode_arr


def _calc_cc_mode_energies(nc, nm, theta):
    """As implemented in lines 141-154 and 158-165 of pah_modes.f"""
    mode_energy_arr = np.zeros((nm))  # erg
    beta = _calc_beta(
        nc, nm
    )  # note that we are intentionally using the same beta for C-C ip and op modes

    for j in range(1, nm + 1):
        mode_energy_arr[j - 1] = (
            _KB_CGS * theta * np.sqrt((1 - beta) * (j - _calc_delta_j(j)) / nm + beta)
        )

    return mode_energy_arr


def _calc_ch_mode_energies(nh, theta):
    """As implemented in lines 170-184 of pah_modes.f"""
    mode_energy_arr = np.zeros((nh))  # erg
    mode_energy_arr[:] = _KB_CGS * theta

    return mode_energy_arr


def _debye_2(x):
    """As implemented in lines 28-29 of debye.f"""
    x2 = x * x
    return (2 / x2) * (
        2 * 1.20206
        - x2 * exp(-x) * (1 + 2 / x + 2 / x2 + exp(-x) * (0.5 + 0.5 / x + 0.25 / x2))
    )


def _check_param(param, unit, iterable=False, force_iterable=False):
    """Check that input parameters have the correct units (and optionally check if they are array-like).

    Parameters
    ----------
    param : astropy.units.Quantity
        Input parameter to check
    unit : astropy.units.Unit
        Expected units of param
    iterable : Boolean, optional
        Whether to check if param is array-like

    Returns
    -------
    None

    Raises
    ------
    TypeError
        If the input is not an astropy.units.Quantity object
    ValueError
        If the astropy.units.Quantity object has incorrect units (or optionally is not array-like)
    """
    if not isinstance(param, u.Quantity):
        raise TypeError(f"expects astropy.units.Quantity objects of unit {unit}")

    if param.unit != unit:
        raise ValueError(f"parameter expects units of {unit}")

    if iterable and (not isinstance(param.value, (list, tuple, np.ndarray))):
        raise TypeError("expects array-like input")

    if force_iterable:
        if not isinstance(param.value, (list, tuple, np.ndarray)):
            param = np.array([param.value]) * param.unit
        return param

    return None


def _planck_function_nu(nu, T):
    _check_param(nu, u.Hz)
    _check_param(T, u.K)
    return (2 * h.cgs * nu**3 / c.cgs**2) * 1 / (np.exp(h.cgs * nu / (k_B.cgs * T)) - 1)


def _planck_function_lambd(lambd, T):
    """expects lambd in units of cm and T in units of K"""
    return (
        (2 * _H_CGS * _C_CGS**2 / lambd**5)
        * 1
        / (np.exp(_H_CGS * _C_CGS / (lambd * _KB_CGS * T)) - 1)
    )


def _generate_photon_wavelengths(dlambda, lambda_min, lambda_max):
    _check_param(lambda_min, u.um)
    _check_param(lambda_max, u.um)

    # define photon/basis vector wavelengths
    lambda_i = lambda_min.value
    photon_wavelengths = [lambda_i]
    while lambda_i < lambda_max.value:
        lambda_i += dlambda * lambda_i
        photon_wavelengths.append(lambda_i)
    photon_wavelengths *= u.um

    return photon_wavelengths


def _read_size_dist(internal_data_dir=None):
    path = build_data_file_path(
        kind=DataKind.INTERNAL_DATA,
        fname="defaults/pahspec_dnda.out_st_std",
        override_path=internal_data_dir,
    )

    df = pd.read_csv(path, sep="\\s+", skiprows=1)
    rad = df["rad"].to_numpy() * u.um
    size_dist_neu = df["dn_{PAH0}"].to_numpy()
    size_dist_ion = df["dn_{PAH+}"].to_numpy()

    return rad, size_dist_neu, size_dist_ion


def _read_radiation_field(internal_data_dir=None):
    path = build_data_file_path(
        kind=DataKind.INTERNAL_DATA,
        fname="defaults/isrf_mmpisrf_0.00",
        override_path=internal_data_dir,
    )

    df = pd.read_csv(path, sep="\\s+", skiprows=6)
    wavelength_arr_u = df["(um)"].to_numpy() * u.um
    u_lambda_arr = df["(erg"].to_numpy() * (u.erg / u.cm**3) / wavelength_arr_u.to(u.cm)

    return wavelength_arr_u, u_lambda_arr


def _size_integrate(grain_sizes, scaled_basis_spectra, size_dist):
    # Assumes size_dist acounts for size distribution and ionization function
    weighted_spectra = np.zeros((len(grain_sizes), len(scaled_basis_spectra[0])))
    for i, _ in enumerate(grain_sizes):
        weighted_spectra[i] = size_dist[i] * scaled_basis_spectra[i].value

    size_integrated_spectrum = (
        np.sum(np.array(weighted_spectra), axis=0) * scaled_basis_spectra.unit
    )

    return size_integrated_spectrum
