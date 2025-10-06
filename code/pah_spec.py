from astropy.constants import c, h, k_B
import astropy.units as u
import numpy as np
from numpy import exp
import os
import pandas as pd
from scipy.integrate import trapezoid


GRAIN_SIZES = [3.5481e-04, 3.7584e-04, 3.9811e-04, 4.2170e-04, 4.4668e-04, 4.7315e-04, 5.0119e-04, 5.3088e-04,
                5.6234e-04, 5.9566e-04, 6.3096e-04, 6.6834e-04, 7.0795e-04, 7.4989e-04, 7.9433e-04, 8.4140e-04,
                8.9125e-04, 9.4406e-04, 1.0000e-03, 1.0593e-03, 1.1220e-03, 1.1885e-03, 1.2589e-03, 1.3335e-03,
                1.4125e-03, 1.4962e-03, 1.5849e-03, 1.6788e-03, 1.7783e-03, 1.8836e-03, 1.9953e-03, 2.1135e-03, 
                2.2387e-03, 2.3714e-03, 2.5119e-03] * u.um
_DELTA_LAMBDA = 0.01


class PahSpec:

    def __init__(self, basis_directory="../data/basis_spectra", c_abs_data_directory="../data/c_abs_data"):

        script_directory = os.path.dirname(os.path.abspath(__file__))

        if not os.path.exists(os.path.join(script_directory, c_abs_data_directory)):
          raise FileNotFoundError(
            f"calc_cabs expects to find qabs_001um.dat and draine21_Table4.dat at {c_abs_data_directory}, but the path does not exist"
        )
        if not os.path.exists(os.path.join(script_directory, basis_directory)):
          raise FileNotFoundError(
            f"Exptected to find basis spectra at {basis_directory}, but the path does not exist"
        )

        self.emission_wavelengths, self.photon_wavelengths, self.basis_spectra_neu, self.basis_spectra_ion = _read_basis_spectra(os.path.join(script_directory, basis_directory), GRAIN_SIZES)

        self.wav_graphite, self.qabs = np.genfromtxt(os.path.join(script_directory, c_abs_data_directory, "qabs_001um.dat"), unpack=True, usecols=[0, 2])
        self.wav_graphite *= u.um

        self.lamj_tab, self.gamj_tab, self.sigj_neu_tab, self.sigj_ion_tab, self.hc_tab = np.genfromtxt(
            os.path.join(c_abs_data_directory, os.path.join(script_directory, c_abs_data_directory, "draine21_Table4.dat")), unpack=True
        )
        self.lamj_tab *= u.um
        self.sigj_neu_tab *= 1.0e-20
        self.sigj_ion_tab *= 1.0e-20
        self.sigj_neu_tab *= u.cm
        self.sigj_ion_tab *= u.cm


    def generate_spectrum(self, wavelength_arr, u_lambda_arr, size_dist_neu, size_dist_ion):
        _check_param(wavelength_arr, u.um)
        _check_param(u_lambda_arr, u.erg / u.cm ** 4)

        spectra_neu_a = np.zeros((len(GRAIN_SIZES), len(self.emission_wavelengths))) * u.erg / (u.cm * u.s)
        spectra_ion_a = np.zeros((len(GRAIN_SIZES), len(self.emission_wavelengths))) * u.erg / (u.cm * u.s)
        for i, grain_size in enumerate(GRAIN_SIZES):
            spectra_neu_a[i] = np.sum(self.scale_basis_spectra(self.photon_wavelengths, _DELTA_LAMBDA, self.emission_wavelengths, grain_size.to(u.AA), wavelength_arr, u_lambda_arr, self.basis_spectra_neu[i], ion=False), axis=0)
            spectra_ion_a[i] = np.sum(self.scale_basis_spectra(self.photon_wavelengths, _DELTA_LAMBDA, self.emission_wavelengths, grain_size.to(u.AA), wavelength_arr, u_lambda_arr, self.basis_spectra_ion[i], ion=True), axis=0)

        spectrum_neu = _size_integrate(spectra_neu_a, size_dist_neu)
        spectrum_ion = _size_integrate(spectra_ion_a, size_dist_ion)

        return spectrum_neu, spectrum_ion
    
    def generate_basis_spectra(self, grain_sizes, lambda_min, lambda_max, emission_wavelengths, output_directory="./", c_abs_neu=None, c_abs_ion=None):

        _check_param(grain_sizes, u.AA)
        _check_param(lambda_min, u.um)
        _check_param(lambda_max, u.um)
        _check_param(emission_wavelengths, u.um, iterable=True)
        if c_abs_neu != None:
            _check_param(c_abs_neu, u.cm ** 2, iterable=True)
        if c_abs_ion != None:
            _check_param(c_abs_ion, u.cm ** 2, iterable=True)

        if not os.path.exists(output_directory):
            print(f"Path {output_directory} does not exist.")
            return None

        photon_wavelengths = _generate_photon_wavelengths(_DELTA_LAMBDA, lambda_min, lambda_max)

        for grain_size in grain_sizes:

            if c_abs_neu == None:
                c_abs_neu = self.calc_c_abs(emission_wavelengths, grain_size)[1][0]
            if c_abs_ion == None:
                c_abs_ion = self.calc_c_abs(emission_wavelengths, grain_size)[0][0]
            
            temp_arr = np.linspace(1, 5e3, 10000) * u.K
            energy_arr = calc_pah_energy(grain_size, temp_arr)

            basis_dict_neu = dict()
            basis_dict_ion = dict()
            # basis_spectra_ion = np.zeros((len(photon_wavelengths, len(emission_wavelengths)))) * u.erg / (u.cm * u.s)

            for i, lambda_abs in enumerate(photon_wavelengths):

                basis_dict_neu[str(lambda_abs.value)] = _basis_spectrum(lambda_abs, grain_size.to(u.AA), emission_wavelengths, c_abs_neu, temp_arr, energy_arr)
                basis_dict_ion[str(lambda_abs.value)] = _basis_spectrum(lambda_abs, grain_size.to(u.AA), emission_wavelengths, c_abs_ion, temp_arr, energy_arr)

            df = pd.DataFrame(basis_dict_neu)
            df.to_csv(os.path.join(output_directory, f"basis_neu_{grain_size.to(u.AA):.3f}"), index=False)
            df = pd.DataFrame(basis_dict_ion)
            df.to_csv(os.path.join(output_directory, f"basis_ion_{grain_size.to(u.AA):.3f}"), index=False)

        photon_wavelengths = [s.value for s in photon_wavelengths]
        emission_wavelengths = [s.value for s in emission_wavelengths]
        
        data = dict()
        data["lambda_em [um]"] = emission_wavelengths
        df = pd.DataFrame(data)
        df.to_csv(os.path.join(output_directory, "lambda_em.csv"), index=False)

        data = dict()
        data["lambda_abs [um]"] = photon_wavelengths
        df = pd.DataFrame(data)
        df.to_csv(os.path.join(output_directory, "lambda_abs.csv"), index=False)


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

        Raises
        ------
        AttributeError
            If the input is not an astropy.units.Quantity object
        TypeError
            If the astropy.units.Quantity object has incorrect units (or optionally is not array-like)
        FileNotFoundError
            If the directory data_path does not exist
        """
        wavelength_unit, radius_unit = u.um, u.AA

        # ensure that wavelength and grain radius variables are array-like with correct units
        wavelength_arr = _check_param(wavelength_arr, wavelength_unit, force_iterable=True)
        radius_arr = _check_param(radius_arr, radius_unit, force_iterable=True)

        if not isinstance(radius_arr.value, (list, tuple, np.ndarray)):
            radius_arr = np.array([radius_arr.value]) * radius_arr.unit

        # C_abs per volume doesn't depend on size in the small-grain limit
        cabs_V_graphite = 3.0 * self.qabs / (4 * 1.0e-3 * u.um)
        cabs_V_graphite_out = np.interp(
            wavelength_arr.to(u.um),
            self.wav_graphite[::-1].to(u.um),
            cabs_V_graphite[::-1].to(u.cm**-1),
        )

        rho_C = 2.0 * u.g / u.cm**3
        nc = _calc_nc(radius_arr)  # number of Carbon atoms, Draine 2021
        HC = 0.5 * np.ones_like(nc)  # Note that this is from Equation 4 of DL07
        HC[(nc >= 25) & (nc <= 100)] = 0.5 * np.sqrt(25 / nc[(nc >= 25) & (nc <= 100)])
        HC[nc > 100] = 0.25
        xi_gra = np.zeros_like(radius_arr.value)
        xi_gra[radius_arr < 50 * u.AA] = 0.01
        xi_gra[radius_arr >= 50 * u.AA] = 0.01 + 0.99 * (1.0 - (50 * u.AA / radius_arr[radius_arr >= 50 * u.AA]) ** 3)

        def S_func(lam, lamj, gamma, sigma):
            num = 2 * gamma * lamj * sigma
            denom = np.pi * (((lam / lamj) - (lamj / lam)) ** 2 + gamma**2)
            return num / denom

        def C_func(y):
            return np.arctan(1.0e3 * (y - 1.0) ** 3 / y) / np.pi + 0.5

        # Cutoff function parameters
        M_ring = 0.3 * nc
        M_ring[nc > 40] = 0.4 * nc[nc > 40]
        lamc_neu = 0.951 * u.um / (1.0 + 3.616 / np.sqrt(M_ring))
        lamc_ion = 1.125 * u.um / (1.0 + 2.567 / np.sqrt(M_ring))

        x = (1.0 * u.um / wavelength_arr.to(u.um)).value
        c_abs_ion_out = np.zeros((len(radius_arr), len(wavelength_arr))) * u.cm**2
        c_abs_neu_out = np.zeros((len(radius_arr), len(wavelength_arr))) * u.cm**2

        for i, _ in enumerate(radius_arr):
            # Graphite contribtuion
            volume = (4.0 / 3.0) * np.pi * radius_arr[i] ** 3
            c_abs_ion_out[i, :] += xi_gra[i] * cabs_V_graphite_out * volume
            c_abs_neu_out[i, :] += xi_gra[i] * cabs_V_graphite_out * volume

            # PAH contribution
            C_fac_neu = C_func(
                (wavelength_arr.to(u.um) / lamc_neu[i]).value ** -1
            )  # Note error in Equation A7, should be C(lambda_c/lambda)
            C_fac_ion = C_func(
                (wavelength_arr.to(u.um) / lamc_ion[i]).value ** -1
            )  # Note error in Equation A7, should be C(lambda_c/lambda)
            S_mat_neu = np.zeros((len(self.lamj_tab), len(wavelength_arr))) * u.cm**2
            S_mat_ion = np.zeros((len(self.lamj_tab), len(wavelength_arr))) * u.cm**2
            for j, _ in enumerate(self.lamj_tab):
                if self.hc_tab[j] == 0:
                    S_mat_neu[j, :] = S_func(wavelength_arr.to(u.um), self.lamj_tab[j], self.gamj_tab[j], self.sigj_neu_tab[j])
                    S_mat_ion[j, :] = S_func(wavelength_arr.to(u.um), self.lamj_tab[j], self.gamj_tab[j], self.sigj_ion_tab[j])
                else:
                    S_mat_neu[j, :] = S_func(wavelength_arr.to(u.um), self.lamj_tab[j], self.gamj_tab[j], self.sigj_neu_tab[j] * HC[i])
                    S_mat_ion[j, :] = S_func(wavelength_arr.to(u.um), self.lamj_tab[j], self.gamj_tab[j], self.sigj_ion_tab[j] * HC[i])

            idx = np.where((10 < x) & (x < 15))
            c_abs_ion_out[i, idx] += (
                (S_mat_ion[0, idx] + (1.35 * x[idx] - 3.0) * 1.0e-18 * u.cm**2) * nc[i] * (1.0 - xi_gra[i])
            )
            c_abs_neu_out[i, idx] += (
                (S_mat_neu[0, idx] + (1.35 * x[idx] - 3.0) * 1.0e-18 * u.cm**2) * nc[i] * (1.0 - xi_gra[i])
            )

            idx = np.where((7.7 < x) & (x <= 10))
            c_abs_ion_out[i, idx] += (
                (66.302 - 24.367 * x[idx] + 2.950 * x[idx] ** 2 - 0.1057 * x[idx] ** 3)
                * 1.0e-18
                * u.cm**2
                * nc[i]
                * (1.0 - xi_gra[i])
            )
            c_abs_neu_out[i, idx] += (
                (66.302 - 24.367 * x[idx] + 2.950 * x[idx] ** 2 - 0.1057 * x[idx] ** 3)
                * 1.0e-18
                * u.cm**2
                * nc[i]
                * (1.0 - xi_gra[i])
            )

            idx = np.where((5.9 < x) & (x <= 7.7))
            c0 = 1.8687e-18 * u.cm**2
            c1 = 1.905e-19 * u.cm**2
            c2 = 4.175e-19 * u.cm**2
            c3 = 4.37e-20 * u.cm**2  # Note Equations A11 and A12 are mislabeled
            c_abs_ion_out[i, idx] += (
                (S_mat_ion[1, idx] + c0 + c1 * x[idx] + c2 * (x[idx] - 5.9) ** 2 + c3 * (x[idx] - 5.9) ** 3)
                * nc[i]
                * (1.0 - xi_gra[i])
            )
            c_abs_neu_out[i, idx] += (
                (S_mat_neu[1, idx] + c0 + c1 * x[idx] + c2 * (x[idx] - 5.9) ** 2 + c3 * (x[idx] - 5.9) ** 3)
                * nc[i]
                * (1.0 - xi_gra[i])
            )

            idx = np.where((3.3 < x) & (x <= 5.9))
            c_abs_ion_out[i, idx] += (S_mat_ion[1, idx] + c0 + c1 * x[idx]) * nc[i] * (1.0 - xi_gra[i])
            c_abs_neu_out[i, idx] += (S_mat_neu[1, idx] + c0 + c1 * x[idx]) * nc[i] * (1.0 - xi_gra[i])

            idx = np.where(x <= 3.3)
            c_abs_ion_out[i, idx] += (
                34.58e-18 * 10 ** (-3.431 / x[idx]) * u.cm**2 * C_fac_ion[idx] * nc[i] * (1.0 - xi_gra[i])
            )
            c_abs_neu_out[i, idx] += (
                34.58e-18 * 10 ** (-3.431 / x[idx]) * u.cm**2 * C_fac_neu[idx] * nc[i] * (1.0 - xi_gra[i])
            )

            for j in range(len(self.lamj_tab)):
                if j > 1:
                    c_abs_ion_out[i, idx] += S_mat_ion[j, idx] * nc[i] * (1.0 - xi_gra[i])
                    c_abs_neu_out[i, idx] += S_mat_neu[j, idx] * nc[i] * (1.0 - xi_gra[i])

            c_abs_ion_out[i, idx] += (
                3.5e-19 * 10 ** (-1.45 / x[idx]) * np.exp(-((0.1 * x[idx]) ** 2)) * nc[i] * (1.0 - xi_gra[i]) * u.cm**2
            )  # Note: does not appear in D21+

        return c_abs_ion_out, c_abs_neu_out
    

    def scale_basis_spectra(self, photon_wavelength_arr, dlambda, wavelength_arr, grain_radius, wavelength_arr_u, u_lambda_arr, p_lambda_arr, ion):
        """Calculate the energy conservation normalization to scale basis vectors to the input radiation field.

        Parameters
        ----------
        photon_wavelength_arr : astropy.units.Quantity (array_like)
            Array of absorbed photon wavelengths
        dlambda : float
            Width of the "monochromatic" radiation field, defined as a percentage of lambda_abs
        wavelength_arr : astropy.units.Quantity (array_like)
            Array of emission wavelengths
        grain_radius : astropy.units.Quantity (float)
            Dust grain radius
        wavelength_arr_u : astropy.units.Quantity (array_like)
            Wavelength array for the radiation field u_lambda
        u_lambda_arr : astropy.units.Quantity (array_like)
            Array of length len(wavelengths_u) with the radiation field
        p_lambda_arr : astropy.units.Quantity (array_like)
            Basis vector array of size (len(photon_wavelength_arr), len(wavelength_arr)) for a given grain
        ion : bool
            Specifies whether or not the PAH is ionized, used to calculate the cross-section

        Returns
        -------
        integrated_spectrum : astropy.units.Quantity (array_like)
            Integrated spectrum for the input grain and radiation field (in units of u.erg / (u.cm * u.s) with length len(wavelength_arr))

        Raises
        ------
        AttributeError
            If the input is not an astropy.units.Quantity object
        TypeError
            If the astropy.units.Quantity object has incorrect units (or optionally is not array-like)
        """
        wavelength_unit = u.um
        radiation_field_unit = u.erg / (u.cm**4)
        basis_vector_unit = u.erg / (u.cm * u.s)
        _check_param(photon_wavelength_arr, wavelength_unit, iterable=True)
        _check_param(photon_wavelength_arr, wavelength_unit, iterable=True) 
        _check_param(wavelength_arr, wavelength_unit, iterable=True)
        _check_param(grain_radius, u.AA)
        _check_param(wavelength_arr_u, wavelength_unit, iterable=True)
        _check_param(u_lambda_arr, radiation_field_unit, iterable=True)
        _check_param(p_lambda_arr, basis_vector_unit, iterable=True)

        # call calc_c_abs once for all wavelength values needed by calc_normalization
        lambda1_max = photon_wavelength_arr[-1] + photon_wavelength_arr[-1] * dlambda
        if ion:
            c_abs_arr_u = self.calc_c_abs(wavelength_arr_u, grain_radius)[0][0]
            c_abs_phot_arr = self.calc_c_abs(photon_wavelength_arr, grain_radius)[0][0]
            c_abs_phot_arr = np.insert(c_abs_phot_arr, len(c_abs_phot_arr), self.calc_c_abs(lambda1_max, grain_radius)[0][0])
        else:
            c_abs_arr_u = self.calc_c_abs(wavelength_arr_u, grain_radius)[1][0]
            c_abs_phot_arr = self.calc_c_abs(photon_wavelength_arr, grain_radius)[1][0] 
            c_abs_phot_arr = np.insert(c_abs_phot_arr, len(c_abs_phot_arr), self.calc_c_abs(lambda1_max, grain_radius)[1][0])

        normalizations = np.zeros(len(photon_wavelength_arr))

        for i, lambda_abs in enumerate(photon_wavelength_arr):
            # define wavelengths of monochromatic radiation field (MRF)
            lambda0, lambda1 = lambda_abs, lambda_abs + lambda_abs * dlambda
            wh_mrf = np.where(np.logical_and(wavelength_arr_u > lambda0, wavelength_arr_u < lambda1))
            wavelength_arr_mrf = wavelength_arr_u[wh_mrf]
            wavelength_arr_mrf = np.insert(wavelength_arr_mrf, [0, len(wavelength_arr_mrf)], [lambda0, lambda1])

            # get cross-section values of MRF
            c_abs_arr_mrf = c_abs_arr_u[wh_mrf]
            wh0, wh1 = np.argmin(np.abs(lambda0 - photon_wavelength_arr)), np.argmin(np.abs(lambda1 - photon_wavelength_arr))
            c_abs_arr_mrf = np.insert(c_abs_arr_mrf, [0, len(c_abs_arr_mrf)], [c_abs_phot_arr[wh0], c_abs_phot_arr[wh1]])

            # interpolate radiation field to the MRF wavelength range
            u_lambda_arr_mrf = np.interp(wavelength_arr_mrf, wavelength_arr_u, u_lambda_arr)

            # integrate to determine the power of the radiation field in this wavelength range
            numerator = trapezoid(u_lambda_arr_mrf * c_abs_arr_mrf * c.cgs, x=wavelength_arr_mrf.to(u.cm))

            # integrate to determine the power radiated by the grain over all wavelengths
            denominator = trapezoid(p_lambda_arr[i], x=wavelength_arr.to(u.cm))

            normalizations[i] = numerator / denominator

        # return np.sum(p_lambda_arr * normalizations[:, np.newaxis], axis=0)
        return p_lambda_arr * normalizations[:, np.newaxis]
    

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
    

def _size_integrate(scaled_basis_spectra, size_dist):

    size_integrated_spectrum = []

    for i, _ in enumerate(GRAIN_SIZES):
        size_integrated_spectrum.append(size_dist[i] * scaled_basis_spectra[i].value)

    size_integrated_spectrum = np.sum(np.array(size_integrated_spectrum), axis=0) * scaled_basis_spectra.unit

    return size_integrated_spectrum


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
    energy_arr : astropy.units.Quantity (array_like)
        Resulting PAH energy array (in u.erg)

    Raises
    ------
    AttributeError
        If the input is not an astropy.units.Quantity object
    TypeError
        If the astropy.units.Quantity object has incorrect units (or optionally is not array-like)
    """
    radius_unit, temp_unit = u.AA, u.K
    _check_param(grain_radius, radius_unit)
    _check_param(temp_arr, temp_unit, iterable=True)

    # 5 types of vibration: out-of-plane C-C modes, in-plane C-C modes, out-of-plane C-H bending, in-plane C-H bending, and C-H stretching.
    nc = _calc_nc(grain_radius)  # number of carbon atoms
    nh = _calc_nh(nc)  # number of hydrogen atoms
    nm_cc_op = nc - 2  # total number of C-C out-of-plane modes
    nm_cc_ip = 2 * (nc - 2)  # total number of C-C in-plane modes

    # TODO: consider moving out of function and naming as constants _THETA_OP_CC
    theta_op_cc = 863 * u.K  # C-C out-of-plane bending mode Debye temperature
    theta_ip_cc = 2500 * u.K  # C-C in-plane bending mode Debye temperature
    theta_op_ch = 1275 * u.K  # C-H out-of-plane beinding mode Debye temperature
    theta_ip_ch = 1670 * u.K  # C-H in-plane bending mode Debye temperature
    theta_str_ch = 4360 * u.K  # C-H stretching mode Debye temperature

    # determines if grain energy will be calculated using Debye spectrum method or discrete mode method
    nc_cutoff = 7360

    # For grains with more carbon atoms than nc_cutoff, calculate C-C energies using the Debye spectrum approximation
    # Eq. 33 of Draine & Li (2001)
    if nc > nc_cutoff:
        energy_arr = _calc_pah_energy_debye(
            temp_arr, nh, nm_cc_ip, nm_cc_op, theta_ip_cc, theta_op_cc, theta_ip_ch, theta_op_ch, theta_str_ch
        )

    # For grains with fewer carbon atoms than nc_cutoff, calculate energies by summing contributions from individual modes
    # Eq. 2 of Draine & Li (2001)
    if nc <= nc_cutoff:
        energy_arr = _calc_pah_energy_modes(
            temp_arr, nc, nh, nm_cc_ip, nm_cc_op, theta_ip_cc, theta_op_cc, theta_ip_ch, theta_op_ch, theta_str_ch
        )

    return energy_arr


def _calc_pah_energy_debye(
    temp_arr, nh, nm_cc_ip, nm_cc_op, theta_ip_cc, theta_op_cc, theta_ip_ch, theta_op_ch, theta_str_ch
):
    """Calculate PAH energy using eq. 33 of Draine & Li (2001).

    Parameters
    ----------
    temp_arr : astropy.units.Quantity (array-like)
        Array of temperature to calculate energies for
    nm_cc_ip : int
        Number of C-C in-plane stretching modes
    nm_cc_op : int
        Number of C-C out-of-plane stretching modes
    nh : int
        Number of hydrogen atoms
    theta_ip_cc : astropy.units.Quantity (float)
        C-C in-plane bending mode Debye temperature
    theta_op_cc : astropy.units.Quantity (float)
        C-C out-of-plane bending mode Debye temperature
    theta_ip_ch : astropy.units.Quantity (float)
        C-H in-plane bending mode Debye temperature
    theta_op_ch : astropy.units.Quantity (float)
        C-H out-of-plane bending mode Debye temperature
    theta_str_ch : astropy.units.Quantity (float)
        C-H stretching mode Debye temperature

    Returns
    -------
    energy_arr : astropy.units.Quantity (array_like)
        Resulting PAH energy array (in u.erg)

    Raises
    ------
    AttributeError
        If the input is not an astropy.units.Quantity object
    TypeError
        If the astropy.units.Quantity object has incorrect units (or optionally is not array-like)
    """
    _check_param(temp_arr, u.K, iterable=True)
    _check_param(theta_ip_cc, u.K)
    _check_param(theta_op_cc, u.K)
    _check_param(theta_ip_ch, u.K)
    _check_param(theta_op_ch, u.K)
    _check_param(theta_str_ch, u.K)

    # contributions from individual C-H modes, as implemented in lines 157-181 of pah_spec_heat.f
    energy_ch = np.zeros(len(temp_arr)) * u.erg

    x = theta_op_ch / temp_arr
    y = exp(x)
    tmin = 32 * u.K
    energy_ch[temp_arr > tmin] += (
        nh * (k_B.cgs * temp_arr[temp_arr > tmin]) * (x[temp_arr > tmin] / (y[temp_arr > tmin] - 1))
    )

    x = theta_ip_ch / temp_arr
    y = exp(x)
    tmin = 42 * u.K
    energy_ch[temp_arr > tmin] += (
        nh * (k_B.cgs * temp_arr[temp_arr > tmin]) * (x[temp_arr > tmin] / (y[temp_arr > tmin] - 1))
    )

    x = theta_str_ch / temp_arr
    y = exp(x)
    tmin = 109 * u.K
    energy_ch[temp_arr > tmin] += (
        nh * (k_B.cgs * temp_arr[temp_arr > tmin]) * (x[temp_arr > tmin] / (y[temp_arr > tmin] - 1))
    )

    # contributions from C-C modes (approximated as Debye spectrum)
    energy_cc_op = np.zeros(len(temp_arr)) * u.erg
    energy_cc_op += _debye_2(theta_op_cc / temp_arr) * k_B.cgs * temp_arr

    energy_cc_ip = np.zeros(len(temp_arr)) * u.erg
    energy_cc_ip += _debye_2(theta_ip_cc / temp_arr) * k_B.cgs * temp_arr

    # total energy
    energy_arr = energy_ch + nm_cc_op * energy_cc_op + nm_cc_ip * energy_cc_ip

    return energy_arr


def _calc_pah_energy_modes(
    temp_arr, nc, nh, nm_cc_ip, nm_cc_op, theta_ip_cc, theta_op_cc, theta_ip_ch, theta_op_ch, theta_str_ch
):
    """Calculate PAH energy using eq. 2 of Draine & Li (2001).

    Parameters
    ----------
    temp_arr : astropy.units.Quantity (array-like)
        Array of temperature to calculate energies for
    nc : int
        Number of carbon atoms
    nh : int
        Number of hydrogen atoms
    nm_cc_ip : int
        Number of C-C in-plane stretching modes
    nm_cc_op : int
        Number of C-C out-of-plane stretching modes
    theta_ip_cc : astropy.units.Quantity (float)
        C-C in-plane bending mode Debye temperature
    theta_op_cc : astropy.units.Quantity (float)
        C-C out-of-plane bending mode Debye temperature
    theta_ip_ch : astropy.units.Quantity (float)
        C-H in-plane bending mode Debye temperature
    theta_op_ch : astropy.units.Quantity (float)
        C-H out-of-plane bending mode Debye temperature
    theta_str_ch : astropy.units.Quantity (float)
        C-H stretching mode Debye temperature

    Returns
    -------
    energy_arr : astropy.units.Quantity (array_like)
        Resulting PAH energy array (in u.erg)

    Raises
    ------
    AttributeError
        If the input is not an astropy.units.Quantity object
    TypeError
        If the astropy.units.Quantity object has incorrect units (or optionally is not array-like)
    """
    _check_param(temp_arr, u.K, iterable=True)
    _check_param(theta_ip_cc, u.K)
    _check_param(theta_op_cc, u.K)
    _check_param(theta_ip_ch, u.K)
    _check_param(theta_op_ch, u.K)
    _check_param(theta_str_ch, u.K)

    mode_arr_cc_op = _calc_cc_mode_energies(nc, nm_cc_op, theta_op_cc)
    mode_arr_cc_ip = _calc_cc_mode_energies(nc, nm_cc_ip, theta_ip_cc)
    mode_arr_ch_op = _calc_ch_mode_energies(nh, theta_op_ch)
    mode_arr_ch_ip = _calc_ch_mode_energies(nh, theta_ip_ch)
    mode_arr_ch_str = _calc_ch_mode_energies(nh, theta_str_ch)

    emode_arr = (
        sorted(
            np.concatenate(
                (mode_arr_cc_op.value, mode_arr_cc_ip.value, mode_arr_ch_op.value, mode_arr_ch_ip.value, mode_arr_ch_str.value)
            )
        )
        * mode_arr_cc_op.unit
    )
    nmodes = len(emode_arr)

    energy_arr = np.zeros(len(temp_arr)) * u.erg

    exp_cutoff = 100  # ignore contributions from exp(x) when x is large

    # emodes array contains mode energies from C-C in-plane and out-of-plane bending and C-H stretching and in-plane
    # in-plane and out-of-plane bending
    for j in range(0, nmodes):
        x = emode_arr[j] / (k_B.cgs * temp_arr)
        y = exp(x)
        energy_arr[x < exp_cutoff] += (x[x < exp_cutoff] / (y[x < exp_cutoff] - 1)) * k_B.cgs * temp_arr[x < exp_cutoff]

    return energy_arr


def _calc_pah_cooling(lambda_abs, grain_radius, wavelength_arr, c_abs_arr, temp_arr, energy_arr):
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

    Returns
    -------
    dt_arr_out : astropy.units.Quantity (array_like)
        Array of time-steps used to solve for T(t) (in u.s)
    time_arr_out : astropy.units.Quantity (array_like)
        Array with time values for T(t) (in u.s)
    temp_arr_out : astropy.units.Quantity (array_like)
        Array with temperature values for T(t) (in u.K)

    Raises
    ------
    AttributeError
        If the input is not an astropy.units.Quantity object
    TypeError
        If the astropy.units.Quantity object has incorrect units (or optionally is not array-like)
    """
    wavelength_unit, radius_unit, c_abs_unit, energy_unit, temp_unit = u.um, u.AA, u.cm**2, u.erg, u.K
    _check_param(lambda_abs, wavelength_unit)
    _check_param(grain_radius, radius_unit)
    _check_param(wavelength_arr, wavelength_unit, iterable=True)
    _check_param(c_abs_arr, c_abs_unit, iterable=True)
    _check_param(energy_arr, energy_unit, iterable=True)
    _check_param(temp_arr, temp_unit, iterable=True)

    nu_arr = c.cgs / wavelength_arr.to(u.cm)

    energy_abs = c.cgs * h.cgs / lambda_abs.to(u.cm)
    temp_abs = np.interp(energy_abs.value, energy_arr.value, temp_arr.value) * u.K

    print(f"Photon wavelength: {lambda_abs:.2f}, initial temperature: {temp_abs:.2f}")

    time_i, energy_i, temp_i = 0 * u.s, energy_abs, temp_abs
    time_arr_out, temp_arr_out, dt_arr_out = [0], [temp_i.value], []
    dt_unit, time_unit, temp_unit = None, None, None
    # ensure timestep is not changing too rapidly
    dE_max = 0.005

    while temp_i.value > 5:
        dE_dt = -trapezoid(4 * np.pi * _planck_function_nu(nu_arr, temp_i) * c_abs_arr, x=nu_arr)

        dt = dE_max * energy_i / dE_dt
        dE = dE_dt * dt

        energy_i -= dE
        time_i += dt

        temp_i = np.interp(energy_i.value, energy_arr.value, temp_arr.value) * u.K

        dt_arr_out.append(dt.value)
        time_arr_out.append(time_i.value)
        temp_arr_out.append(temp_i.value)
        dt_unit, time_unit, temp_unit = dt.unit, time_i.unit, temp_i.unit

    dt_arr_out = np.array(dt_arr_out) * dt_unit
    time_arr_out = np.array(time_arr_out) * time_unit
    temp_arr_out = np.array(temp_arr_out) * temp_unit

    print(f"Final temperature of {temp_arr_out[-1]:.2f} at time {time_arr_out[-1]:.2e}")

    return dt_arr_out, time_arr_out, temp_arr_out


def _calc_basis_vector(wavelength_arr, weighting_arr, temp_arr, c_abs_arr):
    """Calculate the basis vector for a single-photon absorption.

    Parameters
    ----------
    wavelength_arr : astropy.units.Quantity (array_like)
        Array of emission wavelengths
    weighting_arr : array_like
        Array of length len(temp_arr) with values to weight temperatures by
    temp_arr : astropy.units.Quantity (array_like)
        Array with grain temperatures
    c_abs_arr : astropy.units.Quantity (array_like)
        Array of length len(wavelength_arr) with C_abs values for a given grain

    Returns
    -------
    basis_vector : astropy.units.Quantity (array_like)
        Array of floats of length len(wavelength_arr) with basis vector for a given grain (in u.erg / (u.cm * u.s))

    Raises
    ------
    AttributeError
        If the input is not an astropy.units.Quantity object
    TypeError
        If the astropy.units.Quantity object has incorrect units (or optionally is not array-like)
    """
    wavelength_unit, temp_unit, c_abs_unit = u.um, u.K, u.cm**2
    _check_param(wavelength_arr, wavelength_unit, iterable=True)
    _check_param(temp_arr, temp_unit, iterable=True)
    _check_param(c_abs_arr, c_abs_unit, iterable=True)
    _check_param(weighting_arr, u.dimensionless_unscaled)

    unit = None

    basis_vector = np.zeros(len(wavelength_arr))
    for i, lambda_i in enumerate(wavelength_arr.to(u.cm)):
        p_lambda_i = np.sum(4 * np.pi * _planck_function_lambd(lambda_i.to(u.cm), temp_arr) * weighting_arr * c_abs_arr[i])
        unit = p_lambda_i.unit
        basis_vector[i] = p_lambda_i.value

    return basis_vector * unit


def _basis_spectrum(lambda_abs, grain_size, emission_wavelengths, c_abs_arr, temp_arr, energy_arr):
    dt_arr, time_arr, temp_arr_t = _calc_pah_cooling(lambda_abs, grain_size.to(u.AA), emission_wavelengths, c_abs_arr, temp_arr, energy_arr)
    temp_arr_t = temp_arr_t[0:-1]  # make array same length as weighting array
    temp_weights = dt_arr / np.sum(dt_arr)

    return _calc_basis_vector(emission_wavelengths, temp_weights, temp_arr_t, c_abs_arr)


################# Utility functions #################

def _calc_nh(nc):
    """Eq. 8 of Draine & Li (2001)"""
    if nc <= 25:
        return round(0.5 * nc + 0.5)
    if (nc > 25) and (nc <= 100):
        return round(2.5 * np.sqrt(nc) + 0.5)
    if nc > 100:
        return round(0.25 * nc + 0.5)


def _calc_nc(a):
    """Eq. 8 of Draine & Li (2021)"""
    _check_param(a, u.AA)
    nc = 418 * (a / (10 * u.AA)) ** 3
    return nc.value.round(0).astype(int)


def _calc_delta_j(j):
    """Eq. 6 and 5 of Draine & Li (2001). Adjusts location of first 3 modes to bring agreement with mode spectrum of
    coronene C_{24}H_{12}.
    """
    if (j == 2) or (j == 3):
        return 1
    else:
        return 1 / 2


def _calc_beta(nc, nm):
    """Eq. 7 of Draine & Li (2001)"""
    if nc <= 54:
        return 0
    elif (nc > 54) and (nc <= 102):
        return (1 / (2 * nm - 1)) * ((nc - 54) / 52)
    elif nc > 102:
        return (1 / (2 * nm - 1)) * (((nc - 2) / 52) * (102 / nc) ** (2 / 3) - 1)


def _calc_cc_mode_energies(nc, nm, theta):
    """As implemented in lines 141-154 and 158-165 of pah_modes.f"""
    _check_param(theta, u.K)
    mode_energy_arr = np.zeros((nm)) * u.erg
    # beta = calc_beta(nc, nc - 2)  # note that we are intentionally using the same beta for C-C ip and op modes
    beta = _calc_beta(nc, nm)  # note that we are intentionally using the same beta for C-C ip and op modes

    for j in range(1, nm + 1):
        mode_energy_arr[j - 1] = k_B.cgs * theta * np.sqrt((1 - beta) * (j - _calc_delta_j(j)) / nm + beta)

    return mode_energy_arr


def _calc_ch_mode_energies(nh, theta):
    """As implemented in lines 170-184 of pah_modes.f"""
    _check_param(theta, u.K)
    mode_energy_arr = np.zeros((nh)) * u.erg
    mode_energy_arr[:] = k_B.cgs * theta

    return mode_energy_arr


def _debye_2(x):
    """As implemented in lines 28-29 of debye.f"""
    x2 = x * x
    return (2 / x2) * (2 * 1.20206 - x2 * exp(-x) * (1 + 2 / x + 2 / x2 + exp(-x) * (0.5 + 0.5 / x + 0.25 / x2)))


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
    AttributeError
        If the input is not an astropy.units.Quantity object
    TypeError
        If the astropy.units.Quantity object has incorrect units (or optionally is not array-like)
    """
    try:
        param.value
    except AttributeError:
        raise AttributeError(f"expects astropy.units.Quantity objects of unit {unit}")

    if param.unit != unit:
        raise TypeError(f"parameter expects units of {unit}")

    if iterable:
        if not isinstance(param.value, (list, tuple, np.ndarray)):
            raise TypeError("expects array-like input")
    
    if force_iterable:
        if not isinstance(param.value, (list, tuple, np.ndarray)):
            param = np.array([param.value]) * param.unit
        return param


def _planck_function_nu(nu, T):
    _check_param(nu, u.Hz)
    _check_param(T, u.K)
    return (2 * h.cgs * nu**3 / c.cgs**2) * 1 / (np.exp(h.cgs * nu / (k_B.cgs * T)) - 1)


def _planck_function_lambd(lambd, T):
    _check_param(lambd, u.cm)
    _check_param(T, u.K)
    return (2 * h.cgs * c.cgs**2 / lambd**5) * 1 / (np.exp(h.cgs * c.cgs / (lambd * k_B.cgs * T)) - 1)

# TODO: decide if this needs to be in here
def convert_temp_to_prob(dt_arr, temp_arr, temp_cutoff):
    dT = abs(np.diff(temp_arr))

    temp_avg = []
    for i, temp in enumerate(temp_arr):
        if (i + 1) < len(temp_arr):
            temp_avg.append((temp_arr[i].value + temp_arr[i + 1].value) / 2)

    temp_avg *= u.K

    temp_avg_sub = temp_avg[temp_avg.value >= temp_cutoff.value]
    print(temp_avg_sub[0], temp_avg_sub[-1])

    dT_sub = dT[temp_avg.value >= temp_cutoff.value]
    dt_arr_sub = dt_arr[temp_avg.value >= temp_cutoff.value]

    area = -trapezoid(dt_arr_sub / dT_sub, x=temp_avg_sub)

    return temp_avg_sub, (dt_arr_sub / dT_sub) / area


def _read_basis_spectra_for_ionization(filename):
    df = pd.read_csv(filename)
    # TODO: move photon wavelengths and emission wavelengths into their own file
    photon_wavelengths = df.columns[1:].to_numpy(dtype=float) * u.um
    emission_wavelengths = df["emission_wavelengths"].to_numpy(dtype=float) * u.um
    
    basis_spectra_i = np.zeros((len(photon_wavelengths), len(emission_wavelengths)))
    for j, col in enumerate(df):
        if df.columns[j] != "emission_wavelengths":
            basis_spectra_i[j-1] = df[col].to_numpy()
    basis_spectra_i *= u.erg / (u.s * u.cm)

    return basis_spectra_i, photon_wavelengths, emission_wavelengths


def _read_basis_spectra(basis_directory, grain_sizes):

    photon_wavelengths, emission_wavelengths  = None, None
    basis_spectra_ion, basis_spectra_neu = [], []
    for grain_size in grain_sizes:

        filename_neu = os.path.join(basis_directory, f"neu/basis_{grain_size.to(u.AA).value:.2f}_neu.csv")
        filename_ion = os.path.join(basis_directory, f"ion/basis_{grain_size.to(u.AA).value:.2f}.csv")

        basis_spectra_neu_i, photon_wavelengths, emission_wavelengths = _read_basis_spectra_for_ionization(filename_neu)
        basis_spectra_ion_i, _, _ = _read_basis_spectra_for_ionization(filename_ion)
        
        basis_spectra_neu.append(basis_spectra_neu_i)
        basis_spectra_ion.append(basis_spectra_ion_i)

    return emission_wavelengths, photon_wavelengths, basis_spectra_neu, basis_spectra_ion
