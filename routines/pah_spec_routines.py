from astropy.constants import c, h, k_B
import astropy.units as u
from mpmath import polylog, zeta, matrix
import numpy as np
from numpy import exp, log
import os
from scipy.integrate import trapezoid
import time


def calc_cabs(wavelength_arr, radius_arr):
    """Calculate the absorption cross-section, C_abs, for input grain sizes and wavelengths based on method from
    Draine et al. (2021).

    Parameters
    ----------
    wavelength_arr : astropy.units.Quantity (array_like)
        Array of wavelengths to calculate C_abs for
    radius_arr : astropy.units.Quantity (array_like)
        Array of dust grain radii to calculate C_abs for

    Returns
    -------
    cabs_ion_out : astropy.units.Quantity (array_like)
        Array with C_abs values for ionized grains (in u.cm ** 2)
    cabs_neu_out : astropy.units.Quantity (array_like)
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
    check_param(wavelength_arr, wavelength_unit, iterable=True)
    check_param(radius_arr, radius_unit, iterable=True)

    script_directory = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(script_directory, "../data/")

    if not os.path.exists(data_path):
        raise FileNotFoundError(
            f"calc_cabs expects to find qabs_001um.dat and draine21_Table4.dat at {data_path}, but the path does not exist"
        )

    (wav_graphite, qabs) = np.genfromtxt(os.path.join(data_path, "qabs_001um.dat"), unpack=True, usecols=[0, 2])
    wav_graphite *= u.um

    # C_abs per volume doesn't depend on size in the small-grain limit
    cabs_V_graphite = 3.0 * qabs / (4 * 1.0e-3 * u.um)
    cabs_V_graphite_out = np.interp(
        wavelength_arr.to(u.um),
        wav_graphite[::-1].to(u.um),
        cabs_V_graphite[::-1].to(u.cm**-1),
    )

    rho_C = 2.0 * u.g / u.cm**3
    N_C = calc_N_C(radius_arr)  # number of Carbon atoms, Draine 2021
    HC = 0.5 * np.ones_like(N_C)  # Note that this is from Equation 4 of DL07
    HC[(N_C >= 25) & (N_C <= 100)] = 0.5 * np.sqrt(25 / N_C[(N_C >= 25) & (N_C <= 100)])
    HC[N_C > 100] = 0.25
    xi_gra = np.zeros_like(radius_arr.value)
    xi_gra[radius_arr < 50 * u.AA] = 0.01
    xi_gra[radius_arr >= 50 * u.AA] = 0.01 + 0.99 * (1.0 - (50 * u.AA / radius_arr[radius_arr >= 50 * u.AA]) ** 3)

    (lamj_tab, gamj_tab, sigj_neu_tab, sigj_ion_tab, hc_tab) = np.genfromtxt(
        os.path.join(data_path, "draine21_Table4.dat"), unpack=True
    )
    lamj_tab *= u.um
    sigj_neu_tab *= 1.0e-20
    sigj_ion_tab *= 1.0e-20
    sigj_neu_tab *= u.cm
    sigj_ion_tab *= u.cm

    def S_func(lam, lamj, gamma, sigma):
        num = 2 * gamma * lamj * sigma
        denom = np.pi * (((lam / lamj) - (lamj / lam)) ** 2 + gamma**2)
        return num / denom

    def C_func(y):
        return np.arctan(1.0e3 * (y - 1.0) ** 3 / y) / np.pi + 0.5

    # Cutoff function parameters
    M_ring = 0.3 * N_C
    M_ring[N_C > 40] = 0.4 * N_C[N_C > 40]
    lamc_neu = 0.951 * u.um / (1.0 + 3.616 / np.sqrt(M_ring))
    lamc_ion = 1.125 * u.um / (1.0 + 2.567 / np.sqrt(M_ring))

    x = (1.0 * u.um / wavelength_arr.to(u.um)).value
    cabs_ion_out = np.zeros((len(radius_arr), len(wavelength_arr))) * u.cm**2
    cabs_neu_out = np.zeros((len(radius_arr), len(wavelength_arr))) * u.cm**2

    for i, _ in enumerate(radius_arr):
        # Graphite contribtuion
        volume = (4.0 / 3.0) * np.pi * radius_arr[i] ** 3
        cabs_ion_out[i, :] += xi_gra[i] * cabs_V_graphite_out * volume
        cabs_neu_out[i, :] += xi_gra[i] * cabs_V_graphite_out * volume

        # PAH contribution
        C_fac_neu = C_func(
            (wavelength_arr.to(u.um) / lamc_neu[i]).value ** -1
        )  # Note error in Equation A7, should be C(lambda_c/lambda)
        C_fac_ion = C_func(
            (wavelength_arr.to(u.um) / lamc_ion[i]).value ** -1
        )  # Note error in Equation A7, should be C(lambda_c/lambda)
        S_mat_neu = np.zeros((len(lamj_tab), len(wavelength_arr))) * u.cm**2
        S_mat_ion = np.zeros((len(lamj_tab), len(wavelength_arr))) * u.cm**2
        for j, _ in enumerate(lamj_tab):
            if hc_tab[j] == 0:
                S_mat_neu[j, :] = S_func(wavelength_arr.to(u.um), lamj_tab[j], gamj_tab[j], sigj_neu_tab[j])
                S_mat_ion[j, :] = S_func(wavelength_arr.to(u.um), lamj_tab[j], gamj_tab[j], sigj_ion_tab[j])
            else:
                S_mat_neu[j, :] = S_func(wavelength_arr.to(u.um), lamj_tab[j], gamj_tab[j], sigj_neu_tab[j] * HC[i])
                S_mat_ion[j, :] = S_func(wavelength_arr.to(u.um), lamj_tab[j], gamj_tab[j], sigj_ion_tab[j] * HC[i])

        idx = np.where((10 < x) & (x < 15))
        cabs_ion_out[i, idx] += (
            (S_mat_ion[0, idx] + (1.35 * x[idx] - 3.0) * 1.0e-18 * u.cm**2) * N_C[i] * (1.0 - xi_gra[i])
        )
        cabs_neu_out[i, idx] += (
            (S_mat_neu[0, idx] + (1.35 * x[idx] - 3.0) * 1.0e-18 * u.cm**2) * N_C[i] * (1.0 - xi_gra[i])
        )

        idx = np.where((7.7 < x) & (x <= 10))
        cabs_ion_out[i, idx] += (
            (66.302 - 24.367 * x[idx] + 2.950 * x[idx] ** 2 - 0.1057 * x[idx] ** 3)
            * 1.0e-18
            * u.cm**2
            * N_C[i]
            * (1.0 - xi_gra[i])
        )
        cabs_neu_out[i, idx] += (
            (66.302 - 24.367 * x[idx] + 2.950 * x[idx] ** 2 - 0.1057 * x[idx] ** 3)
            * 1.0e-18
            * u.cm**2
            * N_C[i]
            * (1.0 - xi_gra[i])
        )

        idx = np.where((5.9 < x) & (x <= 7.7))
        c0 = 1.8687e-18 * u.cm**2
        c1 = 1.905e-19 * u.cm**2
        c2 = 4.175e-19 * u.cm**2
        c3 = 4.37e-20 * u.cm**2  # Note Equations A11 and A12 are mislabeled
        cabs_ion_out[i, idx] += (
            (S_mat_ion[1, idx] + c0 + c1 * x[idx] + c2 * (x[idx] - 5.9) ** 2 + c3 * (x[idx] - 5.9) ** 3)
            * N_C[i]
            * (1.0 - xi_gra[i])
        )
        cabs_neu_out[i, idx] += (
            (S_mat_neu[1, idx] + c0 + c1 * x[idx] + c2 * (x[idx] - 5.9) ** 2 + c3 * (x[idx] - 5.9) ** 3)
            * N_C[i]
            * (1.0 - xi_gra[i])
        )

        idx = np.where((3.3 < x) & (x <= 5.9))
        cabs_ion_out[i, idx] += (S_mat_ion[1, idx] + c0 + c1 * x[idx]) * N_C[i] * (1.0 - xi_gra[i])
        cabs_neu_out[i, idx] += (S_mat_neu[1, idx] + c0 + c1 * x[idx]) * N_C[i] * (1.0 - xi_gra[i])

        idx = np.where(x <= 3.3)
        cabs_ion_out[i, idx] += (
            34.58e-18 * 10 ** (-3.431 / x[idx]) * u.cm**2 * C_fac_ion[idx] * N_C[i] * (1.0 - xi_gra[i])
        )
        cabs_neu_out[i, idx] += (
            34.58e-18 * 10 ** (-3.431 / x[idx]) * u.cm**2 * C_fac_neu[idx] * N_C[i] * (1.0 - xi_gra[i])
        )

        for j in range(len(lamj_tab)):
            if j > 1:
                cabs_ion_out[i, idx] += S_mat_ion[j, idx] * N_C[i] * (1.0 - xi_gra[i])
                cabs_neu_out[i, idx] += S_mat_neu[j, idx] * N_C[i] * (1.0 - xi_gra[i])

        cabs_ion_out[i, idx] += (
            3.5e-19 * 10 ** (-1.45 / x[idx]) * np.exp(-((0.1 * x[idx]) ** 2)) * N_C[i] * (1.0 - xi_gra[i]) * u.cm**2
        )  # Note: does not appear in D21+

    return cabs_ion_out, cabs_neu_out


def calc_pah_energy(grain_radius, temp_arr):
    """Calculate PAH vibrational energy as a function of temperature according to Eq. 33 of Draine & Li (2001).

    Parameters
    ----------
    grain_radius : astropy.units.Quantity (float)
        The PAH effective radius
    temp_arr : astropy.units.Quantity (array-like)
        Array of temperature to calculate energies for

    Returns
    -------
    energies : astropy.units.Quantity (array_like)
        Resulting PAH energy array (in u.ergs)

    Raises
    ------
    AttributeError
        If the input is not an astropy.units.Quantity object
    TypeError
        If the astropy.units.Quantity object has incorrect units (or optionally is not array-like)
    """
    radius_unit, temp_unit = u.AA, u.K

    check_param(grain_radius, radius_unit)
    check_param(temp_arr, temp_unit, iterable=True)

    # 5 types of vibration: out-of-plane C-C modes, in-plane C-C modes, out-of-plane C-H bending, in-plane C-H bending, and C-H stretching.
    N_C = calc_N_C([grain_radius.value] * radius_unit)[0]
    N_H = calc_N_H(N_C)
    theta_op = 863 * u.K  # C-C out-of-plane mode Debye temperature
    theta_ip = 2504 * u.K  # C-C in-plane mode Debye temperature
    EMCH_OP = 886 * (1 / u.cm)  # cm-1
    EMCH_IP = 1161 * (1 / u.cm)
    EMCH_ST = 3030 * (1 / u.cm)

    # contributions from C-H in-plane bending, out-of-plane bending, and stretching
    energy_CH = N_H * (
        EMCH_IP / (np.expm1(h.cgs * c.cgs * EMCH_IP / (k_B.cgs * temp_arr)) - 1 + 1) * h.cgs * c.cgs
        + EMCH_OP / (np.expm1(h.cgs * c.cgs * EMCH_OP / (k_B.cgs * temp_arr)) - 1 + 1) * h.cgs * c.cgs
        + EMCH_ST / (np.expm1(h.cgs * c.cgs * EMCH_ST / (k_B.cgs * temp_arr)) - 1 + 1) * h.cgs * c.cgs
    )
    # contributions from C-C modes (Debye spectrum)
    energy_CC = (
        (N_C - 2)
        * k_B.cgs
        * (theta_op * f2(temp_arr.value / theta_op.value) + 2 * theta_ip * f2(temp_arr.value / theta_ip.value))
    )

    # total energy
    energies = energy_CH + energy_CC

    return energies


def calc_pah_cooling(lambda_abs, grain_radius, wavelength_arr, cabs_arr, temp_arr, energy_arr):
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
    check_param(lambda_abs, wavelength_unit)
    check_param(grain_radius, radius_unit)
    check_param(wavelength_arr, wavelength_unit, iterable=True)
    check_param(cabs_arr, c_abs_unit, iterable=True)
    check_param(energy_arr, energy_unit, iterable=True)
    check_param(temp_arr, temp_unit, iterable=True)

    nu_arr = c.cgs / wavelength_arr.to(u.cm)

    energy_abs = c.cgs * h.cgs / lambda_abs.to(u.cm)
    temp_abs = np.interp(energy_abs.value, energy_arr.value, temp_arr.value) * u.K

    print(f"Photon wavelength: {lambda_abs:.2f}, initial temperature: {temp_abs:.2f}")

    time_i, energy_i, temp_i = 0 * u.s, energy_abs, temp_abs
    time_arr_out, temp_arr_out, dt_arr_out = [0], [temp_i.value], []
    dt_unit, time_unit, temp_unit = None, None, None
    # ensure timestep is not changing too rapidly
    dE_max = 0.001

    while temp_i.value > 5:

        dE_dt = -trapezoid(4 * np.pi * planck_function_nu(nu_arr, temp_i) * cabs_arr, x=nu_arr)

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


def calc_eigenvector(wavelength_arr, weighting, temp_arr, c_abs_arr):
    """Calculate the eigenvector for a single-photon absorption.

    Parameters
    ----------
    wavelength_arr : astropy.units.Quantity (array_like)
        Array of emission wavelengths
    weighting : array_like
        Array of length len(temp_arr) with values to weight temperatures by
    temp_arr : astropy.units.Quantity (array_like)
        Array with grain temperatures
    c_abs_arr : astropy.units.Quantity (array_like)
        Array of length len(wavelength_arr) with C_abs values for a given grain

    Returns
    -------
    eigenvector : astropy.units.Quantity (array_like)
        Array of floats of length len(wavelength_arr) with eigenvector for a given grain (in u.erg / (u.cm * u.s))

    Raises
    ------
    AttributeError
        If the input is not an astropy.units.Quantity object
    TypeError
        If the astropy.units.Quantity object has incorrect units (or optionally is not array-like)
    """
    wavelength_unit, temp_unit, c_abs_unit = u.um, u.K, u.cm**2
    check_param(wavelength_arr, wavelength_unit, iterable=True)
    check_param(temp_arr, temp_unit, iterable=True)
    check_param(c_abs_arr, c_abs_unit, iterable=True)
    # TODO: add check that weighting array is dimensionless

    unit = None

    eigenvector = np.zeros(len(wavelength_arr))
    for i, lambda_i in enumerate(wavelength_arr.to(u.cm)):
        p_lambda_i = np.sum(4 * np.pi * planck_function_lambd(lambda_i.to(u.cm), temp_arr) * weighting * c_abs_arr[i])
        unit = p_lambda_i.unit
        eigenvector[i] = p_lambda_i.value

    return eigenvector * unit


def calc_normalization(lambda_abs, mrf_width, wavelength_arr, c_abs_arr, wavelengths_u, u_lambda, p_lambda):
    """Calculate the energy conservation normalization to scale an eigenvector to an input radiation field.

    Parameters
    ----------
    lambda_abs : astropy.units.Quantity (float)
        Wavelength of absorbed photon
    mrf_width : float
        Width of the "monochromatic" radiation field, defined as a percentage of lambda_abs
    wavelength_arr : astropy.units.Quantity (array_like)
        Array of emission wavelengths
    c_abs_arr : astropy.units.Quantity (array_like)
        Array of length len(wavelength_arr) with C_abs values for a given grain
    wavelengths_u : astropy.units.Quantity (array_like)
        Wavelength array for the radiation field u_lambda
    u_lambda : astropy.units.Quantity (array_like)
        Array of length len(wavelengths_u) with the radiation field
    p_lambda : astropy.units.Quantity (array_like)
        Eigenvector array of length len(wavelength_arr) for a given grain and lambda_abs

    Returns
    -------
    normalization : float
        Incident photon power / grain's radiated power (dimensionless)

    Raises
    ------
    AttributeError
        If the input is not an astropy.units.Quantity object
    TypeError
        If the astropy.units.Quantity object has incorrect units (or optionally is not array-like)
    """
    wavelength_unit, c_abs_unit, radiation_field_unit, eigenvector_unit = (
        u.um,
        u.cm**2,
        u.erg / (u.cm**4),
        u.erg / (u.cm * u.s),
    )
    check_param(lambda_abs, wavelength_unit)
    check_param(wavelength_arr, wavelength_unit, iterable=True)
    check_param(c_abs_arr, c_abs_unit, iterable=True)
    check_param(wavelengths_u, wavelength_unit, iterable=True)
    check_param(u_lambda, radiation_field_unit, iterable=True)
    check_param(p_lambda, eigenvector_unit, iterable=True)

    # TODO: test whether we should use more than two-wavelengths to integrate the MRF

    # define monochromatic radiation field wavlenegth range
    # wh_wav0, wh_wav1 = np.argmin(np.abs(wavelength_arr - wav0)), np.argmin(np.abs(wavelength_arr - wav1))
    # if wh_wav0 == wh_wav1:
    #     raise ValueError(f"emission wavelength binning is too coarse for input mrf_width of {mrf_width}")

    # interpolate radiation field to have the same wavelength resolution as wavelength_arr
    # also define C_abs array over monochromatic radiation field wavelength range
    # wav_mrf, c_abs_mrf = wavelength_arr[wh_wav0 : wh_wav1 + 1], c_abs_arr[wh_wav0 : wh_wav1 + 1]

    wav0, wav1 = lambda_abs, lambda_abs + lambda_abs * mrf_width
    wav_mrf = np.array([wav0.value, wav1.value]) * lambda_abs.unit

    c_abs_mrf = np.interp([wav0, wav1], wavelength_arr, c_abs_arr)

    u_lambda_mrf = np.interp(wav_mrf, wavelengths_u, u_lambda)

    numerator = trapezoid(u_lambda_mrf * c_abs_mrf * c.cgs, x=wav_mrf.to(u.cm))

    denominator = 4 * np.pi * trapezoid(p_lambda, x=wavelength_arr.to(u.cm))

    return numerator, denominator


################# Utility functions #################


def calc_N_H(N_C):
    """Eq. 8 of Draine & Li (2001)"""

    N_H = 0
    if N_C <= 25:
        N_H = int(0.5 * N_C + 0.5)
    if (N_C > 25) and (N_C <= 100):
        N_H = int(2.5 * np.sqrt(N_C) + 0.5)
    if N_C > 100:
        N_H = int(0.25 * N_C + 0.5)
    return N_H


def calc_N_C(a):
    """Eq. 8 of Draine & Li (2021)"""

    check_param(a, u.AA, iterable=True)
    N_C = 418 * (a / (10 * u.AA)) ** 3
    return N_C.astype(int)


def f2(x):
    # From Eq. 10 of Draine & Li (2001), note corrected prefactor of f_n(x) = n * (integral)
    polylog_vectorized = np.frompyfunc(polylog, 2, 1)

    polylog_mpf_3 = matrix(polylog_vectorized(3, np.expm1(-1 / x) + 1))
    polylog_array_3 = np.array(polylog_mpf_3.tolist(), dtype="float64").flatten()

    polylog_mpf_2 = matrix(polylog_vectorized(2, np.expm1(-1 / x) + 1))
    polylog_array_2 = np.array(polylog_mpf_2.tolist(), dtype="float64").flatten()

    return (2 * x) * (
        log(1 - np.expm1(-1 / x) - 1) - 2 * x * (polylog_array_2 + x * (polylog_array_3 - float(zeta(3))))
    )


def check_param(param, unit, iterable=False):
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
            raise TypeError("calc_cabs expects array-like input for wavelength")


def planck_function_nu(nu, T):
    check_param(nu, u.Hz)
    check_param(T, u.K)
    return (2 * h.cgs * nu**3 / c.cgs**2) * 1 / (np.exp(h.cgs * nu / (k_B.cgs * T)) - 1)


def planck_function_lambd(lambd, T):
    check_param(lambd, u.cm)
    check_param(T, u.K)
    return (2 * h.cgs * c.cgs**2 / lambd**5) * 1 / (np.exp(h.cgs * c.cgs / (lambd * k_B.cgs * T)) - 1)
