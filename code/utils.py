import astropy.units as u
import numpy as np
import os
import pandas as pd
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__))))
from pah_spec_routines import check_param, calc_cabs, calc_pah_energy, calc_pah_cooling, calc_basis_vector


def read_basis_vectors(filename):

    df = pd.read_csv(filename)

    photon_wavelenghts = df.columns[1:].to_numpy(dtype=float) * u.um
    emission_wavelengths = df["emission_wavelengths"].to_numpy(dtype=float) * u.um

    basis_vectors = np.zeros((len(photon_wavelenghts), len(emission_wavelengths)))
    for i, col in enumerate(df):
        if df.columns[i] != "emission_wavelengths":
            basis_vectors[i-1] = df[col].to_numpy()
    basis_vectors *= u.erg / (u.s * u.cm)

    return emission_wavelengths, photon_wavelenghts, basis_vectors


def write_basis_vectors(emission_wavelengths, photon_wavelengths, basis_vectors, data_path, filename):
    check_param(emission_wavelengths, u.um, iterable=True)
    check_param(photon_wavelengths, u.um, iterable=True)
    check_param(basis_vectors, u.erg / (u.cm * u.s), iterable=True)

    photon_wavelengths = [s.value for s in photon_wavelengths]
    emission_wavelengths = [s.value for s in emission_wavelengths]

    data = dict()
    data["emission_wavelengths"] = emission_wavelengths
    for i, vec in enumerate(basis_vectors):
        data[str(photon_wavelengths[i])] = vec
    df = pd.DataFrame(data)
    df.to_csv(os.path.join(data_path, filename), index=False)


def generate_photon_wavelengths(dlambda, lambda_min, lambda_max):
    check_param(lambda_min, u.um)
    check_param(lambda_max, u.um)

    # define photon/basis vector wavelengths
    lambda_i = lambda_min.value
    photon_wavelengths = [lambda_i]
    while lambda_i < lambda_max.value:
        lambda_i += dlambda * lambda_i
        photon_wavelengths.append(lambda_i)
    photon_wavelengths *= u.um

    return photon_wavelengths


def generate_basis_vectors(grain_radius, dlambda, lambda_min, lambda_max, emission_wavelengths, cabs_arr=None):
    check_param(grain_radius, u.AA)
    check_param(lambda_min, u.um)
    check_param(lambda_max, u.um)
    check_param(emission_wavelengths, u.um, iterable=True)

    photon_wavelengths = generate_photon_wavelengths(dlambda, lambda_min, lambda_max)

    if cabs_arr == None:
        cabs_arr = calc_cabs(emission_wavelengths, grain_radius)[0][0]
    
    temp_arr = np.linspace(1, 1e4, 1000) * u.K
    energy_arr = calc_pah_energy(grain_radius[0], temp_arr)

    basis_vectors = np.zeros((len(photon_wavelengths), len(emission_wavelengths))) * u.erg / (u.cm * u.s)

    for i, lambda_abs in enumerate(photon_wavelengths):

        dt_arr, time_arr, temp_arr_t = calc_pah_cooling(lambda_abs, grain_radius.to(u.AA), emission_wavelengths, cabs_arr, temp_arr, energy_arr)
        temp_arr_t = temp_arr_t[0:-1]  # make array same length as weighting array
        temp_weights = dt_arr / np.sum(dt_arr)

        basis_vectors[i] = calc_basis_vector(emission_wavelengths, temp_weights, temp_arr_t, cabs_arr)

    return emission_wavelengths, photon_wavelengths, basis_vectors


def write_basis_vectors(emission_wavelengths, photon_wavelengths, basis_vectors, data_path, filename):
    # TODO: enforce naming style for basis vector file
    check_param(emission_wavelengths, u.um, iterable=True)
    check_param(photon_wavelengths, u.um, iterable=True)
    check_param(basis_vectors, u.erg / (u.cm * u.s), iterable=True)

    photon_wavelengths = [s.value for s in photon_wavelengths]
    emission_wavelengths = [s.value for s in emission_wavelengths]

    data = dict()
    data["emission_wavelengths"] = emission_wavelengths
    for i, vec in enumerate(basis_vectors):
        data[str(photon_wavelengths[i])] = vec
    df = pd.DataFrame(data)
    df.to_csv(os.path.join(data_path, filename), index=False)


def read_basis_vectors(filename):

    df = pd.read_csv(filename)

    photon_wavelenghts = df.columns[1:].to_numpy(dtype=float) * u.um
    emission_wavelengths = df["emission_wavelengths"].to_numpy(dtype=float) * u.um

    basis_vectors = np.zeros((len(photon_wavelenghts), len(emission_wavelengths)))
    for i, col in enumerate(df):
        if df.columns[i] != "emission_wavelengths":
            basis_vectors[i-1] = df[col].to_numpy()
    basis_vectors *= u.erg / (u.s * u.cm)

    return emission_wavelengths, photon_wavelenghts, basis_vectors