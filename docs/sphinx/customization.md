# Custom PAH Physics

## Overview

The ``pah_spec`` package relies on a set of pre-computed basis spectra $\tilde{p}_{\lambda_{\rm em}}(\lambda_{\rm abs})$ (Eq. 3 of [Richie & Hensley 2026](https://ui.adsabs.harvard.edu/abs/2025arXiv251016861R/abstract)) to generate PAH spectra. Defining these $\tilde{p}_{\lambda_{\rm em}}(\lambda_{\rm abs})$ requires the specification of PAH energy and absorption cross-section ($C_{\rm abs}$) models. We provide a set of default, pre-computed $\tilde{p}_{\lambda_{\rm em}}(\lambda_{\rm abs})$ using the PAH physics model described in Section 3 of [Richie & Hensley 2026](https://ui.adsabs.harvard.edu/abs/2025arXiv251016861R/abstract) (based on the [Draine et al. 2021](https://ui.adsabs.harvard.edu/abs/2021ApJ...917....3D/abstract) model). 

However, since the single photon approximation (SPA) is agnostic to the choice of PAH physics model, one may wish to employ a different PAH model to generate SPA spectra. To this end, we provide the capability to run ``pah_spec`` with custom PAH physics models. The steps to do this are:

1. Define Python functions encapsulating the desired energy and/or cross-section models
2. Compute a set of basis spectra using the custom energy and/or cross-section model
3. Employ the {py:class}``pah_spec.PahSpec`` API with the custom basis spectra

The following instructions show how to modify the cross-section and energy functions simultaneously, but note that it is also possible to modify these functions one at a time.

## 1. Defining Compatible PAH Physics Python Functions

Up to three Python functions are required as inputs to {py:class}``pah_spec.PahSpec`` to utilize custom PAH physics.

**Required function for custom cross-sections:**
1. ``c_abs_func(wavelength_arr, radius_arr)``

This function takes an array of PAH radii and wavelengths as inputs and returns the absorption cross-sections $C_{\rm abs}$ for ionized and neutral PAHs. This function **MUST** have the same signature as {py:meth}``pah_spec.PahSpec.calc_c_abs``. 

Specifically, ``c_abs_func()`` must take an ``astropy.units.Quantity`` array of PAH radii (in units of Angstroms) and an ``astropy.units.Quantity`` array of wavelengths (in $\rm \mu m$) as inputs and return two ``astropy.units.Quantity`` arrays of dimensions ``(len(grain_sizes), len(wavelength_arr))`` of cross-sections (in $\rm cm^2$) for ionized and neutral PAHs, respectively.

**Required functions for custom energy:**
1. ``pah_energy_func(grain_size, temp_arr)``
2. ``pah_mode_energy_func(grain_size)``

``pah_energy_func()`` **MUST** have the same signature as {py:func}``pah_spec.calc_pah_energy``, which takes a single ``astropy.units.Quantity`` grain radius (in Angstroms) and an ``astropy.units.Quantity`` array of temperatures (in Kelvin) as inputs and returns an ``astropy.units.Quantity`` array of PAH vibrational energies of length ``len(temp_arr)`` (in ergs).

``pah_mode_energy_func()`` **MUST** have the same signature as {py:func}``pah_spec.calc_pah_mode_energies``, which calculates the normal vibrational mode energies of the PAH that arise from $\rm C$--$\rm C$ and $\rm C$--$\rm H$ in-plane and out-of-plane bending and $\rm C$--$\rm H$ stretching. It takes a single ``astropy.units.Quantity`` PAH radius (in Angstroms) and returns an ``astropy.units.Quantity`` array of vibrational mode energies (in ergs) sorted in order of increasing energy.

## 2. Computing Custom Basis Spectra

Once defined, the above functions can be passed to the {py:class}``pah_spec.PahSpec`` API to compute a set of custom basis spectra. This can be done follows:

```
ps = pah_spec.PahSpec(c_abs_func=custom_c_abs_func, pah_energy_func=custom_pah_energy_func)

ps.generate_basis_spectra(ps.grain_sizes, ps.emission_wavelengths, ion=True, pah_mode_energy_func=custom_pah_mode_energy_func)
ps.generate_basis_spectra(ps.grain_sizes, ps.emission_wavelengths, ion=False, pah_mode_energy_func=custom_pah_mode_energy_func)
```

The operative step is the specification of the ``c_abs_func`` and ``pah_energy_func`` arguments in defining the the {py:class}``pah_spec.PahSpec`` and the ``pah_mode_energy_func`` argument in {py:meth}``pah_spec.PahSpec.generate_basis_spectra``. See our [Examples](https://github.com/helenarichie/pah_spec/tree/main/examples) for instructions on using {py:meth}``pah_spec.PahSpec.generate_basis_spectra``. 

Once this function finishes running (it will take several hours to run), new basis_ion.h5 and basis_neu.h5 data files will be generated.

## 3. Generating PAH Spectra With Custom PAH Physics

Once the new basis_ion.h5 and basis_neu.h5 data files have been created, they can be used to generate PAH spectra. To do this, you must define a new {py:class}``pah_spec.PahSpec`` object using the new basis spectra. This can be done by specifying their location with the ``basis_dir`` argument when deifning the new {py:class}``pah_spec.PahSpec`` object:

```
ps = pah_spec.PahSpec(c_abs_func=custom_c_abs_func, pah_energy_func=custom_pah_energy_func, basis_dir="/path/to/custom/basis/spectra")
```

Note that you **MUST** use the same energy and/or cross-section functions used to compute the basis spectra when defining the new {py:class}``pah_spec.PahSpec`` object. 

After this step, you can use {py:meth}``pah_spec.PahSpec.generate_spectrum`` to compute PAH spectra using your own PAH physics!