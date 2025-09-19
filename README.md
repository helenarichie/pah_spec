# Code Documentation

## Dependencies: astropy, numpy, pandas, scipy
Tested with python 3.13.7, astropy 6.1.2, numpy 2.3.2, pandas 2.3.2, scipy 1.16.1


### `calc_cabs(wavelength_arr, radius_arr)`

Calculates $C_{\rm abs}$ following the method described in [Draine et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021ApJ...917....3D/abstract)
and [Draine et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025ApJ...989..232D/abstract) for neutral and ionized PAHs for arbitrary wavelengths and grain sizes, assuming $a\ll \lambda$.

### `calc_pah_energy(grain_radius, temp_arr)`

Calculates the PAH energy function, $E(T)$, for an input temperature range and grain size following the method described in Section 2 of Richie & Hensley (2025). 
Depending on the input effective radius, $a$, `calc_pah_energy()` calls one of two sub-routines: `calc_pah_energy_modes()` for $N_{\rm C}(a)\leq7360$ (i.e. Equation 10 of Richie & Hensley 2025) and `calc_pah_energy_debye()` for $N_{\rm C}(a)>7360$ (i.e. Equation 14 fo Richie & Hensley 2025).

### `calc_pah_cooling(lambda_abs, grain_radius, wavelength_arr, c_abs_arr, temp_arr, energy_arr)`

Solves for the PAH cooling function, $T(t)$, given an input $\lambda_{\rm abs}$ and grain size. 
The routine does this by first integrating Equation 2 of Richie & Hensley (2025) over the input `wavelength_arr` (which must span a wavelength range sufficient to capture the emission wavelength range of PAHs, i.e., $\sim0.1-3000~{\mu m}$) to determine $dE/dt$, assuming that the initial temperature is given by the temperature corresponding to $E_{\rm abs}(T) = hc/\lambda_{\rm abs}$. 
It then applies this $dE/dt$ to $E_{i+1}=E_\mathrm{i}+\left(\frac{dE}{dt}\right)_idt$, with $dt$ chosen such that $dE$ is 0.001\% of the grain's current energy, ensuring that the energy does not change rapidly. 
The updated grain energy is then used to determine its current temperature, and this loop repeats until the grain has cooled to $5~\textrm{K}$.

### `calc_basis_vectors(wavelength_arr, weighting_arr, temp_arr, c_abs_arr)`
 
Applies the output $T(t)$ from `calc_pah_cooling()` for a given $\lambda_{\rm abs}$ to Equation 3 of Richie & Hensley (2025) to generate the basis spectra.
We apply a weighting array to Equation 3 to account for our use of an adaptive timestep, defined as `weighting = dt_arr / sum(dt_arr)`.

### `calc_normalization(lambda_abs, dlambda, grain_radius, wavelength_arr, p_lambda_arr, wavelength_arr_u, u_lambda_arr)`

Applies the normalization in Equation 5 to the basis spectra according to the power of the input radiation field, `u_lambda`, over the range `[lambda_abs, lambda_abs * (1 + dlambda)]`.

### Defining the basis spectra

Each $\lambda_{\rm abs}$ used in this routine corresponds to its own basis spectrum. 
It would be computationally infeasible to generate basis spectra for every relevant photon wavelength. 
Therefore, we choose to define our $\lambda_{\rm abs}$ wavelengths such that they span the relevant range (i.e., where the energy emitted by $u_\lambda$ contributes significantly to PAH heating) with a fixed spacing of $\lambda_{\rm abs}\Delta \lambda$ between each $\lambda_{\rm abs}$. 
Then, the power absorbed from the radiation field ``photon" is integrated over a small range, $[\lambda_\mathrm{abs}, \lambda_\mathrm{abs}(1+\Delta\lambda)]$. 
We find that the resulting spectrum is converged for $\Delta\lambda\leq0.01$. 
The input basis spectra must be defined for all $\lambda_{\rm abs}$ where the $u_\lambda$ array is defined, or at least where the power from $u_\lambda$ contributes significantly to PAH heating. 
If they are not, then the integrated spectrum given by Equation 6 will not include energy absorbed from photons at those wavelengths. 
It is also imperative that the input `dlambda` and `lambda_abs` values are equivalent to those used to generate the input basis spectra. 
If `dlambda` is larger or smaller than the value used to define the basis spectra, this normalization will either over- or under-count the energy absorbed from the radiation field. 
We provide a utility function for defining the $\lambda_{\rm abs}$ array, `generate_photon_wavelengths()`.
