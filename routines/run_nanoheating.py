from astropy.constants import h, c
from astropy.io import fits
import astropy.units as u
import numpy as np
import os
import subprocess
from scipy.integrate import trapezoid


def main():
    """A wrapper for running the DL01 code."""
    grain_sizes = np.linspace(4, 30, 167)  # Angstroms
    wavelengths = np.logspace(np.log10(0.1), np.log10(2), 60) * u.um  # micron
    input_dir = "/Users/helenarichie/Documents/Grad School/research/projects/kspa/kspa_nanoheating/test_case/"

    for wav_i, wav in enumerate(wavelengths):
        bin_crit = 0.05
        wav_range = wav * bin_crit

        for a_i, size in enumerate(grain_sizes):
            c_ulambd = calc_energy_density(wav, size * u.AA, bin_crit)

            with open(os.path.join(input_dir, "vsg_mrf.par"), "w") as f:
                f.write(f"{wav.value:.5f}d-4 = DELTAMIN (cm)\n")
                f.write(f"{wav.value + wav_range.value:.5f}d-4 = DELTAMAX (cm)\n")
                str_python = str(f"{c_ulambd.value:.5e} = LAMBD_ULAMBD (erg cm-3)")
                str_fortran = str_python.replace("e+", "d+")
                f.write(f"{str_fortran}")
                f.close()
            with open(os.path.join(input_dir, "vsg_size.par"), "w") as f:
                f.write("'graD16emtPAHib' = COMPOSITION (A14)\n")
                f.write("24         = ICASE\n")
                f.write("0         = ISHAP\n")
                f.write("1        = NSIZE = number of sizes (default = 167)\n")
                f.write(f"{size:.3f}d-8 = AMIN (cm) [NC=460*(a/1.e-7)**3]\n")
                f.write("100.7d-8 = AMAX (cm)\n")
            f.close()

            subprocess.run(["mv", "vsg_stat_therm.dpdtlib.out", f"5A/dpdtlib_{size:.2f}_{wav.value:.3f}.out"])


def calc_energy_density(photon_lambd, grain_size, bin_crit):
    """
    photon_lambd : float (length unit)
        photon wavelength, specified as lower limit of monochromatic radiation field delta function
    grain_size : float (length unit)
        dust grain radius
    bin_crit : float
        criterion for bin width relative to photon wavelength, d(lambda) / lambda
    """
    # constants
    freq_day = 1 / (1e5 * u.s)

    bin_width = photon_lambd.to(u.cm) * bin_crit
    lambd0, lambd1 = photon_lambd.to(u.cm), photon_lambd.to(u.cm) + bin_width

    # TODO: implement calc_cabs() function to get grain cross-sections instead of reading in from data file
    hdul_cion = fits.open("../../data/dataverse_files/cion_csecs.fits")

    # read in list of grain sizes and wavelengths that have corresponding C_abs values
    a_cion = hdul_cion[1].data * u.um
    lambd_cion = hdul_cion[2].data * u.um

    # locate array indices corresponding to photon wavelength and dust grain size
    wh_a_cion = np.argmin(np.abs(a_cion - grain_size.to(u.um)))
    wh_lambd0, wh_lambd1 = (
        np.argmin(np.abs(lambd_cion - lambd0.to(u.um))),
        np.argmin(np.abs(lambd_cion - lambd1.to(u.um))),
    )

    # get values of C_abs corresponding to specified wavelength range
    cabs_range = hdul_cion[4].data[:, wh_a_cion][wh_lambd0 : wh_lambd1 + 1] * (u.cm**2)

    lambd_range = lambd_cion[wh_lambd0 : wh_lambd1 + 1]

    # calculate energy density
    ulambda = freq_day / trapezoid(cabs_range.to(u.cm**2) * lambd_range.to(u.cm) / h.cgs, x=lambd_range.to(u.cm))

    print(f"Energy density: {photon_lambd.to(u.cm) * ulambda:.3e}")

    return c.cgs * ulambda


if __name__ == "__main__":
    main()
