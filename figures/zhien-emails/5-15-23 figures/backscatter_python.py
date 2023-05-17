# find the molecular backscatter coefficient using pressure (hPa), temperature (C), and crl wavelength (nm)
def find_molec_backscatter(p, T, crl_wavelength):
    # find the number density first
    k = 1.381e-23 # boltzmann constant, m2 kg / s2 K
    Molecular_number_density = (p*100.) / (k * (T+273.15))
    MND = Molecular_number_density* 10**-32 # convert number density to m^(-3)
    Molecular_backscattering_Coefficient = MND *5.45 *(550./ crl_wavelength)**4 # units: m^(-1)sr^(-1)

    return Molecular_backscattering_Coefficient * 10**3 # units: km^(-1)sr^(-1)
