import numpy as np

# ------------------------------------------------------------------------------------------------
# COMPONENTS
# ------------------------------------------------------------------------------------------------

CH4 = 0
CO2 = 1
H2O = 2
H2 = 3
CO = 4
AR = 5

# ------------------------------------------------------------------------------------------------
# THERMODYNAMIC PROPERTIES
# ------------------------------------------------------------------------------------------------

# Coefficients in heat capacity equation Reid et al (1987)
A = np.array([19.3, 19.8, 32.2, 27.1, 30.9, 20.8])
B = np.array([52.1e-3, 73.4e-3, 1.9e-3, 9.3e-3, -12.9e-3, 0.0e-3])
C = np.array([12.0e-6, -56.0e-6, 10.6e-6, -13.8e-6, 27.9e-6, 51.7e-6])
D = np.array([-11.3e-9, 17.2e-9, -3.6e-9, 7.65e-9, -12.7e-9, 0.0e-9])

# Heat, Gibbs free energy, and Entropy of formation in kJ/kmol
HF298 = np.array([-74.9e3, -393.8e3, -242.0e3, 0.0e3, -110.6e3, 0.0e3])
GF298 = np.array([-50.9e3, -394.6e3, -228.8e3, 0.0e3, -137.4e3, 0.0e3])
SF298 = (np.array(HF298) - np.array(GF298)) / 298.15

def get_delta_coefficient(coef, positions, y):
    coef = np.array(coef)
    positions = np.array(positions)
    components = np.zeros(coef.shape, dtype=int) 
    components[positions] = y
    return coef.dot(components)

def get_delta_multiple_coefficient(coefs, positions, y):
    return np.array([get_delta_coefficient(coef, positions, y) for coef in coefs])

DELTA_DRM = get_delta_multiple_coefficient([A, B, C, D, HF298, GF298, SF298], [CO2, CH4, CO, H2], [-1, -1, 2, 2])
DELTA_SRM1 = get_delta_multiple_coefficient([A, B, C, D, HF298, GF298, SF298], [CH4, H2O, CO, H2], [-1, -1, 1, 3])
DELTA_SRM2 = get_delta_multiple_coefficient([A, B, C, D, HF298, GF298, SF298], [CH4, H2O, CO2, H2], [-1, -2, 1, 4])
DELTA_WGS = get_delta_multiple_coefficient([A, B, C, D, HF298, GF298, SF298], [CO, H2O, CO2, H2], [-1, -1,  1,  1])

def get_deltas(*args):
    mat = np.vstack(args)
    return mat.T

DELTA_A, DELTA_B, DELTA_C, DELTA_D, DELTA_HR298, DELTA_GF298, DELTA_SR298 = \
    get_deltas(DELTA_DRM, DELTA_SRM1, DELTA_SRM2, DELTA_WGS)

# ------------------------------------------------------------------------------------------------
# FLUIDDYNAMIC PROPERTIES
# ------------------------------------------------------------------------------------------------

MM = np.array([16.043, 44.01, 18.015, 2.016, 28.01, 39.948])  # molar mass in kg/kmol
TC = np.array([190.4, 304.2, 647.3, 33.2, 132.9, 150.8])  # critical temperature in K
PC = np.array([46.0, 73.8, 220.5, 13.0, 35.0, 48.7])  # critical pressure in bar
SIGMA = np.array([3.758, 3.941, 2.641, 2.827, 2.69, 3.542])  # Lennard-Jones parameter in Angstron
EK = np.array([148.6, 195.2, 809.1, 59.7, 91.7, 93.3])  # potential parameter in K
DELTA_POT = np.zeros_like(MM)
DELTA_POT[H2O] = 1.0
NU = np.array([21.42, 26.9, 12.7, 7.07, 18.9, 16.1]) # molecular diffusion volume

# For Methane: Thodos Equation
# For Carbon Dioxide, Water, Hydrogen, Carbon Monoxide, Argon: Chapman-Enskog

HC = np.array([CH4])
LOW = np.array([CO2, H2O, H2, CO, AR])