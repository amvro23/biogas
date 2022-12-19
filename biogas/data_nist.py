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
A = np.array([-7.03029e-1,  2.499735e1,  3.0092e1, 3.3066178e1,  2.556759e1, 20.78600])
B = np.array([ 1.084773e2,  5.518696e1,  6.832514,-1.1363417e+01,  6.09613, 2.825911e-07])
C = np.array([-4.252157e1, -3.36913700e1,  6.793435, 1.1432816e+01,  4.054656, -1.464191e-07])
D = np.array([ 5.862788,  7.948387, -2.53448, -2.772874, -2.671301, 1.092131e-08])
E = np.array([ 6.78565e-1, -1.36638e-1,  8.2139e-2, -1.58558e-1,  1.31021e-1, -3.661371e-08])
F = np.array([-7.684376e1, -4.036075e2, -2.50881e2, -9.980797, -1.180089e2, -6.197350])
G = np.array([ 1.587163e2,  2.282431e2,  2.233967e2, 1.72707974e2,  2.273665e2, 179.9990])
H = np.array([-7.48731e1, -3.935224e2, -2.418264e2, 0.0, -1.105271e2, 0.0])

# Heat, Gibbs free energy, and Entropy of formation in kJ/kmol
HF298 = np.array([-74.85, -393.51, -241.826, 0.0, -110.53, 0.0]) # kJ/mol

# ------------------------------------------------------------------------------------------------
# INLET RATIO OF BIOGAS & TEMPERATURE RANGE
# ------------------------------------------------------------------------------------------------

# 50% biogas & 50% argon
inlet_species = np.array([0.3, 0.2, 0.0, 0.0, 0.0, 0.5])
temp_range = np.linspace(573.15, 1213.15)

# ------------------------------------------------------------------------------------------------
# ATOMS
# ------------------------------------------------------------------------------------------------               

C_ = 0
O_ = 1
H_ = 2
AR_ = 3
atoms_number = np.array([C_, O_, H_, AR_])

def atoms_coeffs(atoms, positions, y):
    atoms = np.array(atoms)
    positions = np.array(positions)
    components = np.zeros(atoms.shape, dtype=int)
    components[positions] = y
    return components

a_jk = np.zeros([len(inlet_species), len(atoms_number)])
a_jk[CH4] = atoms_coeffs(atoms_number, [C_, H_], [1, 4])
a_jk[CO2] = atoms_coeffs(atoms_number, [C_, O_], [1, 2])
a_jk[H2O] = atoms_coeffs(atoms_number, [H_, O_], [2, 1])
a_jk[H2] = atoms_coeffs(atoms_number, [H_], [2])
a_jk[CO] = atoms_coeffs(atoms_number, [C_, O_], [1, 1])
a_jk[AR] = atoms_coeffs(atoms_number, [AR_], [1])

# select the inlet molecules which is biogas (i.e., CH4 & CO2) and inert gas argon
a_jk_inlet = np.array([a_jk[CH4], a_jk[CO2], a_jk[AR]])
