import numpy as np

# -----------------------------------------------------------------------------
# ARRHENIUS & VANT HOFF PARAMETERS
# -----------------------------------------------------------------------------
def arrhenius(A, B, T):
    return A*np.exp(B*(1-1123.15/T))

def vant_hoff(A, B, T, R=8.314):
    return A * np.exp(-B*1e3 / (R * T))

# -----------------------------------------------------------------------------
# EQUILIBRIUM
# -----------------------------------------------------------------------------
def exp_term2(A, B, T):
    return np.exp(A/T+B)

# -----------------------------------------------------------------------------
# COMPONENTS
# -----------------------------------------------------------------------------
CH4, CO2, H2O, H2, CO, AR = np.arange(6)

# -----------------------------------------------------------------------------
# ADSORPTION
# -----------------------------------------------------------------------------
def fnum_ad(p, T): 
    K_ch4 = vant_hoff(6.65e-04, -38.280, T)
    K_h2o = vant_hoff(1.77e05, 88.680, T)
    K_co = vant_hoff(8.23e-05, -70.650, T) 
    K_h2 = vant_hoff(6.12e-09, -82.900, T)
    return (1 + K_co*p[CO] + K_h2*p[H2] + K_ch4*p[CH4] + K_h2o*(p[H2O]/p[H2]))

# -----------------------------------------------------------------------------
# CATALYST REACTIONS
# -----------------------------------------------------------------------------
def rxn_srm1(p, T):
    k_SRM1 = arrhenius(416.547, 25.433, T)
    Kpsrm1 = exp_term2(-26200, 29.71, T)
    num = (k_SRM1/p[H2]**2.5)*(p[CH4]*p[H2O] - p[H2]**3*p[CO]/Kpsrm1)
    den = fnum_ad(p, T)
    return num / den

def rxn_srm2(p, T):
    k_SRM2 = arrhenius(57.512, 28.492, T)
    Kpsrm2 = exp_term2(-23400, 27.40, T)
    num = (k_SRM2/p[H2]**3.5)*(p[CH4]*p[H2O]**2 - p[H2]**4*p[CO2]/Kpsrm2)
    den = fnum_ad(p, T)
    return num/den

def rxn_wgs(p, T):
    k_WGS = arrhenius(14.497, 8.716, T)
    KpWGS = exp_term2(-23400, 27.40, T)/exp_term2(-26200, 29.71, T)
    num = (k_WGS/p[H2])*(p[CO]*p[H2O] - p[H2]*p[CO2]/KpWGS)
    den = fnum_ad(p, T)
    return num/den

def rxn_drm(p, T): 
    k_DRM = arrhenius(775.106, 24.160, T)
    KpWGS = exp_term2(-23400, 27.40, T)/exp_term2(-26200, 29.71, T)
    KpDRM = exp_term2(-26200, 29.71, T)/KpWGS # bar2
    num = k_DRM*(p[CH4]*p[CO2] - p[H2]**2*p[CO]**2/KpDRM)
    den = (1 + 0.5*p[CH4] + 9.71*p[CO])*(1 + 26.1*p[CO2])
    return num/den