import numpy as np

# ------------------------------------------------------------------------------------------------
# INDIVIDUAL COMPONENTS
# ------------------------------------------------------------------------------------------------

def fuller_ab_mat(M, nu, T, P):
    """
    Function to get the diffusion coefficients of A in B pure using Fuller.
    Parameters
    ----------
    M : 1d array like
        Molar masses of components in kg/kmol.
    nu : 1d array like
        Coefficients of Fuller of components.
    T : float or int
        Temperature in K.
    P : float or int
        Pressure in bar.
    Returns
    -------
    mat: 2d array
        Diffusion coefficient of each pair of A in B pure in m^2/s.
    """
    
    Ma = np.atleast_2d(M).reshape([-1, 1])
    Mb = np.atleast_2d(M).reshape([1, -1])
    
    nua = np.atleast_2d(nu).reshape([-1, 1])
    nub = np.atleast_2d(nu).reshape([1, -1])
    
    mat = 1e-9 * (T**1.75) * ((1/Ma + 1/Mb) ** 0.5) / \
        (P * ((nua**(1/3) + nub**(1/3))**2))
    
    return mat


# ------------------------------------------------------------------------------------------------
# MIXTURE
# ------------------------------------------------------------------------------------------------


def wilke_mist(y, mat):
    """
    Returns the diffusion coefficients for each component in the mixture using the Wilke equation.
    Parameters
    ----------
    y : list od 1d array
        Contains fractions of each component in the phase of mixture.
    mat : 2d array
        Matrix containing individual diffusion coefficients for each pair AB.
    Returns
    -------
    vect: 1d array
        Contains each component diffusion coefficient in the mixture.
    """
    y = np.array(y)
    vect = (1 - y) / ((y / mat).sum(axis=1) - y / np.diag(mat))
    
    return vect


def sm_mist(y, a_est, mat):
    """
    Returns the diffusion coefficients for each component in the reacting mixture using the Stefan-Maxwell equation.
    Parameters
    ----------
    y : list od 1d array
        Contains fractions of each component in the phase of mixture.
    a_est : 1d array
        Stoichiometric coefficients.
    mat : 2d array
        Matrix containing individual diffusion coefficients for each pair AB.
    Returns
    -------
    vect: 1d array
        Contains each component diffusion coefficient in the mixture.
    """
    return (
        (1 - y * np.sum(a_est.reshape((-1, 1)) / a_est.reshape((-1, 1)).T, axis=0))
        / (np.sum((a_est.reshape((-1, 1)) / a_est.reshape((-1, 1)).T)
                  * (y.reshape((-1, 1)) - y.reshape((-1, 1)).T)
                  * (1 / mat.T), axis=0))
    )


def knudseen(r_pore, T, Mm):
    """
    Returns the Knudsen diffusion coefficient in a porous system.
    Parameters
    ----------
    r_pore : float
        Pore radius in m.
    T : float or int
        Temperature in K.
    Mm : float, int, or 1d array
        Molar mass of each component in kg/kmol.
    Returns
    -------
    Dk: float or 1d array
        Knudsen diffusion coefficient of each component in a porous system.
    """
    return 97 * r_pore * (T / Mm) ** 0.5


def effective_diff(Da, tao, es):
    """
    Returns the effective diffusion coefficients of components through the solid.
    Parameters
    ----------
    Da : float or 1d array
        Contains the diffusion coefficients of each component in the mixture.
    tao : float or int
        Tortuosity factor.
    es : float or int
        Porosity.
    Returns
    -------
    Dae: float or 1d array
        Contains the diffusion coefficients of each component through the solid.
    """
    return Da*es/tao


def effective_diff_pore(Dk, Dm, tao, es):
    """
    Returns the Knudsen diffusion coefficient in a porous system using diffusion coefficients in the mixture and knudseen.
    Parameters
    ----------
    Dk : 1d array
        Knudssen diffusion coefficients of each component.
    Dm : 1d array
        Diffusion coefficients of each component in the mixture.
    tao : float or int
        Tortuosity factor.
    es : float or int
        Porosity.
    Returns
    -------
    Daek : 1d array
        Contains the effective diffusion coefficients of each component through the solid.
    """
    Dk = np.array(Dk)
    Dm = np.array(Dm)
    
    Daek = es/tao * Dk * Dm / (Dk + Dm)
    
    return Daek

# atomic diffusion volume
_nuC = 16.5
_nuH = 1.98
_nuO = 5.48
_nuArom = -20.2

def fnu(C, H, O, Arom):
    """
    Returns the nu coefficients of the Fuller difussion equation.
    Parameters
    ----------
    C : int
        Number of Carbons in substance.
    H : int
        Number of Hidrogens in substance.
    O : int
        Number of Oxigens in substance.
    Arom : int
        Number of Aromatic rings in substance.
    Returns
    -------
    nu: int or float
        nu coefficient for Fuller equation.
    """
    return C*_nuC + H*_nuH + O*_nuO + Arom*_nuArom