import numpy as np

def std_hr(T, hf298, a, b, c, d, e, f, g, h):
    """
    Returns the standard enthalpy of species in kJ/mol.
    Parameters
    ----------
    T : float or int
        Temperature in K.
    a : 1d array, float or int
        Contains Nist coefficient of the species.
    b : 1d array, float or int
        Contains Nist coefficient of the species.
    c : 1d array, float or int
        Contains Nist coefficient of the species.
    d : 1d array, float or int
        Contains Nist coefficient of the species.
    e : 1d array, float or int
        Contains Nist coefficient of the species.
    f : 1d array, float or int
        Contains Nist coefficient of the species.
    g : 1d array, float or int
        Contains Nist coefficient of the species.
    h : 1d array, float or int
        Contains Nist coefficient of the species.
    Returns
    -------
    1d array, float or int
        Contains standard enthalpy of species at T in kJ/mol.
    """
    T = T/1000
    
    return np.array(hf298) + np.array(a)*(T) + np.array(b)/2*(T**2)\
        + np.array(c)/3*(T**3) + np.array(d)/4*(T**4) - np.array(e)/(T)\
        + 1*np.array(f) + 0*np.array(g) - 1*np.array(h)


def std_sr(T, a, b, c, d, e, f, g):
    """
    Returns the standard entropy of species in J/mol/K.
    Parameters
    ----------
    T : float or int
        Temperature in K.
    a : 1d array, float or int
        Contains Nist coefficient of the species.
    b : 1d array, float or int
        Contains Nist coefficient of the species.
    c : 1d array, float or int
        Contains Nist coefficient of the species.
    d : 1d array, float or int
        Contains Nist coefficient of the species.
    e : 1d array, float or int
        Contains Nist coefficient of the species.
    f : 1d array, float or int
        Contains Nist coefficient of the species.
    g : 1d array, float or int
        Contains Nist coefficient of the species.
    Returns
    -------
    1d array, float or int
        Contains standard entropy of species at T in J/mol/K.
    """
    T = T/1000
    
    return np.array(a)*np.log(T) + np.array(b)*(T)\
        + np.array(c)/2*(T**2) + np.array(d)/3*(T**3) - np.array(e)/(2*T**2)\
        + 0*f + 1*g

def std_gibbs(T, standard_hr, standard_sr):
    """
    Returns the gibbs energy of reaction in kJ/mol.
    Parameters
    ----------
    T : float or int
        Temperature in K.
    standard_hr : float, int, or 1d array
        Contains heat of reaction at T in kJ/kmol.
    standard_hr : float, int, or 1d array
        Contains entropy of reaction at T in kJ/mol/K.
    Returns
    -------
    float, int, or 1d array
        Contains gibbs energy of reaction at T in kJ/mol.
    """
    return (np.array(standard_hr) - T*np.array(standard_sr)/1000)
