import numpy as np
import pandas as pd
from biogas.thermodynamics import (calc_hf_temp, calc_sf_temp, calc_gibbs_temp)
from biogas.data import (CH4, CO2, H2O, H2, CO, C_,
                              a_jk, a_jk_inlet, 
                              atoms_coeffs, atoms_number,
                              HF298, SF298, A, B, C, D)
from scipy.optimize import fmin_slsqp
import matplotlib.pyplot as plt


class GasEquilibrium:
    
    def __init__(self, T = 500, P = 1, P0 = 1, R = 8.314):
        """Class for equilibrium compositions for gas mixture.
        Parameters
        ----------
        T  : float or integer, optional
             Temperature of the system [K], by default 500
        P  : float or integer, optional
             Pressure of the system [atm], by default 1
        P0 : float or integer, optional
             Pressure of the standard state [atm], by default 1
        R  : float, optional
             Gas constant [kJ/kmol/K], by default 8.314
        """
        self.T = T
        self.P = P
        self.P0 = P0
        self.R = R
        self.T_range = np.linspace(573.15, 1213.15)


    def set_inlet(self, ych4=0.3, yco2=0.2, yh2o=0.0, yh2=0.0, yco=0.0, yar=0.5):
        """Set inlet conditions.
        ych4 : float, optional
               Methane inlet molar ratio, by default 0.3
        yco2 : float, optional
               Carbon dioxide inlet molar ratio, by default 0.2
        yh2o : float, optional
               Water inlet molar ratio, by default 0.0
        yh2  : float, optional
               Hydrogen inlet molar ratio, by default 0.0
        yco  : float, optional
               Carbon monoxide inlet molar ratio, by default 0.0
        yar  : float, optional
               Argon inlet molar ratio, by default 0.5
        """
        params = np.array([ych4, yco2, yh2o, yh2, yco, yar])
        keys = ['ych4', 'yco2', 'yh2o', 'yh2', 'yco', 'yar']
        self.inlet = dict(zip(keys, params))
        self._inlet_values = params
        self._keys = keys
        
        params_biogas = np.array([ych4, yco2, yar])
        keys_biogas = ['ych4', 'yco2', 'yar']
        self.inlet_biogas = dict(zip(keys_biogas, params_biogas))
        self._inlet_values_biogas = params_biogas
        self._keys_biogas = keys_biogas

           
    def gibbs_min(self, nj):
        """Minimized Gibbs energy of a mixture"""
        hr = calc_hf_temp(self.T, HF298, A, B, C, D)
        sr = calc_sf_temp(self.T, SF298, A, B, C, D)
        gr = calc_gibbs_temp(self.T, hr, sr)
        nj = np.array(nj)
        gj = gr + self.R*self.T*np.log(nj / np.sum(nj) * self.P / self.P0)
        return nj.dot(gj)


    def ec1(self, nj):
        """conservation of atoms constraint - equality constraint"""
        nj = np.array(nj)
        return a_jk.T.dot(nj) - a_jk_inlet.T.dot(self._inlet_values_biogas)

    
    def ic1(self, nj):
        """inequality constraint all n>=0"""
        nj = np.array(nj)
        return nj
    
    @property
    def min_fun(self):
        n0 = np.random.rand(len(self._inlet_values))
        sol = fmin_slsqp(self.gibbs_min, n0, f_eqcons=self.ec1, f_ieqcons=self.ic1)
        yj = sol/np.sum(sol)
        return yj


    def multiple_temperatures(self, temps):
        return np.array([self.min_fun for self.T in temps]).T

    @property
    def conversions(self):
        yj = self.multiple_temperatures(self.T_range)
        Xch4 = ((self.inlet['ych4'] - yj[CH4])/self.inlet['ych4'])*100
        Xco2 = ((self.inlet['yco2'] - yj[CO2])/self.inlet['yco2'])*100    
        return {'Xch4': Xch4, 'Xco2': Xco2}
  
    @property
    def yields(self):
        yj = self.multiple_temperatures(self.T_range)
        Yh2 = ((yj[H2])/(2*self.inlet['ych4']))*100
        Yco = ((yj[CO])/(self.inlet['ych4'] + self.inlet['yco2']))*100 
        return {'Yh2': Yh2, 'Yco': Yco}   
   
    @property
    def ratio_h2_co(self):
        yj = self.multiple_temperatures(self.T_range)
        H2_CO_ratio = yj[H2]/yj[CO]     
        return {'H2/CO': H2_CO_ratio}

    
    def get_dataframe(self):
        yj = self.multiple_temperatures(self.T_range)
        df = pd.DataFrame()
        df[['Xch4', 'Xco2']] = np.array(list(self.conversions.values())).T
        df[['Yh2', 'Yco']] = np.array(list(self.yields.values())).T
        df['H2/CO'] = np.array(list(self.ratio_h2_co.values())).T
        df['ych4'] = yj[CH4]
        df['yco2'] = yj[CO2]
        df['yh2'] = yj[H2]
        df['yco'] = yj[CO]
        return df


    def to_excel(self, filename, **options):
        """Saves the pandas.DataFrame of profiles in an Excel file.
        Parameters
        ----------
        filename : str
            Name of destination file without suffix .xlsx.
        """
        path = filename + ".xlsx"
        with pd.ExcelWriter(path) as writer:
            self.get_dataframe(**options).to_excel(writer,
                                                   sheet_name='GasEquilibrium')

    def plot_conversions_yields(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.T_range - 273.15, self.conversions['Xch4'], label = '$CH_4$')
        ax.plot(self.T_range - 273.15, self.conversions['Xco2'], label = '$CO_2$')
        ax.plot(self.T_range - 273.15, self.yields['Yh2'], label = '$H_2$')
        ax.plot(self.T_range - 273.15, self.yields['Yco'], label = '$CO$')
        ax.set_title('Equilibrium Conversions/Yields', fontsize=10, fontweight='bold')
        ax.set_ylabel('[%]', fontsize=12, fontweight='bold')
        ax.set_xlabel('Temperature $^o$C', fontsize=12, fontweight='bold')
        ax.grid(ls="-", alpha = 0.2)
        ax.legend()
        ax_b = ax.twinx()
        ax_b.plot(self.T_range-273.15, self.ratio_h2_co['H2/CO'], "k:", label = '$H_2$/$CO$', lw = 1)
        ax_b.legend(loc = "lower right")
        fig.show()


    def plot_molar_ratio(self):
        yj = self.multiple_temperatures(self.T_range)
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.T_range - 273.15, yj[CH4], label = '$CH_4$')
        ax.plot(self.T_range - 273.15, yj[CO2], label = '$CO_2$')
        ax.plot(self.T_range - 273.15, yj[H2O], label = '$H_2O$')
        ax.plot(self.T_range - 273.15, yj[H2], label = '$H_2$')
        ax.plot(self.T_range - 273.15, yj[CO], label = '$CO$')
        ax.set_title('Equilibrium Compositions', fontsize=10, fontweight='bold')
        ax.set_ylabel('molar ratio', fontsize=12, fontweight='bold')
        ax.set_xlabel('Temperature $^o$C', fontsize=12, fontweight='bold')
        ax.grid(ls="-", alpha = 0.2)
        ax.legend()
        fig.show()

            
class Equilibrium(GasEquilibrium):
    
    def __init__(self, T = 500, P = 1, P0 = 1, R = 8.314):
        """Class for equilibrium compositions for gas mixture considering carbon.
        Parameters
        ----------
        T  : float or integer, optional
             Temperature of the system [K], by default 500
        P  : float or integer, optional
             Pressure of the system [atm], by default 1
        P0 : float or integer, optional
             Pressure of the standard state [atm], by default 1
        R  : float, optional
             Gas constant [kJ/kmol/K], by default 8.314
        """
        super(Equilibrium, self).__init__(T=T, P=P, P0=P0, R=R)
        
        self.dGf_carbon = 0
        self.a_jk_carbon = np.vstack([a_jk, atoms_coeffs(atoms_number, [C_], [1])])
        
        
    def add_carbon(self):
        self.inlet_species_carbon = np.append([self._inlet_values], 0.0)
        

    def gibbs_min_carbon(self, n):
        """Minimized Gibbs energy of a mixture considering carbon"""
        hr = calc_hf_temp(self.T, HF298, A, B, C, D)
        sr = calc_sf_temp(self.T, SF298, A, B, C, D)
        gr = calc_gibbs_temp(self.T, hr, sr)
        n = np.array(n)
        nj = n[:-1] # moles of gas
        nc = n[-1]  # moles of carbon
        gj = gr + self.R*self.T*np.log(nj / np.sum(nj) * self.P / self.P0)
        return nj.dot(gj) + nc*self.dGf_carbon

    
    def ec2(self, nj):
        """conservation of atoms constraint - equality constraint"""
        nj = np.array(nj)
        return self.a_jk_carbon.T.dot(nj) - a_jk_inlet.T.dot(self._inlet_values_biogas)


    @property
    def min_fun_carbon(self):
        n0 = np.random.rand(len(self.inlet_species_carbon))
        sol = fmin_slsqp(self.gibbs_min_carbon, n0, f_eqcons=self.ec2, f_ieqcons=self.ic1)
        yj = sol/np.sum(sol)
        return yj


    def multiple_temperatures_carbon(self, temps):
        return np.array([self.min_fun_carbon for self.T in temps]).T

    @property
    def conversions_carbon(self):
        yj = self.multiple_temperatures_carbon(self.T_range)
        Xch4 = ((self.inlet['ych4'] - yj[CH4])/self.inlet['ych4'])*100
        Xco2 = ((self.inlet['yco2'] - yj[CO2])/self.inlet['yco2'])*100    
        return {'Xch4': Xch4, 'Xco2': Xco2}
   
    @property
    def yields_carbon(self):
        yj = self.multiple_temperatures_carbon(self.T_range)
        Yh2 = ((yj[H2])/(2*self.inlet['ych4']))*100
        Yco = ((yj[CO])/(self.inlet['ych4'] + self.inlet['yco2']))*100 
        return {'Yh2': Yh2, 'Yco': Yco}   
    
    @property
    def ratio_h2_co_carbon(self):
        yj = self.multiple_temperatures_carbon(self.T_range)
        H2_CO_ratio = yj[H2]/yj[CO]     
        return {'H2/CO': H2_CO_ratio}

    
    def get_dataframe_carbon(self):
        yj = self.multiple_temperatures_carbon(self.T_range)
        df = pd.DataFrame()
        df[['Xch4', 'Xco2']] = np.array(list(self.conversions_carbon.values())).T
        df[['Yh2', 'Yco']] = np.array(list(self.yields_carbon.values())).T
        df['H2/CO'] = np.array(list(self.ratio_h2_co_carbon.values())).T
        df['ych4'] = yj[CH4]
        df['yco2'] = yj[CO2]
        df['yh2'] = yj[H2]
        df['yco'] = yj[CO]
        df['yc'] =  yj[-1]       
        return df
    

    def to_excel_carbon(self, filename, **options):
        """Saves the pandas.DataFrame of profiles in an Excel file.
        Parameters
        ----------
        filename : str
            Name of destination file without suffix .xlsx.
        """
        path = filename + ".xlsx"
        with pd.ExcelWriter(path) as writer:
            self.get_dataframe_carbon(**options).to_excel(writer,
                                                   sheet_name='GasEquilibriumCarbon')    

    def plot_conversions_yields_carbon(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.T_range - 273.15, self.conversions_carbon['Xch4'], label = '$CH_4$')
        ax.plot(self.T_range - 273.15, self.conversions_carbon['Xco2'], label = '$CO_2$')
        ax.plot(self.T_range - 273.15, self.yields_carbon['Yh2'], label = '$H_2$')
        ax.plot(self.T_range - 273.15, self.yields_carbon['Yco'], label = '$CO$')
        ax.set_title('Equilibrium Conversions/Yields', fontsize=10, fontweight='bold')
        ax.set_ylabel('[%]', fontsize=12, fontweight='bold')
        ax.set_xlabel('Temperature $^o$C', fontsize=12, fontweight='bold')
        ax.grid(ls="-", alpha = 0.2)
        ax.legend()
        ax.legend(loc = 'upper left')
        ax_b = ax.twinx()
        ax_b.plot(self.T_range-273.15, self.ratio_h2_co_carbon['H2/CO'], "k:", label = '$H_2$/$CO$', lw = 1)
        ax_b.legend(loc = 'lower right')
        fig.show()
     

    def plot_molar_ratio_carbon(self):
        yj = self.multiple_temperatures_carbon(self.T_range)
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.T_range - 273.15, yj[CH4], label = '$CH_4$')
        ax.plot(self.T_range - 273.15, yj[CO2], label = '$CO_2$')
        ax.plot(self.T_range - 273.15, yj[H2O], label = '$H_2O$')
        ax.plot(self.T_range - 273.15, yj[H2], label = '$H_2$')
        ax.plot(self.T_range - 273.15, yj[CO], label = '$CO$')
        ax.plot(self.T_range - 273.15, yj[-1], label = '$C$')
        ax.set_title('Equilibrium Compositions', fontsize=10, fontweight='bold')
        ax.set_ylabel('molar ratio', fontsize=12, fontweight='bold')
        ax.set_xlabel('Temperature $^o$C', fontsize=12, fontweight='bold')
        ax.grid(ls="-", alpha = 0.2)
        ax.legend()
        fig.show()
