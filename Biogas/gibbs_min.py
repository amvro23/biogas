import numpy as np
from data_nist import (CH4, CO2, H2O, H2, CO, temp_range, C_,
                              a_jk, a_jk_inlet, inlet_species, atoms_coeffs, atoms_number,
                              HF298, A, B, C, D, E, F, G, H)
from thermodynamics_nist import std_hr, std_sr, std_gibbs
from scipy.optimize import fmin_slsqp
import matplotlib.pyplot as plt

class GasEquilibrium:
    
    # class variables
    inlet = inlet_species[inlet_species!=0]
    T_range = temp_range
    
    def __init__(self, T = 500, P = 1, P0 = 1, R = 8.314e-3):
        self.T = T
        self.P = P
        self.P0 = P0
        self.R = R
           
    def gibbs_min(self, nj): 
        hr = std_hr(self.T, HF298, A, B, C, D, E, F, G, H)
        sr = std_sr(self.T, A, B, C, D, E, F, G)
        gr = std_gibbs(self.T, hr, sr)
        nj = np.array(nj)
        gj = gr + self.R*self.T*np.log(nj / np.sum(nj) * self.P / self.P0)
        return nj.dot(gj)   

    def ec1(self, nj):
        nj = np.array(nj)
        '''conservation of atoms constraint - equality constraint'''
        return a_jk.T.dot(nj) - a_jk_inlet.T.dot(self.inlet)
    
    def ic1(self, nj):
        nj = np.array(nj)
        '''inequality constraint all n>=0'''
        return nj
    
    @property
    def min_fun(self):
        n0 = np.random.rand(len(inlet_species))  # initial guesses
        sol = fmin_slsqp(self.gibbs_min, n0, f_eqcons=self.ec1, f_ieqcons=self.ic1)
        yj = sol/np.sum(sol)
        return yj

    def multiple_temperatures(self, temps):
        return np.array([self.min_fun for self.T in temps]).T

    @property
    def conversions(self):
        yj = self.multiple_temperatures(self.T_range)
        Xch4 = ((inlet_species[CH4] - yj[CH4])/inlet_species[CH4])*100
        Xco2 = ((inlet_species[CO2] - yj[CO2])/inlet_species[CO2])*100    
        return {'Xch4': Xch4, 'Xco2': Xco2}
    
    @property
    def yields(self):
        yj = self.multiple_temperatures(self.T_range)
        Yh2 = ((yj[H2])/(2*inlet_species[CH4]))*100
        Yco = ((yj[CO])/(inlet_species[CH4] + inlet_species[CO2]))*100 
        return {'Yh2': Yh2, 'Yco': Yco}   
    
    @property
    def ratio_h2_co(self):
        yj = self.multiple_temperatures(self.T_range)
        H2_CO_ratio = yj[H2]/yj[CO]     
        return {'H2/CO': H2_CO_ratio}
    
    @property
    def plot_conversions_yields(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 150)
        ax.plot(self.T_range - 273.15, self.conversions['Xch4'], label = 'CH$_4$')
        ax.plot(self.T_range - 273.15, self.conversions['Xco2'], label = 'CO$_2$')
        ax.plot(self.T_range - 273.15, self.yields['Yh2'], label = 'H$_2$')
        ax.plot(self.T_range - 273.15, self.yields['Yco'], label = 'CO')
        ax.set_title('Equilibrium Conversions/Yields', fontsize=10, fontweight='bold')
        ax.set_ylabel('[%]', fontsize=12, fontweight='bold')
        ax.set_xlabel('Temperature $^o$C', fontsize=12, fontweight='bold')
        ax.grid(ls="-", alpha = 0.2)
        ax.legend()
        ax_b = ax.twinx()
        ax_b.plot(self.T_range-273.15, self.ratio_h2_co['H2/CO'], "k:", label = 'H$_2$/CO', lw = 1)
        ax_b.legend(loc = "lower right")
        fig.show()
     
    @property
    def plot_molar_ratio(self):
        yj = self.multiple_temperatures(self.T_range)
        fig, ax = plt.subplots(figsize = (6,4), dpi = 150)
        ax.plot(self.T_range - 273.15, yj[CH4], label = 'CH$_4$')
        ax.plot(self.T_range - 273.15, yj[CO2], label = 'CO$_2$')
        ax.plot(self.T_range - 273.15, yj[H2O], label = 'H$_2$O')
        ax.plot(self.T_range - 273.15, yj[H2], label = 'H$_2$')
        ax.plot(self.T_range - 273.15, yj[CO], label = 'CO')
        ax.set_title('Equilibrium Compositions', fontsize=10, fontweight='bold')
        ax.set_ylabel('molar ratio', fontsize=12, fontweight='bold')
        ax.set_xlabel('Temperature $^o$C', fontsize=12, fontweight='bold')
        ax.grid(ls="-", alpha = 0.2)
        ax.legend()
        fig.show()
            
class CarbonEquilibrium(GasEquilibrium):
    
    dGf_carbon = 0
    a_jk_carbon = np.vstack([a_jk, atoms_coeffs(atoms_number, [C_], [1])])
    inlet_species_carbon = np.append([inlet_species], 0.0)
    
    def _init(self, T, P):
        super().__init__(T, P)
    
    def gibbs_min_carbon(self, n):
        hr = std_hr(self.T, HF298, A, B, C, D, E, F, G, H)
        sr = std_sr(self.T, A, B, C, D, E, F, G)
        gr = std_gibbs(self.T, hr, sr)
        n = np.array(n)
        nj = n[:-1] # moles of gas
        nc = n[-1]  # moles of carbon
        gj = gr + self.R*self.T*np.log(nj / np.sum(nj) * self.P / self.P0)
        return nj.dot(gj) + nc*self.dGf_carbon
    
    def ec2(self, nj):
        nj = np.array(nj)
        '''conservation of atoms constraint - equality constraint'''
        return self.a_jk_carbon.T.dot(nj) - a_jk_inlet.T.dot(self.inlet)
    
    def ic2(self, nj):
        nj = np.array(nj)
        '''inequality constraint all n>=0'''
        return nj 

    @property
    def min_fun_carbon(self):
        n0 = np.random.rand(len(self.inlet_species_carbon))  # initial guesses
        sol = fmin_slsqp(self.gibbs_min_carbon, n0, f_eqcons=self.ec2, f_ieqcons=self.ic2)
        yj = sol/np.sum(sol)
        return yj

    def multiple_temperatures_carbon(self, temps):
        return np.array([self.min_fun_carbon for self.T in temps]).T

    @property
    def conversions_carbon(self):
        yj = self.multiple_temperatures_carbon(self.T_range)
        Xch4 = ((inlet_species[CH4] - yj[CH4])/inlet_species[CH4])*100
        Xco2 = ((inlet_species[CO2] - yj[CO2])/inlet_species[CO2])*100    
        return {'Xch4': Xch4, 'Xco2': Xco2}
    
    @property
    def yields_carbon(self):
        yj = self.multiple_temperatures_carbon(self.T_range)
        Yh2 = ((yj[H2])/(2*inlet_species[CH4]))*100
        Yco = ((yj[CO])/(inlet_species[CH4] + inlet_species[CO2]))*100 
        return {'Yh2': Yh2, 'Yco': Yco}   
    
    @property
    def ratio_h2_co_carbon(self):
        yj = self.multiple_temperatures_carbon(self.T_range)
        H2_CO_ratio = yj[H2]/yj[CO]     
        return {'H2/CO': H2_CO_ratio}
    
    @property
    def plot_conversions_yields_carbon(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 150)
        ax.plot(self.T_range - 273.15, self.conversions_carbon['Xch4'], label = 'CH$_4$')
        ax.plot(self.T_range - 273.15, self.conversions_carbon['Xco2'], label = 'CO$_2$')
        ax.plot(self.T_range - 273.15, self.yields_carbon['Yh2'], label = 'H$_2$')
        ax.plot(self.T_range - 273.15, self.yields_carbon['Yco'], label = 'CO')
        ax.set_title('Equilibrium Conversions/Yields', fontsize=10, fontweight='bold')
        ax.set_ylabel('[%]', fontsize=12, fontweight='bold')
        ax.set_xlabel('Temperature $^o$C', fontsize=12, fontweight='bold')
        ax.grid(ls="-", alpha = 0.2)
        ax.legend()
        ax_b = ax.twinx()
        ax_b.plot(self.T_range-273.15, self.ratio_h2_co_carbon['H2/CO'], "k:", label = 'H$_2$/CO', lw = 1)
        ax_b.legend(loc = "lower right")
        fig.show()
     
    @property
    def plot_molar_ratio_carbon(self):
        yj = self.multiple_temperatures_carbon(self.T_range)
        fig, ax = plt.subplots(figsize = (6,4), dpi = 150)
        ax.plot(self.T_range - 273.15, yj[CH4], label = 'CH$_4$')
        ax.plot(self.T_range - 273.15, yj[CO2], label = 'CO$_2$')
        ax.plot(self.T_range - 273.15, yj[H2O], label = 'H$_2$O')
        ax.plot(self.T_range - 273.15, yj[H2], label = 'H$_2$')
        ax.plot(self.T_range - 273.15, yj[CO], label = 'CO')
        ax.plot(self.T_range - 273.15, yj[-1], label = 'C')
        ax.set_title('Equilibrium Compositions', fontsize=10, fontweight='bold')
        ax.set_ylabel('molar ratio', fontsize=12, fontweight='bold')
        ax.set_xlabel('Temperature $^o$C', fontsize=12, fontweight='bold')
        ax.grid(ls="-", alpha = 0.2)
        ax.legend()
        fig.show()