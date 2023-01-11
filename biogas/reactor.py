import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import math

from biogas.thermodynamics import (get_Cp, calc_delta_hr, calc_hf_temp)
from biogas.fluid_dynamics import (calc_mu_mist, calc_pressure_drop)
from biogas.kinetics import (rxn_srm1, rxn_srm2, rxn_wgs, rxn_drm)
from biogas.kinetics import CH4, CO2, H2O, H2, CO, AR
from biogas.data import (A, B, C, D, HF298,
                         DELTA_A, DELTA_B, DELTA_C, DELTA_D,
                         DELTA_HR298)
from biogas.data import (MM, EK, DELTA_POT, PC, TC, SIGMA)

# -----------------------------------------------------------------------------
# BASE CLASS
# -----------------------------------------------------------------------------
class CatalystBed(object):

    def __init__(
        self, z,
        rhos=950.0, rhob=570.0, es=0.4,
        dp=0.425e-3, tao=3.0, inner_R=9.2e-3, Pmin=0.5, Pterm=None,
        terminal=True, ivp_rtol=1e-6
    ):
        """Class for catalyst bed.
        Parameters
        ----------
        z        : float
                   Reactor length [m].
        rhos     : float, optional
                   Catalyst solid density [kg/m3], by default 950.0
        rhob     : float, optional
                   Catalyst bulk density [kg/m3], by default 570.0
        es       : float, optional
                   Catalyst solid void fraction, by default 0.4
        dp       : float, optional
                   Pellet equivalent diameter [m], by default 0.425e-3
        tao      : float, optional
                   Pellet tortuosity, by default 3.0
        inner_R  : float, optional
                   Catalyst bed inner radius [m], by default 9.2e-3
        Pmin     : float, optional
                   Minimum pressure allowed [bar]. Useful in optimization. By default 0.5
        Pterm    : float or None, optional
                   Terminal pressure to interrupt ODE system [bar], by default None
        ivp_rtol : float, optional
                   Relative tolerance for ODE system, by default 1e-6
        eg       : Void fraction of the catalyst bed"
        """

        self.z = z
        self.rhos = rhos
        self.rhob = rhob
        self.es = es
        self.eg = 1 - rhob/rhos
        self.dp = dp
        self.tao = tao
        self.inner_R = inner_R
        self.Pmin = Pmin
        self.terminal = terminal
        self._ivp_rtol = ivp_rtol
        
        
        if Pterm is None:
            self.Pterm = Pmin / 10
        else:
            self.Pterm = Pterm


    def ode_system(self, z, params):
        pass


    def _f_term(self, z, params):
        return params[-1] - self.Pterm

    _f_term.terminal = True

    
    def set_inlet(self, Fch4=30, Fco2=20, Fh2o=0, Fh2=0, Fco=0,
                  Far=50, T=773, P=1.0):
        """Set inlet conditions of catalyst bed.
        Parameters
        ----------
        Fch4 : float, optional
             Methane feed ratio [ml/min], by default 30
        Fco2 : float, optional
             Carbon dioxide feed ratio [ml/min], by default 
        Fh2o : float, optional
             Water feed ratio [ml/min], by default 0.0
        Fh2  : float, optional
             Hydrogen feed ratio [ml/min], by default 0.0
        Fco  : float, optional
             Carbon monoxide feed ratio [ml/min], by default 0.0
        Far  : float, optional
             Argon [ml/min], by default 50
        T    : int, optional
             Temperature [K], by default 773
        P    : float, optional
             Pressure [bar], by default 1.0
        """
        convert_to_mol_sec = 1/22414/60
        Fch4 = Fch4*convert_to_mol_sec
        Fco2 = Fco2*convert_to_mol_sec
        Fh2o = Fh2o*convert_to_mol_sec
        Fh2 = Fh2*convert_to_mol_sec
        Fco = Fco*convert_to_mol_sec
        Far = Far*convert_to_mol_sec
        
        params = np.array([Fch4, Fco2, Fh2o, Fh2, Fco, Far, T, P])
        keys = ['Fch4', 'Fco2', 'Fh2o', 'Fh2', 'Fco', 'Far', 'T', 'P']
        self.inlet = dict(zip(keys, params))
        self._inlet_values = params
        self._keys = keys


    def solve(self, points_eval=None, **kwargs):
        """Solve ODE system.
        Parameters
        ----------
        points_eval : int or None, optional
            Number of points to eval in ivp solution, by default None
        **kwargs    : any
            Additional keyword arguments passed to scipy.integrate solve_ivp.
        """
        t_eval = None
        if not (points_eval is None):
            t_eval = np.linspace(0, self.z, points_eval)

        CatalystBed._f_term.terminal = self.terminal

        self.ivp_solution = solve_ivp(self.ode_system, (0, self.z), self._inlet_values,
                                      events=(self._f_term), t_eval=t_eval,
                                      rtol=self._ivp_rtol, **kwargs)

        self._outlet_values = self.ivp_solution.y[:, -1]
        self.outlet = dict(zip(self._keys, self._outlet_values))


    def get_outlet(self, **options):
        """Obtain reactor outlet dictionary.
        
        **options are passed to scipy solve_ivp
        Returns
        -------
        dict
            Reactor outlet
        """
        try:
            return self.outlet
        except:
            self.solve(**options)
            return self.outlet


    def get_pre_heating(self, T_before):
        T0 = self.inlet['T']
        h_before = calc_hf_temp(T_before, HF298, A, B, C, D)
        h0 = calc_hf_temp(T0, HF298, A, B, C, D)
        self.pre_heat = (h0 - h_before).dot(self._inlet_values[:-2]) / 3.6e6
        return self.pre_heat  # MW

        
    def get_dataframe(self, **options):
        if self.ivp_solution is None:
            self.solve(**options)
        df = pd.DataFrame(self.ivp_solution.y.T, columns=self._keys)
        df['z'] = self.ivp_solution.t
        df.set_index('z', inplace=True)
        return df
    

class AxialBed(CatalystBed):

    def __init__(
        self, z,
        rhos=950.0, rhob=570.0, es=0.4,
        dp=0.425e-3, tao=3.0, inner_R=9.2e-3,
        **options
    ):
        """Class for catalyst bed.
        Parameters
        ----------
        z        : float
                   Reactor length [m].
        rhos     : float, optional
                   Catalyst solid density [kg/m3], by default 950.0
        rhob     : float, optional
                   Catalyst bulk density [kg/m3], by default 570.0
        es       : float, optional
                   Catalyst solid void fraction, by default 0.4
        dp       : float, optional
                   Pellet equivalent diameter [m], by default 0.425e-3
        tao      : float, optional
                   Pellet tortuosity, by default 3.0
        inner_R  : float, optional
                   Catalyst bed inner radius [m], by default 9.2e-3
        Pmin     : float, optional
                   Minimum pressure allowed [bar]. Useful in optimization. By default 0.5
        Pterm    : float or None, optional
                   Terminal pressure to interrupt ODE system [bar], by default None
        ivp_rtol : float, optional
                   Relative tolerance for ODE system, by default 1e-6
        """

        super(AxialBed, self).__init__(
            z, rhos=rhos, rhob=rhob, es=es, dp=dp, tao=tao, 
            inner_R=inner_R, **options
        )

        self.Ac = math.pi*inner_R**2
        self.eff1 = 1.086e-2
        self.eff2 = 0.0
        self.eff3 = 1.0
        self.eff4 = 1.0


    def ode_system(self, z, params):

        F = np.array(params[0:-2])
        T = params[-2]
        P = params[-1]
        Ft = F.sum(axis=-1)
        y = F / Ft
        pp = y * P
        
        """factor to convert dW = dz.rhob.Ac"""
        factor = self.Ac*self.rhob

        rxn_SRM1 = self.eff1*rxn_srm1(pp, T)*factor
        rxn_SRM2 = self.eff2*rxn_srm2(pp, T)*factor
        rxn_WGS = self.eff3*rxn_wgs(pp, T)*factor
        rxn_DRM = self.eff4*rxn_drm(pp, T)*factor
        
        dF = np.zeros(len(F))
        dF[CH4] = -rxn_SRM1 - rxn_SRM2 - rxn_DRM
        dF[CO2] = rxn_SRM2 + rxn_WGS - rxn_DRM
        dF[H2O] = -rxn_SRM1 -2*rxn_SRM2 - rxn_WGS
        dF[H2] = 3*rxn_SRM1 + 4*rxn_SRM2 + rxn_WGS + 2*rxn_DRM
        dF[CO] = rxn_SRM1 - rxn_WGS + 2*rxn_DRM
        dF[AR] = 0

        Cp = get_Cp(T, A, B, C, D)
        delta_hr = calc_delta_hr(T, DELTA_HR298, DELTA_A, DELTA_B, DELTA_C, DELTA_D)
        dT = -factor*(np.array([rxn_SRM1, rxn_SRM2, rxn_WGS, rxn_DRM]).T.dot(delta_hr)) / (Cp.T.dot(F))

        mu = calc_mu_mist(T, y, MM, TC, PC, SIGMA, EK, DELTA_POT)
        rhog = (y.T.dot(np.array(MM))) * P / 8.314e-2 / T
        G = (np.array(F).T.dot(np.array(MM))) / self.Ac
        dP = -calc_pressure_drop(G, rhog, mu, self.Ac, self.dp, self.rhob, self.eg)*factor
        
        if P < self.Pmin:
            dF = np.zeros(len(dF))
            dT = 0
            if P < self.Pterm:
                dP = 0

        return np.append(dF, [dT, dP])
    

class Cooler(CatalystBed):

    def __init__(
        self, z,
        rhos=950.0, rhob=570.0, es=0.4,
        dp=0.425e-3, tao=3.0, inner_R=9.2e-3,
        **options
    ):
        """Class for catalyst bed with a cooloer incorporated.
        Parameters
        ----------
        z        : float
                   Reactor length [m].
        rhos     : float, optional
                   Catalyst solid density [kg/m3], by default 950.0
        rhob     : float, optional
                   Catalyst bulk density [kg/m3], by default 570.0
        es       : float, optional
                   Catalyst solid void fraction, by default 0.4
        dp       : float, optional
                   Pellet equivalent diameter [m], by default 0.425e-3
        tao      : float, optional
                   Pellet tortuosity, by default 3.0
        inner_R  : float, optional
                   Catalyst bed inner radius [m], by default 9.2e-3
        Pmin     : float, optional
                   Minimum pressure allowed [bar]. Useful in optimization. By default 0.5
        Pterm    : float or None, optional
                   Terminal pressure to interrupt ODE system [bar], by default None
        ivp_rtol : float, optional
                   Relative tolerance for ODE system, by default 1e-6
        """

        super(Cooler, self).__init__(
            z, rhos=rhos, rhob=rhob, es=es, dp=dp, tao=tao, 
            inner_R=inner_R, **options
        )

        self.Ac = math.pi*inner_R**2
        """overall heat transfer coefficient [J/m2/s/K]"""
        self.U = 500
        """cooler is at ambient temperature"""
        self.Tcool = 298.15
        self.eff1 = 1.086e-2
        self.eff2 = 0.0
        self.eff3 = 1.0
        self.eff4 = 1.0


    def ode_system(self, z, params):

        F = np.array(params[0:-2])
        T = params[-2]
        P = params[-1]
        Ft = F.sum(axis=-1)
        y = F / Ft
        pp = y * P
        
        """factor to convert dW = rhob.Ac.dz"""
        factor = self.Ac*self.rhob

        rxn_SRM1 = self.eff1*rxn_srm1(pp, T)*factor
        rxn_SRM2 = self.eff2*rxn_srm2(pp, T)*factor
        rxn_WGS = self.eff3*rxn_wgs(pp, T)*factor
        rxn_DRM = self.eff4*rxn_drm(pp, T)*factor
        
        dF = np.zeros(len(F))
        dF[CH4] = -rxn_SRM1 - rxn_SRM2 - rxn_DRM
        dF[CO2] = rxn_SRM2 + rxn_WGS - rxn_DRM
        dF[H2O] = -rxn_SRM1 -2*rxn_SRM2 - rxn_WGS
        dF[H2] = 3*rxn_SRM1 + 4*rxn_SRM2 + rxn_WGS + 2*rxn_DRM
        dF[CO] = rxn_SRM1 - rxn_WGS + 2*rxn_DRM
        dF[AR] = 0

        
        Cp = get_Cp(T, A, B, C, D)
        delta_hr = calc_delta_hr(T, DELTA_HR298, DELTA_A, DELTA_B, DELTA_C, DELTA_D)
        cooler = ((4*self.Ac*self.U/2*self.inner_R)*(T-self.Tcool))/(Cp.T.dot(F))
        dT = -factor*(np.array([rxn_SRM1, rxn_SRM2, rxn_WGS, rxn_DRM]).T.dot(delta_hr))\
            / (Cp.T.dot(F)) - cooler

        mu = calc_mu_mist(T, y, MM, TC, PC, SIGMA, EK, DELTA_POT)
        rhog = (y.T.dot(np.array(MM))) * P / 8.314e-2 / T
        G = (np.array(F).T.dot(np.array(MM))) / self.Ac
        dP = -calc_pressure_drop(G, rhog, mu, self.Ac, self.dp, self.rhob, self.eg)*factor
        
        if P < self.Pmin:
            dF = np.zeros(len(dF))
            dT = 0
            if P < self.Pterm:
                dP = 0

        return np.append(dF, [dT, dP])


class SingleBed(object):
    
    def __init__(self):
        """
        Creates multibed reactor.
        """
        self.n_beds = 0
        self.beds = {}
        self._keys = ['Fch4', 'Fco2', 'Fh2o', 'Fh2', 'Fco', 'Far', 'T', 'P']

    
    def add_axial_bed(
        self, z,
        rhos=950.0, rhob=570.0, es=0.4,
        dp=0.425e-3, tao=3.0, inner_R=9.2e-3,
        **options
    ):
        """Class for catalyst bed.
        Parameters
        ----------
        z        : float
                   Reactor length [m].
        rhos     : float, optional
                   Catalyst solid density [kg/m3], by default 950.0
        rhob     : float, optional
                   Catalyst bulk density [kg/m3], by default 570.0
        es       : float, optional
                   Catalyst solid void fraction, by default 0.4
        dp       : float, optional
                   Pellet equivalent diameter [m], by default 0.425e-3
        tao      : float, optional
                   Pellet tortuosity, by default 3.0
        inner_R  : float, optional
                   Catalyst bed inner radius [m], by default 9.2e-3
        Pmin     : float, optional
                   Minimum pressure allowed [bar]. Useful in optimization. By default 0.5
        Pterm    : float or None, optional
                   Terminal pressure to interrupt ODE system [bar], by default None
        ivp_rtol : float, optional
                   Relative tolerance for ODE system, by default 1e-6
        """
        self.n_beds = self.n_beds + 1
        self.beds[self.n_beds] = AxialBed(
            z, rhos=rhos, rhob=rhob, es=es, dp=dp, tao=tao, inner_R=inner_R, **options
        )

    def add_cooler(
        self, z,
        rhos=950.0, rhob=570.0, es=0.4,
        dp=0.425e-3, tao=3.0, inner_R=9.2e-3,
        **options
    ):
        """Class for catalyst bed with a cooler incorporated.
        Parameters
        ----------
        z        : float
                   Reactor length [m].
        rhos     : float, optional
                   Catalyst solid density [kg/m3], by default 950.0
        rhob     : float, optional
                   Catalyst bulk density [kg/m3], by default 570.0
        es       : float, optional
                   Catalyst solid void fraction, by default 0.4
        dp       : float, optional
                   Pellet equivalent diameter [m], by default 0.425e-3
        tao      : float, optional
                   Pellet tortuosity, by default 3.0
        inner_R  : float, optional
                   Catalyst bed inner radius [m], by default 9.2e-3
        Pmin     : float, optional
                   Minimum pressure allowed [bar]. Useful in optimization. By default 0.5
        Pterm    : float or None, optional
                   Terminal pressure to interrupt ODE system [bar], by default None
        ivp_rtol : float, optional
                   Relative tolerance for ODE system, by default 1e-6
        """
        self.n_beds = self.n_beds + 1
        self.beds[self.n_beds] = Cooler(
            z, rhos=rhos, rhob=rhob, es=es, dp=dp, tao=tao, inner_R=inner_R, **options
        )

    def set_inlet(self, Fch4=30, Fco2=20, Fh2o=1e-2, Fh2=1e-2, Fco=0.0,
                  Far=50, T=773, P=1.0):
        """Set inlet conditions of the first catalyst bed.
        Parameters
        ----------
        Fch4 : float, optional
             Methane feed ratio [ml/min], by default 30
        Fco2 : float, optional
             Carbon dioxide feed ratio [ml/min], by default 
        Fh2o : float, optional
             Water feed ratio [ml/min], by default 0.0
        Fh2  : float, optional
             Hydrogen feed ratio [ml/min], by default 0.0
        Fco  : float, optional
             Carbon monoxide feed ratio [ml/min], by default 0.0
        Far  : float, optional
             Argon [ml/min], by default 50
        T    : int, optional
             Temperature [K], by default 773
        P    : float, optional
             Pressure [bar], by default 1.0
        """

        if self.n_beds == 0:
            raise IndexError('Please add a bed before the inlet')

        self.beds[1].set_inlet(
            Fch4=Fch4,
            Fco2=Fco2,
            Fh2o=Fh2o,
            Fh2=Fh2,
            Fco=Fco,
            Far=Far,
            T=T,
            P=P)

        self.inlet = self.beds[1].inlet.copy()


    def solve(self, **options):
        """Solves sequence of ODE systems with given adjustments. 
        
        **options is parsed to scipy solve_ivp.
        Raises
        ------
        IndexError
            If no bed was included yet
        """

        self.beds[1].solve(**options)

        if self.n_beds == 0:
            raise IndexError('Please add a bed before the inlet')

        elif self.n_beds == 1:
            self.outlet = self.beds[1].outlet
            
        elif self.n_beds > 1:
            raise IndexError('You can only one bed')


    def get_outlet(self, **options):
        """Returns dictionary in reactor outlet.
        Returns
        -------
        dict
            Reactor outlet
        """
        try:
            return self.outlet
        except:
            self.solve(**options)
            return self.outlet


    def get_heat_consumed(self, initial_T=781.00):
        if self.n_beds == 0:
            raise IndexError('Please add a bed before the inlet')
        elif self.n_beds == 1:
            return self.beds[1].get_pre_heating(initial_T)
        
    @property
    def conversion(self):    
        Xch4 = 100*(self.inlet['Fch4'] - self.outlet['Fch4'])/self.inlet['Fch4']
        Xco2 = 100*(self.inlet['Fco2'] - self.outlet['Fco2'])/self.inlet['Fco2']
        Xh2 = 100*(self.outlet['Fh2'])/(2*self.inlet['Fch4'])
        Xco = 100*(self.outlet['Fco'])/(self.inlet['Fco2']+self.inlet['Fch4'])
        return {'Xch4': Xch4, 'Xco2': Xco2, 'Xh2': Xh2, 'Xco': Xco}

        
    def get_bed_dataframe(self, bed_number, **options):
        """Returns the pandas.DataFrame of profiles in a given catalyst bed.
        Parameters
        ----------
        bed_number : int
            Catalyst bed number (starts at 1).
        Returns
        -------
        pandas.DataFrame
            Profiles along reactor length.
        """
        try:
            return self.beds[bed_number].get_dataframe()
        except:
            self.solve(**options)
            return self.beds[bed_number].get_dataframe()


    def get_dataframe(self, **options):
        """Returns the pandas.DataFrame of profiles in the reactor.
        Returns
        -------
        pandas.DataFrame
            Profiles along reactor length.
        """
        df = self.get_bed_dataframe(1, **options)
        z = df.index.values
        df['z'] = z
        df.set_index('z', inplace=True)
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
                                                   sheet_name='SingleBbed')
