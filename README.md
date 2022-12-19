# Biogas_reforming
Thermodynamic and kinetic analysis of Biogas Dry Reforming

# Equilibrium - Usage

```Python
from biogas_gibbs_minimization import CarbonEquilibrium
```
Create an object for the desired temperature and pressure of the system (e.g., 800 K, 1 atm)
```
mix = CarbonEquilibrium(800, 1)
```
Create plot for the gas mixture compositions
```Python
mix.plot_molar_ratio
```
Create plot for the gas mixture conversions, yields and ratios
```Python
mix.plot_conversions_yields_carbon
```
Create plot for the gas mixture compositions considering carbon formation
```Python
mix.plot_molar_ratio_carbon
```
Create plot for the gas mixture conversions, yields and ratios considering carbon formation
```Python
mix.plot_molar_ratio_carbon
```
