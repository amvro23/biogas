# Biogas_reforming
Thermodynamic and kinetic analysis of Biogas Dry Reforming

# Install
```Python
pip install -e git+https://github.com/amvro23/biogas/#egg=biogas
```

Note: It might be useful to write "git+https://github.com/amvro23/biogas/#egg=biogas" if installing directly from a Python interpreter as # can be interpreted as a comment.

# Equilibrium - Usage

```Python
from gibbs_min import CarbonEquilibrium
```
Create an object for the gas mixture of biogas (default values are T = 500K, P = 1atm, P0 = 1atm).
```
mix = CarbonEquilibrium()
```
Define the temperature range for equilibrium analysis (default is 573.15K to 1213.15K)
```Python
mix.T_range = np.linspace(573.15, 1213.15)
```
Define the inlet ratio of biogas mixture which can be diluted or not (default value is 50% biogas & 50% argon inert gas).
```Python
mix.inlet = np.array([0.3, 0.2, 0.5]) # for diluted mixture
```
```Python
mix.inlet = np.array([0.6, 0.4, 0.0]) # for pure biogas
```
Create plot for the gas mixture compositions.
```Python
mix.plot_molar_ratio
```
![molar_ratio_gas](https://user-images.githubusercontent.com/91277572/208469749-f7682117-3dae-471d-bbfe-5d32b41e7533.png)

Create plot for the gas mixture conversions, yields and ratios.
```Python
mix.plot_conversions_yields
```
![conversions_yields_gas](https://user-images.githubusercontent.com/91277572/208469178-6e58a363-3ff3-46bc-80a9-f09ca3b23ffe.png)

Create plot for the gas mixture compositions considering carbon.
```Python
mix.plot_molar_ratio_carbon
```
![molar_ratio_carbon](https://user-images.githubusercontent.com/91277572/208422921-4a0a8ed5-585d-45b9-9cfa-d442563caef5.png)

Create plot for the gas mixture conversions, yields and ratios considering carbon.
```Python
mix.plot_molar_ratio_carbon
```
![conversions_yields_carbon](https://user-images.githubusercontent.com/91277572/208424006-013497e9-451a-496b-b252-3d88f6cbb111.png)

You can also obtain numerical values of conversions, yields, molar ratio compositions etc., for both cases (i.e., with and without carbon) by using the following:
```Python
mix.conversions
mix.conversions_
mix.yields
mix.yields_carbon
mix.ratio_h2_co
mix.ratio_h2_co_carbon
mix.multiple_temperatures(mix.T_range) # molar compositions for gas mixture
mix.multiple_temperatures_carbon(mix.T_range) # molar compositions for gas mixture including carbon
```
