# biogas
Thermodynamic and kinetic analysis of Biogas Dry Reforming

[Install](#Install) / [Usage](#Usage) / [Equilibrium](#Equilibrium) / [Reactor](#Reactor) / [Contact](#Contact)

# Install
```Python
pip install -e git+https://github.com/amvro23/biogas/#egg=biogas
```

Note: It might be useful to write "git+https://github.com/amvro23/biogas/#egg=biogas" if installing directly from a Python interpreter as # can be interpreted as a comment.

# Usage
```Python
import matplotlib.pyplot as plt
import numpy as np
from biogas.reactor import SingleBed
from biogas.gibbs_min import Equilibrium
import warnings
warnings.filterwarnings("ignore")
```
# Equilibrium
Create an object for the gas mixture of biogas (default values are T = 500K, P = 1atm, P0 = 1atm).
```
mix = Equilibrium()
```
Define the temperature range for equilibrium analysis (default is 573.15K to 1213.15K)
```Python
mix.T_range = np.linspace(573.15, 1213.15)
```
Define the inlet ratio of biogas mixture which can either be diluted (default values),
```Python
mix.set_inlet(ych4=0.3, yco2=0.2, yar=0.5)
```

or pure biogas (default value is 50% biogas & 50% argon inert gas).
```Python
mix.set_inlet(ych4=0.6, yco2=0.4, yar=0.0)
```
You can create a plot regarding the desired gas mixture (e.g., diluted mixture).
```Python
mix.plot_molar_ratio()
```
![molar](https://user-images.githubusercontent.com/91277572/211278653-cd7515ce-0b07-49b3-a556-2a51cee4a74c.png)

You can also create a plot for the gas mixture conversions, yields and ratios (e.g., diluted mixture).
```Python
mix.plot_conversions_yields()
```
![conv](https://user-images.githubusercontent.com/91277572/211278532-7adfea83-c938-418a-81ec-f1eca6f4b500.png)

The package allows you to consider carbon during the thermodynamic analysis.

First, you have to add carbon.
```Python
mix.add_carbon()
```
Then, by applying that same logic you can have access to molar composition plots,
```Python
mix.plot_molar_ratio_carbon()
```
![eq_comp](https://user-images.githubusercontent.com/91277572/211210025-6f60cfa8-072c-403b-aee1-3cde791a9f7c.png)

or, you can create a plot for the gas mixture conversions, yields and ratios.
```Python
mix.plot_molar_ratio_carbon()
```
![eqconv](https://user-images.githubusercontent.com/91277572/211210090-1c8d4313-a6c9-493b-a703-afb751728621.png)

You can also obtain a dataframe of all the obtained numerical values of conversions, yields, molar ratio compositions etc., for both cases (i.e., with and without carbon) by using the following.
```Python
mix.get_dataframe()
mix.get_dataframe_carbon()
```
Finally, you can get the excel files with the dataframes to create your own plots or manipulate data, if necessary.
```Python
mix.to_excel('GasEquilibrium')
mix.to_excel_carbon('GasEquilibriumCarbon')
```

# Reactor
First you have to create an object for the test reactor.
```Python
test_reac = SingleBed()
```
Next you have to add the single bed reactor with one required parameter z (i.e., length of the reactor). Other optional parameters are (rhos=950.0, rhob=570.0, es=0.4, dp=0.425e-3, tao=3.0, inner_R=9.2e-3, Pmin=0.5) which are described in greater detail within the package.

You can either choose a single axial bed reactor,
```Python
test_reac.add_axial_bed(5.1e-3)
```
or you can take advantage of the flexibility of the package and incorporate an external cooling fluid of variable temperature. Default value of cooling temperature is the ambient temperature (T=298K)
```Python
test_reac.add_cooler(5.1e-3)
test_reac.beds[1].Tcool = 298
```
Then you have to set inlet conditions (50% biogas & 50% inert gas). Default values are (Fch4=30, Fco2=20, Fh2o=1e-2, Fh2=1e-2, Fco=0.0, Far=50, T=773, P=1.0) with molar rates (Fi) in [ml/min], temperature (T) in K, and pressure (P) in bar. Note that small amounts of hydrogen and water should be inserted too in order to avoid division with zero and maintain numerical stability.
```Python
test_reac.set_inlet()
```
After that, you have to solve the system of differential equations,
```Python
test_reac.solve()
```
and save the results obtained in a variable.
```Python
profiles = test_reac.get_dataframe()
```
From the created profiles your can obtain the plots of molar rates and temperature profiles.
For example, for the axial bed reactor at default values:
```Python
fig, ax = plt.subplots(figsize=[6, 4], dpi=100, sharex=True)
ax.plot(profiles.index, profiles["Fch4"], color="darkgreen", label="$CH_4$")
ax.plot(profiles.index, profiles["Fco2"], color="black", label="$CO_2$")
ax.plot(profiles.index, profiles["Fh2"], color="red", label="$H_2$")
ax.plot(profiles.index, profiles["Fco"], color="orange", label="$CO$")
ax.set_ylabel("$F$ [mol/s]")
ax.set_xlabel("$z$ [m]")
ax.legend()
ax.grid(ls=":")
fig.tight_layout()
plt.show()
```
![773_mr](https://user-images.githubusercontent.com/91277572/211167335-dcc1993a-fccf-423d-9de2-f765b202c21e.png)
![773T_mr](https://user-images.githubusercontent.com/91277572/211167369-d9ce6a15-929e-4e0c-8a56-03c507dfa6e5.png)

or for the axial bed reactor at T=950K,

```Python
test_reac.add_axial_bed(5.1e-3)
test_reac.set_inlet(T=950)
test_reac.solve()
profiles = test_reac.get_dataframe()

fig, ax = plt.subplots(figsize=[6, 4], dpi=100, sharex=True)
ax.plot(profiles.index, profiles["Fch4"], color="darkgreen", label="$CH_4$")
ax.plot(profiles.index, profiles["Fco2"], color="black", label="$CO_2$")
ax.plot(profiles.index, profiles["Fh2"], color="red", label="$H_2$")
ax.plot(profiles.index, profiles["Fco"], color="orange", label="$CO$")
ax.set_ylabel("$F$ [mol/s]")
ax.set_xlabel("$z$ [m]")
ax.legend()
ax.grid(ls=":")
fig.tight_layout()
plt.show()
```
![950T_mr](https://user-images.githubusercontent.com/91277572/211167584-4a1a3396-0e7a-491b-adb7-7feacb7bd531.png)

or for the cooler using the default values.

```Python
test_reac.add_axial_bed(5.1e-3)
test_reac.set_inlet()
test_reac.solve()
profiles = test_reac.get_dataframe()

fig, ax = plt.subplots(figsize=[6, 4], dpi=100, sharex=True)
ax.plot(profiles.index, profiles["T"], color="purple", label="$T$")
ax.set_ylabel("$T$ [K]")
ax.set_xlabel("$z$ [m]")
ax.legend()
ax.grid(ls=":")
fig.tight_layout()
plt.show()
```
![773T_mr](https://user-images.githubusercontent.com/91277572/211167473-ed98254c-d83a-4060-bd2d-c2b15a2bbd46.png)

The effectiveness factors can also be adjusted accordingly.
```Python
test_reac = SingleBed()

test_reac.add_axial_bed(5.1e-3)
test_reac.set_inlet(T=1000)

test_reac.beds[1].eff1 = 0
test_reac.beds[1].eff2 = 0
test_reac.beds[1].eff3 = 1
test_reac.beds[1].eff4 = 1

test_reac.solve()
profiles = test_reac.get_dataframe()


fig, ax = plt.subplots(figsize=[6, 4], dpi=100, sharex=True)
ax.plot(profiles.index, profiles["Fch4"], color="darkgreen", label="$CH_4$")
ax.plot(profiles.index, profiles["Fco2"], color="black", label="$CO_2$")
ax.plot(profiles.index, profiles["Fh2"], color="red", label="$H_2$")
ax.plot(profiles.index, profiles["Fco"], color="orange", label="$CO$")
ax.set_ylabel("$F$ [mol/s]")
ax.set_xlabel("$z$ [m]")
ax.legend()
ax.grid(ls=":")
fig.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=[6, 4], dpi=100, sharex=True)
ax.plot(profiles.index, profiles["T"], color="purple", label="$T$")
ax.set_ylabel("$T$ [K]")
ax.set_xlabel("$z$ [m]")
ax.legend()
ax.grid(ls=":")
fig.tight_layout()
plt.show()
```

![mf](https://user-images.githubusercontent.com/91277572/211168849-c742187d-93e0-421c-ab06-1c5ae15cf8ed.png)
![temp](https://user-images.githubusercontent.com/91277572/211168846-f7c818e5-6e09-4c8c-be4d-772787affac2.png)

The user can also have access to the values of conversions and yields (e.g., at T = 950K)

```Python
conversions = test_reac.conversion
conversions
```

```
Out: 
{'Xch4': 55.10987440137124,
 'Xco2': 88.3423110736619,
 'Xh2': 53.23404124416936,
 'Xco': 68.40284907028743}
```

The package also returns the heat consumption in MW with default value of (initial_T=781) to give an indicative overview of the εnergy requirements.

```Python
test_reac.get_heat_consumed(initial_T=781.00)
```

```
Out
2.792282394738481e-07
```

Finally, you can get an excell file of the molar rates, temperature, and pressure profiles to reate your own plots.

```Python
test_reac.to_excel("SingleBed")
```

The kinetic model, effectiveness factors, reactor characteristics and particl properties were reproduced based on the Charisiou et al. work, while other equations and parameters used herein were adopted from the Leite et al. articles.

# References
[Charisiou, N. D., Siakavelas, G., Papageridis, K. N., Baklavaridis, A., Tzounis, L., Avraam, D. G., & Goula, M. A. (2016). Syngas production via the biogas dry reforming reaction over nickel supported on modified with CeO2 and/or La2O3 alumina catalysts. Journal of Natural Gas Science and Engineering, 31, 164-183.](https://doi.org/10.1016/j.jngse.2016.02.021)

[Leite, B., Costa, A. O. S., Costa, E. F., 2023. Multi-objective optimization of adiabatic styrene reactors using Generalized Differential Evolution 3 (GDE3). Chem. Eng. Sci., Volume 265, Article 118196.](https://doi.org/10.1016/j.ces.2022.118196)

[Leite, B., Costa, A. O. S. & Costa Junior, E. F., 2021. Simulation and optimization of axial-flow and radial-flow reactors for dehydrogenation of ethylbenzene into styrene based on a heterogeneous kinetic model. Chem. Eng. Sci., Volume 244, Article 116805.](https://doi.org/10.1016/j.ces.2021.116805)

# Contact
amvro23@gmail.com
