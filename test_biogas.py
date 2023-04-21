import matplotlib.pyplot as plt
import numpy as np
from biogas.reactor import SingleBed
from biogas.gibbs_min import Equilibrium
import warnings
warnings.filterwarnings("ignore")

test_reac = SingleBed()

test_reac.add_axial_bed(5.1e-3)

test_reac.set_inlet()

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
ax.plot(profiles.index[1:], profiles["T"].iloc[1:], color="purple", label="$T$")
ax.set_ylabel("$T$ [K]")
ax.set_xlabel("$z$ [m]")
ax.legend()
ax.grid(ls=":")
fig.tight_layout()
plt.show()