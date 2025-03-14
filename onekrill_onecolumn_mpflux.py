import numpy as np
import matplotlib.pyplot as plt
from onekrill_onecolumn import (
    calc_clearance_rate,
    calc_krill_mp_consumption,
    calc_mp_fp_production_rate,
    calc_initial_sinking_velocity_ko,
    calc_flux
)

# Parameters
krill_length_mm = 50  # mm
mp_conc = 200  # particles/m3
egestion_rate = 100  # ugCkrill-1d-1 (one krill)
gut_passage_time = 4  # hours

depth = 50  # m (changed to 50m for mass accumulation calculation)
b = 0.36 # Attenuation of velocity with depth

mu = 0.001  # Viscosity of water
rho = 1025  # Density of water
rho_s = 1121  # Density of krill FP (Atkinson et al 2012)
L_ko = 2928 * 10 ** (-6)  # Length of FP (Atkinson et al 2012)
D_ko = 183 * 10 ** (-6)  # Width/diameter of FP (Atkinson et al 2012)

time = np.linspace(0, 500, 500)  # Simulation time in hours

# Compute rates
clearance_rate = calc_clearance_rate(krill_length_mm)
krill_mp_consumption = calc_krill_mp_consumption(clearance_rate, mp_conc)
mp_fp_production_rate = calc_mp_fp_production_rate(krill_mp_consumption, gut_passage_time)
initial_sinking_velocity_ko = calc_initial_sinking_velocity_ko(mu, rho, rho_s, L_ko, D_ko)

# Convert sinking velocity from m/day to m/hour
sinking_velocity_hour_ko = initial_sinking_velocity_ko / 24

# Time for fecal pellets to reach 50m
time_to_50m_hours_ko = depth / sinking_velocity_hour_ko

# Convert fecal pellet production rate to pellets per day
fp_production_rate_day = mp_fp_production_rate * 24  

# Volume and mass of a single fecal pellet
fp_volume = np.pi * (D_ko / 2) ** 2 * L_ko  # m^3
fp_mass = rho_s * fp_volume  # kg

# Compute fecal pellet mass flux at 50m (kg/day)
fp_mass_flux_50m = fp_production_rate_day * fp_mass

# Compute microplastic flux at 50m (particles/day)
mp_flux_50m = fp_production_rate_day  # Assuming each FP carries MPs

z = np.linspace(1, 600, 1000)

fp_flux = calc_flux(fp_mass_flux_50m, z, b)
mp_flux = calc_flux(mp_flux_50m, z, b)

# Plotting with two x-axes
fig, ax1 = plt.subplots(figsize=(7, 5))

color = 'brown'
ax1.set_xlabel("Flux")
ax1.set_ylabel("Depth (m)", color=color)
ax1.plot(fp_flux, (z + 50), label="FP Mass Flux", color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.invert_yaxis()  # Invert y-axis so depth increases downward

ax2 = ax1.twiny()  # Create a second x-axis sharing the y-axis
color = 'blue'
ax2.set_xlabel("Microplastic Flux", color=color)
ax2.plot(mp_flux, (z + 50), label="Microplastic Flux", color=color, linestyle='dashed')
ax2.tick_params(axis='x', labelcolor=color)

fig.suptitle("Fecal Pellet Mass Flux and Microplastic Flux with Depth")
fig.tight_layout()
plt.grid()
plt.show()