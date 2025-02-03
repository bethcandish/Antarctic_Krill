import numpy as np
import matplotlib.pyplot as plt
from onekrill_onecolumn import (
    calc_clearance_rate,
    calc_krill_mp_consumption,
    calc_mp_fp_production_rate,
    calc_initial_sinking_velocity_at,
    calc_initial_sinking_velocity_ko,
)

# Parameters
krill_length_mm = 50  # mm
mp_conc = 200  # particles/m3
egestion_rate = 100  # ugCkrill-1d-1 (one krill)
gut_passage_time = 4  # hours

depth = 650  # m
b = 0.05  # Attenuation of velocity with depth

mu = 0.001  # Viscosity of water
rho = 1025  # Density of water
rho_s = 1121  # Density of krill FP (Atkinson et al 2012)
L_ko = 2928 * 10 ** (-6)  # Length of FP (Atkinson et al 2012)
D_ko = 183 * 10 ** (-6)  # Width/diameter of FP (Atkinson et al 2012)

L_at = 2928
D_at = 183
rho_s = 1121  # all Atkinson et al 2012 values
time = np.linspace(0, 500, 500)  # Simulation time in hours

# Compute rates
clearance_rate = calc_clearance_rate(krill_length_mm)
krill_mp_consumption = calc_krill_mp_consumption(clearance_rate, mp_conc)
mp_fp_production_rate = calc_mp_fp_production_rate(krill_mp_consumption, gut_passage_time)
initial_sinking_velocity_at = calc_initial_sinking_velocity_at(L_at, D_at, rho_s)
initial_sinking_velocity_ko = calc_initial_sinking_velocity_ko(mu, rho, rho_s, L_ko, D_ko)

## TRACKING MP ACCUMULATION AT 650M ###

# Convert sinking velocity from m/day to m/hour
sinking_velocity_hour_at = initial_sinking_velocity_at / 24 
sinking_velocity_hour_ko = initial_sinking_velocity_ko / 24

# Time for fecal pellets to reach 650m
time_to_650m_hours_at = depth / sinking_velocity_hour_at
time_to_650m_hours_ko = depth / sinking_velocity_hour_ko

# Convert MP FP production rate to MPs/day
mp_fp_production_rate_day = mp_fp_production_rate * 24  

# Time array in days
time_days = np.linspace(0, 200, 500)  # 200 days, 500 time points

# Initialize MP accumulation arrays
mp_accumulation_at = np.zeros_like(time_days)
mp_accumulation_ko = np.zeros_like(time_days)

# Accumulate MPs over time
for i, t in enumerate(time_days):
    if t >= time_to_650m_hours_at / 24:  # MPs start arriving after travel time
        mp_accumulation_at[i] = mp_accumulation_at[i - 1] + mp_fp_production_rate_day 
    if t >= time_to_650m_hours_ko / 24:  # MPs start arriving after travel time
        mp_accumulation_ko[i] = mp_accumulation_ko[i - 1] + mp_fp_production_rate_day 

# Calculate the difference in MP accumulation
mp_accumulation_diff = mp_accumulation_at - mp_accumulation_ko

# Plot results
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

# First plot: MP accumulation over time
ax1.plot(time_days, mp_accumulation_at, label="Atkinson et al. (2012) Velocity", color="blue")
ax1.plot(time_days, mp_accumulation_ko, label="Ko Model Velocity", color="red")
ax1.set_ylabel("Accumulated MP at 650m")
ax1.set_xlabel("Time (days)")
ax1.set_title("Accumulation of Microplastics at 650m Due to Sinking Fecal Pellets")
ax1.legend()
ax1.grid()

# Second plot: Difference in MP accumulation over time
ax2.plot(time_days, mp_accumulation_diff, label="Difference in MP Accumulation", color="green", linestyle="dashed")
ax2.set_xlabel("Time (days)")
ax2.set_ylabel("Difference in MP Accumulation")
ax2.set_title("Difference in MP Accumulation Over Time")
ax2.legend()
ax2.grid()

plt.tight_layout()
plt.show()

print(mp_accumulation_diff[-1])
print(mp_accumulation_at[-1])