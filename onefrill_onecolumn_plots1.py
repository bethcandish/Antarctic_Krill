import numpy as np
import matplotlib.pyplot as plt
from onekrill_onecolumn import (
    calc_clearance_rate,
    calc_krill_mp_consumption,
    calc_mp_fp_production_rate,
    calc_sinking_velocity,
    calc_initial_sinking_velocity_ko,
    calc_initial_sinking_velocity_at
)

# Parameters
krill_length_mm = 50  # mm
mp_conc = 500  # particles/m3
depth = 650  # m
time = np.linspace(0, 200, 500)  # Simulation time in hours
b = 0.05  # Attenuation of velocity with depth
mu = 0.001  # Viscosity of water
rho = 1025  # Density of water
rho_s = 1121  # Density of krill FP (Atkinson et al 2012)
L = 2928 * 10 ** (-6)  # Length of FP (Atkinson et al 2012)
D = 183 * 10 ** (-6)  # Width/diameter of FP (Atkinson et al 2012)

# L = 2928
# D = 183
# rho_s = 1121 # all Atkinson et al 2012 values

gut_passage_time = 4  # Gut passage time in hours

# Compute rates
clearance_rate = calc_clearance_rate(krill_length_mm)
krill_mp_consumption = calc_krill_mp_consumption(clearance_rate, mp_conc)
time_produce_one_mp_fp = calc_mp_fp_production_rate(krill_mp_consumption, gut_passage_time)
initial_sinking_velocity = calc_initial_sinking_velocity_ko(mu, rho, rho_s, L, D)

# Generate fecal pellet release times
fp_release_times = np.arange(0, max(time), gut_passage_time)

# Generate MP FP release times (fixed intervals)
mp_fp_release_times = np.arange(time_produce_one_mp_fp, max(time), time_produce_one_mp_fp)

# Create figure
fig, ax = plt.subplots(figsize=(12, 10))

# Simulate sinking for each fecal pellet
for release_time in fp_release_times:
    pellet_time = time[time >= release_time]  # Time after release
    time_since_release = pellet_time - release_time  # Time elapsed since egestion

    # Apply exponential decay to sinking velocity
    sinking_velocity = calc_sinking_velocity(time_since_release, initial_sinking_velocity, b)
    sinking_velocity_per_hour = sinking_velocity / 24
    
    # Compute depth with decayed sinking velocity
    pellet_depth = np.cumsum(sinking_velocity_per_hour) * (time[1] - time[0])  # Numerical integration
    
    # Stop tracking pellets that exceed 650m
    pellet_depth = np.clip(pellet_depth, 0, depth)

    # Determine if this FP is an MP FP
    contains_mp = np.any(np.isclose(release_time, mp_fp_release_times, atol=gut_passage_time / 2))
    color = "red" if contains_mp else "blue"  # Red for MP, Blue for normal

    ax.plot(pellet_time, pellet_depth, color=color)

# Formatting
ax.invert_yaxis()  # Invert y-axis to show depth increasing downward
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Depth (m)")
ax.set_title(f"Sinking Fecal Pellets (Gut Passage Time: {gut_passage_time} hours)")
ax.grid()

# Show plot
plt.show()