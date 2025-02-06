import numpy as np
import matplotlib.pyplot as plt
from onekrill_onecolumn import (
    calc_clearance_rate,
    calc_krill_mp_consumption,
    calc_mp_fp_production_rate,
    calc_sinking_velocity,
    calc_fp_width_um,
    calc_length_decrease
)
import scipy.stats as stats

# Parameters
krill_length_mm = 50  # mm
mp_conc = 500  # particles/m3
depth_limit = 650  # m
time = np.linspace(0, 200, 500)  # Simulation time in hours
b = -0.3  # Attenuation coefficient
mu = 0.001  # Viscosity of water
rho = 1025  # Density of water
rho_s = 1121  # Density of krill FP (Atkinson et al 2012)
L_init = 2928 * 10 ** (-6)  # Initial FP length (m)
D = calc_fp_width_um(krill_length_mm) * 10 ** (-6)  # Width/diameter of FP (m)
gut_passage_time = 2  # Gut passage time in hours

# Compute ingestion and egestion rates
clearance_rate = calc_clearance_rate(krill_length_mm)
krill_mp_consumption = calc_krill_mp_consumption(clearance_rate, mp_conc)
time_produce_one_mp_fp = calc_mp_fp_production_rate(krill_mp_consumption, gut_passage_time)
initial_sinking_velocity = calc_sinking_velocity(mu, rho, rho_s, L_init, D)

# Fecal pellet release times
fp_release_times = np.arange(0, max(time), gut_passage_time)
mp_fp_release_times = np.arange(time_produce_one_mp_fp, max(time), time_produce_one_mp_fp)

# Create figure
fig, ax = plt.subplots(figsize=(12, 10))

# Simulate sinking for each fecal pellet
for release_time in fp_release_times:
    pellet_time = time[time >= release_time]  # Time after release
    time_since_release = pellet_time - release_time  # Time elapsed since egestion

    # Initialize values
    current_depth = 100
    L = L_init
    ws = initial_sinking_velocity
    dt = (time[1] - time[0])  # Time step in hours
    sinking_depths = []

    for t in time_since_release:
        if current_depth >= depth_limit:
            break  # Stop tracking if pellet reaches the max depth

        # Update length based on depth
        L = calc_length_decrease(L_init, b, current_depth)
        
        # Recalculate sinking velocity with updated length
        ws = calc_sinking_velocity(mu, rho, rho_s, L, D)
        ws_per_hour = ws / 24  # Convert m/day to m/hour
        
        # Update depth
        current_depth += ws_per_hour * dt
        sinking_depths.append(current_depth)

    # Stop at max depth
    sinking_depths = np.clip(sinking_depths, 0, depth_limit)

    # Determine if this FP is an MP FP
    contains_mp = np.any(np.isclose(release_time, mp_fp_release_times, atol=gut_passage_time / 2))
    color = "red" if contains_mp else "blue"

    ax.plot(pellet_time[:len(sinking_depths)], sinking_depths, color=color)

# Formatting
ax.invert_yaxis()  # Invert y-axis to show depth increasing downward
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Depth (m)")
ax.set_title(f"Sinking Fecal Pellets with Dynamic Length & Velocity")
ax.grid()

# Show plot
plt.show()
