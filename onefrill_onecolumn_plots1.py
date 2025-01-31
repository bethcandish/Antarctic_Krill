import numpy as np
import matplotlib.pyplot as plt
from onekrill_onecolumn import calc_clearance_rate, calc_krill_mp_consumption, calc_mp_fp_production_rate, calc_sinking_velocity


### Function for plotting different gut passage times 
def simulate_sinking(gut_passage_time, ax):
    # Parameters
    krill_length_mm = 50  # mm
    mp_conc = 200  # particles/m3
    sinking_velocity_initial = 304  # m/day
    depth = 650  # m
    time = np.linspace(0, 500, 500)  # Simulation time in hours
    k = 0.01 # attenuatio$n of the velocity with depth

    # Compute rates
    clearance_rate = calc_clearance_rate(krill_length_mm)
    krill_mp_consumption = calc_krill_mp_consumption(clearance_rate, mp_conc)
    time_produce_one_mp_fp = calc_mp_fp_production_rate(krill_mp_consumption, gut_passage_time)
    


    # Generate fecal pellet release times
    fp_release_times = np.arange(0, max(time), gut_passage_time)

    # Generate MP FP release times (fixed intervals)
    mp_fp_release_times = np.arange(time_produce_one_mp_fp, max(time), time_produce_one_mp_fp)

    # Simulate sinking for each fecal pellet
    for release_time in fp_release_times:
        pellet_time = time[time >= release_time]  # Time after release
        time_since_release = pellet_time - release_time  # Time elapsed since egestion

        # Apply exponential decay to sinking velocity
        sinking_velocity = calc_sinking_velocity(time_since_release, sinking_velocity_initial, k)
        sinking_velocity_per_hour = sinking_velocity / 24
        
        # Compute depth with decayed sinking velocity
        pellet_depth = np.cumsum(sinking_velocity_per_hour) * (time[1] - time[0])  # Numerical integration

        # Stop tracking pellets that exceed 650m
        pellet_depth = np.clip(pellet_depth, 0, depth)

        # Determine if this FP is an MP FP (if it's in the mp_fp_release_times)
        contains_mp = np.any(np.isclose(release_time, mp_fp_release_times, atol=gut_passage_time/2))  
        color = "red" if contains_mp else "blue"  # Red for MP, Blue for normal

        ax.plot(pellet_time, pellet_depth, color=color)

    # Formatting
    ax.invert_yaxis()  # Invert y-axis to show depth increasing downward
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Depth (m)")
    ax.set_title(f"Sinking Fecal Pellets (Gut Passage Time: {gut_passage_time} hours)")
    ax.grid()

# Create figure and subplots
fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

# Simulate for different gut passage times
simulate_sinking(4, axes[0])
simulate_sinking(8, axes[1])

# Show plot
plt.tight_layout()
plt.show()
