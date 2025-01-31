# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 09:34:12 2025

@author: elican27
"""

import numpy as np
import matplotlib.pyplot as plt
from onekrill_onecolumn import calc_clearance_rate, calc_krill_mp_consumption, calc_mp_fp_production_rate

# Parameters
krill_length_mm = 50  # mm
mp_conc = 200  # particles/m3
egestion_rate = 100  # ugCkrill-1d-1 (one krill)
gut_passage_time = 4  # hours
sinking_velocity = 304  # m/day
depth = 650  # m
time = np.linspace(0, 500, 500)  # Simulation time in hours

# Compute rates
clearance_rate = calc_clearance_rate(krill_length_mm)
krill_mp_consumption = calc_krill_mp_consumption(clearance_rate, mp_conc)
mp_fp_production_rate = calc_mp_fp_production_rate(krill_mp_consumption, gut_passage_time)



## TRACKING MP ACCUMULATION AT 650M ###

# Convert sinking velocity from m/day to m/hour
sinking_velocity_hour = sinking_velocity / 24  

# Time for fecal pellets to reach 650m
time_to_650m_hours = depth / sinking_velocity_hour  

# Convert MP FP production rate to MPs/day
mp_fp_production_rate_day = mp_fp_production_rate * 24  

# Time array in days
time_days = np.linspace(0, 200, 500)  # 200 days, 500 time points

# Initialize MP accumulation array
mp_accumulation = np.zeros_like(time_days)

# Accumulate MPs over time
for i, t in enumerate(time_days):
    if t >= time_to_650m_hours / 24:  # MPs start arriving after travel time
        mp_accumulation[i] = mp_accumulation[i-1] + mp_fp_production_rate_day 

# Plot results

plt.plot(time_days, mp_accumulation)
plt.xlabel("Time (days)")
plt.ylabel("Accumulated MP at 650m")
plt.title("Accumulation of Microplastics at 650m Due to Sinking Fecal Pellets")
plt.show()
