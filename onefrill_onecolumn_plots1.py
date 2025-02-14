import numpy as np
import matplotlib.pyplot as plt
from onekrill_onecolumn import (
    calc_clearance_rate,
    calc_krill_mp_consumption,
    calc_mp_fp_production_rate,
    calc_sinking_velocity,
    calc_fp_width_um,
    calc_length_decrease,
    generate_random,
    swdens
)
import scipy.stats as stats
import netCDF4 as nc
import pandas as pd

# Parameters
krill_length_mm = 50  # mm
mp_conc = 500  # particles/m3
depth_limit = 2000  # m
time = np.linspace(0, 200, 500)  # Simulation time in hours
b = -0.3  # Attenuation coefficient
mu = 0.001  # Viscosity of water
rho = 1025  # Density of seawater 
gut_passage_time = 2  # Gut passage time in hours
mass_loss_threshold = 0.4  # Stop sinking when 40% of mass is lost

##importing the temperature and salinity 
temp_data = nc.Dataset('C:/Users/elican27/Documents/Antarctic_krill/Model/Ocean_data/cmems_mod_glo_phy-thetao_anfc_0.083deg_P1M-m_1739443916403.nc')
sal_data = nc.Dataset('C:/Users/elican27/Documents/Antarctic_krill/Model/Ocean_data/cmems_mod_glo_phy-so_anfc_0.083deg_P1M-m_1739443856060.nc')
# Get the temperature variable (adjust the name if different in your file)
temp = temp_data.variables['thetao'][:]  # This will return a numpy array
# Get the depth dimension (adjust the name if necessary)
depth_t = temp_data.variables['depth'][:]  # Depth values as a numpy array
# Get the time dimension (if relevant)
time_t = temp_data.variables['time'][:]  # Time values as a numpy array (months or days, etc.)
# Get the latitude and longitude (if available, adjust names if necessary)
lat_t = temp_data.variables['latitude'][:]  # Latitude values as a numpy array
lon_t = temp_data.variables['longitude'][:]  # Longitude values as a numpy array
# Average over time, latitude, and longitude dimensions (axis 0, 2, and 3)
avg_temp = np.mean(temp, axis=(0, 2, 3))
temp_avg = pd.DataFrame({
    'Depth': depth_t,
    'Average Temperature': avg_temp})
temp_avg = temp_avg.dropna()


# Get the temperature variable (adjust the name if different in your file)
sal = sal_data.variables['so'][:]  # This will return a numpy array
# Get the depth dimension (adjust the name if necessary)
depth_s = sal_data.variables['depth'][:]  # Depth values as a numpy array
# Get the time dimension (if relevant)
time_s = sal_data.variables['time'][:]  # Time values as a numpy array (months or days, etc.)
# Get the latitude and longitude (if available, adjust names if necessary)
lat_s = sal_data.variables['latitude'][:]  # Latitude values as a numpy array
lon_s = sal_data.variables['longitude'][:]  # Longitude values as a numpy array
# Average over time, latitude, and longitude dimensions (axis 0, 2, and 3)
avg_sal = np.mean(sal, axis=(0, 2, 3))
sal_avg = pd.DataFrame({
    'Depth': depth_s,
    'Average Salinity': avg_sal})
sal_avg = sal_avg.dropna()


##merge the two data sets on the depth column 
temp_sal_data = pd.merge(temp_avg, sal_avg, how = 'inner')
temp_sal_data['Density'] = swdens(temp_sal_data['Average Temperature'], temp_sal_data['Average Salinity'])


# Compute ingestion and egestion rates
clearance_rate = calc_clearance_rate(krill_length_mm)
krill_mp_consumption = calc_krill_mp_consumption(clearance_rate, mp_conc)
time_produce_one_mp_fp = calc_mp_fp_production_rate(krill_mp_consumption, gut_passage_time)

# Fecal pellet release times
fp_release_times = np.arange(0, max(time), gut_passage_time)
mp_fp_release_times = np.arange(time_produce_one_mp_fp, max(time), time_produce_one_mp_fp)

# Create figure
fig, ax = plt.subplots(figsize=(12, 10))

# Simulate sinking for each fecal pellet
for release_time in fp_release_times:
    pellet_time = time[time >= release_time]  # Time after release
    time_since_release = pellet_time - release_time  # Time elapsed since egestion

    # **Generate unique length, width, and density for each fecal pellet**
    L_init = generate_random(2927, 2667, 517, 34482) * 10**-6  # Initial FP length (m)
    D = generate_random(183, 178, 80, 600) * 10**-6  # Width/diameter of FP (m)
    rho_s = generate_random(1121, 1116, 1038, 1391)  # Density of krill FP

    # Compute initial sinking velocity for this pellet
    initial_sinking_velocity = calc_sinking_velocity(mu, rho, rho_s, L_init, D)

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
        nearest_depth_index = (temp_sal_data['Depth'] - current_depth).abs().idxmin()
        rho_at_depth = temp_sal_data.loc[nearest_depth_index, 'Density'] *1000

        # Recalculate sinking velocity with updated length
        ws = calc_sinking_velocity(mu, rho_at_depth, rho_s, L, D)
        ws_per_hour = ws / 24  # Convert m/day to m/hour
        
        # Update depth
        current_depth += ws_per_hour * dt
        sinking_depths.append(current_depth)

    # Stop at max depth
    sinking_depths = np.clip(sinking_depths, 0, depth_limit)

    # Determine if this FP contains microplastics
    contains_mp = np.any(np.isclose(release_time, mp_fp_release_times, atol=gut_passage_time / 2))
    color = "red" if contains_mp else "blue"

    ax.plot(pellet_time[:len(sinking_depths)], sinking_depths, color=color)

# Formatting
ax.invert_yaxis()  # Invert y-axis to show depth increasing downward
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Depth (m)")
ax.set_title(f"Sinking Fecal Pellets with Stopping Condition at 40% Mass Loss")
ax.grid()

# Add legend
ax.plot([], [], color='red', label='Contains Microplastics')
ax.plot([], [], color='blue', label='No Microplastics')
ax.legend()

# Show plot
plt.show()


# %% 
import numpy as np
import matplotlib.pyplot as plt
from onekrill_onecolumn import (
    calc_clearance_rate,
    calc_krill_mp_consumption,
    calc_mp_fp_production_rate,
    calc_sinking_velocity,
    calc_fp_width_um,
    calc_length_decrease,
    generate_random
)

# Parameters
krill_length_mm = 50  # mm
depth_limit = 650  # m
time = np.linspace(0, 1000, 2000)  # Simulation time in hours
b = -0.3  # Attenuation coefficient
mu = 0.001  # Viscosity of water
rho = 1025  # Density of water
gut_passage_time = 2  # Gut passage time in hours
mass_loss_threshold = 0.4  # Stop sinking when 40% of mass is lost

def accumulate_microplastics(mp_conc):
    # Compute ingestion and egestion rates
    clearance_rate = calc_clearance_rate(krill_length_mm)
    krill_mp_consumption = calc_krill_mp_consumption(clearance_rate, mp_conc)
    time_produce_one_mp_fp = calc_mp_fp_production_rate(krill_mp_consumption, gut_passage_time)

    # Fecal pellet release times
    fp_release_times = np.arange(0, max(time), gut_passage_time)
    mp_fp_release_times = np.arange(time_produce_one_mp_fp, max(time), time_produce_one_mp_fp)

    # Initialize array to store microplastic accumulation over time
    mp_accumulation = np.zeros_like(time)

    # Simulate sinking for each fecal pellet
    for release_time in fp_release_times:
        pellet_time = time[time >= release_time]  # Time after release
        time_since_release = pellet_time - release_time  # Time elapsed since egestion

        # **Generate unique length, width, and density for each fecal pellet**
        L_init = generate_random(2927, 2667, 517, 34482) * 10 ** (-6)  # Initial FP length (m)
        D = generate_random(183, 178, 80, 600) * 10 ** (-6)  # Width/diameter of FP (m)
        rho_s = generate_random(1121, 1116, 1038, 1391)  # Density of krill FP

        # Compute initial sinking velocity for this pellet
        initial_sinking_velocity = calc_sinking_velocity(mu, rho, rho_s, L_init, D)

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

            # **Stop sinking when 40% of the mass is lost**
            if L <= 0.6 * L_init:
                break  # Stop sinking if 40% of mass is lost

            # Recalculate sinking velocity with updated length
            ws = calc_sinking_velocity(mu, rho, rho_s, L, D)
            ws_per_hour = ws / 24  # Convert m/day to m/hour

            # Update depth
            current_depth += ws_per_hour * dt
            sinking_depths.append(current_depth)

        # Stop at the actual sinking depth (not necessarily 650m)
        final_depth = current_depth  # Depth where sinking stops

        # Determine if this FP contains microplastics
        contains_mp = np.any(np.isclose(release_time, mp_fp_release_times, atol=gut_passage_time / 2))

        # **Accumulate microplastics at the depth where sinking stopped**
        if contains_mp:
            for t in pellet_time[:len(sinking_depths)]:
                if np.isclose(sinking_depths[pellet_time.tolist().index(t)], final_depth):
                    mp_accumulation[time >= t] += mp_conc  # Add microplastics to the accumulation

    return mp_accumulation

# Simulate for two different concentrations
mp_conc_200 = 200  # particles/m³
mp_conc_2000 = 2000  # particles/m³

# Accumulate microplastics for both concentrations
mp_accumulation_200 = accumulate_microplastics(mp_conc_200)
mp_accumulation_2000 = accumulate_microplastics(mp_conc_2000)

# Create figure
fig, ax = plt.subplots(figsize=(12, 10))

# Plot the microplastic accumulation for both concentrations
ax.plot(time, mp_accumulation_200, label="MP Concentration 200 particles/m³", color='blue')
ax.plot(time, mp_accumulation_2000, label="MP Concentration 2000 particles/m³", color='red')

# Formatting
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Microplastics Accumulated (particles/m³)")
ax.set_title("Accumulation of Microplastics at Final Sinking Depth Over Time")
ax.grid()

# Add legend
ax.legend()

# Show plot
plt.show()

#%% This is for breaking at 40% loss against not breaking at all 

import numpy as np
import matplotlib.pyplot as plt
from onekrill_onecolumn import (
    calc_clearance_rate,
    calc_krill_mp_consumption,
    calc_mp_fp_production_rate,
    calc_sinking_velocity,
    calc_fp_width_um,
    calc_length_decrease,
    generate_random
)

# Parameters
krill_length_mm = 50  # mm
depth_limit = 2000  # m
time = np.linspace(0, 1000, 2000)  # Simulation time in hours
b = -0.3  # Attenuation coefficient
mu = 0.001  # Viscosity of water
rho = 1025  # Density of water
gut_passage_time = 2  # Gut passage time in hours
mass_loss_threshold = 0.4  # Stop sinking when 40% of mass is lost
mp_conc = 500  # Microplastic concentration (particles/m³)

# Compute ingestion and egestion rates
clearance_rate = calc_clearance_rate(krill_length_mm)
krill_mp_consumption = calc_krill_mp_consumption(clearance_rate, mp_conc)
time_produce_one_mp_fp = calc_mp_fp_production_rate(krill_mp_consumption, gut_passage_time)

# Fecal pellet release times
fp_release_times = np.arange(0, max(time), gut_passage_time)
mp_fp_release_times = np.arange(time_produce_one_mp_fp, max(time), time_produce_one_mp_fp)

# Initialize array to store microplastic accumulation over time
mp_accumulation_break = np.zeros_like(time)
mp_accumulation = np.zeros_like(time)

# Simulate sinking for each fecal pellet with breaking
for release_time in fp_release_times:
    pellet_time_break = time[time >= release_time]  # Time after release
    time_since_release_break = pellet_time_break - release_time  # Time elapsed since egestion

    # **Generate unique length, width, and density for each fecal pellet**
    L_init_break = generate_random(2927, 2667, 517, 34482) * 10 ** (-6)  # Initial FP length (m)
    D_break = generate_random(183, 178, 80, 600) * 10 ** (-6)  # Width/diameter of FP (m)
    rho_s_break = generate_random(1121, 1116, 1038, 1391)  # Density of krill FP

    # Compute initial sinking velocity for this pellet
    initial_sinking_velocity_break = calc_sinking_velocity(mu, rho, rho_s_break, L_init_break, D_break)

    # Initialize values
    current_depth_break = 100
    L_break = L_init_break
    ws_break = initial_sinking_velocity_break
    dt_break = (time[1] - time[0])  # Time step in hours
    sinking_depths_break = []

    for t in time_since_release_break:
        if current_depth_break >= depth_limit:
            break  # Stop tracking if pellet reaches the max depth

        # Update length based on depth
        L_break = calc_length_decrease(L_init_break, b, current_depth_break)

        # **Stop sinking when 40% of the mass is lost**
        if L_break <= 0.6 * L_init_break:
            break  # Stop sinking if 40% of mass is lost

        # Recalculate sinking velocity with updated length
        ws_break = calc_sinking_velocity(mu, rho, rho_s_break, L_break, D_break)
        ws_per_hour_break = ws_break / 24  # Convert m/day to m/hour

        # Update depth
        current_depth_break += ws_per_hour_break * dt_break
        sinking_depths_break.append(current_depth_break)

    # Stop at the actual sinking depth (not necessarily 650m)
    final_depth_break = current_depth_break  # Depth where sinking stops

    # Determine if this FP contains microplastics
    contains_mp_break = np.any(np.isclose(release_time, mp_fp_release_times, atol=gut_passage_time / 2))

    # **Accumulate microplastics at the depth where sinking stopped**
    if contains_mp_break:
        for t in pellet_time_break[:len(sinking_depths_break)]:
            if np.isclose(sinking_depths_break[pellet_time_break.tolist().index(t)], final_depth_break):
                mp_accumulation_break[time >= t] += mp_conc
  # Add microplastics to the accumulation


# Simulate sinking for each fecal pellet no breaking 
for release_time in fp_release_times:
    pellet_time = time[time >= release_time]  # Time after release
    time_since_release = pellet_time - release_time  # Time elapsed since egestion

    # **Generate unique length, width, and density for each fecal pellet**
    L_init = generate_random(2927, 2667, 517, 34482) * 10 ** (-6)  # Initial FP length (m)
    D = generate_random(183, 178, 80, 600) * 10 ** (-6)  # Width/diameter of FP (m)
    rho_s = generate_random(1121, 1116, 1038, 1391)  # Density of krill FP

    # Compute initial sinking velocity for this pellet
    initial_sinking_velocity = calc_sinking_velocity(mu, rho, rho_s, L_init, D)

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

    # Stop at the actual sinking depth (not necessarily 650m)
    final_depth = current_depth  # Depth where sinking stops

    # Determine if this FP contains microplastics
    contains_mp = np.any(np.isclose(release_time, mp_fp_release_times, atol=gut_passage_time / 2))

    # **Accumulate microplastics at the depth where sinking stopped**
    if contains_mp:
        for t in pellet_time[:len(sinking_depths)]:
            if np.isclose(sinking_depths[pellet_time.tolist().index(t)], final_depth):
                mp_accumulation[time >= t] += mp_conc  # Add microplastics to the accumulation

# Create figure
fig, ax = plt.subplots(figsize=(20, 10))

# Plot the microplastic accumulation for the selected concentration
ax.plot(time, mp_accumulation_break, label="MP Concentration with breakage", color='red')
ax.plot(time, mp_accumulation, label="MP Concentration", color='blue')

# Formatting
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Microplastics Accumulated (particles/m³)")
ax.set_title("Accumulation of Microplastics at Final Sinking Depth Over Time")
ax.grid()

# Add legend
ax.legend()

# Show plot
plt.show()

#%%

import numpy as np
import matplotlib.pyplot as plt
from onekrill_onecolumn import (
    calc_clearance_rate,
    calc_krill_mp_consumption,
    calc_mp_fp_production_rate,
    calc_sinking_velocity,
    calc_fp_width_um,
    calc_length_decrease,
    generate_random,
    swdens
)
import netCDF4 as nc
import pandas as pd


# Parameters
krill_length_mm = 50  # mm
depth_limit = 600  # m
time = np.linspace(0, 10000, 20000)  # Simulation time in hours
b = -0.3  # Attenuation coefficient
mu = 0.001  # Viscosity of water
rho = 1025  # Density of water
gut_passage_time = 2  # Gut passage time in hours
mass_loss_threshold = 0.4  # Stop sinking when 40% of mass is lost
mp_conc = 500  # Microplastic concentration (particles/m³)

##importing the temperature and salinity 
temp_data = nc.Dataset('C:/Users/elican27/Documents/Antarctic_krill/Model/Ocean_data/cmems_mod_glo_phy-thetao_anfc_0.083deg_P1M-m_1739443916403.nc')
sal_data = nc.Dataset('C:/Users/elican27/Documents/Antarctic_krill/Model/Ocean_data/cmems_mod_glo_phy-so_anfc_0.083deg_P1M-m_1739443856060.nc')
# Get the temperature variable (adjust the name if different in your file)
temp = temp_data.variables['thetao'][:]  # This will return a numpy array
# Get the depth dimension (adjust the name if necessary)
depth_t = temp_data.variables['depth'][:]  # Depth values as a numpy array
# Get the time dimension (if relevant)
time_t = temp_data.variables['time'][:]  # Time values as a numpy array (months or days, etc.)
# Get the latitude and longitude (if available, adjust names if necessary)
lat_t = temp_data.variables['latitude'][:]  # Latitude values as a numpy array
lon_t = temp_data.variables['longitude'][:]  # Longitude values as a numpy array
# Average over time, latitude, and longitude dimensions (axis 0, 2, and 3)
avg_temp = np.mean(temp, axis=(0, 2, 3))
temp_avg = pd.DataFrame({
    'Depth': depth_t,
    'Average Temperature': avg_temp})
temp_avg = temp_avg.dropna()


# Get the temperature variable (adjust the name if different in your file)
sal = sal_data.variables['so'][:]  # This will return a numpy array
# Get the depth dimension (adjust the name if necessary)
depth_s = sal_data.variables['depth'][:]  # Depth values as a numpy array
# Get the time dimension (if relevant)
time_s = sal_data.variables['time'][:]  # Time values as a numpy array (months or days, etc.)
# Get the latitude and longitude (if available, adjust names if necessary)
lat_s = sal_data.variables['latitude'][:]  # Latitude values as a numpy array
lon_s = sal_data.variables['longitude'][:]  # Longitude values as a numpy array
# Average over time, latitude, and longitude dimensions (axis 0, 2, and 3)
avg_sal = np.mean(sal, axis=(0, 2, 3))
sal_avg = pd.DataFrame({
    'Depth': depth_s,
    'Average Salinity': avg_sal})
sal_avg = sal_avg.dropna()


##merge the two data sets on the depth column 
temp_sal_data = pd.merge(temp_avg, sal_avg, how = 'inner')
temp_sal_data['Density'] = swdens(temp_sal_data['Average Temperature'], temp_sal_data['Average Salinity'])

# Compute ingestion and egestion rates
clearance_rate = calc_clearance_rate(krill_length_mm)
krill_mp_consumption = calc_krill_mp_consumption(clearance_rate, mp_conc)
time_produce_one_mp_fp = calc_mp_fp_production_rate(krill_mp_consumption, gut_passage_time)

# Fecal pellet release times
fp_release_times = np.arange(0, max(time), gut_passage_time)
mp_fp_release_times = np.arange(time_produce_one_mp_fp, max(time), time_produce_one_mp_fp)

# Initialize arrays to store microplastic accumulation at max depth
mp_accumulation_at_2000m_break = np.zeros_like(time)
mp_accumulation_at_2000m = np.zeros_like(time)
pellets_reaching_2000m_break = 0
pellets_reaching_2000m = 0

# Track when MP reaches 2000m
def find_nearest_index(array, value):
    return (np.abs(array - value)).argmin()

# Simulate sinking for each fecal pellet without breakage
for release_time in fp_release_times:
    pellet_time = time[time >= release_time]
    time_since_release = pellet_time - release_time

    L_init = generate_random(2927, 2667, 517, 34482) * 10 ** (-6)
    D = generate_random(183, 178, 80, 600) * 10 ** (-6)
    rho_s = generate_random(1121, 1116, 1038, 1391)

    initial_sinking_velocity = calc_sinking_velocity(mu, temp_sal_data['Density'].iloc[0], rho_s, L_init, D)
    
    current_depth = 100
    ws = initial_sinking_velocity
    dt = (time[1] - time[0])
    reached_2000m = False
    time_at_2000m = None
    
    for t in time_since_release:
        if current_depth >= depth_limit:
            reached_2000m = True
            time_at_2000m = release_time + t
            pellets_reaching_2000m += 1
            break
        
        # Update length based on depth
        L = calc_length_decrease(L_init, b, current_depth)
        
        #calculate the density at the current depth
        nearest_depth_index = (temp_sal_data['Depth'] - current_depth).abs().idxmin()
        rho_at_depth = temp_sal_data.loc[nearest_depth_index, 'Density'] * 1000
        
        ws = calc_sinking_velocity(mu, rho_at_depth, rho_s, L, D)
        ws_per_hour = ws / 24
        current_depth += ws_per_hour * dt
    
    if reached_2000m and time_at_2000m is not None:
        contains_mp = np.any(np.isclose(release_time, mp_fp_release_times, atol=gut_passage_time / 2))
        if contains_mp:
            index = find_nearest_index(time, time_at_2000m)
            mp_accumulation_at_2000m[index] += mp_conc

mp_accumulation_at_2000m = np.cumsum(mp_accumulation_at_2000m)

# Plot results
fig, ax = plt.subplots(figsize=(20, 10))
ax.plot(time, mp_accumulation_at_2000m, label="MP Concentration at 2000m without breakage", color='blue')
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Microplastics Accumulated at 600m (particles/m³)")
ax.set_title("Accumulation of Microplastics at 600m Over Time")
ax.grid()
ax.legend()
plt.show()

print(f"Total pellets released: {len(fp_release_times)}")
print(f"Pellets reaching 2000m: {pellets_reaching_2000m}")

#%% analysing the FP_length data 

import pandas as pd 
import matplotlib.pyplot as plt 
import scipy.stats as stats
import numpy as np

FP_length = pd.read_csv('FP_length.csv')
print(FP_length.columns)

def truncated_normal_pdf(mean, min_val, max_val, num_points=1000):
    """
    Computes x-values and corresponding PDF values for a truncated normal distribution.

    Parameters:
        mean (float): Mean of the distribution.
        min_val (float): Minimum possible value.
        max_val (float): Maximum possible value.
        num_points (int): Number of points for smooth curve.

    Returns:
        x (numpy array): Range of values from min_val to max_val.
        pdf_y (numpy array): PDF values corresponding to x.
    """
    # Estimate standard deviation
    std_dev = (max_val - min_val) / 6  # Rough estimate assuming normal-like distribution

    # Define bounds in standard normal form
    lower_bound = (min_val - mean) / std_dev
    upper_bound = (max_val - mean) / std_dev

    # Create truncated normal distribution
    distribution = stats.truncnorm(lower_bound, upper_bound, loc=mean, scale=std_dev)

    # Generate x values (range for PDF calculation)
    x = np.linspace(min_val, max_val, num_points)

    # Compute the probability density function (PDF)
    pdf_y = distribution.pdf(x)

    return x, pdf_y  # Return values for plotting

x, pdf_y = truncated_normal_pdf(2.927, 0.517, 4.0, num_points=1000)
pdf_y_scaled = pdf_y * 16 / max(pdf_y)  # Normalize PDF to histogram

mean_length = FP_length['Length '].mean()

plt.hist(FP_length['Length '], bins = 40, label='300m data')
plt.plot(x, pdf_y_scaled, color='red', linewidth=2, label="Atkinson et al 2012 data")
plt.axvline(x=mean_length, color='blue', linestyle='--', linewidth=2, label="Mean at 300m")
plt.axvline(x=2.927, color='red', linestyle='--', linewidth=2, label="Mean at surface")
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))



plt.show()



