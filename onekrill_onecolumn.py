# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 17:06:53 2025

@author: elican27
"""

# functions for ingestion of microplastics

import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats as stats
import pandas as pd

def calc_clearance_rate(krill_length_mm): #Atkinsons et al 2002 I think  
    
    a = 0.00036
    b = 3.277
    
    dry_mass_mg = a * (krill_length_mm**b)
    
    clearance_rate = 0.5 * dry_mass_mg ###units?? - ml/hour I think 
    
    return clearance_rate

def calc_krill_mp_consumption(clearance_rate, mp_conc):
    
    #convert microplastic conc to particles/ml 
    mp_conc_ml= mp_conc / 1000000
    
    krill_mp_consumption = clearance_rate * mp_conc_ml #particles/hour 
    
    return krill_mp_consumption #particles/hour



def calc_mp_fp_production_rate(krill_mp_consumption, gut_passage_time):
    
    one_mp_consumption_time = 1/krill_mp_consumption #hours
    
    time_produce_one_mp_fp = one_mp_consumption_time + gut_passage_time #hour
    
    #mp_fp_production_rate = 1 / time_produce_one_mp_fp 
    
    return time_produce_one_mp_fp #hour

def calc_sinking_velocity(mu, rho, rho_s, L, D): #calculating sinking velcoity as in komar et al 1981
    
    #all the units need to be SI
    
    g= 9.81

    sinking_velocity_s = 0.079 * (1/mu)*(rho_s - rho) * g * L**2 *(L/D)**(-1.664)
    
    #converting the sinking velocity to m/d 
    sinking_velocity = sinking_velocity_s * (60*60*24)
    
    return sinking_velocity

def calc_rising_velocity(D, rho_p, rho, mu ):
    
    g = 9.81
    
    
    w_s = (( D**2 * (rho - rho_p) * g ) / (18 * mu )) * 0.1
    
    #convert to m/day
    w = w_s * (60*60*24)
    
    return w


def calc_fp_width_um(krill_length_mm): 
    
    a = 0.00036
    b = 3.277
    
    dry_mass_mg = a * (krill_length_mm**b)
    dry_mass_ug = dry_mass_mg/1000
    
    fp_width_um = 0.165 * dry_mass_ug +153
    
    return fp_width_um
    

##This is not based on anythong at all##
def calc_length_decrease( L_init, b, z):
    
    #calculate the 
    L = L_init * (z/100) ** b   
    
    return L


##DISTRIBUTION FROM THE DATA IN THE ATKINSON ET AL 2012## 

def generate_random(mean, median, min_val, max_val, size=1):
    # Estimate standard deviation (assume mean-median difference as skew indicator)
    std_dev = (max_val - min_val) / 6  # Rough estimate assuming normal-like distribution

    # Define bounds in standard normal form
    lower_bound = (min_val - mean) / std_dev
    upper_bound = (max_val - mean) / std_dev

    # Create truncated normal distribution
    distribution = stats.truncnorm(lower_bound, upper_bound, loc=mean, scale=std_dev)

    # Sample from the distribution
    random_values = distribution.rvs(size=size)
    
    return random_values if size > 1 else random_values[0]

def swdens(TempC, Sal):
    """
    Seawater Density (kg / L) from Temp (C) and Sal (PSU)
    Chapter 5, Section 4.2 of Dickson, Sabine and Christian
    (2007, http://cdiac.ornl.gov/oceans/Handbook_2007.html)
    Parameters
    ----------
    TempC : array-like
    Temperature in celcius.
    Sal : array-like
    Salinity in PSU
    Returns
    -------
    Density in kg / L
    """
    # convert temperature to IPTS-68
    T68 = (TempC + 0.0002) / 0.99975
    pSMOW = (
    999.842594
    + 6.793952e-2 * T68
    + -9.095290e-3 * T68 ** 2
    + 1.001685e-4 * T68 ** 3
    + -1.120083e-6 * T68 ** 4
    + 6.536332e-9 * T68 ** 5
    )
    A = (
    8.24493e-1
    + -4.0899e-3 * T68
    + 7.6438e-5 * T68 ** 2
    + -8.2467e-7 * T68 ** 3
    + 5.3875e-9 * T68 ** 4
    )
    B = -5.72466e-3 + 1.0227e-4 * T68 + -1.6546e-6 * T68 ** 2
    C = 4.8314e-4
    return (pSMOW + A * Sal + B * Sal ** 1.5 + C * Sal ** 2) / 1000


def assign_mp_size(size=1):
    """
    Assigns microplastic sizes based on predefined Feret diameters and their corresponding normalized frequencies.
    
    Parameters:
    size (int): The number of microplastic sizes to generate. Default is 1.
    
    Returns:
    numpy.ndarray: An array of randomly selected microplastic sizes in meters (m).
    """
    # Feret diameters in micrometers
    feret_diameters = [
        40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 250, 260, 270,
        280, 290, 300, 320, 340, 360, 380, 400, 450, 500
    ]
    
    # Normalized frequencies corresponding to each Feret diameter
    normalized_frequencies = [
        0.00243309, 0.00851582, 0.02433090, 0.04257908, 0.06690998,
        0.09732360, 0.09489051, 0.09124088, 0.08515815, 0.07907543,
        0.07055961, 0.03649635, 0.03892944, 0.04257908, 0.04501217,
        0.04014599, 0.03406326, 0.03041363, 0.02433090, 0.01824818,
        0.01216545, 0.00729927, 0.00486618, 0.00243309
    ]
    
    # Normalize frequencies to sum to 1
    normalized_frequencies = np.array(normalized_frequencies)
    normalized_frequencies /= normalized_frequencies.sum()
    
    # Randomly select Feret diameters based on the distribution
    mp_size_um = np.random.choice(feret_diameters, size=size, p=normalized_frequencies)
    mp_size = mp_size_um * 10**(-6)  # Convert micrometers to meters
    
    return mp_size

def assign_krill_length(size=1):

    df = pd.read_csv('krill_length_clean.csv')

    krill_length = df['length'].values

    frequency = df['frequency'].values

    noralized_frequency = frequency / frequency.sum()

    krill_length_mm = np.random.choice(krill_length, size=size, p=noralized_frequency)

    return krill_length_mm if size > 1 else krill_length_mm.item()
    

def calculate_fp_density(date):
    """
    Calculate the density of krill faecal pellets based on the date provided, which determines the season.
    
    Parameters:
    date (datetime.date): The date for which the pellet density is to be calculated.

    Returns:
        float: Calculated density of the faecal pellets (kg/m^3).
    """
    # Seasonal diatom biomass ratios in decimal form
    seasonal_diatom_ratios = {
        'Spring': 0.90,  # 90% diatoms
        'Summer': 0.16,  # 16% diatoms
        'Autumn': 0.03   # 3% diatoms
    }
    
    # Determine season based on date
    month = date.month
    if month in [9, 10, 11]:  # September, October, November
        season = 'Spring'
    elif month in [12, 1, 2]:  # December, January, February
        season = 'Summer'
    elif month in [3, 4, 5]:  # March, April, May
        season = 'Autumn'
    else:
        return 1050  # Default density if not in specified seasons
    
    density_diatom = 1080  # kg/m^3
    density_non_diatom = 1050  # kg/m^3
    diatom_ratio = seasonal_diatom_ratios[season]
    
    # Calculate the density based on the ratio
    fp_density = (diatom_ratio * density_diatom) + ((1 - diatom_ratio) * density_non_diatom)
    return fp_density







    

    





  
                                                   
                                                   
        