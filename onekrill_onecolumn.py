# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 17:06:53 2025

@author: elican27
"""

# functions for ingestion of microplastics

import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats as stats
    
def calc_clearance_rate(krill_length_mm): #Atkinsons but not 2012 one 
    
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
    
    return time_produce_one_mp_fp #fp/hour

def calc_sinking_velocity(mu, rho, rho_s, L, D): #calculating sinking velcoity as in komar et al 1981
    
    #all the units need to be SI
    
    g= 9.81

    sinking_velocity_s = 0.079 * (1/mu)*(rho_s - rho) * g * L**2 *(L/D)**(-1.664)
    #converting the sinking velocity to m/d 
    sinking_velocity = sinking_velocity_s * (60*60*24)
    
    return sinking_velocity

# def calc_initial_sinking_velocity_at(L, D, rho_s): # Atkinson et al 2012
        
#     #[L] = um
#     #[D] = um
#     #[rho_s] = kg/m^3
    
#     log_sinking_velocity = 1.96 * np.log10(D) + 0.539 * np.log10(L) + 0.00405 * rho_s - 8.33
    
#     initial_sinking_velocity = 10 ** log_sinking_velocity #m/d
    
#     return initial_sinking_velocity

#from atkinson et al 2012 but the R2 is 23??

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

#want to calculate the density based on the type of food that the krill is consuming
#def calc_density(food):
    





    

    





  
                                                   
                                                   
        