# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 17:06:53 2025

@author: elican27
"""

# functions for ingestion of microplastics

import numpy as np 
import matplotlib.pyplot as plt
    
def calc_clearance_rate(krill_length_mm):
    
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

def calc_initial_sinking_velocity(mu, rho, rho_s, L, D): #calculating sinking velcoity as in komar et al 1981
    
    #all the units need to be SI
    
    g= 9.81

    initial_sinking_velocity_s = 0.079 * (1/mu)*(rho_s - rho) * g * L**2 *(L/D)**(-1.664)
    #converting the sinking velocity to m/d 
    initial_sinking_velocity = initial_sinking_velocity_s * (60*60*24)
    
    return initial_sinking_velocity

# def calc_initial_sinking_velocity(L, D, rho_s): # Atkinson et al 2012
        
#     #[L] = um
#     #[D] = um
#     #[rho_s] = kg/m^3
    
#     log_sinking_velocity = 1.96 * np.log10(D) + 0.539 * np.log10(L) + 0.00405 * rho_s - 8.33
    
#     initial_sinking_velocity = 10 ** log_sinking_velocity #m/d
    
#     return initial_sinking_velocity

def calc_sinking_velocity(time_since_release, initial_sinking_velocity, b):
    
    #Calculate varying sinking velocity with depth
    sinking_velocity = initial_sinking_velocity * np.exp(-b * time_since_release)
    
    return sinking_velocity




    
                                                   
                                                   
                                                   
        








    
    



# print(f"Clearance Rate: {clearance_rate} ml/hour"):