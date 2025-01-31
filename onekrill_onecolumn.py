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

def calc_sinking_velocity(time_since_release, sinking_velocity_initial, k):
    
    #Calculate varying sinking velocity with depth
    sinking_velocity = sinking_velocity_initial * np.exp(-k * time_since_release)
    
    return sinking_velocity




    
                                                   
                                                   
                                                   
        








    
    



# print(f"Clearance Rate: {clearance_rate} ml/hour"):