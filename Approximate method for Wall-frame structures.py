# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:36:36 2024

@author: ljp70
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 09:22:26 2023

@author: ljp70
"""

from openseespy.opensees import *
# import openseespy.postprocessing.Get_Rendering as opsplt
import numpy as np
import math
import matplotlib.pyplot as plt
from math import asin, sqrt, pi
import opsvis as opsv
import vfo.vfo as vfo
import pandas as pd

wipe()

model('basic', '-ndm', 2, '-ndf', 3)
###################################################################################################
"Define nodes and cordinates"
###################################################################################################

# Load data from Excel file
file_path = r"Saved Location\Wall-frame dataset from GitHUB.xlsx"
df = pd.read_excel(file_path)

# Define section properties and elastic beam column elements
E = 25 * 10**6  # kN / m^2


build_ID = 1
for z in range(build_ID, build_ID + 1):   
    wipe()
    
    # Extract parameters from the DataFrame for each row

    NStories = int(df.iloc[z, 0])  
    Wall_stories = int(NStories)
    HStoryTyp = float(df.iloc[z, 1])
    NBays = int(df.iloc[z, 2])
    WBay = float(df.iloc[z, 7])/1000            # unit m
    dc = float(df.iloc[z, 4])/1000              # unit m
    wb = float(df.iloc[z, 5])/1000              # unit m
    db = float(df.iloc[z, 6])/1000              # unit m
    Lw = float(df.iloc[z, 8])/1000              # unit m
    tw = float(df.iloc[z, 9])/1000              # unit m
    
    H = NStories * HStoryTyp  # height of building
    HStoryTyp = HStoryTyp
    HBuilding = HStoryTyp + (NStories - 1) * HStoryTyp  # height of building
    
    Ic = dc**4 / 12 /3  # m^4
    Ac = dc**2 # m^2
    
    Ib = wb * db**3 / 12/3   # m^4
    Ab = wb*db  # m^2
    
    Iw = Lw**3 * tw / 12 /3  # m^4
    Aw = Lw*tw # m^2
    
    frame_area = WBay**2 * (NBays)
    SW = 8  # kPa
    mass_frame = frame_area * SW / 9.81                         # in Tonnes
    
        
    """Estimating storey stiffness of the wall"""
    
    # Arbitrary linear force distribution in kN
    Fwi = np.arange(1, Wall_stories + 1) * 100 
    
    # Calculating shear demand at each floor 
    Vwi = np.array([np.sum(Fwi[i:]) for i in range(Wall_stories)])
    
    # Calculating moment demand at each floor 
    Mwi = np.array([np.sum(Vwi[i:]*HStoryTyp) for i in range(Wall_stories)])
    
    # Calculating curvature
    phi_i = Mwi / (E * Iw)
    
    # Calculating change in slope
    delta_theta = phi_i * HStoryTyp
    delta_theta[0] = phi_i[0] * HStoryTyp /2
    
    # Calculating slope at each floor
    theta = np.cumsum(delta_theta)
    
    # Calculating storey distortion
    delta_w = theta * HStoryTyp
    
    "Storey stiffness of the wall "
    Kw = np.zeros(NStories)
    for i in range(Wall_stories):
        Kw[i] = Vwi[i] / delta_w[i]
           
    "Calculating storey stiffness of the frame"
    kc = Ic / HStoryTyp     # Relative stiffness of the column
    kb = Ib / WBay          # Relative stiffness of the beam
    
    nc = NBays   # No. of column rows in the structure
    nb = NBays   # Number of beams on each floor
    
    # First floor storey stiffness
    K1 = 24 * E / HStoryTyp**2 * 1 / ((2 / (nc * kc)) + (1 / (nb * kb))) 
    # Typical floor storey stiffness
    Ktyp = 24 * E / HStoryTyp**2 * 1 / ((2 / (nc * kc)) + (1 / (nb * kb)) + (1 / (nb * kb))) 
    
    # Creatign storey stiffness array
    kf = np.array([K1] + [Ktyp] * (NStories - 1))
    
    "Combining storey stiffness of the wall and the frame"
    K_tot = kf + Kw
    
    
    "Applying Rayleigh's Principle to estimate the period and mode-shape"
    # Calculating shear demand and deflection
    Fi = np.arange(1, NStories + 1) * 100  # Arbitrary linear force distribution in kN
    Vi = np.array([np.sum(Fi[i:]) for i in range(NStories)])
    delta = Vi / K_tot
    disp = np.cumsum(delta)
    
    mode_shape = disp/(disp[-1])                            # normalised the mode shape (sum equal to 1)
    mode_shape = np.round(mode_shape, 3)
    mode_shape = np.insert(mode_shape, 0, 0)
    
    # Potential energy and kinetic energy
    PE = np.sum(0.5 * Fi * disp)
    KE = np.sum(0.5 * mass_frame * disp**2)
    
    "Fundamental period of the structure "
    omega = np.sqrt(PE / KE)
    est_T1 = 2 * np.pi / omega
    print("Period estimate = " + str(np.round(est_T1, 2)))   
    

    

