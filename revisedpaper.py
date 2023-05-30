import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
import seaborn as sns
from scipy.stats import special_ortho_group as sog

class Generate_SIR_Neff_fixed_dur:
"""
    Generate_SIR_Neff_fixed_dur
    
    Given an infected data time series and duration of infection generates
    a SIR (N*, beta) (For more details consult the paper) with the given
    duration of infection. To start off least square optimisation, initial
    estimates for (N*, beta) are required. 
"""
    def __init__(self, inf_data, param, duration):
        self.I = inf_data
        self.dur = duration
        self.initialval = param
    
    def obj_func(self, params):
        N = params[0]
        beta = params[1]
        gamma = self.dur
        t = len(self.I)
        SIR = solve_ivp(SIR_model, (0, 2*t), y0=[self.I[0], 0, beta, N, gamma], t_eval=range(t))
        return SIR.y[0]-self.I

    def generate_parameters(self):
        N_ = self.initialval[0]
        beta_ = self.initialval[1]
        optim_params = least_squares(self.obj_func, (N_, beta_))
        return optim_params

class Generate_SIR:
"""
    Generate_SIR
    
    Given an infected data time series and pop size N generates
    a SIR (beta, gamma) (For more details consult the paper) with the given
    pop size N. To start off least square optimisation, initial
    estimates for (beta, gamma) are required. 
"""
    def __init__(self, inf_data, param, N):
        self.I = inf_data
        self.initialval = param
        self.N = N
    
    def obj_func(self, params):
        beta = params[0]
        gamma = params[1]
        t = len(self.I)
        SIR = solve_ivp(SIR_model, (0, t), y0=[self.I[0], 0, beta, self.N, gamma], t_eval=range(t))
        return SIR.y[0]-self.I

    def generate_parameters(self):
        beta_ = self.initialval[0]
        gamma_ = self.initialval[1]
        optim_params = least_squares(self.obj_func, (beta_, gamma_))
        return optim_params      
      
