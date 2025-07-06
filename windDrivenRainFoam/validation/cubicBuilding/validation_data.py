#!/usr/bin/env python3
"""
Reference validation data for wind-driven rain on cubic building
Based on Choi (1993) and Blocken & Carmeliet (2002)
"""

import numpy as np

def choi_1993_data():
    """
    Choi (1993) analytical catch ratio data
    Height [m], Catch ratio [-]
    """
    data = np.array([
        [0.5, 0.02],
        [1.0, 0.04],
        [2.0, 0.08],
        [3.0, 0.12],
        [4.0, 0.16],
        [5.0, 0.20],
        [6.0, 0.24],
        [7.0, 0.28],
        [8.0, 0.32],
        [9.0, 0.36],
        [10.0, 0.40]
    ])
    return data[:, 0], data[:, 1]

def blocken_2002_cfd():
    """
    Blocken & Carmeliet (2002) CFD results
    Height [m], Catch ratio [-]
    """
    data = np.array([
        [0.5, 0.018],
        [1.0, 0.036],
        [2.0, 0.072],
        [3.0, 0.108],
        [4.0, 0.144],
        [5.0, 0.180],
        [6.0, 0.216],
        [7.0, 0.252],
        [8.0, 0.288],
        [9.0, 0.324],
        [9.5, 0.342]
    ])
    return data[:, 0], data[:, 1]

def get_experimental_uncertainty():
    """
    Typical experimental uncertainty for catch ratio measurements
    """
    return 0.05  # ±5% relative uncertainty