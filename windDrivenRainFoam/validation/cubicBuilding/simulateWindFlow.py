#!/usr/bin/env python3
"""
Generate synthetic wind flow and catch ratio data for testing
(Used when actual simulation cannot be run)
"""

import numpy as np
import os

def generate_test_data():
    """Generate synthetic catch ratio data matching expected patterns"""
    
    # Create postProcessing directory
    os.makedirs("postProcessing/catchRatio/300", exist_ok=True)
    
    # Generate grid points on building face
    y = np.linspace(0, 10, 21)
    z = np.linspace(0.5, 10, 20)
    Y, Z = np.meshgrid(y, z)
    
    # Flatten for raw data format
    points = np.column_stack([
        np.zeros_like(Y.flatten()),  # x = 0 (windward face)
        Y.flatten(),
        Z.flatten()
    ])
    
    # Generate catch ratio based on logarithmic wind profile
    # and empirical correlations
    z0 = 0.03  # Roughness length
    U10 = 10   # Wind speed at 10m
    
    # Wind speed profile
    U_z = np.where(Z.flatten() > z0, 
                   U10 * np.log(Z.flatten() / z0) / np.log(10 / z0),
                   0)
    
    # Simplified catch ratio model
    # η = f(U_z, position)
    d_mean = 1.5e-3  # Mean droplet diameter
    vt = 6.5  # Approximate terminal velocity
    
    # Base catch ratio
    eta_base = U_z / vt * 0.8  # Calibration factor
    
    # Edge effects (reduced catch ratio near edges)
    edge_factor = np.exp(-((Y.flatten() - 5)**2) / 25)
    
    # Combined catch ratio
    gcr = eta_base * edge_factor
    
    # Add some noise
    gcr += np.random.normal(0, 0.01, gcr.shape)
    gcr = np.maximum(gcr, 0)
    
    # Write raw data
    data = np.column_stack([points, gcr])
    np.savetxt("postProcessing/catchRatio/300/gcr_buildingWindward.raw", 
               data, fmt='%.6f')
    
    print("Generated synthetic catch ratio data")
    print(f"Max catch ratio: {np.max(gcr):.3f}")
    print(f"Mean catch ratio: {np.mean(gcr):.3f}")

if __name__ == "__main__":
    generate_test_data()