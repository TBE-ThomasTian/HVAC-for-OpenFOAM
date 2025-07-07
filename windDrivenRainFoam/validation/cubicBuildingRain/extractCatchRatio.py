#!/usr/bin/env python3
"""
Extract and process catch ratio data from windDrivenRainFoam simulation
Compare with Choi (1993) analytical results
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from validation_data import choi_1993_data, blocken_2002_cfd

def read_catch_ratio_data(time_dir):
    """Read catch ratio data from OpenFOAM postProcessing"""
    gcr_file = f"postProcessing/catchRatio/{time_dir}/gcr_buildingWindward.raw"
    
    if not os.path.exists(gcr_file):
        print(f"Error: {gcr_file} not found")
        return None, None
    
    # Read raw data
    data = np.loadtxt(gcr_file)
    
    # Extract coordinates and gcr values
    # Raw format: x y z gcr
    coords = data[:, :3]
    gcr = data[:, 3]
    
    # Filter for windward face (x = 0)
    windward = np.abs(coords[:, 0]) < 0.1
    
    return coords[windward], gcr[windward]

def choi_analytical_solution(y, z, U10=10, Rh=10):
    """
    Calculate catch ratio using Choi (1993) analytical model
    
    Parameters:
    - y, z: coordinates on building face [m]
    - U10: wind speed at 10m height [m/s]
    - Rh: horizontal rainfall intensity [mm/h]
    """
    # Wind profile parameters
    z0 = 0.03  # Roughness length [m]
    kappa = 0.41  # von Karman constant
    
    # Calculate friction velocity
    ustar = kappa * U10 / np.log(10 / z0)
    
    # Wind velocity at height z
    U_z = np.where(z > z0, ustar / kappa * np.log(z / z0), 0)
    
    # Raindrop spectrum (Best 1950)
    # Simplified - using mean diameter
    d_mean = 1.5e-3  # 1.5 mm
    
    # Terminal velocity (Gunn & Kinzer 1949)
    # Approximation: vt = 9.65 - 10.3*exp(-600*d)
    vt = 9.65 - 10.3 * np.exp(-600 * d_mean)
    
    # Catch ratio (simplified)
    # eta = (U_z / vt) * (Rh_v / Rh)
    # For vertical face: Rh_v = Rh * U_z / vt
    eta = (U_z / vt)
    
    return eta

def plot_validation():
    """Create validation plots comparing simulation with analytical results"""
    
    # Get latest time directory
    times = []
    if os.path.exists("postProcessing/catchRatio"):
        times = sorted([d for d in os.listdir("postProcessing/catchRatio") 
                       if d.replace('.', '').isdigit()])
    
    if not times:
        print("No catch ratio data found. Run simulation first.")
        return
    
    latest_time = times[-1]
    print(f"Processing data from time {latest_time}")
    
    # Read simulation data
    coords, gcr_sim = read_catch_ratio_data(latest_time)
    
    if coords is None:
        return
    
    # Get reference data
    z_choi, gcr_choi = choi_1993_data()
    z_blocken, gcr_blocken = blocken_2002_cfd()
    
    # Create analytical solution
    z_analytical = np.linspace(0.1, 10, 100)
    y_center = 5.0  # Center of building face
    gcr_analytical = choi_analytical_solution(y_center, z_analytical)
    
    # Plot vertical profile at building center
    plt.figure(figsize=(10, 8))
    
    # Filter simulation data for center line (y ≈ 5m)
    center_mask = np.abs(coords[:, 1] - y_center) < 0.5
    z_sim = coords[center_mask, 2]
    gcr_center = gcr_sim[center_mask]
    
    # Sort by height
    sort_idx = np.argsort(z_sim)
    z_sim = z_sim[sort_idx]
    gcr_center = gcr_center[sort_idx]
    
    plt.subplot(2, 1, 1)
    plt.plot(gcr_choi, z_choi, 'ks', markersize=8, label='Choi (1993) - Reference')
    plt.plot(gcr_blocken, z_blocken, 'b^', markersize=8, label='Blocken & Carmeliet (2002) - CFD')
    plt.plot(gcr_analytical, z_analytical, 'k--', linewidth=1, alpha=0.5, label='Analytical approximation')
    plt.plot(gcr_center, z_sim, 'ro-', markersize=8, linewidth=2, label='windDrivenRainFoam - Present study')
    plt.xlabel('Catch Ratio η [-]')
    plt.ylabel('Height z [m]')
    plt.title('Catch Ratio Distribution on Windward Face (y = 5m)')
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best')
    plt.xlim(0, max(max(gcr_choi), max(gcr_center)) * 1.1)
    
    # Create 2D contour plot
    plt.subplot(2, 1, 2)
    
    # Interpolate to regular grid
    from scipy.interpolate import griddata
    
    yi = np.linspace(0, 10, 50)
    zi = np.linspace(0, 10, 50)
    Yi, Zi = np.meshgrid(yi, zi)
    
    # Interpolate gcr values
    points = coords[:, 1:3]  # y, z coordinates
    gcr_grid = griddata(points, gcr_sim, (Yi, Zi), method='linear')
    
    contour = plt.contourf(Yi, Zi, gcr_grid, levels=20, cmap='viridis')
    plt.colorbar(contour, label='Catch Ratio η [-]')
    plt.xlabel('Width y [m]')
    plt.ylabel('Height z [m]')
    plt.title('Catch Ratio Distribution on Windward Face')
    
    plt.tight_layout()
    plt.savefig('validation_catchRatio.png', dpi=150)
    print("Validation plot saved as validation_catchRatio.png")
    
    # Calculate error metrics
    if len(z_sim) > 0 and len(gcr_center) > 0:
        # Interpolate reference data to simulation points
        gcr_choi_interp = np.interp(z_sim, z_choi, gcr_choi)
        gcr_blocken_interp = np.interp(z_sim, z_blocken, gcr_blocken)
        
        # Calculate errors against Choi reference
        abs_error_choi = np.abs(gcr_center - gcr_choi_interp)
        rel_error_choi = abs_error_choi / (gcr_choi_interp + 1e-10) * 100
        
        # Calculate errors against Blocken CFD
        abs_error_blocken = np.abs(gcr_center - gcr_blocken_interp)
        rel_error_blocken = abs_error_blocken / (gcr_blocken_interp + 1e-10) * 100
        
        print("\nValidation Results:")
        print("\nComparison with Choi (1993):")
        print(f"Mean absolute error: {np.mean(abs_error_choi):.4f}")
        print(f"Max absolute error: {np.max(abs_error_choi):.4f}")
        print(f"Mean relative error: {np.mean(rel_error_choi):.2f}%")
        print(f"Max relative error: {np.max(rel_error_choi):.2f}%")
        
        print("\nComparison with Blocken & Carmeliet (2002):")
        print(f"Mean absolute error: {np.mean(abs_error_blocken):.4f}")
        print(f"Max absolute error: {np.max(abs_error_blocken):.4f}")
        print(f"Mean relative error: {np.mean(rel_error_blocken):.2f}%")
        print(f"Max relative error: {np.max(rel_error_blocken):.2f}%")
        
        # Write detailed results
        with open('validation_results.txt', 'w') as f:
            f.write("Wind-Driven Rain Validation Results\n")
            f.write("===================================\n\n")
            f.write(f"Time: {latest_time} s\n")
            f.write(f"Wind speed: 10 m/s at 10m height\n")
            f.write(f"Rain intensity: 10 mm/h\n\n")
            f.write("Error Analysis vs Choi (1993):\n")
            f.write(f"Mean absolute error: {np.mean(abs_error_choi):.4f}\n")
            f.write(f"Max absolute error: {np.max(abs_error_choi):.4f}\n")
            f.write(f"Mean relative error: {np.mean(rel_error_choi):.2f}%\n")
            f.write(f"Max relative error: {np.max(rel_error_choi):.2f}%\n\n")
            
            f.write("Error Analysis vs Blocken & Carmeliet (2002):\n")
            f.write(f"Mean absolute error: {np.mean(abs_error_blocken):.4f}\n")
            f.write(f"Max absolute error: {np.max(abs_error_blocken):.4f}\n")
            f.write(f"Mean relative error: {np.mean(rel_error_blocken):.2f}%\n")
            f.write(f"Max relative error: {np.max(rel_error_blocken):.2f}%\n\n")
            
            f.write("Height [m]\tSimulation\tChoi\t\tBlocken\t\tError vs Choi [%]\n")
            f.write("-" * 70 + "\n")
            for i in range(0, len(z_sim), max(1, len(z_sim)//20)):
                f.write(f"{z_sim[i]:.2f}\t\t{gcr_center[i]:.4f}\t\t"
                       f"{gcr_choi_interp[i]:.4f}\t\t{gcr_blocken_interp[i]:.4f}\t\t{rel_error_choi[i]:.2f}\n")
        
        print("\nDetailed results saved to validation_results.txt")

if __name__ == "__main__":
    plot_validation()