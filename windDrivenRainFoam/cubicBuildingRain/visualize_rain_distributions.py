#!/usr/bin/env python3
"""
Visualize rain droplet distributions for different climate zones
"""

import matplotlib.pyplot as plt
import numpy as np

# City configurations with their characteristics
cities_data = {
    'London': {
        'diameters': [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5],
        'fractions': [0.001, 0.011, 0.046, 0.116, 0.213, 0.364, 0.205, 0.045],
        'Rh': 15,
        'color': 'gray',
        'description': 'Maritime drizzle'
    },
    'Dubai': {
        'diameters': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        'fractions': [0.004, 0.030, 0.101, 0.240, 0.372, 0.254],
        'Rh': 50,
        'color': 'orange',
        'description': 'Desert storms'
    },
    'Mumbai': {
        'diameters': [0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0],
        'fractions': [0.000, 0.000, 0.001, 0.002, 0.006, 0.011, 0.020, 0.031, 0.046, 0.066, 0.091, 0.121, 0.157, 0.199, 0.248],
        'Rh': 200,
        'color': 'darkblue',
        'description': 'Extreme monsoon'
    },
    'Vancouver': {
        'diameters': [0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0],
        'fractions': [0.000, 0.002, 0.008, 0.029, 0.077, 0.146, 0.293, 0.323, 0.121],
        'Rh': 25,
        'color': 'green',
        'description': 'Pacific coast'
    },
    'Miami': {
        'diameters': [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5],
        'fractions': [0.001, 0.006, 0.021, 0.050, 0.097, 0.167, 0.200, 0.195, 0.155, 0.085, 0.024],
        'Rh': 100,
        'color': 'red',
        'description': 'Hurricane'
    },
    'Bergen': {
        'diameters': [0.3, 0.5, 0.7, 0.9, 1.2, 1.5, 1.8, 2.2, 2.6, 3.0],
        'fractions': [0.001, 0.004, 0.014, 0.029, 0.070, 0.136, 0.233, 0.266, 0.183, 0.065],
        'Rh': 30,
        'color': 'purple',
        'description': "Europe's rainiest"
    }
}

# Create figure with subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Plot 1: Droplet size distributions
for city, data in cities_data.items():
    ax1.bar(data['diameters'], data['fractions'], 
            width=0.15, alpha=0.7, color=data['color'], 
            label=f"{city} (Rh={data['Rh']} mm/h)")

ax1.set_xlabel('Droplet Diameter (mm)', fontsize=12)
ax1.set_ylabel('Volume Fraction', fontsize=12)
ax1.set_title('Rain Droplet Size Distributions for Different Climate Zones', fontsize=14)
ax1.legend(loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 7.5)

# Plot 2: Rainfall intensity comparison with droplet range
cities = list(cities_data.keys())
intensities = [cities_data[city]['Rh'] for city in cities]
min_diams = [min(cities_data[city]['diameters']) for city in cities]
max_diams = [max(cities_data[city]['diameters']) for city in cities]

x = np.arange(len(cities))
width = 0.35

bars = ax2.bar(x, intensities, width, color=[cities_data[city]['color'] for city in cities], alpha=0.7)

# Add droplet size range as error bars
yerr_lower = np.zeros(len(cities))
yerr_upper = np.zeros(len(cities))
ax2.errorbar(x, intensities, yerr=[yerr_lower, yerr_upper], 
            fmt='none', color='black', capsize=5, capthick=2)

# Add text annotations for droplet ranges
for i, city in enumerate(cities):
    ax2.text(i, intensities[i] + 5, f'{min_diams[i]:.1f}-{max_diams[i]:.1f} mm', 
            ha='center', va='bottom', fontsize=9)

ax2.set_xlabel('City', fontsize=12)
ax2.set_ylabel('Rainfall Intensity (mm/h)', fontsize=12)
ax2.set_title('Rainfall Intensity and Droplet Size Ranges by Climate Zone', fontsize=14)
ax2.set_xticks(x)
ax2.set_xticklabels([f"{city}\n({cities_data[city]['description']})" for city in cities])
ax2.grid(True, axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('climate_zones_rain_distributions.png', dpi=300, bbox_inches='tight')
print("Visualization saved as 'climate_zones_rain_distributions.png'")

# Create a table summary
print("\n" + "="*80)
print("SUMMARY OF RAIN DISTRIBUTIONS")
print("="*80)
print(f"{'City':<12} {'Rh (mm/h)':<12} {'Phases':<8} {'D_min (mm)':<12} {'D_max (mm)':<12} {'Peak (mm)':<12}")
print("-"*80)

for city, data in cities_data.items():
    peak_idx = np.argmax(data['fractions'])
    peak_diam = data['diameters'][peak_idx]
    print(f"{city:<12} {data['Rh']:<12} {len(data['diameters']):<8} "
          f"{min(data['diameters']):<12.1f} {max(data['diameters']):<12.1f} {peak_diam:<12.1f}")

print("\nNote: Peak diameter is the droplet size with highest volume fraction")