# Wind-Driven Rain Validation Case: Cubic Building

## Overview

This validation case simulates wind-driven rain on a cubic building, comparing against experimental data from Choi (1993) and computational results from Blocken & Carmeliet (2002).

## Case Description

- **Geometry**: 10m × 10m × 10m cubic building
- **Domain**: 100m × 60m × 60m (10H × 6H × 6H)
- **Wind Speed**: U₁₀ = 10 m/s at 10m height
- **Wind Profile**: Logarithmic with z₀ = 0.03m
- **Rain Intensity**: 10 mm/h horizontal rainfall
- **Raindrop Spectrum**: Based on Best (1950) distribution

## Validation Data

The catch ratio distribution on the windward facade is compared against:
1. Choi (1993) - Analytical model based on particle trajectories
2. Blocken & Carmeliet (2002) - CFD validation study

## Key Parameters

- Reynolds number (building): Re = 6.7×10⁶
- Raindrop diameters: 0.5mm to 3.0mm (5 phases)
- Terminal velocities calculated using Gunn & Kinzer (1949) data

## References

1. Choi, E.C.C. (1993). "Simulation of wind-driven-rain around a building." Journal of Wind Engineering and Industrial Aerodynamics, 46, 721-729.
2. Blocken, B., & Carmeliet, J. (2002). "Spatial and temporal distribution of driving rain on a low-rise building." Wind and Structures, 5(5), 441-462.
3. Best, A.C. (1950). "The size distribution of raindrops." Quarterly Journal of the Royal Meteorological Society, 76(327), 16-36.