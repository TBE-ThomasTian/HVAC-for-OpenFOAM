# UTCIFoam - Universal Thermal Climate Index Calculator

## Overview

UTCIFoam is an OpenFOAM post-processing utility that calculates the Universal Thermal Climate Index (UTCI) from CFD simulation results. It's designed to work seamlessly with buoyantHumiditySimpleFoam and other OpenFOAM solvers.

## Features

### Version 2.0 Improvements

1. **Sky View Factor (SVF) Calculation**
   - Automatic ray-tracing based SVF calculation for outdoor environments
   - Simplified height-based SVF for large meshes
   - Optional SVF field input

2. **Enhanced Mean Radiant Temperature (Tmrt)**
   - Improved calculation considering SVF
   - Support for direct Tmrt input
   - Outdoor environment considerations

3. **UTCI Polynomial Approximation**
   - Implementation of 6th-order polynomial structure
   - Proper input validation for UTCI ranges
   - Wind speed adjustment to 10m height

4. **Physical Accuracy**
   - Validated input ranges for all parameters
   - Proper unit conversions
   - Clear warnings for out-of-range values

## Required Fields

- `T`: Air temperature [K]
- `U`: Velocity field [m/s]
- `thermo:relHum` or `relHum`: Relative humidity [%]
  - buoyantHumiditySimpleFoam writes as `thermo:relHum`
  - Other solvers may use `relHum`
- `qr`: Radiation flux density [W/m²] (or provide Tmrt directly)

## Optional Fields

- `SVF`: Sky View Factor [-] (default: calculated or 0.7)
- `Tmrt`: Mean radiant temperature [K or °C] (if not provided, calculated from qr)

## Usage

After running your CFD simulation (e.g., with buoyantHumiditySimpleFoam):

### Basic usage (requires qr field):
```bash
UTCIFoam
```

### With EPW weather data (recommended for outdoor):
```bash
UTCIFoam -epw weather.epw
```

### With specific date/time:
```bash
UTCIFoam -epw weather.epw -month 7 -day 21 -hour 15
```

### For specific CFD time steps:
```bash
UTCIFoam -time 1000 -epw weather.epw
```

## EPW Weather Data

When using EPW files, UTCIFoam:
- Reads location (latitude, longitude, timezone, elevation)
- Extracts solar radiation data (direct, diffuse)
- Calculates mean radiant temperature using:
  - Direct normal irradiance (DNI)
  - Diffuse horizontal irradiance (DHI)
  - Sky view factor (SVF) from geometry
  - Ground surface temperature
- Uses actual weather conditions for the specified date/time

EPW files can be downloaded from:
- [EnergyPlus Weather Data](https://energyplus.net/weather)
- [Climate.OneBuilding.Org](http://climate.onebuilding.org/)

## Output

The utility creates a `UTCI` field containing the Universal Thermal Climate Index values in °C.

## Valid Ranges

UTCI calculations are valid for:
- Air temperature: -50 to +50°C
- Wind speed (10m): 0.5 to 30 m/s
- Vapor pressure: 0 to 50 hPa
- Tmrt - Ta: -30 to +70°C

Values outside these ranges will be marked as -9999.

## Implementation Notes

### Sky View Factor Calculation

The SVF is calculated using:
- Ray tracing for meshes < 100,000 cells (50 rays per cell)
- Simplified height-based estimation for larger meshes

### Mean Radiant Temperature

Tmrt calculation considers:
- Direct radiation flux (qr)
- Sky View Factor
- Human body emissivity (0.97)
- Stefan-Boltzmann law

### Wind Speed Adjustment

Wind speeds are adjusted from measurement height (assumed 1.5m) to standard 10m height using:
```
v10 = v_h * (10/h)^0.25
```

## Limitations

The current implementation uses example polynomial coefficients. For production use, the full UTCI polynomial coefficients should be obtained from:
- Bröde et al. (2012) Int J Biometeorol 56:481-494
- Official UTCI website (www.utci.org)

## References

1. Bröde, P., Fiala, D., Błażejczyk, K., et al. (2012). "Deriving the operational procedure for the Universal Thermal Climate Index (UTCI)". International Journal of Biometeorology, 56(3), 481-494.

2. Fiala, D., Havenith, G., Bröde, P., et al. (2012). "UTCI-Fiala multi-node model of human heat transfer and temperature regulation". International Journal of Biometeorology, 56(3), 429-441.

3. Lindberg, F., Holmer, B., & Thorsson, S. (2008). "SOLWEIG 1.0 – Modelling spatial variations of 3D radiant fluxes and mean radiant temperature in complex urban settings". International Journal of Biometeorology, 52(7), 697-713.