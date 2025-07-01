# ASHRAE55Foam

## Overview

ASHRAE55Foam is an OpenFOAM post-processing utility that calculates thermal comfort metrics according to the ANSI/ASHRAE Standard 55-2020 Adaptive Comfort Model. It evaluates indoor and outdoor thermal comfort conditions in CFD simulations by computing operative temperature and determining compliance with ASHRAE 55 acceptability criteria. The tool works globally by automatically extracting location data from EPW weather files.

## Features

- **Global Location Support**
  - Automatically extracts latitude, longitude, and time zone from EPW files
  - Works for any location worldwide without manual coordinate input
  - Supports command-line override of location parameters

- **Thermal Comfort Assessment**
  - Calculates Operative Temperature (TOp) - the average of air temperature and mean radiant temperature
  - Evaluates compliance with ASHRAE 55 adaptive comfort model at 80% and 90% acceptability levels
  - Accounts for cooling effects from elevated air speeds (>0.6 m/s)

- **Solar Radiation Integration**
  - Optimized for OpenFOAM 2412+ with native solar radiation models
  - Automatically reads radiation field 'G' from OpenFOAM radiation models (solarLoad, P1, fvDOM)
  - Includes fallback EPW-based solar calculations with proper time zone corrections

- **Climate Data Processing**
  - Parses EPW (EnergyPlus Weather) files to calculate running mean outdoor temperature
  - Uses 30-day exponentially weighted running mean per ASHRAE 55 requirements
  - Extracts hourly radiation data for accurate MRT calculations

## Installation

### Prerequisites

- OpenFOAM (v2412 or later, optimized for v2412+)
- Standard OpenFOAM development environment

### Build Instructions

```bash
# Navigate to the ASHRAE55 directory
cd $FOAM_RUN/HVAC-for-OpenFOAM/ASHRAE55

# Clean previous builds (optional)
wclean

# Compile the utility
wmake
```

## Usage

### Basic Usage

```bash
ASHRAE55Foam [OPTIONS]
```

### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-epw <fileName>` | EPW weather file for running mean calculation and location data | - |
| `-dayOfYear <scalar>` | Day of year (1-365) for EPW calculation | 180 |
| `-hour <scalar>` | Hour of day for solar calculations (0-24) | 12.0 |
| `-solarData` | Use detailed solar calculations from EPW (only if no OpenFOAM solar model) | false |
| `-latitude <scalar>` | Site latitude in degrees (overrides EPW value) | EPW value |
| `-longitude <scalar>` | Site longitude in degrees (overrides EPW value) | EPW value |
| `-runningMean <scalar>` | Directly specify running mean outdoor temperature in C | - |

### Examples

1. **Using EPW file with automatic location detection:**
   ```bash
   # Automatically uses location from EPW file
   ASHRAE55Foam -epw weather.epw -dayOfYear 200
   ```

2. **Bangkok example with solar calculations:**
   ```bash
   # Automatically detects Bangkok coordinates from EPW
   ASHRAE55Foam -epw THA_Bangkok.484560_IWEC.epw -dayOfYear 180 -hour 14 -solarData
   ```

3. **Override EPW coordinates if needed:**
   ```bash
   ASHRAE55Foam -epw weather.epw -solarData -latitude 40.7 -longitude -74.0 -hour 15
   ```

4. **Direct specification of running mean temperature:**
   ```bash
   ASHRAE55Foam -runningMean 25.5
   ```

## Required Input Fields

The following fields must be present in your OpenFOAM case:

- **T** (Temperature) - MUST_READ
- **U** (Velocity) - MUST_READ
- **G** (Radiation) - READ_IF_PRESENT (automatically used if available)

## Output Fields

The utility creates three fields in each time directory:

1. **TOp** - Operative temperature [K]
2. **ASHRAELevel80** - Binary field (0/1) indicating 80% acceptability compliance
3. **ASHRAELevel90** - Binary field (0/1) indicating 90% acceptability compliance

## Technical Details

### Adaptive Comfort Model

The ASHRAE 55-2020 adaptive comfort model is implemented as:

- Neutral temperature: `t_cmf = 0.31 × T_rm_out + 17.8` (C)
- 80% acceptability: ±3.5 C from neutral temperature
- 90% acceptability: ±2.5 C from neutral temperature
- Valid for naturally ventilated spaces with running mean outdoor temperature between 10 C and 33.5 C

### Cooling Effect from Air Speed

When operative temperature > 25 C and air speed ≥ 0.6 m/s:
- 0.6-0.9 m/s: 1.2K cooling effect
- 0.9-1.2 m/s: 1.8K cooling effect
- >1.2 m/s: 2.2K cooling effect

### Mean Radiant Temperature Calculation

- **With OpenFOAM radiation model**: Converts irradiance G to MRT using empirical correlations
- **With EPW solar data**: Calculates solar position based on location and time zone
- **Without radiation data**: Uses area-weighted wall temperature average
- Solar effects are limited to realistic ranges (5-20 C temperature rise)
- Time zone corrections ensure accurate solar calculations globally

### Location Data from EPW

The tool automatically extracts from EPW files:
- Latitude and longitude for solar position calculations
- Time zone (GMT offset) for accurate local solar time
- Location name for reference
- Hourly weather and radiation data

## Applications

- HVAC system design and evaluation
- Building thermal comfort assessment
- Urban microclimate studies
- Natural ventilation design
- Compliance checking with ASHRAE 55-2020 standards

## Notes

- The tool is designed to work globally with any standard EPW weather file
- MRT formulas are universally valid but not specifically optimized for tropical regions
- For most accurate results, use OpenFOAM v2412+ with native solar radiation models

## Version Information

- **Version**: 1.3
- **Author**: Thomas Tian
- **License**: GPL-3.0
- **Updates**: Added automatic global location support from EPW files

## References

1. ANSI/ASHRAE Standard 55-2020: Thermal Environmental Conditions for Human Occupancy
2. de Dear, R. J., & Brager, G. S. (1998). Developing an adaptive model of thermal comfort and preference

## Support

For issues, questions, or contributions, please visit the [GitHub repository](https://github.com/your-repo/HVAC-for-OpenFOAM).
