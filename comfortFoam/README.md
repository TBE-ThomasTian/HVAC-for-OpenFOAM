# comfortFoam

## Overview

comfortFoam is an OpenFOAM post-processing utility that calculates thermal comfort indices according to ISO 7730:2005. It computes Predicted Mean Vote (PMV), Predicted Percentage of Dissatisfied (PPD), and Draft Rating (DR) from CFD simulation results.

## Features

- **ISO 7730 Compliance**: Implements the standard PMV/PPD model by P.O. Fanger
- **Draft Rating Calculation**: Assesses local thermal discomfort due to air drafts
- **Regional Analysis**: Support for cellSet/cellZone-based comfort evaluation
- **Volume-Weighted Averaging**: Accurate spatial averaging for room-level metrics
- **Flexible Humidity Input**: Compatible with various OpenFOAM solver outputs
- **Multiple Turbulence Models**: Supports k-epsilon, k-omega SST, and other models
- **Advanced Radiation Handling**: Compatible with various radiation models (P1, fvDOM, etc.)
- **Validation Mode**: Built-in ISO 7730 validation capability

## Thermal Comfort Indices

### PMV (Predicted Mean Vote)
Seven-point thermal sensation scale:
- +3: Hot
- +2: Warm
- +1: Slightly warm
- 0: Neutral
- -1: Slightly cool
- -2: Cool
- -3: Cold

### PPD (Predicted Percentage of Dissatisfied)
Percentage of people predicted to be dissatisfied with thermal conditions:
- PPD = 100 - 95 × exp(-0.03353×PMV⁴ - 0.2179×PMV²)
- Minimum PPD ≈ 5% at PMV = 0

### DR (Draft Rating)
Percentage of people dissatisfied due to draft:
- DR = (34 - Ta) × (v - 0.05)^0.62 × (0.37 × v × Tu + 3.14)
- Where Ta is air temperature, v is velocity, Tu is turbulence intensity

## Installation

```bash
cd $FOAM_RUN/../HVAC-for-OpenFOAM/comfortFoam
wmake
```

## Usage

### Basic Usage

```bash
comfortFoam
```

### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-region <name>` | Specify mesh region for multi-region cases | - |
| `-cellSet <name>` | Calculate comfort only for specified cellSet | - |
| `-cellZone <name>` | Calculate comfort only for specified cellZone | - |
| `-noWrite` | Calculate without writing fields | false |
| `-validate` | Run ISO 7730 validation mode | false |
| `-validateAirTemp <T>` | Air temperature for validation [°C] | - |
| `-validateRadTemp <T>` | Radiant temperature for validation [°C] | - |
| `-validateVelocity <v>` | Air velocity for validation [m/s] | - |
| `-validateRH <RH>` | Relative humidity for validation [%] | - |
| `-validateMet <met>` | Metabolic rate for validation [met] | - |
| `-validateClo <clo>` | Clothing insulation for validation [clo] | - |

### Configuration File

Create a `comfortFoamDict` file in the `system` directory:

```c
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      comfortFoamDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Clothing insulation [clo]
// 1 clo = 0.155 m²·K/W
// Typical values:
// - Summer clothing: 0.5 clo
// - Light working clothes: 0.7 clo
// - Typical indoor clothing: 1.0 clo
// - Heavy suit: 1.5 clo
clo         1.0;

// Metabolic rate [met]
// 1 met = 58.15 W/m²
// Typical values:
// - Seated, quiet: 1.0 met
// - Seated, light work: 1.2 met
// - Standing, light work: 1.6 met
// - Walking slowly: 2.0 met
met         1.2;

// External work [met]
// Usually 0 for most activities
wme         0;

// Relative humidity [%]
// Used only if humidity field is not available in the case
// Range: 0-100
RH          50;

// ************************************************************************* //
```

## Required Input Fields

The following fields must be present in your OpenFOAM case:

- **T** (Temperature) [K] - MUST_READ
- **U** (Velocity) [m/s] - MUST_READ
- **p** (Pressure) [Pa] - READ_IF_PRESENT

### Optional Humidity Fields

comfortFoam searches for humidity fields in this order:
1. `RH` - Relative humidity [%]
2. `relativeHumidity` - Relative humidity [%]
3. `thermo:relHum` - From buoyantHumiditySimpleFoam
4. `relHum` - Alternative naming
5. Falls back to value in `comfortFoamDict`

### Optional Turbulence Fields

For draft rating calculation:
- `k` - Turbulent kinetic energy [m²/s²]
- `epsilon` - Turbulent dissipation rate [m²/s³] (k-epsilon models)
- `omega` - Specific dissipation rate [1/s] (k-omega SST models)
- If not present, assumes Tu = 40% (typical for mixed convection)

### Optional Radiation Fields

For accurate mean radiant temperature calculation (in order of preference):
- `G` - Incident radiation field [W/m²] from radiation models
- `qr` - Radiative heat flux field [W/m²] from some radiation models
- `IDefault` - Default radiation intensity [W/m²/sr] from DOM models
- If not present, uses area-weighted wall temperature

## Output Fields

comfortFoam creates four fields in each time directory:

1. **PMV** - Predicted Mean Vote [-3 to +3]
2. **PPD** - Predicted Percentage of Dissatisfied [%]
3. **DR** - Draft Rating [%]
4. **TOp** - Operative Temperature [K]

## Examples

### Example 1: Basic comfort analysis
```bash
# Run CFD simulation first
buoyantHumiditySimpleFoam

# Calculate comfort indices
comfortFoam

# Visualize in ParaView
paraFoam
```

### Example 2: Regional analysis
```bash
# Calculate comfort for a specific zone
comfortFoam -cellZone occupiedZone
```

### Example 3: Multi-region case
```bash
# Calculate comfort for a specific region
comfortFoam -region indoorAir
```

## Physical Parameters

### Environmental Parameters
- Air temperature: 10-30°C (50-86°F)
- Mean radiant temperature: Assumed equal to air temperature
- Air velocity: 0-1 m/s
- Relative humidity: 0-100%

### Personal Parameters
- Clothing insulation: 0-2 clo
- Metabolic rate: 0.8-4 met
- External work: Usually 0 met

## Limitations

1. **Steady-State Model**: PMV/PPD assumes thermal equilibrium
2. **Uniform Conditions**: Best suited for spaces with relatively uniform conditions
3. **Activity Level**: Assumes constant metabolic rate
4. **Comfort Categories**: For EN 15251 comfort categories, use ASHRAE55Foam

## Validation

The implementation has been validated against:
- ISO 7730:2005 reference tables
- ASHRAE Fundamentals Handbook
- Published PMV/PPD calculation tools

### Validation Example

To verify ISO 7730 compliance:
```bash
comfortFoam -validate -validateAirTemp 22 -validateRadTemp 22 \
            -validateVelocity 0.1 -validateRH 60 \
            -validateMet 1.2 -validateClo 0.5
# Expected: PMV ≈ -0.75, PPD ≈ 17%
```

## References

1. Fanger, P.O. (1970). "Thermal Comfort: Analysis and Applications in Environmental Engineering"
2. ISO 7730:2005. "Ergonomics of the thermal environment"
3. ASHRAE Standard 55-2020. "Thermal Environmental Conditions for Human Occupancy"

## Troubleshooting

### Common Issues

1. **Missing humidity field**: Ensure your solver outputs humidity or set RH in comfortFoamDict
2. **Out-of-range PMV**: Check input parameters (clothing, activity, temperature)
3. **High DR values**: May indicate excessive air velocities or temperature gradients

### Debug Mode

For detailed output, modify Info statements in the source code and recompile.

## Version History

- v1.0: Initial implementation for OpenFOAM v2312
- v1.1: Added cellSet/cellZone support
- v1.2: Improved humidity field detection
- v2.0: Updated for OpenFOAM v2412+
- v3.0: Added validation mode and improved ISO 7730 compliance
- v3.1: Enhanced turbulence model support (k-epsilon, k-omega SST)
- v3.2: Advanced radiation field handling (G, qr, IDefault)

## Support

For issues or questions:
- Check the [GitHub repository](https://github.com/tian-pvam/HVAC-for-OpenFOAM)
- Review tutorial cases in `tutorials/` directory
- Open an issue for bugs or feature requests