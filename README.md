# HVAC-for-OpenFOAM

[![OpenFOAM v2412+](https://img.shields.io/badge/OpenFOAM-v2412+-blue.svg)](https://www.openfoam.com/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Advanced HVAC simulation toolkit for OpenFOAM, providing specialized solvers and utilities for humidity modeling, thermal comfort analysis, and indoor climate assessment according to international standards (ISO 7730, ASHRAE-55, UTCI).

## Features

### Core Capabilities
- **Humidity Transport Modeling**: Full integration of humidity physics into OpenFOAM's thermophysical framework
- **Thermal Comfort Analysis**: Comprehensive comfort metrics (PMV, PPD, UTCI, ASHRAE-55)
- **Solar Radiation Integration**: Support for EPW weather data and OpenFOAM's native solar models
- **Buoyancy-Driven Flows**: Specialized solvers for natural and mixed convection with humidity effects
- **Building Physics**: Advanced boundary conditions for multi-layer walls with thermal mass
- **Wind-Driven Rain**: Lagrangian particle tracking for rain droplet impingement

### Components

#### CFD Solvers
- **buoyantHumidityPimpleFoam**: Transient solver for buoyant, turbulent flow with humidity transport
- **buoyantHumiditySimpleFoam**: Steady-state solver using SIMPLE algorithm
- **buoyantBoussinesqPimpleDyMFoam**: Dynamic mesh solver for Boussinesq approximation
- **windDrivenRainFoam**: Wind-driven rain simulation with multiple droplet phases

#### Comfort Analysis Tools
- **comfortFoam**: ISO 7730 thermal comfort (PMV, PPD, Draft Rating)
- **UTCIFoam**: Universal Thermal Climate Index calculator
- **ASHRAE55Foam**: ASHRAE Standard 55 compliance analysis
- **ISO7243Foam**: ISO 7243 WBGT heat stress assessment
- **ISO7933Foam**: ISO 7933 heat strain screening
- **AoAFoam**: Age of Air calculation for ventilation assessment
- **NEN8100Foam**: Multi-direction pedestrian wind comfort and danger assessment

#### fvOptions
- **thermostatSource**: Thermostat-controlled heating source with proportional regulation, sensor placement, and supply temperature limiting (Vorlauf)

#### Boundary Conditions
- **buildingElementBC**: Advanced thermal boundary condition for building walls with multi-layer support
- **heatFluxRadiation**: Lumped wall-temperature boundary condition with prescribed total power, optional `qr` coupling, and temperature limits for radiative HVAC loads

#### Libraries
- **humidityRhoThermo**: Thermophysical model extension for humidity calculations
- **solarCalculator**: Solar position and radiation calculations

## Documentation

### Solver Documentation
- [buoyantHumidityPimpleFoam](buoyantHumidityPimpleFoam/README.md) - Transient humidity solver
- [buoyantHumiditySimpleFoam](buoyantHumiditySimpleFoam/README.md) - Steady-state humidity solver
- [humidityRhoThermo](humidityRhoThermo/README.md) - Humidity thermophysical library
- [windDrivenRainFoam](windDrivenRainFoam/README.md) - Wind-driven rain simulation

### Comfort Analysis Tools
- [comfortFoam](comfortFoam/README.md) - ISO 7730 comfort metrics
- [UTCIFoam](UTCIFoam/README.md) - UTCI calculation with weather data support
- [ASHRAE55Foam](ASHRAE55/README.md) - ASHRAE-55 adaptive comfort model
- [ISO7243Foam](ISO7243/README.md) - ISO 7243 WBGT post-processing
- [ISO7933Foam](ISO7933/README.md) - ISO 7933 heat strain screening
- [NEN8100Foam](NEN8100/README.md) - NEN 8100 pedestrian wind comfort and danger

### fvOptions
- [thermostatSource](thermostatSource/README.md) - Thermostat-controlled heating with sensor and Vorlauf limit

### Boundary Conditions
- [buildingElementBC](buildingElementBC/README.md) - Multi-layer wall boundary condition with radiation

## Additional ThermoTools Patch Set

This repository also contains an additional `thermoTools` patch tree at:

- `src/thermoTools/`

The patch currently adds:

- `heatFluxRadiation`: a patch-temperature boundary condition for cases that need a prescribed total power (for example occupants or equipment) while allowing the surface temperature to react to convection and radiation
- temperature guards for `externalWallHeatFluxTemperature` via the optional entries `minTemperature` and `maxTemperature`

### `heatFluxRadiation`

`heatFluxRadiation` is intended for single-region cases where `fixedValue` is too restrictive, but `externalWallHeatFluxTemperature` in `mode flux/power` with `kappaMethod fluidThermo` can drive the patch temperature to unphysical values.

The model advances a single patch temperature from a lumped energy balance:

`m Cp dT/dt = Q + integral(q dA) + Qrad - Qconv`

with:

- `Q`: prescribed total power input in W
- `q`: optional prescribed heat flux in W/m^2
- `Qrad`: integrated contribution from the optional `qr` field
- `Qconv`: heat loss to the adjacent fluid based on the local patch gradient

Example:

```foam
person
{
    type            heatFluxRadiation;
    kappaMethod     fluidThermo;
    Q               constant 60;
    mass            8;
    Cp              3500;
    qr              qr;
    qrRelaxation    0.7;
    relaxation      1.0;
    minTemperature  280;
    maxTemperature  320;
    value           uniform 304.15;
}
```

For area-based cooling or heating loads, `q` can be prescribed directly:

```foam
floor
{
    type            heatFluxRadiation;
    kappaMethod     fluidThermo;
    q               uniform -20;
    mass            500;
    Cp              900;
    qr              qr;
    value           uniform 295.15;
}
```

Validation cases for the boundary condition are provided under:

- `src/thermoTools/derivedFvPatchFields/heatFluxRadiation/validation/prescribedQr`
- `src/thermoTools/derivedFvPatchFields/heatFluxRadiation/validation/fvDOMParallelPlates`

`prescribedQr` validates the pure `qr` source term against an exact
first-step balance. `fvDOMParallelPlates` validates the coupled
`fvDOM + heatFluxRadiation` response against the black parallel-plates
solution, with the radiative flux appearing in time `1` and the patch
temperature update in time `2`.

### Temperature Guards for `externalWallHeatFluxTemperature`

`externalWallHeatFluxTemperature` now also supports:

- `minTemperature`
- `maxTemperature`

These bounds are intended as a numerical safeguard when a strongly negative `qr` or a very small effective conductivity would otherwise drive the patch temperature below physically meaningful values.

Example:

```foam
wall
{
    type            externalWallHeatFluxTemperature;
    mode            flux;
    qr              qr;
    q               uniform 0;
    qrRelaxation    0.7;
    relaxation      1.0;
    kappaMethod     function;
    kappaValue      constant 0.2;
    alphaValue      constant 1e-7;
    minTemperature  250;
    maxTemperature  350;
    value           uniform 299.15;
}
```

## OpenFOAM-v2512 Patch Set

This repository contains a ready-to-copy patch tree for OpenFOAM v2512 at:

- `OpenFOAM-v2512/src/parallel/distributed/distributedTriSurfaceMesh/`
- `OpenFOAM-v2512/src/thermophysicalModels/radiation/`

The patch extends solar radiation handling for partially transparent surfaces (`0 < t < 1`) by:

- enabling transmissivity from `transparent` boundary radiation properties
- propagating ray transmissivity through `solarLoad` face shading
- scaling direct and reflected solar load contributions consistently
- clipping distributed ray-search segments to processor-local bounds

To apply in an OpenFOAM-v2512 source tree, copy the matching files from this folder into the OpenFOAM source and rebuild the affected OpenFOAM libraries with `wmake`.

## Technical Details

### Humidity Modeling
The humidity transport is integrated into OpenFOAM's thermophysical framework through:
- Extended equation of state with humidity effects
- Humidity-dependent thermophysical properties
- Custom boundary conditions (e.g., `fixedHumidity`)
- Fields: relative humidity, specific humidity, water vapor partial pressure

### Comfort Calculations
All comfort tools follow established standards:
- **ISO 7730**: Fanger's PMV/PPD model with local draft rating
- **ISO 7243**: WBGT-based heat stress assessment
- **ISO 7933**: Heat strain screening based on required evaporation
- **UTCI**: 6th-order polynomial approximation with radiation effects
- **ASHRAE-55**: Adaptive comfort model for naturally ventilated spaces

## Applications

- HVAC system design and optimization
- Indoor air quality assessment
- Natural ventilation design
- Building energy simulation with dynamic thermal mass
- Urban microclimate studies
- Thermal comfort evaluation
- Compliance checking with international standards
- Building facade performance analysis
- Rain penetration and moisture risk assessment

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.
