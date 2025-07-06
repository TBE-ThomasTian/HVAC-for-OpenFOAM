# windDrivenRainFoam

## Description

windDrivenRainFoam is an advanced OpenFOAM solver for simulating wind-driven rain (WDR) using an Eulerian multiphase approach. The solver calculates how raindrops of different sizes are transported through a wind field and determines rain loads on building facades.

**Original Author**: Aytac Kubilay, ETH Zurich/Empa (March 2012)  
**Enhanced Version**: 2024  
**OpenFOAM Version**: v2412+

## Features

### Core Features
- **Multiphase rain modeling**: Simulates multiple raindrop size classes simultaneously
- **Drag force modeling**: Reynolds number-dependent drag coefficient based on Gunn & Kinzer experimental data
- **Catch ratio calculation**: Computes specific and global catch ratios for building facades
- **Temperature-dependent properties**: Accounts for temperature effects on air and water properties

### Enhanced Features (2024)
- **Generic turbulence modeling**: Supports all OpenFOAM turbulence models (RAS/LES), not limited to k-epsilon
- **Adaptive time stepping**: CFL-based time step control for rain phases ensures numerical stability
- **Custom boundary conditions**: 
  - `windDrivenRainInlet`: Specialized inlet BC for rain with wind angle effects
  - `catchRatio`: Direct calculation of catch ratio at building surfaces
- **Modern C++ implementation**: Uses modern OpenFOAM constructs and memory management
- **Parallel-ready architecture**: Structured for efficient MPI parallelization
- **Extended droplet physics** (available in library):
  - Evaporation modeling with Ranz-Marshall correlation
  - Breakup criteria (bag, shear breakup) based on Weber number
  - Collision and coalescence efficiency
  - Drag enhancement due to droplet deformation

## Physics

The solver implements:
1. Transport equations for rain volume fraction (alpharain) and momentum (Urain) for each droplet phase
2. Drag force between air and raindrops with Reynolds-dependent drag coefficient
3. Gravitational settling of raindrops
4. Optional turbulent dispersion based on fluid and particle time scales
5. Calculation of specific catch ratio (scr) and global catch ratio (gcr)

## Compilation

```bash
cd $FOAM_RUN/../applications/solvers/windDrivenRainFoam
wmake
```

## Required Files

### Case Directory Structure
```
case/
├── 0/
│   ├── U           # Wind velocity field (pre-computed steady-state)
│   ├── p           # Pressure field
│   ├── k           # Turbulent kinetic energy (if using turbulence)
│   └── epsilon     # Turbulent dissipation rate (if using turbulence)
├── constant/
│   ├── g           # Gravitational acceleration
│   └── transportProperties  # Rain phase properties
└── system/
    ├── controlDict
    ├── fvSchemes
    └── fvSolution
```

### transportProperties Dictionary

Example configuration:
```c++
phases
(
    // Phase 1: Small droplets (0.5 mm diameter)
    (0.5e-3  0.05)  // (diameter [m], volume fraction [-])
    
    // Phase 2: Medium droplets (1.0 mm diameter)  
    (1.0e-3  0.10)
    
    // Phase 3: Large droplets (2.0 mm diameter)
    (2.0e-3  0.15)
);

// Reference values
Rh              10.0;    // Horizontal rainfall intensity [mm/h]
temp            283.15;  // Temperature [K] (updated: use 'temp' not 'temperature')
scalingFactor   1.0;     // Wind field scaling factor [-] (updated: use 'scalingFactor')

// Turbulent dispersion
solveTD         true;    // Enable turbulent dispersion (updated: boolean flag)
```

### controlDict Settings (for adaptive time stepping)

```c++
// Add to system/controlDict for adaptive time stepping
adjustTimeStep  yes;
maxCoRain       0.5;     // Maximum Courant number for rain phases
maxDeltaT       1.0;     // Maximum time step [s]
```

## Usage

1. Prepare a steady-state wind field using simpleFoam or another appropriate solver
2. Configure rain phases in `constant/transportProperties`
3. Set appropriate boundary conditions for alpharain and Urain fields
4. Select turbulence model in `constant/turbulenceProperties` (if using turbulent dispersion)
5. Run the solver:
   ```bash
   windDrivenRainFoam
   ```

### Boundary Conditions Example

```c++
// 0/U1 - Rain velocity for phase 1
boundaryField
{
    inlet
    {
        type            windDrivenRainInlet;
        terminalVelocity 4.0;  // Terminal velocity [m/s]
        windAngle       30;    // Wind angle [degrees]
        value           uniform (0 0 0);
    }
    
    building
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
    
    outlet
    {
        type            inletOutlet;
        inletValue      uniform (0 0 0);
    }
}

// 0/gcr - Global catch ratio
boundaryField
{
    building
    {
        type            catchRatio;
        Rh              10;    // Must match transportProperties
        phi             phi1;  // Flux field name
        value           uniform 0;
    }
}
```

## Output

The solver produces:
- **alpharain**: Volume fraction fields for each rain phase
- **Urain**: Velocity fields for each rain phase
- **scr**: Specific catch ratio for each phase
- **gcr**: Global catch ratio (sum of all phases)

## Applications

- Building physics: Assessment of rain loads on facades
- Urban planning: WDR analysis for pedestrian comfort
- Building envelope design: Moisture ingress studies
- Heritage conservation: Rain impact on historical buildings

## Notes

- The solver assumes a pre-computed steady-state wind field
- Multiple raindrop sizes can be simulated to represent realistic rain spectra
- Catch ratios help determine rain deposition patterns on building surfaces
- Temperature effects on fluid properties are included for accurate modeling

## References

Kubilay, A., Derome, D., Blocken, B., & Carmeliet, J. (2013). CFD simulation and validation of wind-driven rain on a building facade with an Eulerian multiphase model. Building and Environment, 61, 69-81.