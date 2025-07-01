# windDrivenRainFoam

## Description

windDrivenRainFoam is an OpenFOAM solver for simulating wind-driven rain (WDR) using an Eulerian multiphase approach. The solver calculates how raindrops of different sizes are transported through a wind field and determines rain loads on building facades.

**Author**: Aytac Kubilay, ETH Zurich/Empa (March 2012)  
**OpenFOAM Version**: v2412+

## Features

- **Multiphase rain modeling**: Simulates multiple raindrop size classes simultaneously
- **Drag force modeling**: Reynolds number-dependent drag coefficient based on Gunn & Kinzer experimental data
- **Turbulent dispersion**: Optional modeling of turbulence effects on raindrop trajectories using k-epsilon model
- **Catch ratio calculation**: Computes specific and global catch ratios for building facades
- **Temperature-dependent properties**: Accounts for temperature effects on air and water properties

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
```
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
Uref            10.0;    // Reference wind velocity [m/s]
Rh              10.0;    // Horizontal rainfall intensity [mm/h]
temperature     283.15;  // Temperature [K]
windScaling     1.0;     // Wind field scaling factor [-]

// Turbulent dispersion
Ct              0.0;     // Turbulent dispersion switch (0=off, 1=on)
```

## Usage

1. Prepare a steady-state wind field using simpleFoam or another appropriate solver
2. Configure rain phases in `constant/transportProperties`
3. Set appropriate boundary conditions for alpharain and Urain fields
4. Run the solver:
   ```bash
   windDrivenRainFoam
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