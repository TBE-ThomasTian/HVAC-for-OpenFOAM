# AoAFoam - Age of Air Solver for OpenFOAM

## Overview

AoAFoam is a transient solver for calculating the Age of Air (AoA) in ventilated spaces. It is designed for HVAC applications to evaluate ventilation effectiveness and air quality distribution in buildings.

## Features

- **Transient Age of Air calculation** with unit source term
- **PIMPLE algorithm** for time-stepping with adjustable time step control
- **Turbulence model support** for enhanced diffusion modeling
- **Ventilation efficiency calculation** (η = τ̄/(2·τ))
- **Parallel execution support** via MPI
- **fvOptions framework** for additional source terms
- **Automatic calculation** of mean age and local ventilation efficiency

## Theory

The solver solves the transport equation for Age of Air:

```
∂τ/∂t + ∇·(Uτ) - ∇·((ν_t + D_T)∇τ) = 1
```

Where:
- τ = Age of Air [s]
- U = Velocity field [m/s]
- ν_t = Turbulent kinematic viscosity [m²/s]
- D_T = Molecular diffusivity [m²/s]

## Installation

### Prerequisites
- OpenFOAM v2412 or later
- Compiled HVAC-for-OpenFOAM humidity library (if using humidity solvers)

### Compilation
```bash
cd $FOAM_RUN/../applications/solvers/HVAC/AoAFoam
wmake
```

## Usage

### Basic Usage
```bash
AoAFoam
```

### Parallel Execution
```bash
decomposePar
mpirun -np 4 AoAFoam -parallel
reconstructPar
```

### Required Files

#### Initial Conditions (0/)
- `U`: Velocity field (from CFD simulation or prescribed)
- `AoA`: Age of Air field (typically initialized to 0)

#### Constant Directory
- `transportProperties`: Contains diffusivity `DT`
- `turbulenceProperties`: Turbulence model selection

#### System Directory
- Standard OpenFOAM dictionaries: `controlDict`, `fvSchemes`, `fvSolution`

### Boundary Conditions

Typical boundary conditions for AoA:
- **Inlets**: `fixedValue` with `value uniform 0` (fresh air)
- **Outlets**: `zeroGradient`
- **Walls**: `zeroGradient`

## Output Fields

The solver automatically calculates and writes:
- `AoA`: Age of Air distribution [s]
- `localMeanAge`: Local mean age field [s]
- `ventilationEfficiency`: Local ventilation efficiency [-]

Console output includes:
- Mean Age of Air for the entire domain
- Min/Max AoA values
- Execution time statistics

## Example Case

A tutorial case is provided in `tutorials/roomVentilation/`:
```bash
cd tutorials/roomVentilation
./Allrun
paraFoam
```

This simulates a simple room (4×3×2.5m) with ceiling inlet and floor outlet.

## Transport Properties

Example `transportProperties` file:
```
DT              DT [ 0 2 -1 0 0 0 0 ] 2.5e-5;  // Molecular diffusivity [m²/s]
nu              nu [ 0 2 -1 0 0 0 0 ] 1.5e-5;  // Kinematic viscosity [m²/s]
```

## Ventilation Efficiency

The solver calculates ventilation efficiency as:
```
η = τ̄/(2·τ_local)
```

Where:
- η > 1: Better than perfect mixing
- η = 1: Perfect mixing
- η < 1: Worse than perfect mixing
- η = 0.5: Piston flow

## Tips for Usage

1. **Steady-State Flow**: First run a steady-state flow solver (simpleFoam, buoyantSimpleFoam) to obtain the velocity field
2. **Time Step**: Use adjustable time stepping with `maxCo` around 0.5 for stability
3. **Simulation Time**: Run until AoA field reaches steady state (typically 3-5 air changes)
4. **Mesh Quality**: Ensure good mesh quality for accurate diffusion calculations

## Common Applications

- Office ventilation assessment
- Cleanroom air age distribution
- Hospital room air quality
- Industrial ventilation effectiveness
- Displacement ventilation analysis
- Natural ventilation studies

## Troubleshooting

1. **Negative AoA values**: Check initial conditions and boundary conditions
2. **Divergence**: Reduce time step or relaxation factors
3. **Unrealistic values**: Verify velocity field and diffusivity settings

## References

- Sandberg, M. (1981). "What is ventilation efficiency?" Building and Environment, 16(2), 123-135.
- Mundt, E., et al. (2004). "Ventilation Effectiveness." REHVA Guidebook No. 2.

## License

This solver is part of HVAC-for-OpenFOAM and follows the GNU General Public License v3.0.