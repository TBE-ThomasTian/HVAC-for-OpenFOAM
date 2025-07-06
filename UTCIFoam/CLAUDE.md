# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is an OpenFOAM-based project focused on HVAC simulations and thermal comfort analysis. It provides specialized solvers and utilities for modeling humidity, buoyancy-driven flows, and calculating various thermal comfort indices (UTCI, ASHRAE-55, PMV/PPD).

## Key Commands

### Building Components

```bash
# Build individual solvers (from solver directory)
wmake

# Build libraries (e.g., from humidityRhoThermo/)
wmake libso

# Clean build
wclean

# Build all components in order:
# 1. First build the humidity library
cd ../humidityRhoThermo && wmake libso

# 2. Then build solvers that depend on it
cd ../buoyantHumidityPimpleFoam && wmake
cd ../buoyantHumiditySimpleFoam && wmake

# 3. Build comfort calculation utilities
cd ../UTCIFoam && wmake
cd ../ASHRAE55 && wmake
cd ../comfortFoam && wmake
```

### Running Tests

```bash
# Navigate to tutorial case (example)
cd ../humidityRhoThermo/tutorials/laminar/fixedHumidityBC

# Clean previous results
./clean

# Run the case (includes mesh conversion and solver execution)
./run

# View results with ParaView
paraFoam
```

## Architecture Overview

### Project Structure

The repository contains several interconnected components:

1. **Core Library**: `humidityRhoThermo` - Provides thermophysical models for humidity calculations that other solvers depend on

2. **CFD Solvers with Humidity**:
   - `buoyantHumidityPimpleFoam`: Transient solver using PIMPLE algorithm
   - `buoyantHumiditySimpleFoam`: Steady-state solver using SIMPLE algorithm
   - Both extend standard OpenFOAM buoyant solvers with humidity transport

3. **Thermal Comfort Calculators**:
   - `UTCIFoam`: Calculates Universal Thermal Climate Index
   - `ASHRAE55`: Implements ASHRAE-55 thermal comfort standards
   - `comfortFoam`: Calculates PMV/PPD according to DIN EN ISO 7730
   - These are post-processing tools that read CFD results and calculate comfort indices

4. **Additional Utilities**:
   - `solarCalculator`: Calculates solar radiation and sun position
   - `AoAFoam`: Age of Air calculation solver

### Key Design Patterns

1. **Solver Structure**: All solvers follow standard OpenFOAM patterns with:
   - `createFields.H` for field initialization
   - Separate equation files (UEqn.H, TEqn.H, EEqn.H, pEqn.H)
   - Standard time loop with PIMPLE/SIMPLE iterations

2. **Humidity Integration**: The humidity functionality is implemented as:
   - Custom thermophysical model extending standard OpenFOAM models
   - Additional transport equation for humidity ratio
   - Custom boundary conditions (e.g., `fixedHumidityFvPatchScalarField`)

3. **Comfort Calculations**: Post-processing utilities that:
   - Read existing temperature, velocity, humidity fields
   - Calculate comfort indices based on established standards
   - Write results as new fields for visualization

### Build Dependencies

1. Libraries must be built before solvers that depend on them
2. Installation locations:
   - System solvers: `$FOAM_APPBIN` (UTCIFoam, ASHRAE55, comfortFoam)
   - User solvers: `$FOAM_USER_APPBIN` (buoyantHumidity solvers)
   - Libraries: `$FOAM_LIBBIN` or `$FOAM_USER_LIBBIN`

### OpenFOAM Version Compatibility

The codebase is maintained for OpenFOAM v2312, with recent updates ensuring compatibility. The Make directories show cross-platform support including MinGW builds.