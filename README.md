# HVAC-for-OpenFOAM

[![OpenFOAM v2412+](https://img.shields.io/badge/OpenFOAM-v2412+-blue.svg)](https://www.openfoam.com/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Advanced HVAC simulation toolkit for OpenFOAM, providing specialized solvers and utilities for humidity modeling, thermal comfort analysis, and indoor climate assessment according to international standards (ISO 7730, ASHRAE-55, UTCI).

## 🌟 Features

### Core Capabilities
- **Humidity Transport Modeling**: Full integration of humidity physics into OpenFOAM's thermophysical framework
- **Thermal Comfort Analysis**: Comprehensive comfort metrics (PMV, PPD, UTCI, ASHRAE-55)
- **Solar Radiation Integration**: Support for EPW weather data and OpenFOAM's native solar models
- **Buoyancy-Driven Flows**: Specialized solvers for natural and mixed convection with humidity effects
- **Building Physics**: Advanced boundary conditions for multi-layer walls with thermal mass
- **Wind-Driven Rain**: Lagrangian particle tracking for rain droplet impingement

### Components

#### 🔧 CFD Solvers
- **buoyantHumidityPimpleFoam**: Transient solver for buoyant, turbulent flow with humidity transport
- **buoyantHumiditySimpleFoam**: Steady-state solver using SIMPLE algorithm
- **buoyantBoussinesqPimpleDyMFoam**: Dynamic mesh solver for Boussinesq approximation
- **windDrivenRainFoam**: Wind-driven rain simulation with multiple droplet phases

#### 📊 Comfort Analysis Tools
- **comfortFoam**: ISO 7730 thermal comfort (PMV, PPD, Draft Rating)
- **UTCIFoam**: Universal Thermal Climate Index calculator
- **ASHRAE55Foam**: ASHRAE Standard 55 compliance analysis
- **AoAFoam**: Age of Air calculation for ventilation assessment

#### 🏗️ Boundary Conditions
- **buildingElementBC**: Advanced thermal boundary condition for building walls with multi-layer support

#### 📚 Libraries
- **humidityRhoThermo**: Thermophysical model extension for humidity calculations
- **solarCalculator**: Solar position and radiation calculations

## 🚀 Quick Start

### Prerequisites
- OpenFOAM v2412 or later
- Standard OpenFOAM development environment
- C++ compiler (GCC/Clang)

### Installation

```bash
# Clone the repository
git clone https://github.com/tian-pvam/HVAC-for-OpenFOAM.git
cd HVAC-for-OpenFOAM

# IMPORTANT: Build order matters due to dependencies

# 1. First, install the humidity library into OpenFOAM source
cd $FOAM_SRC/thermophysicalModels/basic
# Add these lines to Make/files after liquidThermo.C:
# humidityRhoThermo/humidityRhoThermo.C
# humidityRhoThermo/humidityRhoThermos.C 
# humidityRhoThermo/derivedFvPatchFields/fixedHumidity/fixedHumidityFvPatchScalarField.C
wmake libso

# 2. Build CFD solvers with humidity support
cd $FOAM_RUN/../HVAC-for-OpenFOAM/buoyantHumidityPimpleFoam && wmake
cd ../buoyantHumiditySimpleFoam && wmake

# 3. Build comfort analysis utilities
cd ../comfortFoam && wmake
cd ../UTCIFoam && wmake
cd ../ASHRAE55 && wmake

# 4. Build boundary conditions
cd ../buildingElementBC && wmake

# 5. Build additional utilities as needed
cd ../AoAFoam && wmake
cd ../solarCalculator && wmake libso
cd ../windDrivenRainFoam && wmake
```

### Basic Usage Example

```bash
# Run a humidity test case
cd humidityRhoThermo/tutorials/laminar/fixedHumidityBC
./clean  # Clean previous results
./run    # Run mesh conversion and solver
paraFoam # Visualize results

# Post-process with comfort analysis
comfortFoam  # Calculates PMV, PPD, DR indices
UTCIFoam     # Calculates UTCI index
ASHRAE55Foam # ASHRAE-55 compliance check
```

## 📖 Documentation

### Solver Documentation
- [buoyantHumidityPimpleFoam](buoyantHumidityPimpleFoam/README.md) - Transient humidity solver
- [buoyantHumiditySimpleFoam](buoyantHumiditySimpleFoam/README.md) - Steady-state humidity solver
- [humidityRhoThermo](humidityRhoThermo/README.md) - Humidity thermophysical library
- [windDrivenRainFoam](windDrivenRainFoam/README.md) - Wind-driven rain simulation

### Comfort Analysis Tools
- [comfortFoam](comfortFoam/README.md) - ISO 7730 comfort metrics
- [UTCIFoam](UTCIFoam/README.md) - UTCI calculation with weather data support
- [ASHRAE55Foam](ASHRAE55/README.md) - ASHRAE-55 adaptive comfort model

### Boundary Conditions
- [buildingElementBC](buildingElementBC/README.md) - Multi-layer wall boundary condition with radiation

### Tutorials
Each component includes tutorial cases in its `tutorials/` directory demonstrating typical usage scenarios.

## 🔬 Technical Details

### Humidity Modeling
The humidity transport is integrated into OpenFOAM's thermophysical framework through:
- Extended equation of state with humidity effects
- Humidity-dependent thermophysical properties
- Custom boundary conditions (e.g., `fixedHumidity`)
- Fields: relative humidity, specific humidity, water vapor partial pressure

### Comfort Calculations
All comfort tools follow established standards:
- **ISO 7730**: Fanger's PMV/PPD model with local draft rating
- **UTCI**: 6th-order polynomial approximation with radiation effects
- **ASHRAE-55**: Adaptive comfort model for naturally ventilated spaces

### Workflow Integration

```
CFD Simulation → Humidity Calculation → Comfort Analysis → Visualization
     ↓                    ↓                    ↓                ↓
buoyantHumidity*   humidityRhoThermo    comfortFoam      ParaView
                   ::correct()          UTCIFoam
                                       ASHRAE55Foam
```

## 🛠️ Development

### Coding Standards
- Follow OpenFOAM coding conventions
- Use existing patterns for new implementations
- Maintain backward compatibility where possible

### Contributing
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📚 Applications

- HVAC system design and optimization
- Indoor air quality assessment
- Natural ventilation design
- Building energy simulation with dynamic thermal mass
- Urban microclimate studies
- Thermal comfort evaluation
- Compliance checking with international standards
- Building facade performance analysis
- Rain penetration and moisture risk assessment

## 🤝 Support

For issues, questions, or contributions:
- Open an issue on [GitHub](https://github.com/tian-pvam/HVAC-for-OpenFOAM/issues)
- Check the [Wiki](https://github.com/tian-pvam/HVAC-for-OpenFOAM/wiki) for detailed guides
- Contact the maintainers

## 📄 License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- OpenFOAM Foundation and ESI-OpenCFD for the OpenFOAM framework
- ASHRAE for thermal comfort standards
- UTCI consortium for the Universal Thermal Climate Index
- All contributors and users of this toolkit

## 📊 Citation

If you use HVAC-for-OpenFOAM in your research, please cite:

```bibtex
@software{hvac_for_openfoam,
  author = {Thomas Tian},
  title = {HVAC-for-OpenFOAM: Advanced HVAC Simulation Toolkit},
  year = {2024},
  url = {https://github.com/tian-pvam/HVAC-for-OpenFOAM}
}
```

---

**Note**: This project is under active development. Check the release notes for the latest features and improvements.