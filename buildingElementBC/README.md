# BuildingElement Boundary Condition

A comprehensive thermal boundary condition for building walls in OpenFOAM that accounts for:
- Conduction through wall (U-value or multi-layer)
- Internal and external convection
- Solar radiation gains
- Long-wave radiation to sky
- Integration with OpenFOAM radiation models

## Features

### Simple Mode (U-value)
- Specify overall thermal transmittance (U-value)
- Quick setup for standard walls

### Multi-layer Mode
- Define individual layer properties (conductivity, density, specific heat)
- Transient heat storage in wall layers
- Automatic temperature distribution calculation

### Radiation Support
- Solar gains through g-value
- Sky radiation with view factor
- Integration with OpenFOAM radiation models (e.g., P1, viewFactor)

### Debug Options
- Temperature and heat flux monitoring
- Face-specific or patch-wide statistics
- Configurable output intervals

## Usage

### Simple Mode Example
```
wall
{
    type            buildingElement;
    
    // Thermal properties
    uValue          0.24;        // W/(m²K)
    gValue          0.6;         // Solar heat gain coefficient
    
    // External conditions
    tempExt         283;         // External air temperature [K]
    tempSky         273;         // Sky temperature [K]
    solarExt        500;         // Solar radiation [W/m²]
    hExt            25;          // External convection coefficient [W/(m²K)]
    
    // Radiation properties
    emissivity      0.9;         // Surface emissivity
    viewSky         0.5;         // View factor to sky
    radField        qr;          // Internal radiation field name (optional)
    
    // Initial value
    value           uniform 293;
}
```

### Multi-layer Mode Example
```
wall
{
    type            buildingElement;
    
    // Layer definition
    layers
    (
        // Inner layer (e.g., plaster)
        {
            lambda      0.7;         // Thermal conductivity [W/(mK)]
            thickness   0.02;        // Layer thickness [m]
            rho         1200;        // Density [kg/m³]
            cp          1000;        // Specific heat [J/(kgK)]
            nNodes      5;           // Discretization nodes
        }
        // Insulation layer
        {
            lambda      0.04;
            thickness   0.2;
            rho         30;
            cp          1000;
            nNodes      10;
        }
        // Outer layer (e.g., concrete)
        {
            lambda      1.4;
            thickness   0.1;
            rho         2200;
            cp          900;
            nNodes      5;
        }
    );
    
    gValue          0.6;
    tempExt         283;
    tempSky         273;
    solarExt        500;
    hExt            25;
    emissivity      0.9;
    viewSky         0.5;
    radField        qr;
    
    value           uniform 293;
}
```

### Debug Options
```
    // Debug options (optional)
    debugTemp       true;        // Output temperatures
    debugFlux       true;        // Output heat fluxes
    debugFace       0;           // Face index to monitor (-1 for none)
    debugInterval   60;          // Output interval [s]
```

## Compilation

```bash
cd buildingElementBC
wmake
```

## Integration with Solvers

Add to your solver's controlDict:
```
libs ("libbuildingElementBC.so");
```

Or for Windows:
```
libs ("libbuildingElementBC.dll");
```

## Heat Flux Convention

- Positive flux = heat entering the wall from inside
- Negative flux = heat leaving the wall to inside

## Output Variables

When debug is enabled:
- `q_conv_int`: Interior convective heat flux
- `q_rad_int`: Interior radiative heat flux (from OpenFOAM radiation model)
- `q_conduction`: Conduction through wall
- `q_conv_ext`: Exterior convective heat flux
- `q_rad_ext`: Exterior radiative heat flux to sky/environment
- `q_solar`: Solar heat gains
- `q_total`: Net heat flux balance

## Notes

- The boundary condition automatically handles the interaction with OpenFOAM radiation models
- For multi-layer mode, ensure time step satisfies stability criteria (Fourier number < 0.5)
- External surface resistance (0.04 m²K/W) and internal surface resistance (0.13 m²K/W) are included in U-value calculation