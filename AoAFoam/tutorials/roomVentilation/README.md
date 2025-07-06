# Room Ventilation Test Case for AoAFoam

## Description
Simple room ventilation test case for Age of Air (AoA) calculation.

### Geometry
- Room dimensions: 4m x 3m x 2.5m (L x W x H)
- Inlet: Top wall (ceiling)
- Outlet: Bottom wall (floor)
- Mesh: 40 x 30 x 25 cells

### Physics
- Inlet velocity: 0.5 m/s downward
- Laminar flow
- Diffusivity DT: 2.5e-5 m2/s

### Running the case
```bash
./Allrun
```

### Expected results
- Mean Age of Air: ~150-200 seconds
- Ventilation efficiency: 0.5-0.7
- Fresh air enters from ceiling, old air exits at floor

### Parallel execution
```bash
# Decompose domain
decomposePar

# Run parallel
mpirun -np 4 AoAFoam -parallel

# Reconstruct
reconstructPar
```