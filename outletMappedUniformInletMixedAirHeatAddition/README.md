# outletMappedUniformInletMixedAirHeatAddition

Runtime-selectable inlet temperature boundary condition for OpenFOAM that maps
an outlet temperature, mixes additional air with its own mass flow and
temperature, and applies optional heat addition.

## Build

Source your OpenFOAM environment first, then:

```bash
cd /home/flowsimpro/OpenFOAM/HVAC-for-OpenFOAM/outletMappedUniformInletMixedAirHeatAddition
wmake libso
```

The library is built to:

```text
$FOAM_USER_LIBBIN/liboutletMappedUniformInletMixedAirHeatAddition.so
```

Add to `system/controlDict`:

```foam
libs ("liboutletMappedUniformInletMixedAirHeatAddition.so");
```

## Usage

Example in `0/T`:

```foam
inlet
{
    type                    outletMappedUniformInletMixedAirHeatAddition;
    outletPatch             outlet;
    phi                     phi;

    Q                       1500;       // heat addition [W]
    additionalMassFlowRate  0.10;       // added air mass flow [kg/s]
    additionalTemperature   293.15;     // added air temperature [K]

    TMin                    250;
    TMax                    350;
    debug                   1;

    value                   uniform 293.15;
}
```

The aliases `freshMassFlowRate` and `freshTemperature` are also accepted.
`Q`, `additionalMassFlowRate`, and `additionalTemperature` can be constants or
OpenFOAM `Function1<scalar>` entries evaluated over time.

`TMin` and `TMax` limit the final imposed inlet temperature after outlet
mapping, additional-air mixing, and heat addition. Their defaults are:

```foam
TMin 0;
TMax 5000;
```

For the velocity field, set the inlet flow consistently with the total flow
rate:

```text
mDot_in = mDot_recirc + additionalMassFlowRate
```

## Formula

```text
T_in =
(
    mDot_recirc * Cp_recirc * T_out
  + mDot_add    * Cp_add    * T_add
  + Q
)
/
(
    mDot_recirc * Cp_recirc
  + mDot_add    * Cp_add
)
```

`Cp_add` is evaluated from the thermophysical model at
`additionalTemperature`. You can override it with:

```foam
additionalCp 1005;
```
