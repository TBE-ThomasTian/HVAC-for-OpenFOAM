# NEN8100Foam

`NEN8100Foam` is a serial OpenFOAM post-processing utility for pedestrian
wind-comfort assessment using the public criteria of NEN 8100:2006. It reads
multiple directional CFD cases, calculates full-volume velocity amplification
and writes ordinary OpenFOAM `volScalarField` results for direct use in
ParaView.

The utility does not validate the CFD model, atmospheric boundary layer, mesh,
or meteorological source. Its output is an engineering assessment and is not by
itself a certificate of standards compliance.

## Current scope

- multiple steady directional OpenFOAM cases;
- meteorological `windFrom` convention;
- directional Weibull wind climate;
- NEN comfort probability at 5 m/s and classes A to E;
- NEN danger probability at 15 m/s;
- standard OpenFOAM fields in a directional result case;
- configurable pedestrian-height mask for ParaView;
- optional amplification fields for every direction.

All directional cases currently need the same cell count, cell order, and cell
coordinates. The tool checks this before combining values. Terrain-following
pedestrian surfaces and mesh-to-mesh mapping are planned extensions.

## Build

The current source targets OpenCFD OpenFOAM-v2512:

```bash
source /path/to/OpenFOAM-v2512/etc/bashrc
cd /path/to/HVAC-for-OpenFOAM/NEN8100
wmake -j
```

The executable is written to `$FOAM_USER_APPBIN/NEN8100Foam`.

A port for OpenCFD OpenFOAM-v2606 with stricter mesh and input checking is kept
in [OpenFOAM-v2606](OpenFOAM-v2606/README.md).

## Case arrangement

Run the utility from an OpenFOAM driver case containing
`system/NEN8100Dict`. The driver may also be one of the directional cases.
Paths in `directions` are resolved relative to the driver:

```text
project/
├── analysis/system/NEN8100Dict
├── wind_000/{constant,system,1000/U}
├── wind_030/{constant,system,1000/U}
└── ...
```

Copy the bundled [`NEN8100Dict`](NEN8100Dict) into the driver case and adjust
the case paths, reference speeds, pedestrian plane, and climate data.

```bash
cd project/analysis
NEN8100Foam
```

Every directional entry contains:

```c
{
    windFrom        270;            // meteorological direction [degrees]
    case            "../wind_270";
    Uref            10;             // mean CFD reference speed [m/s]
    referenceHeight 60;             // CFD reference height [m]
    time            latest;         // optional, default: latest
    // region        region0;        // optional
}
```

Only one CFD speed is needed per direction. The utility uses
`gamma = mag(Ulocal)/Uref` in every cell and applies the full climate
distribution during post-processing.

## Pedestrian plane and mask

NEN fields are calculated throughout the common volume mesh. The configured
pedestrian plane is additionally written as the field
`NEN8100PedestrianMask`: cells cut by the plane have value `1`, all other cells
have value `0`.

For flat ground at `z=0`:

```c
evaluationSurface
{
    type                horizontalPlane;
    groundPoint         (0 0 0);
    verticalDirection   (0 0 1);
    pedestrianHeight    1.75;
    triangulate         true;

    // Optional restriction:
    // bounds (-50 -50 0) (50 50 5);
    // zone pedestrianZone;
}
```

The plane point is calculated as:

```text
groundPoint + pedestrianHeight * unit(verticalDirection)
```

In ParaView either create an exact geometric `Slice` at 1.75 m or apply a
`Threshold` to `NEN8100PedestrianMask` with value `1` to inspect the intersected
cell layer.

## Climate input

The initial climate input is a directional Weibull distribution:

```c
climate
{
    type                directionalWeibull;
    referenceHeight     60;
    directionTolerance  0.1;

    sectors
    (
        { windFrom 0;  frequency 0.08; scale 6.2; shape 2.1; }
        // one sector for every CFD direction
    );
}
```

`frequency` is the joint fraction of all annual hours belonging to the sector,
not a percentage and not a conditional fraction. Its total may be below one;
the remainder is treated as calm/non-exceeding time. `scale` is in m/s and
`shape` is dimensionless. Climate and CFD reference heights must already agree.

For Dutch projects the distribution can be derived from NPR 6097 data. For
other locations it can be fitted to suitable local long-term hourly wind data.
A later reader can add raw hourly or binned climate input without changing the
directional CFD interface.

## Output fields

By default the fields are written into the selected time directory of the first
case in `directions`. `outputCase`, `outputTime`, and `outputRegion` can override
the destination, but its mesh must match the directional cases:

```c
fieldPrefix            NEN8100;
// outputCase          "../wind_000";
outputTime             latest;
// outputRegion        region0;
writeAmplification     false;
```

| OpenFOAM field | Meaning |
| --- | --- |
| `NEN8100PedestrianMask` | 1 for cells cut by the configured pedestrian plane |
| `NEN8100P5Percent` | Annual probability of hourly mean local speed above 5 m/s [%] |
| `NEN8100ComfortClass` | 1=A, 2=B, 3=C, 4=D, 5=E |
| `NEN8100TraversingRating` | 0=good, 1=moderate, 2=poor |
| `NEN8100StrollingRating` | 0=good, 1=moderate, 2=poor |
| `NEN8100SittingRating` | 0=good, 1=moderate, 2=poor |
| `NEN8100P15Percent` | Annual probability of hourly mean local speed above 15 m/s [%] |
| `NEN8100DangerClass` | 0=no danger, 1=limited risk, 2=dangerous |
| `NEN8100Gamma_<direction>` | Optional local amplification per CFD direction |

The comfort A/B/C/D/E boundaries for `P(U > 5 m/s)` are 2.5, 5, 10, and
20 percent. For `P(U > 15 m/s)`, values below 0.05 percent mean no danger,
values below 0.30 percent mean limited risk, and values from 0.30 percent are
classified as dangerous.

