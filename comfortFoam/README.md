# comfortFoam

`comfortFoam` is an OpenFOAM post-processing utility for assessing moderate
thermal environments according to ISO 7730:2025. It calculates predicted mean
vote (PMV), predicted percentage dissatisfied (PPD), draught rate (DR), and
operative temperature from CFD results. It can also evaluate the informative
thermal-environment categories I to IV from Annex A when the required local
discomfort inputs are available. This implementation is version 5.0.

## Scope and interpretation

PMV and PPD describe whole-body thermal sensation and dissatisfaction. DR and
the optional local criteria describe local discomfort caused by draught,
vertical air-temperature difference, floor temperature, and radiant-temperature
asymmetry.

ISO 7730:2025 Annex A is informative. A category reported by `comfortFoam` is
therefore an assessment result, not a certification. A volume-weighted regional
average can also hide locally uncomfortable cells; inspect the written fields,
validity information, and local-criterion results as well as the summary.

## Build

```bash
cd $FOAM_RUN/../HVAC-for-OpenFOAM/comfortFoam
wmake
```

## Run

Run the utility in an OpenFOAM case:

```bash
comfortFoam
```

Select a region or time in the usual way:

```bash
comfortFoam -cellZone occupiedZone
comfortFoam -cellSet occupiedCells
comfortFoam -region indoorAir -latestTime
mpirun -np 4 comfortFoam -parallel
```

The region can alternatively be selected with `cellSet`, its legacy alias
`setFields`, or `cellZone` in `constant/comfortFoamDict`. Specify at most one.
Command-line region selection takes precedence over the dictionary.

### Relevant command-line options

| Option | Description |
| --- | --- |
| `-region <name>` | Select an OpenFOAM mesh region. |
| `-cellSet <name>` | Assess only the named cell set. |
| `-setFields <name>` | Alias for `-cellSet`. |
| `-cellZone <name>` | Assess only the named cell zone. |
| `-validate` | Run the ISO 7730:2025 Annex D regression suite. |

Standard OpenFOAM time-selection and parallel options are also supported. The
utility writes its result fields; it does not provide a `-noWrite` option.

## Required and optional fields

The selected time directory must contain:

- `T`: air temperature [K]
- `U`: air velocity [m/s]

A pressure field is not required for the comfort calculation.

### Relative humidity

The first available relative-humidity field is used in this order:

1. `thermo:relHum`
2. `thermoRelHum`
3. `relHum`
4. `RH`
5. `relativeHumidity`

If none is present, the scalar `RH` from `constant/comfortFoamDict` is used.
Relative humidity is expressed in percent. For the PMV applicability check it is
converted to water-vapour partial pressure.

Alternatively, provide water-vapour partial pressure directly [Pa] as a
`waterVapourPressure` field or as a non-negative dictionary scalar. The field
overrides the scalar cell by cell. With the default scalar value `-1`, relative
humidity is used to calculate the pressure.

### Turbulence and draught

DR uses the local mean air speed and turbulence intensity. If `k` is available,
the turbulence intensity is derived from the CFD result; `epsilon` or `omega`
may also be read to identify the turbulence model. If `k` is absent, the
utility assumes `Tu = 40 %`.

In the ISO DR equation, `Tu` is a percentage, not a fraction: for example, a
turbulence intensity of 0.4 is supplied to the equation as `40`. The DR model is
applicable only when all of these conditions hold:

- air temperature: 20 degrees C to 26 degrees C
- local mean air speed: less than 0.5 m/s
- turbulence intensity: 10 % to 60 %
- metabolic rate: not greater than 1.2 met

For the equation, speeds below 0.05 m/s are evaluated as 0.05 m/s. DR is limited
to the range 0 % to 100 %. The model represents people performing light, mainly
sedentary activity with nearly neutral whole-body thermal sensation and predicts
draught discomfort at the neck. DR always uses the local room-air speed
`mag(U)`; a motion-adjusted `relativeAirVelocity` override affects PMV and the
clothing correction, not DR.

### Mean radiant temperature

An optional `MRT` field supplies mean radiant temperature directly [K] and has
precedence over a derived radiation estimate. Otherwise, `comfortFoam` can
derive the radiant temperature from supported radiation fields (`G`, `qr`, or
`IDefault`) and ultimately fall back to an area-weighted wall-temperature
estimate. Check the run log to see which source was selected; the fallback is
only an approximation of the occupant's angle-factor-weighted mean radiant
temperature.

### Optional local-discomfort fields

The following optional `volScalarField` names override their scalar dictionary
counterparts cell by cell:

- `verticalAirTemperatureDifference`: air-temperature difference between head
  and ankle level [K], conventionally at 1.1 m and 0.1 m above the floor
- `floorTemperature`: absolute floor surface temperature [K]
- `radiantTemperatureAsymmetry`: radiant-temperature asymmetry [K]

An optional `relativeAirVelocity` field [m/s] similarly overrides the constant
dictionary value. Field temperatures follow the usual OpenFOAM convention:
absolute temperatures are in kelvin, while temperature differences are in
kelvin and have the same numerical value as differences in degrees Celsius. In
the dictionary, `floorTemperature` is specified in degrees Celsius.

`radiantAsymmetryType` in the dictionary selects how the asymmetry is assessed:
`warmCeiling`, `coolWall`, `coolCeiling`, or `warmWall`. A negative dictionary
value disables the corresponding local criterion when no field is present.

The vertical-difference PD model applies to a rising head-to-ankle difference
below 8 K. The radiant-asymmetry PD model applies below 23 K for a warm ceiling,
15 K for a cool wall or cool ceiling, and 35 K for a warm wall. Inputs outside
these model ranges leave the corresponding local criterion unevaluated. The
floor-temperature model assumes light indoor footwear and is not applicable to
long exposures over electrically heated floors.

## Configuration

Create `constant/comfortFoamDict`, for example:

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      comfortFoamDict;
}

clo     1.0;    // Clothing insulation [clo]
met     1.2;    // Metabolic rate [met]
wme     0;      // External work [met]
RH      50;     // Fallback relative humidity [%]

// Motion and clothing correction
applyDynamicClothing       true;
walkingSpeed               -1;     // [m/s]; -1 = derive from met
clothingAlreadyResultant   false;
relativeAirVelocity        -1;     // [m/s]; -1 = derive from U and motion
waterVapourPressure        -1;     // [Pa]; -1 = calculate from RH

// Optional local-discomfort inputs; negative = unavailable
verticalAirTemperatureDifference  -1;  // [K]
floorTemperature                  -1;  // [degrees C]
radiantTemperatureAsymmetry       -1;  // [K]
radiantAsymmetryType              warmCeiling;

// Optional region selection; choose at most one
// cellSet   occupiedCells;
// cellZone  occupiedZone;
```

### Relative air velocity and dynamic clothing

ISO 7730 uses air velocity relative to the person (`v_ar`) and resultant
clothing insulation (`I_cl,r`) for PMV. The following settings control the
motion correction:

- `applyDynamicClothing true` enables the ISO 7730:2025 Annex C / ISO 9920
  clothing correction. In this mode `clo` is the static clothing insulation.
  `walkingSpeed` is interpreted in m/s; an explicit value is limited to
  1.2 m/s. A value of `-1` derives the equivalent walking speed from metabolic
  rate as `max(0, min(0.7, 0.0052*(M - 58)))`, with `M` in W/m2.
- `clothingAlreadyResultant true` declares that `clo` already represents
  resultant/dynamic insulation and prevents a second clothing correction.
- Exactly one of `applyDynamicClothing` and `clothingAlreadyResultant` must be
  true for an ISO-valid clothing input. Setting both to false preserves legacy
  behavior but is not marked ISO-valid. Setting both to true is an input
  error because it would request a second clothing correction.
- For backward compatibility, an omitted `applyDynamicClothing` entry defaults
  to `false`; `clothingAlreadyResultant` also defaults to `false`. Existing
  dictionaries therefore still run, but must explicitly select one clothing
  mode before their result is marked ISO-valid. The shipped example selects
  dynamic correction explicitly.
- `relativeAirVelocity` values greater than or equal to zero provide a constant
  override in m/s. It is the final `v_ar`, already including person motion;
  `walkingSpeed` is not added to an explicit override. A same-named
  `volScalarField` takes precedence cell by cell. With the default value `-1`,
  the utility uses the Annex C scalar approximation
  `v_ar = mag(U) + walkingSpeed`. Walking speed is zero unless supplied
  explicitly or derived while dynamic clothing correction is enabled. This is a
  conservative scalar assumption; it does not resolve the person's walking
  direction relative to the CFD velocity vector.

For the clothing correction, `v_ar` is limited to 0.15 m/s to 3.5 m/s. PMV
itself uses the actual, unlimited `v_ar`, and its ISO applicability limit remains
1 m/s. The normal/light clothing correction covers `0.6 < I_cl < 1.4 clo`; the
implementation interpolates the ISO 9920 clothed and nude endpoints from 0 clo
through 0.6 clo. Above 1.4 clo, use the applicable ISO 9920 cold-weather model
outside `comfortFoam` and provide the result with
`clothingAlreadyResultant true`.

The supplied converted standard document prints `0.44` as the quadratic
coefficient in Annex C.1 and presents C.4 in a dimensionally inconsistent form.
`comfortFoam` intentionally follows ISO 9920:2007, Equations (32) to (34): it
uses `0.044` and interpolates the endpoint total-insulation ratios. This avoids
encoding those apparent C.1/C.4 defects into the calculation.

The clothing values in the Annex D regression cases are treated as the values
required directly by the reference program; the optional dynamic correction is
not applied to those cases.

## ISO 7730 PMV applicability

The PMV model is intended for input combinations within all of the following
ranges:

| Quantity | Applicability range |
| --- | ---: |
| Metabolic rate | 0.8 met to 4.0 met (46 W/m2 to 232 W/m2) |
| Clothing insulation | 0 clo to 2 clo (0 m2 K/W to 0.310 m2 K/W) |
| Air temperature | 10 degrees C to 30 degrees C |
| Mean radiant temperature | 10 degrees C to 40 degrees C |
| Relative air velocity | 0 m/s to 1 m/s |
| Water-vapour partial pressure | 0 Pa to 2700 Pa |
| Calculated PMV | -2 to +2 |

Relative humidity itself is accepted from 0 % to 100 %, but ISO applicability
is determined by the resulting water-vapour partial pressure. Results outside
the applicability ranges must not be interpreted as ISO 7730 PMV results.

## Informative categories I to IV

The whole-body and draught limits from ISO 7730:2025 Annex A are:

| Category | PMV | PPD | DR |
| --- | --- | ---: | ---: |
| I | -0.2 < PMV < +0.2 | < 6 % | < 10 % |
| II | -0.5 < PMV < +0.5 | < 10 % | < 20 % |
| III | -0.7 < PMV < +0.7 | < 15 % | < 30 % |
| IV | -1.0 < PMV < +1.0 | < 25 % | no local limit |

For categories I to III, the local criteria are evaluated as well:

| Criterion | I | II | III |
| --- | ---: | ---: | ---: |
| Vertical air-temperature difference | < 2 K | < 3 K | < 4 K |
| Floor temperature | 19 to 29 degrees C | 19 to 29 degrees C | 17 to 31 degrees C |
| Warm-ceiling asymmetry | < 5 K | < 5 K | < 7 K |
| Cool-wall asymmetry | < 10 K | < 10 K | < 13 K |
| Cool-ceiling asymmetry | < 14 K | < 14 K | < 18 K |
| Warm-wall asymmetry | < 23 K | < 23 K | < 35 K |

A complete category I, II, or III can be determined only when DR is applicable
and all three additional local inputs can be evaluated. If one is missing or
outside its model's range, `comfortFoam` reports only the whole-body result and
the criteria that could be evaluated; it must not present that partial result as
a complete ISO category. Category IV has no local limits in Annex A.

## Output

For every processed time, `comfortFoam` writes these ten fields:

- `PMV`: predicted mean vote [dimensionless]
- `PPD`: predicted percentage dissatisfied [%]
- `DR`: draught rate [%]
- `TOp`: operative temperature [K]
- `PDVertical`: percentage dissatisfied due to the vertical air-temperature
  difference [%]
- `PDFloor`: percentage dissatisfied due to floor temperature [%]
- `PDRadiantAsymmetry`: percentage dissatisfied due to radiant-temperature
  asymmetry [%]
- `ISO7730Valid`: PMV applicability mask (`1` = valid, `0` = invalid,
  `-1` = outside the selected region)
- `ISO7730WholeBodyCategory`: category based on PMV and PPD only
- `ISO7730Category`: complete category including the local criteria

The two category fields use these dimensionless codes:

| Code | Meaning |
| ---: | --- |
| `1` | Category I |
| `2` | Category II |
| `3` | Category III |
| `4` | Category IV |
| `0` | No category |
| `-1` | Full category incomplete because a required local criterion cannot be evaluated, or outside the selected region |

The three local PD fields use `-1` where their criterion cannot be evaluated or
the cell is outside the selected region. `ISO7730WholeBodyCategory` remains
useful when `ISO7730Category` is incomplete, but it is not a complete category I
to III. In the whole-body field, `-1` only denotes a cell outside the selected
region.

The terminal summary identifies the analysis region, input sources,
volume-weighted results, applicability, and category completeness. Use the
cell-level fields when checking an occupied zone; the regional averages alone
are not a compliance assessment.

## Annex D validation

Run the built-in regression suite:

```bash
comfortFoam -validate
```

The repository helper runs the same check from its validation case and
propagates the utility's exit status:

```bash
cd comfortFoam
./validation/Allrun
```

The suite evaluates all 13 ISO 7730:2025 Annex D reference cases through the
same PMV implementation used during post-processing. It checks the calculated
PMV against the tabulated reference values with a tolerance of `0.001` PMV and
`0.05` percentage points PPD, and returns a non-zero exit status if a case
fails. The literal `273.0` offset from the Annex D program remains internal to
that calculation; physical OpenFOAM `T`, `MRT`, and floor-temperature fields
are converted from kelvin with `273.15`. This is a numerical regression test of
the normative reference implementation; it does not validate the environmental
inputs of a particular CFD case.

## Limitations

- PMV/PPD describes steady or lightly varying conditions and assumes thermal
  equilibrium. Time-weighted assessment of varying conditions is not generated
  automatically from unrelated time directories.
- A CFD cell velocity is not necessarily the velocity at the occupied location
  or height required by the standard. Define the analysis region accordingly.
- Radiation-field and wall-temperature fallbacks are approximations of the
  occupant's mean radiant temperature.
- The categories in Annex A are informative. Building-level acceptance may
  require additional project, national, or contractual criteria.

## References

1. ISO 7730:2025, *Ergonomics of the thermal environment - Analytical
   determination and interpretation of thermal comfort using calculation of the
   PMV and PPD indices and local thermal comfort criteria*.
2. DIN EN ISO 7730:2025-12, German adoption of EN ISO 7730:2025.
3. ISO 9920:2007, *Ergonomics of the thermal environment - Estimation of thermal
   insulation and water vapour resistance of a clothing ensemble*.
4. Fanger, P. O. (1970), *Thermal Comfort: Analysis and Applications in
   Environmental Engineering*.
