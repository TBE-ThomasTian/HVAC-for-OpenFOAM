# OpenCFD OpenFOAM-v2606 NEN 8100 port

This self-contained port builds `NEN8100Foam` against the OpenFOAM-v2606 source
tree without modifying it.

```sh
source /path/to/OpenFOAM-v2606/etc/bashrc
cd /path/to/HVAC-for-OpenFOAM/NEN8100/OpenFOAM-v2606
wmake -j
```

Unlike the v2512 source in the parent directory, the build writes the executable
to `FOAM_APPBIN` rather than `FOAM_USER_APPBIN`.

The assessment model is unchanged: velocity amplification per direction, a
directional Weibull wind climate, NEN comfort probability at 5 m/s with classes
A to E, and danger probability at 15 m/s. The dictionary layout of
`NEN8100Dict` is the same, and the case arrangement, pedestrian mask, climate
input, and output fields are documented in the [parent README](../README.md).

## Differences from the v2512 source

**Stricter mesh consistency.** The reference mesh is captured as a signature of
point, face, internal-face, and cell counts together with cell centres and cell
volumes. Every subsequent directional case is checked against it, and a
mismatch in counts, in cell centres beyond `geometryTolerance`, or in cell
volumes beyond the new `relativeVolumeTolerance` is a fatal error naming the
offending case and the measured deviation. The v2512 source compares cell
centres only.

**Non-finite value guards.** Dictionary inputs, Weibull exceedance values,
cell-centre velocities, and every field about to be written are tested for
`NaN` and infinity, and the run aborts at the point of origin instead of
writing corrupt fields.

**Explicit calm fraction.** The climate dictionary accepts `calmFraction`, the
annual fraction excluded from the directional Weibull fits. Sector frequencies
plus `calmFraction` must sum to one within the new `frequencyTolerance`
(default `1e-8`); the v2512 source requires the sector frequencies alone to sum
to one.

**Duplicate direction detection.** Two directions that would produce the same
amplification field name are rejected at read time.

**`strictDirectionCount` defaults to true.** A direction count differing from
`expectedDirectionCount` is a fatal error rather than a warning.

## Additional dictionary entries

```
climate
{
    frequencyTolerance      1e-8;   // tolerance on the annual frequency sum
    calmFraction            0.05;   // annual fraction excluded from the fits
}

relativeVolumeTolerance     1e-10;  // max relative cell-volume difference
strictDirectionCount        true;   // now fatal rather than a warning
```
