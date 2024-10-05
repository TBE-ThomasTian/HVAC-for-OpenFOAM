#include "IOobject.H"
#include "volFields.H"

// Temperaturfeld (T)
IOobject THeader
(
    "T",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);

volScalarField T
(
    THeader,
    mesh
);

// Geschwindigkeitsfeld (U)
IOobject UHeader
(
    "U",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);

volVectorField U
(
    UHeader,
    mesh
);

// Relative Luftfeuchtigkeit (RH)
IOobject RHHeader
(
    "relHum",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);

volScalarField RH
(
    RHHeader,
    mesh
);

// Strahlungsflussdichte (qr)
IOobject qrHeader
(
    "qr",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);

volScalarField qr
(
    qrHeader,
    mesh
);

// UTCI-Feld zur Speicherung der Ergebnisse
volScalarField UTCI
(
    IOobject
    (
        "UTCI",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("UTCI", dimensionSet(0,0,0,1,0,0,0), 0.0)
);
