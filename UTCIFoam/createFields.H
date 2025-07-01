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
// Try to read thermo:relHum first (from buoyantHumiditySimpleFoam)
IOobject RHHeader
(
    "relHum",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT
);

volScalarField RH
(
    RHHeader,
    mesh
);

// Create radiation model
autoPtr<radiation::radiationModel> radiation
(
    radiation::radiationModel::New(T)
);

// Strahlungsflussdichte (qr) - now obtained from radiation model
volScalarField qr = radiation->Rp();

// Sky View Factor (SVF) - optional field
volScalarField SVF
(
    IOobject
    (
        "SVF",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("SVF", dimless, 0.7)  // Default value 0.7 for partially obstructed outdoor environments
);

// Mean Radiant Temperature field - optional, can be provided directly
volScalarField Tmrt
(
    IOobject
    (
        "Tmrt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Tmrt", dimensionSet(0,0,0,1,0,0,0), 0.0)
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
