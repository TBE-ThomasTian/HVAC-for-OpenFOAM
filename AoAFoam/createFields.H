Info<< "Reading field AoA\n" << endl;

volScalarField AoA
(
    IOobject
    (
        "AoA",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p\n" << endl;

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("p", dimPressure, 0)
);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading diffusivity DT\n" << endl;

dimensionedScalar DT
(
    transportProperties.lookupOrDefault<dimensionedScalar>
    (
        "DT",
        dimensionedScalar("DT", dimArea/dimTime, 0.0)
    )
);

#include "createPhi.H"

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createFvOptions.H"

// Check and print solver settings
Info<< "\nAge of Air solver settings:" << endl;
Info<< "    Molecular/laminar diffusivity DT = " << DT.value() << " m2/s" << endl;

// Validate initial conditions
if (min(AoA).value() < 0)
{
    FatalErrorInFunction
        << "Negative Age of Air values detected. AoA must be >= 0."
        << exit(FatalError);
}


// Find outlet patches for reference
wordList outletPatches;
forAll(mesh.boundary(), patchi)
{
    const fvPatch& patch = mesh.boundary()[patchi];
    const fvPatchScalarField& AoApatch = AoA.boundaryField()[patchi];
    
    if (AoApatch.type() == "zeroGradient" || AoApatch.type() == "inletOutlet")
    {
        outletPatches.append(patch.name());
    }
}

if (outletPatches.size() > 0)
{
    Info<< "    Detected outlet patches: " << outletPatches << endl;
}

Info<< endl;
