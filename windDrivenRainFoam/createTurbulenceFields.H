// Generic turbulence model support for windDrivenRainFoam

#include "turbulentTransportModel.H"
#include "singlePhaseTransportModel.H"

Info<< "Creating turbulence model\n" << endl;

// Create a dummy laminar transport model for turbulence initialization
singlePhaseTransportModel laminarTransport(U, phiwind);

// Create turbulence model selector
turbulence = incompressible::turbulenceModel::New(U, phiwind, laminarTransport);

// Scale turbulence fields if needed
if (scalingFactor != 1.0)
{
    volScalarField& k = const_cast<volScalarField&>(turbulence->k()());
    k *= sqr(scalingFactor);
    k.correctBoundaryConditions();
    
    volScalarField& epsilon = const_cast<volScalarField&>(turbulence->epsilon()());
    epsilon *= pow3(scalingFactor);
    epsilon.correctBoundaryConditions();
    
    volScalarField& nut = const_cast<volScalarField&>(turbulence->nut()());
    nut *= scalingFactor;
    nut.correctBoundaryConditions();
}
