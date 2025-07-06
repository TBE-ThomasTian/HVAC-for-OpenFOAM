/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "rainPhase.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcAverage.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::volScalarField Foam::rainPhase::calculateCdRe(const volScalarField& Re)
{
    volScalarField CdRe(Re);
    
    forAll(Re, celli)
    {
        scalar Re_val = Re[celli];
        scalar CdRe_val;
        
        if (Re_val < 0.1)
        {
            CdRe_val = 24.0;
        }
        else if (Re_val < 2.0)
        {
            CdRe_val = 24.0 * (1.0 + 0.1315*pow(Re_val, 0.82 - 0.05*log10(Re_val)));
        }
        else if (Re_val < 20.0)
        {
            CdRe_val = 24.0 * (1.0 + 0.1935*pow(Re_val, 0.6305));
        }
        else if (Re_val < 260.0)
        {
            CdRe_val = 10.0 * pow(Re_val, 0.45);
        }
        else
        {
            CdRe_val = 0.44 * Re_val;
        }
        
        CdRe[celli] = CdRe_val;
    }
    
    return CdRe;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rainPhase::rainPhase
(
    const label index,
    const scalar diameter,
    const scalar volumeFraction,
    const fvMesh& mesh,
    const word& phaseName
)
:
    index_(index),
    diameter_(diameter),
    volumeFraction_(volumeFraction),
    mesh_(mesh),
    Urain_
    (
        IOobject
        (
            "U" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    phirain_
    (
        IOobject
        (
            "phi" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux(Urain_)
    ),
    alpharain_
    (
        IOobject
        (
            "alpha" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Ctrain_
    (
        IOobject
        (
            "Ct" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Ct", dimless, 0)
    ),
    Re_
    (
        IOobject
        (
            "Re" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Re", dimless, 500)
    ),
    CdRe_
    (
        IOobject
        (
            "CdRe" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("CdRe", dimless, 250)
    ),
    scr_
    (
        IOobject
        (
            "scr" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("scr", dimless, 0)
    ),
    tp_
    (
        IOobject
        (
            "tp" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("tp", dimTime, 0)
    ),
    tfl_
    (
        IOobject
        (
            "tfl" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("tfl", dimTime, 0)
    ),
    nutrain_
    (
        IOobject
        (
            "nutrain" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("nutrain", dimensionSet(0,2,-1,0,0,0,0), 1)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rainPhase::solveAlpha()
{
    // Create temporary phi field for boundary conditions
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phirain_
    );
    
    // Volume fraction equation
    fvScalarMatrix alphaEqn
    (
        fvm::ddt(alpharain_)
      + fvm::div(phirain_, alpharain_, "div(phirain,alpharain)")
      - fvm::Sp(fvc::div(phirain_), alpharain_)
    );
    
    alphaEqn.relax();
    alphaEqn.solve();
    
    alpharain_.correctBoundaryConditions();
}


void Foam::rainPhase::solveMomentum
(
    const volVectorField& U,
    const volVectorField& g,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& rhop,
    const dimensionedScalar& mua,
    const bool solveTD,
    const autoPtr<incompressible::turbulenceModel>& turbulence
)
{
    // Update properties first
    updateProperties(U, rhoa, mua);
    
    // Turbulent dispersion
    if (solveTD && turbulence.valid())
    {
        tmp<volScalarField> tk = turbulence->k();
        tmp<volScalarField> tepsilon = turbulence->epsilon();
        tfl_ = 0.2*(tk()/tepsilon());
        tp_ = (4*rhop*diameter_*diameter_)/(3*mua*CdRe_);
        Ctrain_ = sqrt(tfl_/(tfl_+tp_));
        nutrain_ = turbulence->nut()()*sqr(Ctrain_);
    }
    else
    {
        nutrain_ = dimensionedScalar("zero", nutrain_.dimensions(), 0.0);
    }
    
    // Reynolds stress for rain
    volTensorField Rca(-nutrain_*(T(fvc::grad(Urain_))));
    if (solveTD && turbulence.valid())
    {
        tmp<volScalarField> tk = turbulence->k();
        Rca = Rca + (2.0/3.0)*sqr(Ctrain_)*I*tk() - (2.0/3.0)*I*tr(Rca);
    }
    Rca.correctBoundaryConditions();
    
    surfaceScalarField phiRa
    (
        - fvc::interpolate(nutrain_)
        *mesh_.magSf()*fvc::snGrad(alpharain_)/fvc::interpolate(alpharain_ + scalar(0.001))
    );
    
    // Momentum equation
    fvVectorMatrix UrainEqn
    (
        fvm::ddt(Urain_)
      + fvm::div(phirain_, Urain_, "div(phirain,Urain)")
      - fvm::Sp(fvc::div(phirain_), Urain_)
      - fvm::laplacian(nutrain_, Urain_)
      + fvc::div(Rca)
      + fvm::div(phiRa, Urain_, "div(phirain,Urain)")
      - fvm::Sp(fvc::div(phiRa), Urain_)
      + (fvc::grad(alpharain_)/(fvc::average(alpharain_) + scalar(0.001)) & Rca)
     ==
        g
      + ((3*mua*CdRe_)/(4*rhop*diameter_*diameter_))*U
      - fvm::Sp(((3*mua*CdRe_)/(4*rhop*diameter_*diameter_)), Urain_)
    );
    
    UrainEqn.relax();
    UrainEqn.solve();
    
    Urain_.correctBoundaryConditions();
    phirain_ = fvc::flux(Urain_);
}


void Foam::rainPhase::calculateCatchRatio(const dimensionedScalar& Rh)
{
    tmp<surfaceScalarField> tnormalvel
    (
        new surfaceScalarField
        (
            IOobject
            (
                "normalvel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mag((mesh_.Sf()/mesh_.magSf()) & fvc::interpolate(Urain_))
        )
    );
    const surfaceScalarField& normalvel = tnormalvel();
    
    scr_ = (fvc::average(normalvel) * alpharain_) * ((3600*1E3)/(Rh*volumeFraction_));
}


void Foam::rainPhase::updateProperties
(
    const volVectorField& U,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& mua
)
{
    volVectorField Ur = U - Urain_;
    volScalarField magUr = mag(Ur);
    Re_ = (magUr*diameter_*rhoa)/mua;
    CdRe_ = calculateCdRe(Re_);
    CdRe_.correctBoundaryConditions();
}


Foam::scalar Foam::rainPhase::CourantNumber() const
{
    scalar CoNum = 0.5*gMax
    (
        fvc::surfaceSum(mag(phirain_))().primitiveField()
       /mesh_.V().field()
    )*mesh_.time().deltaTValue();
    
    return CoNum;
}


void Foam::rainPhase::write() const
{
    Urain_.write();
    alpharain_.write();
    scr_.write();
}


// ************************************************************************* //