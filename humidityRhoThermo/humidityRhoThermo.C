/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "humidityRhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(humidityRhoThermo, 0);
    defineRunTimeSelectionTable(humidityRhoThermo, fvMesh);
    defineRunTimeSelectionTable(humidityRhoThermo, fvMeshDictPhase);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::humidityRhoThermo::humidityRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    fluidThermo(mesh, phaseName),
    rho_
    (
        IOobject
        (
            phasePropertyName("thermo:rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimDensity
    ),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    ),

    relHum_
    (
        IOobject
        (
            phasePropertyName("thermo:relHum"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    ),

    waterContent_
    (
        IOobject
        (
            phasePropertyName("thermo:waterContent"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    maxWaterContent_
    (
        IOobject
        (
            phasePropertyName("maxWaterContent"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    specificHumidity_
    (
        IOobject
        (
            phasePropertyName("thermo:specificHumidity"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    maxSpecificHumidity_
    (
        IOobject
        (
            phasePropertyName("thermo:maxSpecificHumidity"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    ),

    pSatH2O_
    (
        IOobject
        (
            phasePropertyName("thermo:pSatH2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimPressure
    ),

    partialPressureH2O_
    (
        IOobject
        (
            phasePropertyName("thermo:partialPressureH2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimPressure
    ),

    muEff_
    (
        IOobject
        (
            phasePropertyName("thermo:muEff"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1,-1,-1,0,0,0,0)
    )
{
    method_ = readMethod();
}


Foam::humidityRhoThermo::humidityRhoThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    fluidThermo(mesh, dict, phaseName),
    rho_
    (
        IOobject
        (
            phasePropertyName("thermo:rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimDensity
    ),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    ),

    relHum_
    (
        IOobject
        (
            phasePropertyName("thermo:relHum"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    ),

    waterContent_
    (
        IOobject
        (
            phasePropertyName("thermo:waterContent"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    maxWaterContent_
    (
        IOobject
        (
            phasePropertyName("thermo:maxWaterContent"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    specificHumidity_
    (
        IOobject
        (
            phasePropertyName("thermo:specificHumidity"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    ),

    maxSpecificHumidity_
    (
        IOobject
        (
            phasePropertyName("thermo:maxSpecificHumidity"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    ),

    pSatH2O_
    (
        IOobject
        (
            phasePropertyName("thermo:pSatH2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimPressure
    ),

    partialPressureH2O_
    (
        IOobject
        (
            phasePropertyName("thermo:partialPressureH2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimPressure
    ),

    muEff_
    (
        IOobject
        (
            phasePropertyName("thermo:muEff"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1,-1,-1,0,0,0,0)
    )
{
    method_ = readMethod();
}


Foam::humidityRhoThermo::humidityRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
    fluidThermo(mesh, phaseName, dictionaryName),
    rho_
    (
        IOobject
        (
            phasePropertyName("thermo:rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    ),

    relHum_
    (
        IOobject
        (
            phasePropertyName("thermo:relHum"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    ),

    waterContent_
    (
        IOobject
        (
            phasePropertyName("thermo:waterContent"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    maxWaterContent_
    (
        IOobject
        (
            phasePropertyName("maxWaterContent"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    specificHumidity_
    (
        IOobject
        (
            phasePropertyName("thermo:specificHumidity"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    maxSpecificHumidity_
    (
        IOobject
        (
            phasePropertyName("thermo:maxSpecificHumidity"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    ),

    pSatH2O_
    (
        IOobject
        (
            phasePropertyName("thermo:pSatH2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimPressure
    ),

    partialPressureH2O_
    (
        IOobject
        (
            phasePropertyName("thermo:partialPressureH2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimPressure
    ),

    muEff_
    (
        IOobject
        (
            phasePropertyName("thermo:muEff"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1,-1,-1,0,0,0,0)
    )
{
    method_ = readMethod();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::humidityRhoThermo> Foam::humidityRhoThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<humidityRhoThermo>(mesh, phaseName);
}


Foam::autoPtr<Foam::humidityRhoThermo> Foam::humidityRhoThermo::New
(
     const fvMesh& mesh,
     const word& phaseName,
     const word& dictName
)
{
    return basicThermo::New<humidityRhoThermo>(mesh, phaseName, dictName);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::humidityRhoThermo::~humidityRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::humidityRhoThermo::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::humidityRhoThermo::rho(const label patchi) const
{
    return rho_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::humidityRhoThermo::rho()
{
    return rho_;
}


void Foam::humidityRhoThermo::correctRho
(
    const Foam::volScalarField& deltaRho,
    const dimensionedScalar& rhoMin,
    const dimensionedScalar& rhoMax
)
{
    rho_ += deltaRho;
    rho_ = max(rho_, rhoMin);
    rho_ = min(rho_, rhoMax);
}

void Foam::humidityRhoThermo::correctRho(const Foam::volScalarField& deltaRho)
{
    rho_ += deltaRho;
}


const Foam::volScalarField& Foam::humidityRhoThermo::psi() const
{
    return psi_;
}


Foam::tmp<Foam::volScalarField> Foam::humidityRhoThermo::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField> Foam::humidityRhoThermo::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}


const Foam::word Foam::humidityRhoThermo::readMethod() const
{
    const wordList& bTsH = this->specificHumidity_.boundaryField().types();

    //- Search if fixedHumidity bc is used
    bool found = false;
    forAll(bTsH, patchI)
    {
        if (bTsH[patchI] == "fixedHumidity")
        {
            found = true;
            break;
        }
    }

    //- Default method
    word method = "buck";

    if (found)
    {
        if
        (
            this->specificHumidity_.mesh().foundObject<IOList<word>>
            (
                "methodName"
            )
        )
        {
            Info<< "inside " << endl;
            method =
                this->specificHumidity_.mesh().lookupObject<IOList<word>>
                (
                    "methodName"
                )[0];
        }
    }
    else
    {
        WarningInFunction
            << "No fixedHumidity boundary condition found. Using the default "
            << "method to calculate the saturation pressure\n" << endl;
    }

    Info<< "Saturation pressure calculation based on "
        << method << "\n" << endl;

    return method;
}


// ************************************************************************* //
