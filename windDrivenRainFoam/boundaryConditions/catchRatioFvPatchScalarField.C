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

#include "catchRatioFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::catchRatioFvPatchScalarField::catchRatioFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Rh_(0),
    phiName_("phi")
{}


Foam::catchRatioFvPatchScalarField::catchRatioFvPatchScalarField
(
    const catchRatioFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Rh_(ptf.Rh_),
    phiName_(ptf.phiName_)
{}


Foam::catchRatioFvPatchScalarField::catchRatioFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Rh_(dict.get<scalar>("Rh")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


Foam::catchRatioFvPatchScalarField::catchRatioFvPatchScalarField
(
    const catchRatioFvPatchScalarField& crpsf
)
:
    fixedValueFvPatchScalarField(crpsf),
    Rh_(crpsf.Rh_),
    phiName_(crpsf.phiName_)
{}


Foam::catchRatioFvPatchScalarField::catchRatioFvPatchScalarField
(
    const catchRatioFvPatchScalarField& crpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(crpsf, iF),
    Rh_(crpsf.Rh_),
    phiName_(crpsf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::catchRatioFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the flux field
    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);
    
    const fvsPatchScalarField& phip = phi.boundaryField()[patch().index()];
    
    // Get patch areas
    const scalarField& magSf = patch().magSf();
    
    // Calculate catch ratio
    // CR = rain_flux / (Rh * area)
    // Convert Rh from mm/h to m/s: mm/h * 1e-3 / 3600
    scalar RhSI = Rh_ * 1e-3 / 3600.0;
    
    if (RhSI > SMALL)
    {
        scalarField catchRatio = mag(phip) / (RhSI * magSf);
        operator==(catchRatio);
    }
    else
    {
        operator==(0.0);
    }
    
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::catchRatioFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("Rh", Rh_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        catchRatioFvPatchScalarField
    );
}

// ************************************************************************* //