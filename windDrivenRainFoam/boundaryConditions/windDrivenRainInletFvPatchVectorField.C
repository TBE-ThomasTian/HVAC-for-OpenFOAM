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

#include "windDrivenRainInletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::windDrivenRainInletFvPatchVectorField::
windDrivenRainInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    UwindName_("U"),
    terminalVelocity_(0),
    windAngle_(0)
{}


Foam::windDrivenRainInletFvPatchVectorField::
windDrivenRainInletFvPatchVectorField
(
    const windDrivenRainInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    UwindName_(ptf.UwindName_),
    terminalVelocity_(ptf.terminalVelocity_),
    windAngle_(ptf.windAngle_)
{}


Foam::windDrivenRainInletFvPatchVectorField::
windDrivenRainInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    UwindName_(dict.lookupOrDefault<word>("Uwind", "U")),
    terminalVelocity_(dict.get<scalar>("terminalVelocity")),
    windAngle_(dict.lookupOrDefault<scalar>("windAngle", 0))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


Foam::windDrivenRainInletFvPatchVectorField::
windDrivenRainInletFvPatchVectorField
(
    const windDrivenRainInletFvPatchVectorField& wdripvf
)
:
    fixedValueFvPatchVectorField(wdripvf),
    UwindName_(wdripvf.UwindName_),
    terminalVelocity_(wdripvf.terminalVelocity_),
    windAngle_(wdripvf.windAngle_)
{}


Foam::windDrivenRainInletFvPatchVectorField::
windDrivenRainInletFvPatchVectorField
(
    const windDrivenRainInletFvPatchVectorField& wdripvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(wdripvf, iF),
    UwindName_(wdripvf.UwindName_),
    terminalVelocity_(wdripvf.terminalVelocity_),
    windAngle_(wdripvf.windAngle_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::windDrivenRainInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get wind velocity field
    const volVectorField& Uwind =
        db().lookupObject<volVectorField>(UwindName_);
    
    const fvPatchVectorField& Uwindp = Uwind.boundaryField()[patch().index()];
    
    // Get patch normal (pointing into domain)
    const vectorField n(patch().nf());
    
    // Calculate rain velocity
    // Horizontal component from wind
    vectorField Urain_horizontal = Uwindp - (Uwindp & n)*n;
    
    // Add vertical terminal velocity component
    vectorField Urain = Urain_horizontal - terminalVelocity_*vector(0, 0, 1);
    
    // Apply wind angle correction if specified
    if (mag(windAngle_) > SMALL)
    {
        scalar angleRad = windAngle_ * constant::mathematical::pi / 180.0;
        // Rotate around patch normal
        // This is simplified - a full implementation would need proper rotation
        Urain = cos(angleRad)*Urain + sin(angleRad)*(n ^ Urain);
    }
    
    // Set the rain velocity
    operator==(Urain);
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::windDrivenRainInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntryIfDifferent<word>("Uwind", "U", UwindName_);
    os.writeEntry("terminalVelocity", terminalVelocity_);
    os.writeEntryIfDifferent<scalar>("windAngle", 0, windAngle_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        windDrivenRainInletFvPatchVectorField
    );
}

// ************************************************************************* //