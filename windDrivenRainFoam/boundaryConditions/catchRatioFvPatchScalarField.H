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

Class
    Foam::catchRatioFvPatchScalarField

Description
    Catch ratio boundary condition for building surfaces.
    
    Calculates the catch ratio based on:
    - Rain flux through the surface
    - Horizontal rainfall intensity
    - Surface orientation

Usage
    \table
        Property     | Description                    | Required | Default
        Rh          | Horizontal rainfall intensity  | yes      |
        phi         | Flux field name               | no       | phi
    \endtable

    Example usage:
    \verbatim
    building_wall
    {
        type            catchRatio;
        Rh              10;  // mm/h
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    catchRatioFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef catchRatioFvPatchScalarField_H
#define catchRatioFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class catchRatioFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class catchRatioFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Horizontal rainfall intensity [mm/h]
        scalar Rh_;

        //- Name of flux field
        word phiName_;


public:

    //- Runtime type information
    TypeName("catchRatio");


    // Constructors

        //- Construct from patch and internal field
        catchRatioFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        catchRatioFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given catchRatioFvPatchScalarField
        //  onto a new patch
        catchRatioFvPatchScalarField
        (
            const catchRatioFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        catchRatioFvPatchScalarField
        (
            const catchRatioFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new catchRatioFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        catchRatioFvPatchScalarField
        (
            const catchRatioFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new catchRatioFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //