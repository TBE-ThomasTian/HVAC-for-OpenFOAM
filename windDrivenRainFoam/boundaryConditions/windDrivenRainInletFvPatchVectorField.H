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
    Foam::windDrivenRainInletFvPatchVectorField

Description
    Wind-driven rain inlet boundary condition for rain velocity.
    
    The rain velocity at the inlet is calculated based on:
    - Wind velocity at the boundary
    - Terminal velocity of raindrops
    - Wind angle relative to the boundary normal

Usage
    \table
        Property     | Description                         | Required | Default
        Uwind        | Wind velocity field name            | no       | U
        terminalVelocity | Terminal velocity magnitude     | yes      |
        windAngle    | Wind angle relative to normal [deg] | no       | 0
    \endtable

    Example usage:
    \verbatim
    inlet
    {
        type            windDrivenRainInlet;
        terminalVelocity 5.0;  // m/s
        windAngle       30;    // degrees
        value           uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    windDrivenRainInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef windDrivenRainInletFvPatchVectorField_H
#define windDrivenRainInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class windDrivenRainInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class windDrivenRainInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Name of the wind velocity field
        word UwindName_;

        //- Terminal velocity magnitude
        scalar terminalVelocity_;

        //- Wind angle relative to boundary normal (degrees)
        scalar windAngle_;


public:

    //- Runtime type information
    TypeName("windDrivenRainInlet");


    // Constructors

        //- Construct from patch and internal field
        windDrivenRainInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        windDrivenRainInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given windDrivenRainInletFvPatchVectorField
        //  onto a new patch
        windDrivenRainInletFvPatchVectorField
        (
            const windDrivenRainInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        windDrivenRainInletFvPatchVectorField
        (
            const windDrivenRainInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new windDrivenRainInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        windDrivenRainInletFvPatchVectorField
        (
            const windDrivenRainInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new windDrivenRainInletFvPatchVectorField(*this, iF)
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