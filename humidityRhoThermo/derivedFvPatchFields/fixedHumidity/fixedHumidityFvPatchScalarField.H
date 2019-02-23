/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::fixedHumidityFvPatchScalarField

Description
    This boundary condition provides a fixed condition for humidity

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedHumidity;
        mode            relative;
        value           uniform 0.1;
        method          buck;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    fixedHumidityFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedHumidityFvPatchScalarField_H
#define fixedHumidityFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class fixedHumidityFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedHumidityFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

private:

    // Private data

        //- Mode used | Absolute | Relative
        word mode_;

        //- Calculation method | magnus | buck
        const word method_;

        //- Value of humidity
        scalar value_;

        //- Exchange value for thermodynamic database
        IOList<word> methodName_;

public:

    //- Runtime type information
    TypeName("fixedHumidity");


    // Constructors

        //- Construct from patch and internal field
        fixedHumidityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedHumidityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedHumidityFvPatchScalarField
        // onto a new patch
        fixedHumidityFvPatchScalarField
        (
            const fixedHumidityFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedHumidityFvPatchScalarField
        (
            const fixedHumidityFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedHumidityFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedHumidityFvPatchScalarField
        (
            const fixedHumidityFvPatchScalarField&,
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
                new fixedHumidityFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Calculate the specific humidity
            virtual const scalarField calcSpecificHumidity
            (
                const basicThermo&,
                const label
            );

    //- Write
    virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
