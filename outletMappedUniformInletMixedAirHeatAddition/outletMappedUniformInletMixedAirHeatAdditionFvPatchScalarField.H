/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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
    Foam::outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    Inlet temperature boundary condition that maps an outlet temperature,
    mixes it with an additional air mass flow at a prescribed temperature,
    and applies an optional heat addition.

    The imposed inlet temperature is

    \f[
        T_{in} =
        \frac
        {
            \dot{m}_{recirc} C_{p,recirc} T_{out}
          + \dot{m}_{add} C_{p,add} T_{add}
          + Q
        }
        {
            \dot{m}_{recirc} C_{p,recirc}
          + \dot{m}_{add} C_{p,add}
        }
    \f]

    For a mass-based phi, the outlet mass flow is sum(phi). For a
    volumetric phi, the outlet mass flow is sum(rho*phi) using the local
    thermodynamic density.

Usage
    \table
        Property     | Description                    | Required | Default
        outletPatch  | name of outlet patch to map    | yes      |
        additionalMassFlowRate | added air mass flow [kg/s] | yes |
        additionalTemperature  | added air temperature [K]   | yes |
        Q            | heat addition [W]              | no       | 0
        phi          | flux field name                | no       | phi
        additionalCp | added air Cp [J/kg/K]          | no       | thermo Cp
        TMin         | min temperature limit [K]      | no       | 0
        TMax         | max temperature limit [K]      | no       | 5000
        debug        | print BC diagnostics           | no       | false
    \endtable

    The entries Q, additionalMassFlowRate and additionalTemperature are
    Function1<scalar> values evaluated as functions of time. The aliases
    freshMassFlowRate and freshTemperature are also accepted.

    Example:
    \verbatim
    inlet
    {
        type                    outletMappedUniformInletMixedAirHeatAddition;
        outletPatch             outlet;
        phi                     phi;
        Q                       1500;
        additionalMassFlowRate  0.1;
        additionalTemperature   293.15;
        value                   uniform 293.15;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField
    Foam::outletMappedUniformInletHeatAdditionFvPatchField

SourceFiles
    outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField_H
#define outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
  Class outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField
  Declaration
\*---------------------------------------------------------------------------*/

class outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Optional heat addition
        autoPtr<Function1<scalar>> Qptr_;

        //- Additional air mass flow rate [kg/s]
        autoPtr<Function1<scalar>> additionalMassFlowRatePtr_;

        //- Additional air temperature [K]
        autoPtr<Function1<scalar>> additionalTemperaturePtr_;

        //- Name of the outlet patch to be mapped
        word outletPatchName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Optional fixed Cp for the additional air [J/kg/K]
        scalar additionalCp_;

        //- Minimum temperature limit
        scalar TMin_;

        //- Maximum temperature limit
        scalar TMax_;

        //- Debug output switch
        Switch debug_;

        //- Last time index for debug output
        label debugTimeIndex_;


public:

    //- Runtime type information
    TypeName("outletMappedUniformInletMixedAirHeatAddition");


    // Constructors

        //- Construct from patch and internal field
        outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField
        (
            const outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField
        (
            const outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField
        (
            const outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member functions

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
