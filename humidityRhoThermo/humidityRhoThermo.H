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

Class
    Foam::humidityRhoThermo

Group
    grpRhoThermo

Description
    Basic thermodynamic properties based on density

SourceFiles
    humidityRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef humidityRhoThermo_H
#define humidityRhoThermo_H

#include "fluidThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class humidityRhoThermo Declaration
\*---------------------------------------------------------------------------*/

class humidityRhoThermo
:
    public fluidThermo
{

protected:

    // Protected data

        //- Density field [kg/m^3]
        //  Named 'humidityRhoThermo' to avoid (potential) conflict with solver density
        volScalarField rho_;

        //- Compressibility [s^2/m^2]
        volScalarField psi_;

        //- Dynamic viscosity [kg/m/s]
        volScalarField mu_;

        //- Relative humidity [-]
        volScalarField relHum_;

        //- Water content [kg(water)/m^3(Air)]
        volScalarField waterContent_;

        //- Max water content [kg(water)/m^3(Air)] for stabilization
        volScalarField maxWaterContent_;

        //- Specific humidity [kg/kg] (water/air)
        volScalarField specificHumidity_;

        //- Max specific humidity [kg/kg] (water/air)
        volScalarField maxSpecificHumidity_;

        //- Saturation pressure of water [Pa]
        volScalarField pSatH2O_;

        //- Partial pressure of water [Pa]
        volScalarField partialPressureH2O_;

        //- Partial pressure calculation method
        mutable word method_;


    // Protected Member Functions

        //- Construct as copy (not implemented)
        humidityRhoThermo(const humidityRhoThermo&);


    // Turbulence fields

        volScalarField muEff_;


public:

    //- Runtime type information
    TypeName("humidityRhoThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        humidityRhoThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    //- Declare run-time constructor selection table for dictionary based
    declareRunTimeSelectionTable
    (
        autoPtr,
        humidityRhoThermo,
        fvMeshDictPhase,
        (const fvMesh& mesh, const word& phaseName, const word& dictName),
        (mesh, phaseName, dictName)
    );


    // Constructors

        //- Construct from mesh and phase name
        humidityRhoThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh, dictionary and phase name
        humidityRhoThermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );

        //- Construct from mesh,dictionary,phase name with a single temperature
        humidityRhoThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictionaryName
        );


    //- Selector
    static autoPtr<humidityRhoThermo> New
    (
        const fvMesh&,
        const word& phaseName=word::null
    );


    //- Selector
    static autoPtr<humidityRhoThermo> New
    (
        const fvMesh&,
        const word& phaseName,
        const word& dictName
    );


    //- Destructor
    virtual ~humidityRhoThermo();


    // Member functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            //  Limit thermo rho between rhoMin and rhoMax
            virtual void correctRho
            (
                const volScalarField& deltaRho,
                const dimensionedScalar& rhoMin,
                const dimensionedScalar& rhoMax
            );

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho);

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> mu(const label patchi) const;

        // Read

            virtual const word readMethod() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
