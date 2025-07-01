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
    Foam::rainPhase

Description
    Class representing a single rain phase with specific droplet size
    in wind-driven rain simulations.

SourceFiles
    rainPhase.C

\*---------------------------------------------------------------------------*/

#ifndef rainPhase_H
#define rainPhase_H

#include "fvCFD.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class rainPhase Declaration
\*---------------------------------------------------------------------------*/

class rainPhase
{
    // Private Data

        //- Phase index
        label index_;
        
        //- Droplet diameter
        scalar diameter_;
        
        //- Volume fraction of this phase in total rainfall
        scalar volumeFraction_;
        
        //- Reference to mesh
        const fvMesh& mesh_;
        
        //- Rain velocity field
        volVectorField Urain_;
        
        //- Rain flux field
        surfaceScalarField phirain_;
        
        //- Rain volume fraction field
        volScalarField alpharain_;
        
        //- Turbulent dispersion coefficient
        volScalarField Ctrain_;
        
        //- Reynolds number field
        volScalarField Re_;
        
        //- Drag coefficient times Reynolds number
        volScalarField CdRe_;
        
        //- Specific catch ratio
        volScalarField scr_;
        
        //- Particle time scale
        volScalarField tp_;
        
        //- Fluid time scale
        volScalarField tfl_;
        
        //- Turbulent viscosity for rain
        volScalarField nutrain_;


    // Private Member Functions

        //- Calculate drag coefficient
        volScalarField calculateCdRe(const volScalarField& Re);
        
        //- Disallow default bitwise copy construction
        rainPhase(const rainPhase&) = delete;
        
        //- Disallow default bitwise assignment
        void operator=(const rainPhase&) = delete;


public:

    // Constructors

        //- Construct from components
        rainPhase
        (
            const label index,
            const scalar diameter,
            const scalar volumeFraction,
            const fvMesh& mesh,
            const word& phaseName
        );


    //- Destructor
    ~rainPhase() = default;


    // Member Functions

        // Access

            //- Return phase index
            inline label index() const
            {
                return index_;
            }
            
            //- Return droplet diameter
            inline scalar diameter() const
            {
                return diameter_;
            }
            
            //- Return volume fraction
            inline scalar volumeFraction() const
            {
                return volumeFraction_;
            }
            
            //- Return rain velocity
            inline volVectorField& Urain()
            {
                return Urain_;
            }
            
            //- Return rain velocity (const)
            inline const volVectorField& Urain() const
            {
                return Urain_;
            }
            
            //- Return rain flux
            inline surfaceScalarField& phirain()
            {
                return phirain_;
            }
            
            //- Return rain volume fraction
            inline volScalarField& alpharain()
            {
                return alpharain_;
            }
            
            //- Return specific catch ratio
            inline const volScalarField& scr() const
            {
                return scr_;
            }


        // Evolution

            //- Solve volume fraction equation
            void solveAlpha();
            
            //- Solve momentum equation
            void solveMomentum
            (
                const volVectorField& U,
                const volVectorField& g,
                const dimensionedScalar& rhoa,
                const dimensionedScalar& rhop,
                const dimensionedScalar& mua,
                const bool solveTD,
                const autoPtr<incompressible::turbulenceModel>& turbulence
            );
            
            //- Calculate catch ratio
            void calculateCatchRatio(const dimensionedScalar& Rh);
            
            //- Update phase properties
            void updateProperties
            (
                const volVectorField& U,
                const dimensionedScalar& rhoa,
                const dimensionedScalar& mua
            );
            
            //- Calculate Courant number
            scalar CourantNumber() const;


        // I/O

            //- Write fields
            void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //