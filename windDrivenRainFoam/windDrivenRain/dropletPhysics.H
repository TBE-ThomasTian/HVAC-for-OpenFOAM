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
    Foam::dropletPhysics

Description
    Physics models for raindrops including:
    - Evaporation
    - Breakup (bag breakup, shear breakup)
    - Collision and coalescence
    - Deformation

SourceFiles
    dropletPhysics.C

\*---------------------------------------------------------------------------*/

#ifndef dropletPhysics_H
#define dropletPhysics_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class dropletPhysics Declaration
\*---------------------------------------------------------------------------*/

class dropletPhysics
{
public:

    // Static Functions

        //- Calculate evaporation rate
        static volScalarField evaporationRate
        (
            const volScalarField& diameter,
            const volScalarField& temperature,
            const volScalarField& relativeHumidity,
            const volScalarField& pressure,
            const volVectorField& Urel
        );
        
        //- Calculate breakup time scale
        static volScalarField breakupTimeScale
        (
            const volScalarField& diameter,
            const volScalarField& Urel,
            const volScalarField& rhoa,
            const volScalarField& rhop,
            const volScalarField& sigma
        );
        
        //- Weber number
        static volScalarField WeberNumber
        (
            const volScalarField& diameter,
            const volScalarField& Urel,
            const volScalarField& rhoa,
            const volScalarField& sigma
        );
        
        //- Critical Weber number for breakup
        static scalar criticalWeber(const scalar Ohnesorge);
        
        //- Ohnesorge number
        static volScalarField OhnesorgeNumber
        (
            const volScalarField& diameter,
            const volScalarField& mup,
            const volScalarField& rhop,
            const volScalarField& sigma
        );
        
        //- Terminal velocity (Gunn & Kinzer)
        static scalar terminalVelocity(const scalar diameter);
        
        //- Drag coefficient enhancement due to deformation
        static volScalarField deformationDragFactor
        (
            const volScalarField& WeberNumber
        );
        
        //- Collision efficiency
        static scalar collisionEfficiency
        (
            const scalar d1,
            const scalar d2,
            const scalar U1,
            const scalar U2
        );
        
        //- Coalescence efficiency
        static scalar coalescenceEfficiency
        (
            const scalar d1,
            const scalar d2,
            const scalar Weber
        );
        
        //- Rain evaporation source term for energy equation
        static volScalarField evaporationHeatSource
        (
            const volScalarField& evapRate,
            const volScalarField& temperature
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //