/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    AoAFoam

Group
    grpHVACSolvers

Description
    Transient solver for Age of Air (AoA) calculation in HVAC applications.
    
    Solves the transport equation for Age of Air with a unit source term:
    
    dAoA/dt + div(phi, AoA) - laplacian(DT, AoA) = 1
    
    where:
    - AoA is the age of air [s]
    - phi is the face flux field [m3/s]
    - DT is the turbulent diffusivity [m2/s]
    
    The solver uses the PIMPLE algorithm for time-stepping with 
    adjustable time step control.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "turbulentTransportModel.H"
#include "singlePhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for Age of Air (AoA) calculation in HVAC applications"
    );

    #include "postProcess.H"
    
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop for Age of Air calculation\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- PIMPLE loop
        while (pimple.loop())
        {
            fvScalarMatrix AoAEqn
            (
                fvm::ddt(AoA)
              + fvm::div(phi, AoA)
              - fvm::laplacian(turbulence->nut() + DT, AoA)
              - fvOptions(AoA)
             ==
                dimensionedScalar("ageSource", AoA.dimensions()/dimTime, 1.0)
            );

            AoAEqn.relax();
            
            fvOptions.constrain(AoAEqn);
            
            AoAEqn.solve();
            
            fvOptions.correct(AoA);
        }

        // Calculate local mean age of air
        if (runTime.writeTime())
        {
            volScalarField localMeanAge
            (
                IOobject
                (
                    "localMeanAge",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                AoA
            );
            
            // Calculate ventilation efficiency
            dimensionedScalar volumeTotal = gSum(mesh.V());
            dimensionedScalar meanAoA = fvc::domainIntegrate(AoA)/volumeTotal;
            
            volScalarField ventilationEfficiency
            (
                IOobject
                (
                    "ventilationEfficiency",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                meanAoA/(2.0*AoA + dimensionedScalar("small", dimTime, SMALL))
            );
            
            Info<< "Mean Age of Air = " << meanAoA.value() << " s" << endl;
            Info<< "Min/Max AoA = " << min(AoA).value() << "/" << max(AoA).value() << " s" << endl;
            
            localMeanAge.write();
            ventilationEfficiency.write();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //