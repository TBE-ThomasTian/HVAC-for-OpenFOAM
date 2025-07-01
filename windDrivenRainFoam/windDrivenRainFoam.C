/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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
    windDrivenRainFoam

Description
    Solves for wind-driven rain with an Eulerian multiphase model
    Written by Aytac Kubilay, March 2012, ETH Zurich/Empa
	
	Latest Update: 17.02.2015

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "turbulentTransportModel.H"
#include "singlePhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "updateValues.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readGravitationalAcceleration.H"

    simpleControl simple(mesh);
    
    // Read adaptive time step controls
    bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
    scalar maxCoRain = runTime.controlDict().lookupOrDefault<scalar>("maxCoRain", 0.5);
    scalar maxDeltaT = runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	#include "createRainFields.H"
	
	// Create turbulence model if turbulent dispersion is enabled
	autoPtr<incompressible::turbulenceModel> turbulence;
	if (solveTD)
	{
		Info <<nl<< "Solving for the turbulent dispersion of raindrops" << endl;
		#include "createTurbulenceFields.H"
	}
	else 
	{ 
		Info <<nl<< "Turbulent dispersion of raindrops is neglected" << endl; 
	}
	
	GET_parameters(temp,rhoa,mua,rhop);
	
	Info <<nl<< "Temperature: " << temp.value() << " K" << endl;
	Info << "Air density: " << rhoa.value() << " kg/m3" << endl;
	Info << "Air dynamic viscosity: " << mua.value() << " kg/m-s" << endl;
	Info << "Water density: " << rhop.value() << " kg/m3" << endl;
	
	while (simple.loop())
	{
	    Info<< nl << "Time = " << runTime.timeName() << nl;
		
	    // Adaptive time stepping for rain phases
	    #include "adaptiveTimeStep.H"
		
	    for (int nonOrth=0; nonOrth<=simple.nNonOrthCorr(); nonOrth++)
	    {
		    for (int phase_no = 0; phase_no < phases.size(); phase_no++)
			{
				// phi is used by the inletOutlet boundary condition and courant number calculation
				surfaceScalarField phi
				(
					IOobject
					(
						"phi",
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					phirain[phase_no]
				);
				//////////////////////////////////////////////////

                #include "CourantNo.H"
				
				#include "alphaEqns.H"
				
				dimensionedScalar dp ("dp", dimensionSet(0,1,0,0,0,0,0), phases[phase_no][0]); 		
						    
				volScalarField magUr = mag(U - Urain[phase_no]);
				
        		Re = (magUr*dp*rhoa)/mua;
				CdRe = GET_CdRe(Re);
				CdRe.correctBoundaryConditions();

				if (solveTD && turbulence.valid())
				{ 
					// Calculate turbulent time scales
					tmp<volScalarField> tk = turbulence->k();
					tmp<volScalarField> tepsilon = turbulence->epsilon();
					tfl = 0.2*(tk()/tepsilon());
					tp = (4*rhop*dp*dp)/(3*mua*CdRe);
					Ctrain[phase_no] = sqrt( tfl/(tfl+tp) );
					nutrain = turbulence->nut()()*sqr(Ctrain[phase_no]);
				}
				else
				{
					nutrain = dimensionedScalar("zero", nutrain.dimensions(), 0.0);
				}
				
				#include "UEqns.H"

			}
	    }

		if (runTime.outputTime())
		{
			for (int phase_no = 0; phase_no < phases.size(); phase_no++)
			{
				Urain[phase_no].write();				
				alpharain[phase_no].write();
				//Ctrain[phase_no].write();
			}
		}
	    runTime.write();
		
		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << " ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
	}
	
	Info<< "Writing final output\n" << endl;
	for (int phase_no = 0; phase_no < phases.size(); phase_no++)
	{
		Urain[phase_no].write();				
		alpharain[phase_no].write();
		//Ctrain[phase_no].write();
	}
	
	#include "calculateCatchRatio.H"
	
	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
