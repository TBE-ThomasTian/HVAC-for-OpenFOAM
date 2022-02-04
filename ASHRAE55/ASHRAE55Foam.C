/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    ASHRAE55

Description
    This tool calculates the following values:
    - Operative room air temperature (TOp)

Background
    ANSI/ASHRAE Standard 55-2020

Autor
    Thomas Tian

Version
    1.0

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "wallFvPatch.H"
#include "cellSet.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"

// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * //

// If no Radiation flux is found, than calculate the surface temperatur as average Ts
Foam::scalar radiationTemperature
(
    const Foam::fvMesh& mesh_,
    const Foam::fvPatchList& Patches_
)
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    scalar PatchSumTemp(0), area(0), sumArea(0);

    forAll (Patches_, patchI)
    {
        const label curPatch = Patches_[patchI].index();

        if (isType<wallFvPatch>( Patches_[patchI] ))
        {
            area = gSum(mesh_.magSf().boundaryField()[curPatch]);

            if (area > 0)
            {
                PatchSumTemp +=
                    gSum
                    (
                        mesh_.magSf().boundaryField()[curPatch]
                      * T.boundaryField()[curPatch]
                    );

                sumArea += area;
            }
        }
   }

    return (PatchSumTemp / sumArea) - 273.15;
}


// * * * * * * * * * * * * * * * * * Program * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    Foam::timeSelector::addOptions();
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    //- Get times list
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"


    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.timeName() << endl;

        #include "createFields.H"

        volScalarField T(THeader, mesh);
        const fvPatchList& Patches = T.mesh().boundary();

        scalar STemp(20);

        volVectorField U(UHeader, mesh);

		scalar t_cmf(0), ce(0);

		//- Radiation Model not available? Lets use classic
        if (G.headerOk()!=1)
		{
            STemp = radiationTemperature(mesh, Patches);
			Info << "Radiation Temperature: " << STemp << " °C" << endl;
		}
		
        forAll (mesh.cells(), cellI)
        {
			ASHRAELevel80[cellI] = 0;
			ASHRAELevel90[cellI] = 0;
			
			//- Lets check weather we can use the radiationModel
			if (G.headerOk() == 1)
			{
				if (G[cellI] < 0)
					G[cellI] = 0;

				//- Because of some bad cells we need to cut values
				if ( G[cellI] > 50000 )
					G[cellI] = 50000;

				STemp = Foam::pow( G[cellI] / ( 4.0 * 0.0000000567), 0.25) - 273.15;
			}

			//- See if the running mean temperature is between 10 °C and 33.5 °C
			if ( (T[cellI]> 283.15) && (T[cellI]< 306.65) )
			{			
				//- Calculate the operative room air temperature
				TOp[cellI] = T[cellI] + (STemp + 273.15) / 2.0;
				
				ce = 0;

				//- Calculate cooling effect (ce) of elevated air speed when Top > 25 degC.		
				if ( (mag(U[cellI]) >= 0.6) && (TOp[cellI]> 298.15) )
				{
					if (mag(U[cellI]) < 0.9) 
					{
						ce 	= 1.2;
					}
					else if (mag(U[cellI]) < 1.9)
					{	
						ce 	= 1.8;
					}
					else
					{					
						ce 	= 2.2;
					}

				}

				//- Figure out the relation between comfort and outdoor
				//- temperature depending on the level of conditioning.
				t_cmf = (0.31 * (T[cellI]-273.15)) + 17.8;
				
				if ( (T[cellI]-273.15 >= ( t_cmf - 3.5 )) && (T[cellI]-273.15 <= ( t_cmf + 3.5 + ce)) )
					ASHRAELevel80[cellI] = 1;
				
				if ( (T[cellI]-273.15 >= ( t_cmf - 2.5 )) && (T[cellI]-273.15 <= ( t_cmf + 2.5 + ce)) )
					ASHRAELevel90[cellI] = 1;

			}

			//- End forAll loop
        }

        ASHRAELevel80.write();
		ASHRAELevel90.write();


        Info << "Done" << endl;

}
    return 0;
}

// ************************************************************************* //

