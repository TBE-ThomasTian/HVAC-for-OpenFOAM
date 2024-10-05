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
    UTCI Calculator

Description
    Dieses Programm berechnet den Universal Thermal Climate Index (UTCI)
    auf Basis der Lufttemperatur, Windgeschwindigkeit, relativen Luftfeuchtigkeit
    und der mittleren Strahlungstemperatur.

Author
    Thomas Tian

Version
    1.0

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volFields.H"
#include "surfaceFields.H"

// Definition der Konstanten
const scalar SIGMA = 5.670374419e-8;  // Stefan-Boltzmann-Konstante in W/m2*K4

// * * * * * * * * * * * * * * * * * Funktionen * * * * * * * * * * * * * * * * //

// Funktion zur Berechnung der mittleren Strahlungstemperatur (Tmrt)
scalarField computeTmrt(const volScalarField& qr)
{
    scalarField Tmrt(qr.size());
    forAll(qr, i)
    {
        // Berechnung von Tmrt in Celsius
        scalar absMrt = Foam::pow(qr[i] / SIGMA, 0.25); // in Kelvin
        Tmrt[i] = absMrt - 273.15; // Umrechnung in Celsius
    }
    return Tmrt;
}

// Funktion zur Berechnung des Wasserdampfdrucks Pa in hPa
scalarField computePa(const volScalarField& Ta, const volScalarField& RH)
{
    scalarField Pa(Ta.size());
    forAll(Ta, i)
    {
        scalar Es = 6.105 * Foam::exp((17.27 * Ta[i]) / (237.7 + Ta[i])); // Saettigungsdampfdruck in hPa
        Pa[i] = (RH[i] / 100.0) * Es; // Wasserdampfdruck in hPa
    }
    return Pa;
}

// Funktion zur Berechnung des UTCI
volScalarField computeUTCI(const volScalarField& Ta, const volVectorField& U, const scalarField& Tmrt, const scalarField& Pa)
{
    volScalarField utci
    (
        IOobject
        (
            "utciTemp",
            Ta.time().timeName(),
            Ta.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ta.mesh(),
        dimensionedScalar("zero", dimTemperature, 0.0)
    );

    forAll(Ta, i)
    {
        scalar windSpeed = mag(U[i]);
        windSpeed = Foam::max(windSpeed, 0.5);
        scalar D_Tmrt = Tmrt[i] - Ta[i];

        // Berechnung des UTCI-Wertes
        utci.internalFieldRef()[i] =  Ta[i] +
                                  (windSpeed) * (-0.281) +
                                  (Pa[i]) * 0.0038 +
                                  (D_Tmrt) * 0.018;
    }

    return utci;
}

// * * * * * * * * * * * * * * * * * Programm * * * * * * * * * * * * * * * * //

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
        Info << "Time = " << runTime.timeName() << endl;

        // Lade die erforderlichen Felder
        #include "createFields.H"

        // Hier die Umrechnung von T von Kelvin zu Grad Celsius**
        // T ist in Kelvin**
	dimensionedScalar T0("T0", dimensionSet(0, 0, 0, 1, 0), 273.15);
	T.internalFieldRef() -= T0;

        // Berechne Tmrt
        scalarField Tmrt = computeTmrt(qr);

        // Berechne Pa
        scalarField Pa = computePa(T, RH);

	// Berechne UTCI
	volScalarField utci = computeUTCI(T, U, Tmrt, Pa);

	// Speichern der Ergebnisse im UTCI-Feld
	UTCI.internalFieldRef() = utci.internalFieldRef();

        // Schreibe das UTCI-Feld in das Verzeichnis
        UTCI.write();
    }

    Info << "UTCI-Berechnung abgeschlossen." << endl;
    return 0;
}

// ************************************************************************* //
