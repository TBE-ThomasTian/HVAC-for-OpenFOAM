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
    
    Required fields:
    - T: Air temperature [K]
    - U: Velocity field [m/s]
    - relHum: Relative humidity [%]
    - qr: Radiation flux density [W/m2] (or provide Tmrt directly)
    
    Optional fields:
    - SVF: Sky View Factor [-] (default: 0.7)
    - Tmrt: Mean radiant temperature [K or °C] (if not provided, calculated from qr)
    
    Note: This implementation currently uses a simplified linear approximation.
    For accurate UTCI values, the full 6th-order polynomial with ~200 coefficients
    is required (Bröde et al. 2012, Int J Biometeorol 56:481-494).
    
    Valid ranges:
    - Air temperature: -50 to +50°C
    - Wind speed (10m): 0.5 to 30 m/s
    - Vapor pressure: 0 to 50 hPa
    - Tmrt - Ta: -30 to +70°C

Author
    Thomas Tian

Version
    2.0 - Improved with SVF, Tmrt options, and wind adjustment

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "utciPolynomial.H"
#include "meanRadiantTemperature.H"
#include "skyViewFactor.H"
#include "epwReader.H"
#include "atmBoundaryLayerInletVelocityFvPatchVectorField.H"
#include "radiationModel.H"
#include "externalWallHeatFluxTemperatureFvPatchScalarField.H"

// Definition der Konstanten
const scalar SIGMA = 5.670374419e-8;  // Stefan-Boltzmann-Konstante in W/m2*K4

// * * * * * * * * * * * * * * * * * Funktionen * * * * * * * * * * * * * * * * //

// Funktion zur Berechnung der mittleren Strahlungstemperatur (Tmrt)
volScalarField computeTmrt(const volScalarField& T, const volScalarField& qr, const volScalarField& SVF)
{
    // Use the improved mean radiant temperature calculation
    return calculateMeanRadiantTemperature(T, qr, SVF);
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

// Funktion zur Anpassung der Windgeschwindigkeit auf 10m Höhe
scalar adjustWindTo10m(scalar windSpeed, scalar height = 1.5)
{
    // Power law wind profile adjustment
    // v10 = v_h * (10/h)^alpha, where alpha ≈ 0.25 for urban areas
    const scalar alpha = 0.25;
    const scalar refHeight = 10.0;
    
    if (height > 0 && height != refHeight)
    {
        return windSpeed * Foam::pow(refHeight/height, alpha);
    }
    return windSpeed;
}

// Funktion zur Berechnung des UTCI
volScalarField computeUTCI(const volScalarField& Ta, const volVectorField& U, const volScalarField& Tmrt, const scalarField& Pa)
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
        // Adjust wind speed to 10m height (assuming measurement at 1.5m)
        scalar v10 = adjustWindTo10m(windSpeed, 1.5);
        v10 = Foam::max(v10, 0.5);  // Minimum wind speed
        v10 = Foam::min(v10, 30.0);  // Maximum valid wind speed for UTCI
        
        scalar D_Tmrt = Tmrt[i] - Ta[i];
        
        // Use the polynomial approximation (it includes validation)
        utci.internalFieldRef()[i] = calculateUTCIPolynomial(Ta[i], v10, D_Tmrt, Pa[i]);
    }

    return utci;
}

// * * * * * * * * * * * * * * * * * Programm * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    Foam::timeSelector::addOptions();
    
    argList::addOption
    (
        "epw",
        "file",
        "EPW weather file for radiation calculation"
    );
    
    argList::addOption
    (
        "dateTime",
        "YYYY-MM-DD HH:MM",
        "Date and time for weather data (default: July 15, 14:00)"
    );
    
    argList::addOption
    (
        "month",
        "value",
        "Month for calculation (1-12, default: 7)"
    );
    
    argList::addOption
    (
        "day",
        "value", 
        "Day of month (1-31, default: 15)"
    );
    
    argList::addOption
    (
        "hour",
        "value",
        "Hour of day (0-23, default: 14)"
    );
    
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    //- Get times list
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"
    
    // Check for EPW file
    autoPtr<epwReader> epwReaderPtr;
    bool useEPW = false;
    scalar month = 7;    // Default: July
    scalar day = 15;     // Default: 15th
    scalar hour = 14;    // Default: 2 PM
    
    if (args.found("epw"))
    {
        fileName epwFile = args.get<fileName>("epw");
        Info<< "\nReading EPW weather file: " << epwFile << endl;
        
        epwReaderPtr.reset(new epwReader(epwFile));
        useEPW = true;
        
        // Get date/time options
        if (args.found("month"))
        {
            month = args.get<scalar>("month");
        }
        if (args.found("day"))
        {
            day = args.get<scalar>("day");
        }
        if (args.found("hour"))
        {
            hour = args.get<scalar>("hour");
        }
        
        Info<< "Using weather data for: "
            << "Month=" << month 
            << " Day=" << day 
            << " Hour=" << hour << ":00" << endl;
        Info<< "Location: Lat=" << epwReaderPtr->latitude() 
            << " Lon=" << epwReaderPtr->longitude() << endl;
    }

    forAll(timeDirs, timei)
    {

        runTime.setTime(timeDirs[timei], timei);
        Info << "Time = " << runTime.timeName() << endl;

        // Lade die erforderlichen Felder
        #include "createFields.H"

        // Calculate radiation if radiation model is active
        radiation->correct();

        // Hier die Umrechnung von T von Kelvin zu Grad Celsius**
        // T ist in Kelvin**
	dimensionedScalar T0("T0", dimensionSet(0, 0, 0, 1, 0), 273.15);
	T.internalFieldRef() -= T0;

        // Check if SVF was provided, otherwise calculate it
        bool SVFProvided = (SVF.headerOk() && max(SVF).value() > SMALL);
        
        if (!SVFProvided)
        {
            Info << "\nSky View Factor not provided. Calculating..." << endl;
            
            // Check if detailed calculation is feasible
            if (mesh.nCells() < 100000)  // For smaller meshes, use ray tracing
            {
                SVF = calculateSkyViewFactor(mesh, 50);  // 50 rays per cell
            }
            else  // For larger meshes, use simplified method
            {
                Info << "Using simplified SVF calculation for large mesh..." << endl;
                SVF = calculateSimplifiedSVF(mesh, 0.7);
            }
        }
        else
        {
            Info << "Using provided Sky View Factor field..." << endl;
            Info << "SVF range: " << min(SVF).value() << " to " << max(SVF).value() << endl;
        }

        // Check if Tmrt was provided directly, otherwise compute it
        bool TmrtProvided = (Tmrt.headerOk() && max(mag(Tmrt)).value() > SMALL);
        
        if (!TmrtProvided)
        {
            if (useEPW)
            {
                Info << "\nComputing mean radiant temperature from EPW weather data..." << endl;
                
                // Calculate time value for EPW interpolation
                scalar dayOfYear = 30*(month-1) + day;  // Simplified
                scalar timeValue = (dayOfYear - 1) * 24 * 3600 + hour * 3600;
                
                // Get weather data
                epwData weather = epwReaderPtr->interpolateData(timeValue);
                
                Info << "Weather conditions:" << endl;
                Info << "  Dry bulb temperature: " << weather.dryBulbTemp << " C" << endl;
                Info << "  Relative humidity: " << weather.relHumidity << " %" << endl;
                Info << "  Direct normal radiation: " << weather.directNormalRadiation << " W/m2" << endl;
                Info << "  Diffuse horizontal radiation: " << weather.diffuseHorizontalRadiation << " W/m2" << endl;
                Info << "  Wind speed: " << weather.windSpeed << " m/s" << endl;
                
                // Calculate Tmrt for each cell using EPW data
                forAll(Tmrt, cellI)
                {
                    // Use EPW data with local SVF
                    scalar groundTemp = T[cellI];  // Approximate ground temp as air temp
                    Tmrt[cellI] = epwReaderPtr->calculateTmrt(
                        weather,
                        T[cellI],
                        SVF[cellI],
                        groundTemp
                    );
                }
            }
            else
            {
                Info << "\nComputing mean radiant temperature from radiation flux..." << endl;
                Tmrt = computeTmrt(T, qr, SVF);
            }
        }
        else
        {
            Info << "Using provided mean radiant temperature field..." << endl;
            // Convert from Kelvin to Celsius if needed
            if (max(Tmrt).value() > 100)  // Likely in Kelvin
            {
                Tmrt -= 273.15;
            }
        }

        // Berechne Pa
        scalarField Pa = computePa(T, RH);

	// Berechne UTCI
	volScalarField utci = computeUTCI(T, U, Tmrt, Pa);

	// Speichern der Ergebnisse im UTCI-Feld
	UTCI.internalFieldRef() = utci.internalFieldRef();

        // Schreibe das UTCI-Feld in das Verzeichnis
        UTCI.write();
    }

    return 0;
}

// ************************************************************************* //
