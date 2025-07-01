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
    This tool calculates thermal comfort according to ASHRAE 55-2020:
    - Operative room air temperature (TOp)
    - ASHRAE 55 adaptive comfort compliance (80% and 90% acceptability)
    - Designed for OpenFOAM 2412+ with solar radiation models
    - Fallback support for older versions using EPW solar calculations

Background
    ANSI/ASHRAE Standard 55-2020 - Adaptive Comfort Model
    Optimized for OpenFOAM 2412+ solar radiation capabilities

Author
    Thomas Tian

Version
    1.2 (Enhanced with improved EPW parser and detailed statistics)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "wallFvPatch.H"
#include "cellSet.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include <fstream>
#include <sstream>
#include <vector>

// Structure to hold EPW weather data
struct EPWData
{
    Foam::List<Foam::scalar> temperature;     // Hourly temperature data
    Foam::List<Foam::scalar> globalRadiation; // Global horizontal radiation (Wh/m2)
    Foam::List<Foam::scalar> directRadiation; // Direct normal radiation (Wh/m2)
    Foam::List<Foam::scalar> diffuseRadiation; // Diffuse horizontal radiation (Wh/m2)
    Foam::scalar latitude;                    // Site latitude from EPW header
    Foam::scalar longitude;                   // Site longitude from EPW header
    Foam::scalar timeZone;                    // Time zone (GMT offset) from EPW header
    Foam::string locationName;                // Location name from EPW header
};

// * * * * * * * * * * * * * * * * * * * Functions  * * * * * * * * * * * * * //

// If no Radiation flux is found, calculate the surface temperature as average Ts
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

// Enhanced solar calculator - ONLY use if no OpenFOAM solar model is active
// Calculate solar MRT using EPW radiation data
Foam::scalar calculateSolarMRT
(
    const EPWData& epwData,
    const Foam::scalar& dayOfYear,
    const Foam::scalar& hour,
    const Foam::scalar& latitude,
    const Foam::scalar& longitude,
    const Foam::scalar& airTemp,
    const bool& showDebug = false
)
{
    // Solar position calculation with proper longitude correction
    scalar dayAngle = 2.0 * M_PI * (dayOfYear - 1) / 365.0;
    scalar declination = 23.45 * Foam::sin(dayAngle + 2.0 * M_PI * (284.0/365.0)) * M_PI/180.0;
    
    // Equation of time correction (in minutes) - standard formula
    scalar eqTimeMinutes = 229.183 * (0.000075 + 0.001868 * Foam::cos(dayAngle) - 
                          0.032077 * Foam::sin(dayAngle) - 0.014615 * Foam::cos(2*dayAngle) - 
                          0.040849 * Foam::sin(2*dayAngle));
    
    // Use time zone from EPW data if available, otherwise estimate from longitude
    scalar timeZoneOffset = epwData.timeZone;  // GMT offset in hours
    scalar standardMeridian = timeZoneOffset * 15.0;  // Convert to degrees
    
    // Time correction from standard meridian to actual longitude
    // 4 minutes per degree of longitude
    scalar longitudeCorrection = 4.0 * (longitude - standardMeridian);
    
    // Solar time = local standard time + equation of time + longitude correction
    scalar solarTime = hour + (eqTimeMinutes + longitudeCorrection) / 60.0;
    
    // Hour angle from solar noon (negative = morning, positive = afternoon)
    scalar hourAngle = (solarTime - 12.0) * 15.0 * M_PI/180.0;
    scalar latRad = latitude * M_PI/180.0;
    
    // Solar elevation angle
    scalar solarElevation = Foam::asin(Foam::sin(latRad) * Foam::sin(declination) + 
                                Foam::cos(latRad) * Foam::cos(declination) * Foam::cos(hourAngle));
    
    if (showDebug)
    {
        Info << "Solar calculation debug:" << endl;
        Info << "  Location: Lat " << latitude << ", Lon " << longitude << endl;
        Info << "  Time zone: GMT" << (timeZoneOffset >= 0 ? "+" : "") << timeZoneOffset << endl;
        Info << "  Standard meridian: " << standardMeridian << " degrees" << endl;
        Info << "  Equation of time: " << eqTimeMinutes << " minutes" << endl;
        Info << "  Longitude correction: " << longitudeCorrection << " minutes" << endl;
        Info << "  Hour: " << hour << ", Solar time: " << solarTime << endl;
        Info << "  Hour angle: " << (hourAngle * 180.0/M_PI) << " degrees" << endl;
        Info << "  Solar elevation: " << (solarElevation * 180.0/M_PI) << " degrees" << endl;
    }
    
    if (solarElevation <= 0) 
    {
        if (showDebug) Info << "Solar elevation <= 0, sun below horizon. Returning air temp." << endl;
        return airTemp; // No sun
    }
    
    // Solar azimuth (needed for directional effects)
    // Note: solarAzimuth calculation commented out as it's not currently used
    // scalar solarAzimuth = Foam::atan2(Foam::sin(hourAngle), 
    //                            Foam::cos(hourAngle) * Foam::sin(latRad) - Foam::tan(declination) * Foam::cos(latRad));
    
    // Get radiation data from EPW for the specific hour
    scalar directNormalIrradiance = 0.0;
    scalar diffuseHorizontalIrradiance = 0.0;
    scalar globalHorizontalIrradiance = 0.0;
    
    // Calculate hour index in EPW data (0-8759)
    // Note: hour should be 0-23, not 1-24
    label hourOfYear = label((dayOfYear - 1) * 24 + hour);
    if (showDebug)
    {
        Info << "\n=== Solar MRT Calculation Debug ===" << endl;
        Info << "Hour of year index: " << hourOfYear << " (day " << dayOfYear << ", hour " << hour << ")" << endl;
    }
    if (hourOfYear >= 0 && hourOfYear < 8760)
    {
        // EPW files store radiation in Wh/m² (energy), not W/m² (power)
        // Since these are hourly integrated values, they are already in W/m² (1 hour average)
        globalHorizontalIrradiance = epwData.globalRadiation[hourOfYear];
        directNormalIrradiance = epwData.directRadiation[hourOfYear];
        diffuseHorizontalIrradiance = epwData.diffuseRadiation[hourOfYear];
        
        // Debug output
        if (showDebug)
        {
            Info << "EPW radiation data for hour " << hour << " on day " << dayOfYear << ":" << endl;
            Info << "  Global horizontal: " << globalHorizontalIrradiance << " W/m2" << endl;
            Info << "  Direct normal: " << directNormalIrradiance << " W/m2" << endl;
            Info << "  Diffuse horizontal: " << diffuseHorizontalIrradiance << " W/m2" << endl;
        }
    }
    else
    {
        WarningInFunction
            << "Hour index " << hourOfYear << " out of range for EPW data" << endl;
    }
    
    scalar skyViewFactor = 0.5;
    scalar groundViewFactor = 0.5;
    
    scalar skyTemp = airTemp - 10.0;
    scalar groundTemp = airTemp + (directNormalIrradiance > 0 ? 5.0 : 0.0);
    
    scalar solarAbsorptivity = 0.7;
    scalar effectiveArea = 0.77;
    
    // Calculate solar heat gain on the person
    scalar solarHeatGain = solarAbsorptivity * effectiveArea * 
                          (directNormalIrradiance * Foam::sin(solarElevation) + 
                           diffuseHorizontalIrradiance * skyViewFactor);
    
    // Debug output for solar heat gain
    if (showDebug)
    {
        Info << "  Solar heat gain: " << solarHeatGain << " W/m²" << endl;
        Info << "  Solar elevation: " << (solarElevation * 180.0/M_PI) << " degrees" << endl;
    }
    
    // Calculate equivalent temperature rise from solar radiation
    // Using empirical correlation for outdoor thermal comfort
    // Typical values: 100 W/m² solar gain ≈ 2-3C MRT increase
    scalar solarTempRise = 0.0;
    if (solarHeatGain > 0)
    {
        // Empirical correlation for outdoor thermal comfort
        // Based on Thorsson et al. (2007) and similar studies
        // MRT rise approximately 0.25°C per 10 W/m² absorbed solar radiation
        solarTempRise = 0.025 * solarHeatGain;
        
        // Limit solar temperature rise to realistic values
        solarTempRise = min(solarTempRise, 25.0);  // Max 25C rise from solar
    }
    
    if (showDebug) Info << "  Solar temperature rise: " << solarTempRise << " C" << endl;
    
    // Calculate mean radiant temperature including solar effects
    scalar meanRadiantTemp = skyViewFactor * skyTemp + 
                            groundViewFactor * groundTemp + 
                            solarTempRise;
    
    if (showDebug) Info << "  Final MRT: " << meanRadiantTemp << " C (Air temp: " << airTemp << " C)" << endl;
    
    return meanRadiantTemp;
}

// Read EPW file and return complete weather data
EPWData readEPWFile(const Foam::fileName& epwFile)
{
    EPWData data;
    data.temperature.setSize(8760);    // 365 days * 24 hours
    data.globalRadiation.setSize(8760);
    data.directRadiation.setSize(8760);
    data.diffuseRadiation.setSize(8760);
    
    std::ifstream file(epwFile);
    if (!file.is_open())
    {
        FatalErrorInFunction
            << "Cannot open EPW file: " << epwFile
            << exit(FatalError);
    }
    
    // Parse EPW header to extract location data
    std::string line;
    
    // Line 1: LOCATION header contains location data
    if (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> tokens;
        
        // Split by comma
        while (std::getline(iss, token, ','))
        {
            tokens.push_back(token);
        }
        
        // EPW format: LOCATION,City,State,Country,Source,WMO,Latitude,Longitude,TimeZone,Elevation
        if (tokens.size() >= 10 && tokens[0] == "LOCATION")
        {
            data.locationName = tokens[1] + ", " + tokens[3];  // City, Country
            data.latitude = std::stod(tokens[6]);
            data.longitude = std::stod(tokens[7]);
            data.timeZone = std::stod(tokens[8]);
            
            Info << "EPW Location Data:" << endl;
            Info << "  Location: " << data.locationName << endl;
            Info << "  Latitude: " << data.latitude << " degrees" << endl;
            Info << "  Longitude: " << data.longitude << " degrees" << endl;
            Info << "  Time Zone: GMT" << (data.timeZone >= 0 ? "+" : "") << data.timeZone << endl;
        }
        else
        {
            WarningInFunction
                << "Could not parse location data from EPW header" << endl;
            // Set default values
            data.latitude = 50.0;
            data.longitude = 8.0;
            data.timeZone = 1.0;
            data.locationName = "Unknown Location";
        }
    }
    
    // Skip remaining header lines (2-8)
    for (int i = 1; i < 8; i++)
    {
        std::getline(file, line);
    }
    
    label hourIndex = 0;
    while (std::getline(file, line) && hourIndex < 8760)
    {
        std::istringstream iss(line);
        std::string token;
        
        // Skip to temperature (column 7)
        for (int i = 0; i < 6; i++)
        {
            std::getline(iss, token, ',');
        }
        
        // Column 7: Dry bulb temperature
        std::getline(iss, token, ',');
        data.temperature[hourIndex] = std::stod(token);
        
        // Skip to radiation data (columns 13-15)
        for (int i = 0; i < 5; i++)
        {
            std::getline(iss, token, ',');
        }
        
        // Column 13: Global horizontal radiation
        std::getline(iss, token, ',');
        data.globalRadiation[hourIndex] = std::stod(token);
        
        // Column 14: Direct normal radiation
        std::getline(iss, token, ',');
        data.directRadiation[hourIndex] = std::stod(token);
        
        // Column 15: Diffuse horizontal radiation
        std::getline(iss, token, ',');
        data.diffuseRadiation[hourIndex] = std::stod(token);
        
        hourIndex++;
    }
    
    file.close();
    return data;
}

// Improved EPW parser with robust day-of-year calculation
Foam::scalar calculateRunningMeanFromEPW
(
    const Foam::fileName& epwFile,
    const Foam::scalar& dayOfYear
)
{
    std::ifstream file(epwFile);
    if (!file.is_open())
    {
        FatalErrorInFunction
            << "Cannot open EPW file: " << epwFile
            << exit(FatalError);
    }

    // Skip EPW header (8 lines)
    std::string line;
    for (int i = 0; i < 8; i++)
    {
        std::getline(file, line);
    }

    // Read hourly data and calculate daily averages
    List<scalar> dailyTemps(365, 0.0);
    List<label> hourCount(365, 0);
    
    // Days in each month (non-leap year)
    int daysInMonth[] = {31,28,31,30,31,30,31,31,30,31,30,31};

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string token;
        
        // Parse CSV: year,month,day,hour,minute,datasource,temp,dewpoint,...
        std::getline(iss, token, ','); // year
        label year = std::stoi(token);
        std::getline(iss, token, ','); // month
        label month = std::stoi(token);
        std::getline(iss, token, ','); // day  
        label day = std::stoi(token);
        std::getline(iss, token, ','); // hour
        std::getline(iss, token, ','); // minute
        std::getline(iss, token, ','); // datasource
        std::getline(iss, token, ','); // dry bulb temperature
        scalar temp = std::stod(token);
        
        // Calculate day of year properly
        label doy = 0;
        
        // Check for leap year
        bool isLeapYear = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
        if (isLeapYear && month > 2) doy += 1; // Add leap day if after February
        
        // Sum days of previous months
        for (int m = 0; m < month-1; m++) {
            doy += daysInMonth[m];
        }
        doy += day - 1; // Add current day (0-indexed)
        
        if (doy < 365)
        {
            dailyTemps[doy] += temp;
            hourCount[doy]++;
        }
    }

    // Average daily temperatures
    label validDays = 0;
    forAll(dailyTemps, dayI)
    {
        if (hourCount[dayI] > 0)
        {
            dailyTemps[dayI] /= hourCount[dayI];
            validDays++;
        }
        else
        {
            // Handle missing days by interpolation
            if (dayI > 0 && dayI < 364)
            {
                // Find nearest valid days
                label prevDay = dayI - 1;
                label nextDay = dayI + 1;
                while (prevDay > 0 && hourCount[prevDay] == 0) prevDay--;
                while (nextDay < 365 && hourCount[nextDay] == 0) nextDay++;
                
                if (hourCount[prevDay] > 0 && hourCount[nextDay] > 0)
                {
                    // Linear interpolation
                    scalar weight = scalar(dayI - prevDay) / scalar(nextDay - prevDay);
                    dailyTemps[dayI] = dailyTemps[prevDay] * (1.0 - weight) + 
                                       dailyTemps[nextDay] * weight;
                }
            }
        }
    }
    
    Info << "EPW file processed: " << validDays << " valid days found" << endl;

    // Calculate running mean (ASHRAE 55 formula)
    scalar alpha = 0.8;
    scalar runningMean = 0.0;
    scalar totalWeight = 0.0;

    label targetDay = label(dayOfYear) - 1;

    for (int i = 1; i <= 30; i++)
    {
        label pastDay = (targetDay - i + 365) % 365;
        scalar weight = alpha * Foam::pow(1.0 - alpha, i - 1);
        runningMean += weight * dailyTemps[pastDay];
        totalWeight += weight;
    }

    scalar T_rm = runningMean / totalWeight;
    
    // Sanity check
    if (T_rm < -50 || T_rm > 50)
    {
        WarningInFunction
            << "Calculated running mean temperature (" << T_rm 
            << "C) seems unrealistic!" << endl;
    }

    return T_rm;
}

// Function to calculate and display comfort statistics
void displayComfortStatistics
(
    const volScalarField& ASHRAELevel80,
    const volScalarField& ASHRAELevel90,
    const volScalarField& TOp,
    const volScalarField& T,
    const volVectorField& U,
    const scalar& T_rm_out
)
{
    // Calculate statistics
    label cells80 = 0;
    label cells90 = 0;
    label totalCells = 0;
    
    scalar minTOp = 1000.0;
    scalar maxTOp = -1000.0;
    scalar avgTOp = 0.0;
    scalar avgAirTemp = 0.0;
    scalar avgVelocity = 0.0;
    
    forAll(ASHRAELevel80, cellI)
    {
        if (T[cellI] > 283.15 && T[cellI] < 306.65) // Valid temperature range
        {
            totalCells++;
            avgTOp += TOp[cellI];
            avgAirTemp += T[cellI];
            avgVelocity += mag(U[cellI]);
            
            if (TOp[cellI] < minTOp) minTOp = TOp[cellI];
            if (TOp[cellI] > maxTOp) maxTOp = TOp[cellI];
            
            if (ASHRAELevel80[cellI] > 0.5) cells80++;
            if (ASHRAELevel90[cellI] > 0.5) cells90++;
        }
    }
    
    if (totalCells > 0)
    {
        avgTOp /= totalCells;
        avgAirTemp /= totalCells;
        avgVelocity /= totalCells;
        
        Info<< nl << "================== ASHRAE 55 Comfort Analysis ==================" << endl;
        Info<< "Running mean outdoor temperature: " << T_rm_out << " C" << endl;
        Info<< "Neutral temperature (t_cmf): " << (0.31 * T_rm_out + 17.8) << " C" << endl;
        Info<< nl << "Domain Statistics:" << endl;
        Info<< "  Total cells evaluated: " << totalCells << endl;
        Info<< "  Average air temperature: " << (avgAirTemp - 273.15) << " C" << endl;
        Info<< "  Average operative temperature: " << (avgTOp - 273.15) << " C" << endl;
        Info<< "  Operative temperature range: " << (minTOp - 273.15) 
            << " to " << (maxTOp - 273.15) << " C" << endl;
        Info<< "  Average air velocity: " << avgVelocity << " m/s" << endl;
        Info<< nl << "Comfort Compliance:" << endl;
        Info<< "  80% acceptability: " << cells80 << " cells (" 
            << (scalar(cells80)/totalCells*100.0) << "%)" << endl;
        Info<< "  90% acceptability: " << cells90 << " cells (" 
            << (scalar(cells90)/totalCells*100.0) << "%)" << endl;
        
        // Comfort zone ranges
        scalar t_cmf = 0.31 * T_rm_out + 17.8;
        Info<< nl << "Comfort Zone Ranges (without cooling effect):" << endl;
        Info<< "  80% acceptability: " << (t_cmf - 3.5) << " to " 
            << (t_cmf + 3.5) << " C" << endl;
        Info<< "  90% acceptability: " << (t_cmf - 2.5) << " to " 
            << (t_cmf + 2.5) << " C" << endl;
        
        // Check if cooling effect is active
        label cellsWithCooling = 0;
        forAll(U, cellI)
        {
            if (mag(U[cellI]) >= 0.6 && TOp[cellI] > 298.15) cellsWithCooling++;
        }
        
        if (cellsWithCooling > 0)
        {
            Info<< nl << "Cooling effect active in " << cellsWithCooling 
                << " cells (" << (scalar(cellsWithCooling)/totalCells*100.0) 
                << "%)" << endl;
        }
        
        Info<< "================================================================" << nl << endl;
    }
    else
    {
        Info<< nl << "Warning: No cells in valid temperature range for comfort analysis!" << endl;
    }
}

// * * * * * * * * * * * * * * * * * * Program  * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    Foam::timeSelector::addOptions();

    // Add command line options
    argList::addOption
    (
        "epw",
        "fileName",
        "EPW weather file to calculate running mean outdoor temperature"
    );

    argList::addOption
    (
        "dayOfYear", 
        "scalar",
        "Day of year (1-365) for EPW calculation (default: 180 = June 29)"
    );

    argList::addBoolOption
    (
        "solarData",
        "Use detailed solar calculations from EPW file for accurate MRT (requires -epw option)"
    );

    argList::addOption
    (
        "latitude",
        "scalar", 
        "Site latitude in degrees (overrides EPW file value if provided)"
    );

    argList::addOption
    (
        "longitude", 
        "scalar",
        "Site longitude in degrees (overrides EPW file value if provided)"
    );

    argList::addOption
    (
        "hour",
        "scalar",
        "Hour of day for solar calculations (0-24, default: 12.0)"
    );

    argList::addOption
    (
        "runningMean",
        "scalar", 
        "Directly specify running mean outdoor temperature in C"
    );

    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    //- Get times list
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    // Determine running mean outdoor temperature and solar calculation settings
    scalar T_rm_out = 20.0; // Default value
    bool useSolarCalculations = args.found("solarData");
    scalar latitude = 50.0;  // Default: Central Europe
    scalar longitude = 8.0;  // Default: Central Europe
    scalar hourOfDay = 12.0; // Default: noon
    fileName epwFile;
    scalar dayOfYear = 180;
    EPWData epwData;  // Store complete EPW weather data
    bool useEPWCoordinates = true; // Flag to track coordinate source

    if (args.found("hour"))
        hourOfDay = args.get<scalar>("hour");

    if (args.found("runningMean"))
    {
        T_rm_out = args.get<scalar>("runningMean");
        Info << "Using directly specified running mean: " << T_rm_out << " C" << endl;
        
        // If no EPW file, must use command line coordinates
        if (args.found("latitude"))
            latitude = args.get<scalar>("latitude");
        if (args.found("longitude"))
            longitude = args.get<scalar>("longitude");
        useEPWCoordinates = false;
    }
    else if (args.found("epw"))
    {
        epwFile = args.get<fileName>("epw");
        
        if (args.found("dayOfYear"))
            dayOfYear = args.get<scalar>("dayOfYear");
            
        Info << "Reading EPW file: " << epwFile << endl;
        Info << "Day of year: " << dayOfYear << endl;
        
        // Read complete EPW data including radiation and location
        epwData = readEPWFile(epwFile);
        
        // Use EPW coordinates by default, but allow command line override
        latitude = epwData.latitude;
        longitude = epwData.longitude;
        
        // Check for command line overrides
        if (args.found("latitude"))
        {
            latitude = args.get<scalar>("latitude");
            useEPWCoordinates = false;
            Info << "Overriding EPW latitude with command line value: " << latitude << endl;
        }
        if (args.found("longitude"))
        {
            longitude = args.get<scalar>("longitude");
            useEPWCoordinates = false;
            Info << "Overriding EPW longitude with command line value: " << longitude << endl;
        }
        
        T_rm_out = calculateRunningMeanFromEPW(epwFile, dayOfYear);
        Info << "Calculated running mean: " << T_rm_out << " C" << endl;
        
        if (useSolarCalculations)
        {
            Info << nl << "Solar calculations enabled" << endl;
            Info << "Using coordinates from " << (useEPWCoordinates ? "EPW file" : "command line") << ":" << endl;
            Info << "  Latitude: " << latitude << " degrees" << endl;
            Info << "  Longitude: " << longitude << " degrees" << endl;
            if (useEPWCoordinates)
            {
                Info << "  Location: " << epwData.locationName << endl;
            }
            Info << "Note: For best results, use OpenFOAM 2412+ native solar radiation model instead" << endl;
        }
    }
    else
    {
        Info << "Warning: No EPW file or running mean specified." << endl;
        Info << "Using default value: " << T_rm_out << " C" << endl;
        Info << "Usage examples:" << endl;
        Info << "  # With OpenFOAM solar model (G contains solar radiation):" << endl;
        Info << "  ASHRAE55 -epw weather.epw -dayOfYear 180" << endl;
        Info << "  ASHRAE55 -runningMean 22.5" << endl;
        Info << "  # With EPW solar calculations (auto-detects location from EPW):" << endl;
        Info << "  ASHRAE55 -epw weather.epw -solarData -hour 14" << endl;
        Info << "  # Override EPW coordinates if needed:" << endl; 
        Info << "  ASHRAE55 -epw weather.epw -solarData -latitude 52.5 -longitude 13.4" << endl;
        
        if (useSolarCalculations)
        {
            FatalErrorInFunction
                << "Solar calculations require EPW file (-epw option)"
                << exit(FatalError);
        }
    }

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

        //- Check if radiation field G is available
        Info << "G field check: G.headerOk() = " << G.headerOk() 
             << ", useSolarCalculations = " << useSolarCalculations 
             << ", epwFile.size() = " << epwFile.size() << endl;
        
        //- Radiation Model not available? Use area-weighted wall temperature
        if (G.headerOk()!=1 && !useSolarCalculations)
        {
            STemp = radiationTemperature(mesh, Patches);
            Info << "Using area-weighted wall temperature for MRT: " << STemp << " C" << endl;
        }

        forAll (mesh.cells(), cellI)
        {
            ASHRAELevel80[cellI] = 0;
            ASHRAELevel90[cellI] = 0;
            
            //- Use radiation model if available
            if (G.headerOk() == 1)
            {
                if (G[cellI] < 0)
                    G[cellI] = 0;

                //- Limit maximum radiation to avoid numerical issues
                if ( G[cellI] > 50000 )
                    G[cellI] = 50000;

                //- For solar radiation: Convert to realistic mean radiant temperature
                //- CRITICAL: Stefan-Boltzmann inversion gives unrealistic values for solar irradiance
                //- Use empirical correlations based on outdoor thermal comfort measurements
                
                scalar T_radiant_equiv = Foam::pow( G[cellI] / ( 4.0 * 5.67e-8), 0.25) - 273.15;
                
                if (T_radiant_equiv > 60.0)  // Unrealistic for outdoor MRT
                {
                    //- Use empirical solar effect correlation
                    //- Based on measured outdoor MRT vs solar irradiance data
                    //- Typical outdoor MRT = air temperature + solar effect (5-20C)
                    scalar solar_effect = 5.0 + Foam::sqrt(G[cellI] / 100.0);  // Empirical: 5-20C rise
                    solar_effect = min(solar_effect, 20.0);  // Max 20C solar heating
                    
                    STemp = (T[cellI] - 273.15) + solar_effect;
                }
                else
                {
                    //- Use Stefan-Boltzmann result for low radiation
                    STemp = T_radiant_equiv;
                }
                
                //- Final safety limits for extreme conditions
                STemp = min(STemp, 65.0);   // Max MRT = 65C (hot asphalt limit)
                STemp = max(STemp, -30.0);  // Min MRT = -30C (extreme cold)
            }
            //- Use EPW solar calculations if radiation model not available and -solarData specified
            else if (G.headerOk() != 1 && useSolarCalculations && epwFile.size() > 0)
            {
                if (cellI == 0)
                {
                    Info << "Entering EPW solar calculation branch" << endl;
                }
                //- Use the specified hour of day with actual EPW radiation data
                bool showDebugInfo = (cellI == 0);  // Only show debug for first cell
                STemp = calculateSolarMRT(
                    epwData,
                    dayOfYear,
                    hourOfDay,  // Use command line specified hour
                    latitude,
                    longitude,
                    T[cellI] - 273.15,  // Air temperature in Celsius
                    showDebugInfo
                );
                
                if (cellI == 0)
                {
                    Info<< "Solar MRT calculation: Air temp = " << (T[cellI] - 273.15) 
                        << " C, Solar MRT = " << STemp << " C" << endl;
                }
            }
            else
            {
                // No radiation field and no solar calculations - use air temperature as MRT
                STemp = T[cellI] - 273.15;
                if (cellI == 0)
                {
                    Info << "No radiation data available, using air temperature as MRT" << endl;
                }
            }

            //- Only evaluate comfort if temperature is in reasonable range (10C to 33.5C)
            if ( (T[cellI] > 283.15) && (T[cellI] < 306.65) )
            {           
                //- Calculate operative temperature: average of air and mean radiant temperature
                TOp[cellI] = ( T[cellI] + (STemp + 273.15)) / 2.0;
                
                ce = 0;

                //- Calculate cooling effect of elevated air speed when Top > 25C
                if ( (mag(U[cellI]) >= 0.6) && (TOp[cellI] > 298.15) )
                {
                    if (mag(U[cellI]) < 0.9) 
                    {
                        ce = 1.2;
                    }
                    else if (mag(U[cellI]) < 1.2)  // CORRECTED: was 1.9, should be 1.2 per ASHRAE 55
                    {   
                        ce = 1.8;
                    }
                    else
                    {                   
                        ce = 2.2;
                    }
                }

                //- Calculate neutral temperature for adaptive comfort model
                //- CORRECTED: Use running mean outdoor temperature, not indoor temperature
                t_cmf = (0.31 * T_rm_out) + 17.8;
                
                //- Convert to Kelvin for comparison with TOp
                scalar t_cmf_K = t_cmf + 273.15;
                
                //- CORRECTED: Proper range checking for comfort zones
                //- 80% acceptability: ±3.5K around neutral temperature + cooling effect
                if ((TOp[cellI] >= (t_cmf_K - 3.5)) && (TOp[cellI] <= (t_cmf_K + 3.5 + ce)))
                    ASHRAELevel80[cellI] = 1;
                
                //- 90% acceptability: ±2.5K around neutral temperature + cooling effect
                if ((TOp[cellI] >= (t_cmf_K - 2.5)) && (TOp[cellI] <= (t_cmf_K + 2.5 + ce)))
                    ASHRAELevel90[cellI] = 1;
            }
        }

        // Display detailed comfort statistics
        displayComfortStatistics(ASHRAELevel80, ASHRAELevel90, TOp, T, U, T_rm_out);

        ASHRAELevel80.write();
        ASHRAELevel90.write();

        Info << "Done" << endl;
    }

    return 0;
}


// ************************************************************************* //
