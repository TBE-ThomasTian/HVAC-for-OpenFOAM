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

Description
    EPW (EnergyPlus Weather) file reader for UTCI calculations
    
    Reads weather data including:
    - Direct normal irradiance (DNI)
    - Diffuse horizontal irradiance (DHI)
    - Dry bulb temperature
    - Relative humidity
    - Wind speed and direction
    
    EPW format specification:
    https://designbuilder.co.uk/cahelp/Content/EnergyPlusWeatherFileFormat.htm

\*---------------------------------------------------------------------------*/

#ifndef epwReader_H
#define epwReader_H

#include "Time.H"
#include "IFstream.H"
#include "DynamicList.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// EPW data structure
struct epwData
{
    scalar year;
    scalar month;
    scalar day;
    scalar hour;
    scalar minute;
    
    // Weather data
    scalar dryBulbTemp;        // °C
    scalar dewPointTemp;       // °C
    scalar relHumidity;        // %
    scalar atmPressure;        // Pa
    scalar extraterrestrialRadiation;     // Wh/m²
    scalar extraterrestrialDirectRadiation; // Wh/m²
    scalar horizontalInfraredRadiation;    // Wh/m²
    scalar globalHorizontalRadiation;      // Wh/m²
    scalar directNormalRadiation;          // Wh/m²
    scalar diffuseHorizontalRadiation;     // Wh/m²
    scalar windDirection;      // degrees
    scalar windSpeed;          // m/s
    
    // Solar position (calculated)
    scalar solarAltitude;      // degrees
    scalar solarAzimuth;       // degrees
};

class epwReader
{
private:
    // EPW file data
    DynamicList<epwData> weatherData_;
    
    // Location data from EPW header
    scalar latitude_;
    scalar longitude_;
    scalar timeZone_;
    scalar elevation_;
    
    // Parse EPW line
    epwData parseLine(const string& line) const;
    
    // Calculate solar position
    void calculateSolarPosition(epwData& data) const;
    
public:
    // Constructor
    epwReader(const fileName& epwFile);
    
    // Destructor
    ~epwReader() = default;
    
    // Member functions
    
    // Get weather data for specific time
    epwData interpolateData(const scalar timeValue) const;
    
    // Get location data
    scalar latitude() const { return latitude_; }
    scalar longitude() const { return longitude_; }
    scalar timeZone() const { return timeZone_; }
    scalar elevation() const { return elevation_; }
    
    // Calculate mean radiant temperature from weather data
    scalar calculateTmrt
    (
        const epwData& weather,
        const scalar Ta,           // Air temperature [°C]
        const scalar SVF,          // Sky view factor [-]
        const scalar groundTemp    // Ground surface temperature [°C]
    ) const;
};

// Implementation
epwReader::epwReader(const fileName& epwFile)
{
    IFstream file(epwFile);
    
    if (!file.good())
    {
        FatalErrorInFunction
            << "Cannot read EPW file " << epwFile
            << exit(FatalError);
    }
    
    string line;
    label lineNo = 0;
    
    // Read header (8 lines)
    while (lineNo < 8 && file.getLine(line))
    {
        if (lineNo == 0)  // Location header
        {
            // Parse location data
            // Format: LOCATION,City,State,Country,DataSource,WMO,lat,lon,tz,elev
            IStringStream iss(line);
            DynamicList<string> tokens;
            
            while (iss.good())
            {
                string token;
                iss.getLine(token, ',');
                tokens.append(token);
            }
            
            if (tokens.size() >= 10)
            {
                latitude_ = readScalar(IStringStream(tokens[6])());
                longitude_ = readScalar(IStringStream(tokens[7])());
                timeZone_ = readScalar(IStringStream(tokens[8])());
                elevation_ = readScalar(IStringStream(tokens[9])());
            }
        }
        lineNo++;
    }
    
    // Read weather data
    while (file.getLine(line))
    {
        if (!line.empty())
        {
            epwData data = parseLine(line);
            calculateSolarPosition(data);
            weatherData_.append(data);
        }
    }
    
    Info<< "Read " << weatherData_.size() << " hours of weather data from " 
        << epwFile << endl;
    Info<< "Location: lat=" << latitude_ << " lon=" << longitude_ 
        << " tz=" << timeZone_ << " elev=" << elevation_ << endl;
}

epwData epwReader::parseLine(const string& line) const
{
    epwData data;
    IStringStream iss(line);
    DynamicList<string> tokens;
    
    // Split by comma
    while (iss.good())
    {
        string token;
        iss.getLine(token, ',');
        tokens.append(token);
    }
    
    // EPW format has specific column positions
    if (tokens.size() >= 31)
    {
        data.year = readScalar(IStringStream(tokens[0])());
        data.month = readScalar(IStringStream(tokens[1])());
        data.day = readScalar(IStringStream(tokens[2])());
        data.hour = readScalar(IStringStream(tokens[3])());
        data.minute = readScalar(IStringStream(tokens[4])());
        
        data.dryBulbTemp = readScalar(IStringStream(tokens[6])());
        data.dewPointTemp = readScalar(IStringStream(tokens[7])());
        data.relHumidity = readScalar(IStringStream(tokens[8])());
        data.atmPressure = readScalar(IStringStream(tokens[9])());
        
        data.extraterrestrialRadiation = readScalar(IStringStream(tokens[10])());
        data.extraterrestrialDirectRadiation = readScalar(IStringStream(tokens[11])());
        data.horizontalInfraredRadiation = readScalar(IStringStream(tokens[12])());
        data.globalHorizontalRadiation = readScalar(IStringStream(tokens[13])());
        data.directNormalRadiation = readScalar(IStringStream(tokens[14])());
        data.diffuseHorizontalRadiation = readScalar(IStringStream(tokens[15])());
        
        data.windDirection = readScalar(IStringStream(tokens[20])());
        data.windSpeed = readScalar(IStringStream(tokens[21])());
    }
    
    return data;
}

void epwReader::calculateSolarPosition(epwData& data) const
{
    // Simplified solar position calculation
    // For accurate UTCI, use proper solar position algorithms
    
    // Day of year
    scalar dayOfYear = data.day;  // Simplified - should account for month
    
    // Solar declination
    scalar declination = 23.45 * Foam::sin(2.0 * constant::mathematical::pi * 
                                          (284.0 + dayOfYear) / 365.0);
    
    // Hour angle
    scalar hourAngle = 15.0 * (data.hour - 12.0);
    
    // Solar altitude
    data.solarAltitude = Foam::asin(
        Foam::sin(declination * constant::mathematical::pi / 180.0) * 
        Foam::sin(latitude_ * constant::mathematical::pi / 180.0) + 
        Foam::cos(declination * constant::mathematical::pi / 180.0) * 
        Foam::cos(latitude_ * constant::mathematical::pi / 180.0) * 
        Foam::cos(hourAngle * constant::mathematical::pi / 180.0)
    ) * 180.0 / constant::mathematical::pi;
    
    // Solar azimuth (simplified)
    data.solarAzimuth = 0.0;  // Would need proper calculation
}

epwData epwReader::interpolateData(const scalar timeValue) const
{
    // For now, return closest hour
    // TODO: Implement proper interpolation
    
    label hour = label(timeValue / 3600.0) % 8760;  // Hour of year
    
    if (hour < weatherData_.size())
    {
        return weatherData_[hour];
    }
    else
    {
        return weatherData_[0];  // Default to first hour
    }
}

scalar epwReader::calculateTmrt
(
    const epwData& weather,
    const scalar Ta,
    const scalar SVF,
    const scalar groundTemp
) const
{
    // Calculate mean radiant temperature from weather data
    
    // Human body parameters
    const scalar EMISSIVITY_HUMAN = 0.97;
    const scalar ALBEDO_HUMAN = 0.37;
    const scalar STEFAN_BOLTZMANN = 5.670374419e-8;
    
    // Convert radiation from Wh/m² to W/m²
    scalar Sdir = weather.directNormalRadiation;
    scalar Sdiff = weather.diffuseHorizontalRadiation;
    scalar Ldown = weather.horizontalInfraredRadiation;
    
    // Simplified ground reflected radiation
    scalar Srefl = (Sdir + Sdiff) * 0.2;  // Assume 20% ground albedo
    
    // Longwave radiation from ground
    scalar Lup = EMISSIVITY_HUMAN * STEFAN_BOLTZMANN * 
                 Foam::pow(groundTemp + 273.15, 4);
    
    // Absorbed shortwave radiation
    scalar Sabs = (1 - ALBEDO_HUMAN) * (Sdir + Sdiff + Srefl);
    
    // Net longwave radiation
    scalar Lnet = EMISSIVITY_HUMAN * (SVF * Ldown + (1 - SVF) * Lup);
    
    // Total absorbed radiation
    scalar Rabs = Sabs + Lnet;
    
    // Calculate Tmrt
    scalar Tmrt_K = Foam::pow(Rabs / (EMISSIVITY_HUMAN * STEFAN_BOLTZMANN), 0.25);
    
    return Tmrt_K - 273.15;  // Return in °C
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //