/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd
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

\*---------------------------------------------------------------------------*/

#include "solarCalculator.H"
#include "Time.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solarCalculator, 0);

    template<>
    const char* NamedEnum
    <
        solarCalculator::sunDirModel,
        2
    >::names[] =
    {
        "sunDirConstant",
        "sunDirTraking"
    };

    template<>
    const char* NamedEnum
    <
        solarCalculator::sunLModel,
        3
    >::names[] =
    {
        "sunLoadConstant",
        "sunLoadFairWeatherConditions",
        "sunLoadTheoreticalMaximum"
    };
}

const Foam::NamedEnum<Foam::solarCalculator::sunDirModel, 2>
  Foam::solarCalculator::sunDirectionModelTypeNames_;

const Foam::NamedEnum<Foam::solarCalculator::sunLModel, 3>
   Foam::solarCalculator::sunLoadModelTypeNames_;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solarCalculator::calculateBetaTetha()
{
    scalar runTime = 0.0;
    switch (sunDirectionModel_)
    {
        case mSunDirTraking:
        {
            runTime = mesh_.time().value();
            break;
        }
        case mSunDirConstant:
        {
            break;
        }
    }
/*
    scalar LSM = 15.0*(readScalar(dict_.lookup("localStandardMeridian")));

    scalar D = readScalar(dict_.lookup("startDay")) + runTime/86400.0;
    scalar M = 6.24004 + 0.0172*D;
    scalar EOT = -7.659*sin(M) + 9.863*sin(2*M + 3.5932);

    startTime_ = readScalar(dict_.lookup("startTime"));
    scalar LST =  startTime_ + runTime/3600.0;

    scalar LON = readScalar(dict_.lookup("longitude"));

    scalar AST = LST + EOT/60.0 + (LON - LSM)/15;

    scalar delta = 23.45*sin(degToRad((360*(284 + D))/365));

    scalar H = degToRad(15*(AST - 12));

    scalar L = degToRad(readScalar(dict_.lookup("latitude")));

    scalar deltaRad = degToRad(delta);
    beta_ = max(asin(cos(L)*cos(deltaRad)*cos(H) + sin(L)*sin(deltaRad)), 1e-3);
    tetha_ = acos((sin(beta_)*sin(L) - sin(deltaRad))/(cos(beta_)*cos(L)));
*/

	    scalar Day = readScalar(dict_.lookup("startDay")) + runTime/86400.0;
	    scalar Month = readScalar(dict_.lookup("startMonth"));
	    scalar Year =  readScalar(dict_.lookup("startYear"));
            scalar LON = readScalar(dict_.lookup("longitude"));
            scalar L   = readScalar(dict_.lookup("latitude"));
	    scalar LST = startTime_ + runTime/3600.0;
            startTime_ = readScalar(dict_.lookup("startTime"));


		Info << " Run Time : " << runTime << endl;
		Info << "LST : " << LST << endl;

            //System.DateTime.ToOADate()
            //jd=1461*(y+4800+(m-14)/12)/4+367*(m-2-12*(m-14))

            // Anzahl der Tage ab dem Jahr 2000 JD0:
            double julianDate = (367.0 * Year -
                ((7.0 / 4.0) * (Year +
                ((Month + 9.0) / 12.0))) +
                ((275.0 * Month) / 9.0) +
                Day - 730531.5);

            // Ein Julian-Jahrhundert seit Jahr 2000 ist JD/36525 Tage T:
            double julianCenturies = julianDate / 36525.0;

	if (debug)
	Info << " julianCenturies : " << julianCenturies << endl;

            // Sternzeit T0: (+ T^2*0,000025862?)
            double siderealTimeHours = 6.697374558 + 2400.051336 * julianCenturies;

            // Ein Sterntag ist etwa 1/365, Länge des Jahres 365,2422 Tage:
            // H * 1.00273790935
            double siderealTimeUT = siderealTimeHours +
                (366.242190402 / 365.242190402) * LST;

	if (debug)
	Info << "  siderealTimeUT : " <<  siderealTimeUT << endl;

            // Frühlingspunkt wandert um 1 Stunde (15°, west of meridian)
            // LST:
            double siderealTime = siderealTimeUT * 15.0 + LON;

	if (debug)
	Info << " siderealTime : " << siderealTime << endl;

            // Refine to number of days (fractional) to specific time.
            julianDate += LST / 24.0;
            julianCenturies = julianDate / 36525.0;

            // Solar Coordinates L (Mean Longtitude of the Sun):
            double  meanLongitude = correctAngle ( degToRad(280.466 + 36000.77 * julianCenturies) );

            double meanAnomaly = correctAngle ( degToRad(357.529 + 35999.05 * julianCenturies) );

	if (debug)
	Info << "meanAnomaly : " << meanAnomaly << endl;

            double equationOfCenter = degToRad( (1.915 - 0.005 * julianCenturies) *
                sin(meanAnomaly) + 0.02 * sin(2.0 * meanAnomaly) );

	if (debug)
	Info << "equationOfCenter : " << equationOfCenter << endl;

            double elipticalLongitude =
               correctAngle(meanLongitude + equationOfCenter);

	if (debug)
	Info << " elipticalLongitude : " << elipticalLongitude << endl;

            double obliquity = degToRad( 23.439 - 0.013 * julianCenturies);

	if (debug)
	Info << " obliquity : " << obliquity << endl;


            double rightAscension = atan2( cos(obliquity) * sin(elipticalLongitude), cos(elipticalLongitude));

	if (debug)
	Info << " Right Ascension " <<  rightAscension << endl;

            double declinationAngle = asin( sin(rightAscension) * sin(obliquity));

	if (debug)
	 Info << " Declination Angle : " << declinationAngle << endl;;

            double hourAngle = correctAngle( degToRad( siderealTime) ) - rightAscension;

	if (debug)
	Info << " hourAngle 2 : " << hourAngle  << endl;

            if (hourAngle > constant::mathematical::pi)
            {
                hourAngle -= 2.0 * constant::mathematical::pi;
            }

            beta_ = max( 
		asin( sin( degToRad(L) ) *
                sin(declinationAngle) + cos( degToRad(L) ) *
                cos(declinationAngle) * cos(hourAngle))
		, 1e-3);

            double aziNom = -sin( hourAngle );
            double aziDenom =
                tan(declinationAngle) * cos( degToRad(L) ) -
                sin( degToRad(L) ) * cos(hourAngle);

            tetha_ = atan(aziNom / aziDenom);      //source http://www.jgiesen.de/SunCompass/azimuth/

            if (aziDenom < 0) // In 2nd or 3rd quadrant
            {
                tetha_ += constant::mathematical::pi;
            }
            else if (aziNom < 0) // In 4th quadrant
            {
                tetha_ += 2.0 * constant::mathematical::pi;
            }

	if (debug)
	{
	Info << tab << "hourAngle : " << radToDeg( hourAngle ) << endl;
        Info << tab << "altitude : " << radToDeg(beta_) << endl;
        Info << tab << "azimuth  : " << radToDeg(tetha_) << endl;
	}
}


void Foam::solarCalculator::calculateSunDirection()
{

    dict_.lookup("gridUp") >> gridUp_;
    gridUp_ /= mag(gridUp_);

    dict_.lookup("gridEast") >> eastDir_;
    eastDir_ /= mag(eastDir_);

    coord_.reset
    (
        new coordinateSystem("grid", vector::zero, gridUp_, eastDir_)
    );

    direction_.z() = -sin(beta_);
    direction_.y() =  cos(beta_)*cos(tetha_); //North
    direction_.x() =  cos(beta_)*sin(tetha_); //East

    direction_ /= mag(direction_);

    if (debug)
    {
        Info<< "Sun direction in absolute coordinates : " << direction_ <<endl;
    }

    direction_ = coord_->R().transform(direction_);

    if (debug)
    {
        Info<< "Sun direction in the Grid coordinates : " << direction_ <<endl;
    }
}


void Foam::solarCalculator::init()
{
    switch (sunDirectionModel_)
    {
        case mSunDirConstant:
        {
            if (dict_.found("sunDirection"))
            {
                dict_.lookup("sunDirection") >> direction_;
                direction_ /= mag(direction_);
            }
            else
            {
                calculateBetaTetha();
                calculateSunDirection();
            }

            break;
        }
        case mSunDirTraking:
        {
            if (word(mesh_.ddtScheme("default")) == "steadyState")
            {
                FatalErrorInFunction
                    << " Sun direction model can not be sunDirtracking if the "
                    << " case is steady " << nl << exit(FatalError);
            }

            dict_.lookup("sunTrackingUpdateInterval") >>
                sunTrackingUpdateInterval_;

            calculateBetaTetha();
            calculateSunDirection();
            break;
        }
    }

    switch (sunLoadModel_)
    {
        case mSunLoadConstant:
        {
            dict_.lookup("directSolarRad") >> directSolarRad_;
            dict_.lookup("diffuseSolarRad") >> diffuseSolarRad_;
            break;
        }
        case mSunLoadFairWeatherConditions:
        {
            A_ = readScalar(dict_.lookup("A"));
            B_ = readScalar(dict_.lookup("B"));

            if (dict_.found("beta"))
            {
                dict_.lookup("beta") >> beta_;
            }
            else
            {
                calculateBetaTetha();
            }

            directSolarRad_ = A_/exp(B_/sin(beta_));

            groundReflectivity_ =
                readScalar(dict_.lookup("groundReflectivity"));

            break;
        }
        case mSunLoadTheoreticalMaximum:
        {
            Setrn_ = readScalar(dict_.lookup("Setrn"));
            SunPrime_ = readScalar(dict_.lookup("SunPrime"));
            directSolarRad_ = Setrn_*SunPrime_;

            groundReflectivity_ =
                readScalar(dict_.lookup("groundReflectivity"));
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarCalculator::solarCalculator
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_(dict),
    direction_(vector::zero),
    directSolarRad_(0.0),
    diffuseSolarRad_(0.0),
    groundReflectivity_(0.0),
    A_(0.0),
    B_(0.0),
    beta_(0.0),
    tetha_(0.0),
    Setrn_(0.0),
    SunPrime_(0.0),
    C_(readScalar(dict.lookup("C"))),
    sunDirectionModel_
    (
        sunDirectionModelTypeNames_.read(dict.lookup("sunDirectionModel"))
    ),
    sunLoadModel_
    (
        sunLoadModelTypeNames_.read(dict.lookup("sunLoadModel"))
    ),
    coord_()
{
    init();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solarCalculator::~solarCalculator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solarCalculator::correctSunDirection()
{
    switch (sunDirectionModel_)
    {
        case mSunDirConstant:
        {
            break;
        }
        case mSunDirTraking:
        {
            calculateBetaTetha();
            calculateSunDirection();
            directSolarRad_ = A_/exp(B_/sin(max(beta_, ROOTVSMALL)));
            break;
        }
    }
}

double Foam::solarCalculator::correctAngle(double angleInRadians)
{
Info << " angleInRadians : " << angleInRadians << endl;


            if (angleInRadians < 0)
            {
                // return 2.0 * constant::mathematical::pi - ( abs(angleInRadians) % int (2.0 * constant::mathematical::pi ) );
                int temp = int( abs(angleInRadians) / (2.0 * constant::mathematical::pi) );
                double temp1 = temp * (2.0 * constant::mathematical::pi);
                return 2.0 *  constant::mathematical::pi - angleInRadians - temp1;

            }
            else if (angleInRadians > (2.0 * constant::mathematical::pi) )
            {
                // return int (angleInRadians) %  int(2.0 * constant::mathematical::pi) );
		int temp = int( angleInRadians / (2.0 * constant::mathematical::pi) ); 
		double temp1 = temp * (2.0 * constant::mathematical::pi);
		return angleInRadians - temp1;
            }
            else
            {
                return angleInRadians;
            }
}


// ************************************************************************* //
