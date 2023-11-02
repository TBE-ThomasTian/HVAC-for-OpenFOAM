/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::fv::SemiImplicitTimeDependentSource

Group
    grpFvOptionsSources

Description
    Semi-implicit source, described using an input dictionary.  The injection
    rate coefficients are specified as pairs of Su-Sp coefficients, i.e.

        \f[
            S(x) = S_u + S_p x
        \f]

    where
    \vartable
        S(x)    | net source for field 'x'
        S_u     | explicit source contribution
        S_p     | linearised implicit contribution
    \endvartable

    Example of the source specification:

    \verbatim
    volumeMode      absolute; 
    injectionRateSuSp
    {
        h           (1 0);
    }
    \endverbatim

    Valid options for the \c volumeMode entry include:
    - absolute: values are given as \<quantity\> - JN: always use this
    - specific: values are given as \<quantity\>/m3 - JN: never use this
	
	JN: We multiply the explicit source S_u from fvOptions with HRR
	JN: Always use 1 for S_u and 0 for S_p if you define a time dependent HRR
	JN: If you want to model the soot as well use a constant value instead of one
	JN: You calcaulate it out of the soot yield and the heat of combustion 
	JN: i.e. Yield_{soot}/Heat_{combution} ---> (kg/kg) / (J/kg) = kg/J 

See also
    Foam::fvOption

SourceFiles
    SemiImplicitTimeDependentSource.C

\*---------------------------------------------------------------------------*/

#ifndef SemiImplicitTimeDependentSource_H
#define SemiImplicitTimeDependentSource_H

#include "Tuple2.H"
#include "cellSetOption.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// Forward declarations
template<class Type> class SemiImplicitTimeDependentSource;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const SemiImplicitTimeDependentSource<Type>&
);


/*---------------------------------------------------------------------------*\
                     Class SemiImplicitTimeDependentSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SemiImplicitTimeDependentSource
:
    public cellSetOption
{
public:

    // Public data

        //- Enumeration for volume types
        enum volumeModeType
        {
            vmAbsolute,
            vmSpecific
        };

        //- Word list of volume mode type names
        static const wordList volumeModeTypeNames_;


protected:

    // Protected data

        //- Volume mode
        volumeModeType volumeMode_;

        //- Volume normalisation
        scalar VDash_;

        //- Source field values
        List<Tuple2<Type, scalar>> injectionRate_;
		
		//- JN: HRR table to be read in
		interpolationTable<scalar> timeDependence_;


    // Protected functions

        //- Helper function to convert from a word to a volumeModeType
        volumeModeType wordToVolumeModeType(const word& vtName) const;

        //- Helper function to convert from a volumeModeType to a word
        word volumeModeTypeToWord(const volumeModeType& vtType) const;

        //- Set the local field data
        void setFieldData(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("SemiImplicitTimeDependentSource");


    // Constructors

        //- Construct from components
        SemiImplicitTimeDependentSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Access

            //- Return const access to the volume mode
            inline const volumeModeType& volumeMode() const;

            //- Return const access to the source field values
            inline const List<Tuple2<Type, scalar>>& injectionRate() const;


        // Edit

            //- Return access to the volume mode
            inline volumeModeType& volumeMode();

            //- Return access to the source field values
            inline List<Tuple2<Type, scalar>>& injectionRate();


        // Evaluation

            //- Add explicit contribution to equation
            virtual void addSup
            (
                fvMatrix<Type>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<Type>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SemiImplicitTimeDependentSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SemiImplicitTimeDependentSourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
