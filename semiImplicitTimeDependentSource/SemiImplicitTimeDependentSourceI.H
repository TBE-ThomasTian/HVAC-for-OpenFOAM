/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "SemiImplicitTimeDependentSource.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline const typename Foam::fv::SemiImplicitTimeDependentSource<Type>::volumeModeType&
Foam::fv::SemiImplicitTimeDependentSource<Type>::volumeMode() const
{
    return volumeMode_;
}


template<class Type>
inline const Foam::List<Foam::Tuple2<Type, Foam::scalar>>&
Foam::fv::SemiImplicitTimeDependentSource<Type>::injectionRate() const
{
    return injectionRate_;
}


template<class Type>
inline typename Foam::fv::SemiImplicitTimeDependentSource<Type>::volumeModeType&
Foam::fv::SemiImplicitTimeDependentSource<Type>::volumeMode()
{
    return volumeMode_;
}


template<class Type>
inline Foam::List<Foam::Tuple2<Type,
Foam::scalar>>& Foam::fv::SemiImplicitTimeDependentSource<Type>::injectionRate()
{
    return injectionRate_;
}


// ************************************************************************* //
