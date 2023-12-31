/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "NoInjection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::NoInjection<CloudType>::parcelsToInject
(
    const scalar,
    const scalar
) const
{
    return 0;
}


template<class CloudType>
Foam::scalar Foam::NoInjection<CloudType>::volumeToInject
(
    const scalar,
    const scalar
) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoInjection<CloudType>::NoInjection
(
    const dictionary&,
    CloudType& owner
)
:
    InjectionModel<CloudType>(owner)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoInjection<CloudType>::~NoInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::NoInjection<CloudType>::active() const
{
    return false;
}


template<class CloudType>
Foam::scalar Foam::NoInjection<CloudType>::timeEnd() const
{
    return 0.0;
}


template<class CloudType>
void Foam::NoInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector&,
    label&
)
{}


template<class CloudType>
void Foam::NoInjection<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    parcel.U() = vector::zero;

    // set particle diameter
    parcel.d() = 0.0;
}


template<class CloudType>
bool Foam::NoInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::NoInjection<CloudType>::validInjection(const label)
{
    return false;
}


// ************************************************************************* //
