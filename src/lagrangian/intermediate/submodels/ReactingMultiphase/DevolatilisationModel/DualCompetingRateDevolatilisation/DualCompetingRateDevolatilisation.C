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

#include "DualCompetingRateDevolatilisation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DualCompetingRateDevolatilisation<CloudType>::
DualCompetingRateDevolatilisation
(
    const dictionary& dict,
    CloudType& owner
)
:
    DevolatilisationModel<CloudType>(dict, owner, typeName),
    a1_(readScalar(this->coeffDict().lookup("a1"))),
    A1_(readScalar(this->coeffDict().lookup("A1"))),
    E1_(readScalar(this->coeffDict().lookup("E1"))),
    a2_(readScalar(this->coeffDict().lookup("a2"))),
    A2_(readScalar(this->coeffDict().lookup("A2"))),
    E2_(readScalar(this->coeffDict().lookup("E2"))),
    volatileResidualCoeff_
    (
        readScalar(this->coeffDict().lookup("volatileResidualCoeff"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::DualCompetingRateDevolatilisation<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::DualCompetingRateDevolatilisation<CloudType>::calculate
(
    const scalar dt,
    const scalar mass0,
    const scalar mass,
    const scalar T,
    const scalar YVolatile0,
    const scalar YVolatile,
    bool& canCombust
) const
{
    const scalar massVolatile0 = YVolatile0*mass0;
    const scalar massVolatile  = YVolatile*mass;

    if (massVolatile <= volatileResidualCoeff_*massVolatile0)
    {
        canCombust = true;
    }

    // Kinetic Rates
    const scalar k1 = A1_*exp(-E1_/(specie::RR()*T));
    const scalar k2 = A2_*exp(-E2_/(specie::RR()*T));

    // Volatile devolatilisation from particle to carrier gas phase
    const scalar dMass = min((a1_*k1 + a2_*k2)*massVolatile*dt, massVolatile);

    return dMass;
}


// ************************************************************************* //
