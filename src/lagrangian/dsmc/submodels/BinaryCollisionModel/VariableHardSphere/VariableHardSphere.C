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

#include "VariableHardSphere.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::VariableHardSphere<CloudType>::VariableHardSphere
(
    const dictionary& dict,
    CloudType& cloud
)
:
    BinaryCollisionModel<CloudType>(dict, cloud, typeName),
    Tref_(readScalar(this->coeffDict().lookup("Tref")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::VariableHardSphere<CloudType>::~VariableHardSphere()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template <class CloudType>
Foam::scalar Foam::VariableHardSphere<CloudType>::sigmaTcR
(
    label typeIdP,
    label typeIdQ,
    const vector& UP,
    const vector& UQ
) const
{
    const CloudType& cloud(this->owner());

    scalar dPQ =
        0.5
       *(
            cloud.constProps(typeIdP).d()
          + cloud.constProps(typeIdQ).d()
        );

    scalar omegaPQ =
        0.5
       *(
            cloud.constProps(typeIdP).omega()
          + cloud.constProps(typeIdQ).omega()
        );

    scalar cR = mag(UP - UQ);

    if (cR < VSMALL)
    {
        return 0;
    }

    scalar mP = cloud.constProps(typeIdP).mass();

    scalar mQ = cloud.constProps(typeIdQ).mass();

    scalar mR = mP*mQ/(mP + mQ);

    // calculating cross section = pi*dPQ^2, where dPQ is from Bird, eq. 4.79
    scalar sigmaTPQ =
        mathematicalConstant::pi*dPQ*dPQ
       *pow(2.0*CloudType::kb*Tref_/(mR*cR*cR), omegaPQ - 0.5)
       /exp(Foam::lgamma(2.5 - omegaPQ));

    return sigmaTPQ*cR;
}


template <class CloudType>
void Foam::VariableHardSphere<CloudType>::collide
(
    label typeIdP,
    label typeIdQ,
    vector& UP,
    vector& UQ,
    scalar& EiP,
    scalar& EiQ
)
{
    CloudType& cloud(this->owner());

    Random& rndGen(cloud.rndGen());

    scalar mP = cloud.constProps(typeIdP).mass();

    scalar mQ = cloud.constProps(typeIdQ).mass();

    vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

    scalar cR = mag(UP - UQ);

    scalar cosTheta = 2.0*rndGen.scalar01() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = 2.0*mathematicalConstant::pi*rndGen.scalar01();

    vector postCollisionRelU =
        cR
       *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );

    UP = Ucm + postCollisionRelU*mQ/(mP + mQ);

    UQ = Ucm - postCollisionRelU*mP/(mP + mQ);
}


// ************************************************************************* //
