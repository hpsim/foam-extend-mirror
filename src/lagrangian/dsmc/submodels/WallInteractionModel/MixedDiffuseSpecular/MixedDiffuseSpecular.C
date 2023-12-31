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

#include "MixedDiffuseSpecular.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::MixedDiffuseSpecular<CloudType>::MixedDiffuseSpecular
(
    const dictionary& dict,
    CloudType& cloud
)
:
    WallInteractionModel<CloudType>(dict, cloud, typeName),
    diffuseFraction_(readScalar(this->coeffDict().lookup("diffuseFraction")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::MixedDiffuseSpecular<CloudType>::~MixedDiffuseSpecular()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::MixedDiffuseSpecular<CloudType>::correct
(
    const wallPolyPatch& wpp,
    const label faceId,
    vector& U,
    scalar& Ei,
    label typeId
)
{
    label wppIndex = wpp.index();

    label wppLocalFace = wpp.whichFace(faceId);

    vector nw = wpp.faceAreas()[wppLocalFace];

    // Normal unit vector
    nw /= mag(nw);

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    CloudType& cloud(this->owner());

    Random& rndGen(cloud.rndGen());

    if (diffuseFraction_ > rndGen.scalar01())
    {
        // Diffuse reflection

        // Wall tangential velocity (flow direction)
        vector Ut = U - U_dot_nw*nw;

        while (mag(Ut) < SMALL)
        {
            // If the incident velocity is parallel to the face normal, no
            // tangential direction can be chosen.  Add a perturbation to the
            // incoming velocity and recalculate.

            U = vector
            (
                U.x()*(0.8 + 0.2*rndGen.scalar01()),
                U.y()*(0.8 + 0.2*rndGen.scalar01()),
                U.z()*(0.8 + 0.2*rndGen.scalar01())
            );

            U_dot_nw = U & nw;

            Ut = U - U_dot_nw*nw;
        }

        // Wall tangential unit vector
        vector tw1 = Ut/mag(Ut);

        // Other tangential unit vector
        vector tw2 = nw^tw1;

        scalar T = cloud.boundaryT().boundaryField()[wppIndex][wppLocalFace];

        scalar mass = cloud.constProps(typeId).mass();

        scalar iDof = cloud.constProps(typeId).internalDegreesOfFreedom();

        U =
            sqrt(CloudType::kb*T/mass)
           *(
                rndGen.GaussNormal()*tw1
              + rndGen.GaussNormal()*tw2
              - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
            );

        U += cloud.boundaryU().boundaryField()[wppIndex][wppLocalFace];

        Ei = cloud.equipartitionInternalEnergy(T, iDof);
    }
    else
    {
        // Specular reflection

        if (U_dot_nw > 0.0)
        {
            U -= 2.0*U_dot_nw*nw;
        }
    }

}


// ************************************************************************* //
