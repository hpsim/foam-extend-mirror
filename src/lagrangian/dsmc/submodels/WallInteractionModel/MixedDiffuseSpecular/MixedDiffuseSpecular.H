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

Class
    Foam::MixedDiffuseSpecular

Description
    Wall interaction setting microscopic velocity to a random one drawn from a
    Maxwellian distribution corresponding to a specified temperature

\*---------------------------------------------------------------------------*/

#ifndef MixedDiffuseSpecular_H
#define MixedDiffuseSpecular_H

#include "WallInteractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class MixedDiffuseSpecular Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MixedDiffuseSpecular
:
    public WallInteractionModel<CloudType>
{
    // Private data

        //- Fraction of wall interactions that are diffuse
        scalar diffuseFraction_;


public:

    //- Runtime type information
    TypeName("MixedDiffuseSpecular");


    // Constructors

        //- Construct from dictionary
        MixedDiffuseSpecular
        (
            const dictionary& dict,
            CloudType& cloud
        );


    // Destructor
    virtual ~MixedDiffuseSpecular();


    // Member Functions

        //- Apply wall correction
        virtual void correct
        (
            const wallPolyPatch& wpp,
            const label faceId,
            vector& U,
            scalar& Ei,
            label typeId
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MixedDiffuseSpecular.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
