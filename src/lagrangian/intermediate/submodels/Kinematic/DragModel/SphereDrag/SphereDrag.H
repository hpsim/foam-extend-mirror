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
    Foam::SphereDrag

Description
    Drag model based on assumption of solid spheres

\*---------------------------------------------------------------------------*/

#ifndef SphereDrag_H
#define SphereDrag_H

#include "DragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class SphereDrag Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SphereDrag
:
    public DragModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("SphereDrag");


    // Constructors

        //- Construct from dictionary
        SphereDrag
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~SphereDrag();


    // Member Functions

        //- Flag to indicate whether model activates drag model
        bool active() const;

        //- Return drag coefficient
        scalar Cd(const scalar Re) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SphereDrag.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
