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
    Foam::LocalInteraction

Description
    Patch interaction specified on a patch-by-patch basis

\*---------------------------------------------------------------------------*/

#ifndef LocalInteraction_H
#define LocalInteraction_H

#include "PatchInteractionModel.H"
#include "patchInteractionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class LocalInteraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LocalInteraction
:
    public PatchInteractionModel<CloudType>
{
    // Private data

        //- List of participating patches
        const List<patchInteractionData> patchData_;

        //- List of participating patch ids
        labelList patchIds_;


    // Private member functions

        //- Returns local patchI if patch is in patchIds_ list
        label applyToPatch(const label globalPatchI) const;



public:

    //- Runtime type information
    TypeName("LocalInteraction");


    // Constructors

        //- Construct from dictionary
        LocalInteraction(const dictionary& dict, CloudType& cloud);


    //- Destructor
    virtual ~LocalInteraction();


    // Member Functions

        //- Flag to indicate whether model activates patch interaction model
        virtual bool active() const;

        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            const polyPatch& pp,
            const label faceId,
            bool& keepParticle,
            vector& U
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LocalInteraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
