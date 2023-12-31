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
    Foam::potential

Description

SourceFiles
    potentialI.H
    potential.C

\*---------------------------------------------------------------------------*/

#ifndef potential_H
#define potential_H

#include "polyMesh.H"
#include "IOdictionary.H"
#include "foamTime.H"
#include "pairPotentialList.H"
#include "electrostaticPotential.H"
#include "tetherPotentialList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class potential Declaration
\*---------------------------------------------------------------------------*/

class potential
{
    // Private data

        const polyMesh& mesh_;

        wordList idList_;

        wordList siteIdList_;

        label nPairPotIds_;

        scalar potentialEnergyLimit_;

        labelList removalOrder_;

        pairPotentialList pairPotentials_;

        tetherPotentialList tetherPotentials_;

        vector gravity_;


    // Private Member Functions

        void setSiteIdList(const IOdictionary& moleculePropertiesDict);

        void readPotentialDict();

        void readMdInitialiseDict
        (
            const IOdictionary& mdInitialiseDict,
            IOdictionary& idListDict
        );

        //- Disallow default bitwise copy construct
        potential(const potential&);

        //- Disallow default bitwise assignment
        void operator=(const potential&);


public:

    // Constructors

        //- Construct from mesh reference
        potential(const polyMesh& mesh);

        //- Construct from mdInitialiseDict
        potential
        (
            const polyMesh& mesh,
            const IOdictionary& mdInitialiseDict,
            IOdictionary& idListDict
        );


    // Destructor

        ~potential();


    // Member Functions

        // Access

            inline label nIds() const;

            inline const wordList& idList() const;

            inline const wordList& siteIdList() const;

            inline scalar potentialEnergyLimit() const;

            inline label nPairPotentials() const;

            inline const labelList& removalOrder() const;

            inline const pairPotentialList& pairPotentials() const;

            inline const tetherPotentialList& tetherPotentials() const;

            inline const vector& gravity() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "potentialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
