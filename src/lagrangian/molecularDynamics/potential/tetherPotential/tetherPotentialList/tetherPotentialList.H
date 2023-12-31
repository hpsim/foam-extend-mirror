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
    Foam::tetherPotentialList

Description

SourceFiles
    tetherPotentialList.C

\*---------------------------------------------------------------------------*/

#ifndef tetherPotentialList_H
#define tetherPotentialList_H

#include "ListOps.H"
#include "PtrList.H"
#include "word.H"
#include "tetherPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class tetherPotentialList Declaration
\*---------------------------------------------------------------------------*/

class tetherPotentialList
:
    public PtrList<tetherPotential>
{
    // Private data

        labelList idMap_;


    // Private Member Functions

        inline label tetherPotentialIndex
        (
            const label a
        ) const;

        void readTetherPotentialDict
        (
            const wordList& siteIdList,
            const dictionary& tetherPotentialDict,
            const wordList& tetherSiteIdList
        );

        //- Disallow default bitwise assignment
        void operator=(const tetherPotentialList&);

        //- Disallow default bitwise copy construct
        tetherPotentialList(const tetherPotentialList&);


public:

    // Constructors

        tetherPotentialList();

        //- Construct from siteIdList and potental dictionaries
        tetherPotentialList
        (
            const wordList& siteIdList,
            const dictionary& tetherPotentialDict,
            const wordList& tetherSiteIdList
        );


    // Destructor

        ~tetherPotentialList();


    // Member Functions

        void buildPotentials
        (
            const wordList& siteIdList,
            const dictionary& tetherPotentialDict,
            const wordList& tetherSiteIdList
        );

        // Access

            inline const wordList& idMap() const;

            const tetherPotential& tetherPotentialFunction(const label a) const;

            vector force(const label a, const vector rIT) const;

            scalar energy (const label a, const vector rIT) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetherPotentialListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
