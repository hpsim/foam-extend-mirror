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
    Foam::directInteractionList

Description

SourceFiles
    directInteractionListI.H
    directInteractionList.C

\*---------------------------------------------------------------------------*/

#ifndef directInteractionList_H
#define directInteractionList_H

#include "polyMesh.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class interactionLists;

/*---------------------------------------------------------------------------*\
                           Class directInteractionList Declaration
\*---------------------------------------------------------------------------*/

class directInteractionList
:
    public labelListList
{
    // Private data

        const interactionLists& il_;


    // Private Member Functions

        void buildDirectInteractionList
        (
            bool pointPointListBuild
        );

        //- Disallow default bitwise copy construct
        directInteractionList(const directInteractionList&);

        //- Disallow default bitwise assignment
        void operator=(const directInteractionList&);


public:

    // Constructors

        //- Construct lists by searching the mesh
        directInteractionList
        (
            const interactionLists& il,
            bool pointPointListBuild
        );

        //- Construct from file
        directInteractionList
        (
            const interactionLists& il
        );


    // Destructor

        ~directInteractionList();


    // Member Functions

        // Access

        inline const interactionLists& il() const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, directInteractionList&);
        friend Ostream& operator<<(Ostream&, const directInteractionList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "directInteractionListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
