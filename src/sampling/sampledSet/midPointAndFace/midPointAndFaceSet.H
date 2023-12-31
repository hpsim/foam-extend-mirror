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
    Foam::midPointAndFaceSet

Description

SourceFiles
    midPointAndFaceSet.C

\*---------------------------------------------------------------------------*/

#ifndef midPointAndFaceSet_H
#define midPointAndFaceSet_H

#include "faceOnlySet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class passiveParticle;
template<class Type> class particle;
class meshSearch;

/*---------------------------------------------------------------------------*\
                           Class midPointAndFaceSet Declaration
\*---------------------------------------------------------------------------*/

class midPointAndFaceSet
:
    public faceOnlySet
{

    // Private Member Functions

        void genSamples();

public:

    //- Runtime type information
    TypeName("midPointAndFace");


    // Constructors

        //- Construct from components
        midPointAndFaceSet
        (
            const word& name,
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const word& axis,
            const point& start,
            const point& end
        );

        //- Construct from dictionary
        midPointAndFaceSet
        (
            const word& name,
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const dictionary& dict
        );


    // Destructor

        virtual ~midPointAndFaceSet();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
