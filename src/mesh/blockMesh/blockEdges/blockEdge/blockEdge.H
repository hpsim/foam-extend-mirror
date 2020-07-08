/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::blockEdge

Description
    Define a curved edge that is parameterized for 0<lambda<1
    between the start and end point.

SourceFiles
    blockEdge.C

Author
    Gaavin Tabor.
    Rewrite by Hrvoje Jasak.

\*---------------------------------------------------------------------------*/

#ifndef blockEdge_H
#define blockEdge_H

#include "edge.H"
#include "pointField.H"
#include "typeInfo.H"
#include "HashTable.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class blockEdge;

Ostream& operator<<(Ostream&, const blockEdge&);


/*---------------------------------------------------------------------------*\
                         Class blockEdge Declaration
\*---------------------------------------------------------------------------*/

class blockEdge
{
protected:

    // Protected data

        //- Reference to point list
        const pointField& points_;

        //- Start point label
        const label start_;

        //- End point label
        const label end_;


    // Protected Member Functions

        //- Return a complete point field by appending the start/end points
        //  to the given list
        static tmp<pointField> appendEndPoints
        (
            const pointField&,
            const label start,
            const label end,
            const pointField& otherKnots
        );


    // Private member functions

        //- Disallow construct as copy
        blockEdge(const blockEdge&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const blockEdge&) = delete;


public:

    //- Runtime type information
    TypeName("blockEdge");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            blockEdge,
            Istream,
            (
                const pointField& points,
                Istream& is
            ),
            (points, is)
        );


    // Constructors

        //- Construct from components
        blockEdge
        (
            const pointField& points,
            const label start,
            const label end
        );

        //- Construct from Istream setting pointsList
        blockEdge(const pointField&, Istream&);

        //- Copy construct
        blockEdge(const blockEdge&);

        //- Clone function
        virtual autoPtr<blockEdge> clone() const = 0;

        //- New function which constructs and returns pointer to a blockEdge
        static autoPtr<blockEdge> New(const pointField&, Istream&);


    //- Destructor
    virtual ~blockEdge() = default;


    // Member Functions

        //- Return label of start point
        inline label start() const;

        //- Return label of end point
        inline label end() const;

        //- Compare the given start and end points with this curve
        //  Return:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        inline int compare(const blockEdge&) const;

        //- Compare the given start and end points with this curve
        //  Return:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        inline int compare(const edge&) const;

        //- Compare the given start and end points with this curve
        //  Return:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        inline int compare(const label start, const label end) const;

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        virtual point position(const scalar) const = 0;

        //- Return the length of the curve
        virtual scalar length() const = 0;


    // Member operators

        void operator=(const blockEdge&);


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const blockEdge&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blockEdgeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
