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
    Foam::block

Description
    Creates a single block of cells from point coordinates, numbers of
    cells in each direction and an expansion ratio.

Note
    The vertices and cells for filling the block are demand-driven.

SourceFiles
    block.C
    blockCreate.C

\*---------------------------------------------------------------------------*/

#ifndef block_H
#define block_H

#include "pointField.H"
#include "labelList.H"

#include "blockDescriptor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of class
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class block;

Ostream& operator<<(Ostream&, const block&);


/*---------------------------------------------------------------------------*\
                           Class block Declaration
\*---------------------------------------------------------------------------*/

class block
:
    public blockDescriptor
{
    // Private data

        //- List of vertices
        mutable pointField vertices_;

        //- List of cells
        mutable labelListList cells_;

        //- Boundary patches
        mutable labelListListList boundaryPatches_;


    // Private Member Functions

        //- Creates vertices for cells filling the block
        void createPoints() const;

        //- Creates cells for filling the block
        void createCells() const;

        //- Creates boundary patch faces for the block
        void createBoundary() const;

        //- Disallow default bitwise copy construct
        block(const block&);

        //- Disallow default bitwise assignment
        void operator=(const block&);

public:

    // Constructors

        //- Construct from components with Istream
        block
        (
            const pointField& blockPointField,
            const curvedEdgeList&,
            Istream&
        );

        //- Construct from a block definition
        block(const blockDescriptor&);

        //- Clone
        autoPtr<block> clone() const
        {
            notImplemented("block::clone()");
            return autoPtr<block>(nullptr);
        }


    //- Destructor
    ~block();


    // Member Functions

        // Access

            //- Return the block definition
            inline const blockDescriptor& blockDef() const
            {
                return *this;
            }

            //- Vertex label offset for a particular i,j,k position
            label vtxLabel(label i, label j, label k) const;

            //- Return the points for filling the block
            const pointField& points() const;

            //- Return the cells for filling the block
            const labelListList& cells() const;

            //- Return the boundary patch faces for the block
            const labelListListList& boundaryPatches() const;


        // Edit

            //- Clear geometry (internal points, cells, boundaryPatches)
            void clearGeom();


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const block&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
