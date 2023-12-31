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
    Foam::refineImmersedBoundaryMesh

Description
    Refine a mesh with immersed boundary

SourceFiles
    refineImmersedBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef refineImmersedBoundaryMesh_H
#define refineImmersedBoundaryMesh_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class refineImmersedBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class refineImmersedBoundaryMesh
{
public:

    // Public enumerations

        //- Cell collection method
        enum ibCellCollection
        {
            UNDEFINED_COLLECTION = 0,
            IB_CELLS,
            IB_CELL_CELLS
        };

        //- Projection method names
        static const NamedEnum<ibCellCollection, 3> ibCellCollectionNames_;


private:

    // Private data

        //- Reference to mesh
        polyMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        refineImmersedBoundaryMesh(const refineImmersedBoundaryMesh&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const refineImmersedBoundaryMesh&) = delete;


        // Mesh refinement

            //- Add ibCells for refinement
            void addIbCells(labelHashSet& refCellSet) const;

            //- Add ibCellCells for refinement
            void addIbCellCells(labelHashSet& refCellSet) const;

            //- Add ibCellCellFaces for refinement
            void addIbCellCellFaces(labelHashSet& refCellSet) const;


            //- From refineMesh application
            label axis(const vector& normal) const;

            //- From refineMesh application
            label twoDNess() const;


public:

    // Constructors

        //- Construct from mesh
        refineImmersedBoundaryMesh(polyMesh& mesh);


    //- Destructor
    ~refineImmersedBoundaryMesh() = default;


    // Member Functions

        //- Return list of cells for refinement based on specifie collection
        //  type (see enumeration)
        labelList refinementCells
        (
            const ibCellCollection& collectionType
        ) const;

        // Refine mesh given a list of cells for refinement
        void refineMesh(const labelList& refCells) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
