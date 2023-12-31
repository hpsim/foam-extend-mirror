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
    Foam::polyDualMesh

Description
    Creates dual of polyMesh.

    Every polyMesh point becomes a cell on the dual mesh

    Every polyMesh cell and patchFace becomes a point on the dual mesh.


SourceFiles
    polyDualMesh.C

\*---------------------------------------------------------------------------*/

#ifndef polyDualMesh_H
#define polyDualMesh_H

#include "polyMesh.H"
#include "labelIOList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyDualMesh Declaration
\*---------------------------------------------------------------------------*/

class polyDualMesh
:
    public polyMesh
{
    // Private data

        //- From polyMesh cell to my point
        labelIOList cellPoint_;

        //- From polyMesh boundary face (face-mesh.nInternalFaces()) to my point
        labelIOList boundaryFacePoint_;


    // Private Member Functions

        static labelList getFaceOrder
        (
            const labelList& faceOwner,
            const labelList& faceNeighbour,
            const cellList& cells,
            label& nInternalFaces
        );

        static void getPointEdges
        (
            const primitivePatch& patch,
            const label faceI,
            const label pointI,
            label& e0,
            label& e1
        );

        static labelList collectPatchSideFace
        (
            const polyPatch& patch,
            const label patchToDualOffset,
            const labelList& edgeToDualPoint,
            const labelList& pointToDualPoint,
            const label pointI,

            label& edgeI
        );

        static void collectPatchInternalFace
        (
            const polyPatch& patch,
            const label patchToDualOffset,
            const labelList& edgeToDualPoint,

            const label pointI,
            const label startEdgeI,

            labelList& dualFace2,
            labelList& featEdgeIndices2
        );

        static void splitFace
        (
            const polyPatch& patch,
            const labelList& pointToDualPoint,

            const label pointI,
            const labelList& dualFace,
            const labelList& featEdgeIndices,

            DynamicList<face>& dualFaces,
            dynamicLabelList& dualOwner,
            dynamicLabelList& dualNeighbour,
            dynamicLabelList& dualRegion
        );

        static void dualPatch
        (
            const polyPatch& patch,
            const label patchToDualOffset,
            const labelList& edgeToDualPoint,
            const labelList& pointToDualPoint,

            const pointField& dualPoints,

            DynamicList<face>& dualFaces,
            dynamicLabelList& dualOwner,
            dynamicLabelList& dualNeighbour,
            dynamicLabelList& dualRegion
        );

        void calcDual
        (
            const polyMesh& mesh,
            const labelList& featureEdges,
            const labelList& featurePoints
        );


        //- Disallow default bitwise copy construct
        polyDualMesh(const polyDualMesh&);

        //- Disallow default bitwise assignment
        void operator=(const polyDualMesh&);


public:

    //- Runtime type information
    ClassName("polyDualMesh");


    // Constructors

        //- Construct from IOobject
        polyDualMesh(const IOobject&);

        //- Construct from polyMesh and list of edges and points to represent.
        //  Feature edge and point labels are in local addressing of a patch
        //  over all boundary faces.
        polyDualMesh
        (
            const polyMesh&,
            const labelList& featureEdges,
            const labelList& featurePoints
        );

        //- Construct from polyMesh and feature edge angle. Uses calcFeatures
        //  below to determine feature edges and points.
        polyDualMesh
        (
            const polyMesh&,
            const scalar featureCos
        );

        //- Helper function to create feature edges and points based on
        //  feature angle and patches.
        static void calcFeatures
        (
            const polyMesh&,
            const scalar featureCos,
            labelList& featureEdges,
            labelList& featurePoints
        );


    // Destructor

        ~polyDualMesh();


    // Member Functions

        // Access

            //- From polyMesh cell to dual point
            const labelIOList& cellPoint() const
            {
                return cellPoint_;
            }

            //- From polyMesh patch face to dual point
            const labelIOList& boundaryFacePoint() const
            {
                return boundaryFacePoint_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
