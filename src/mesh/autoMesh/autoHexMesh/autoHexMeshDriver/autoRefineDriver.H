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
    Foam::autoRefineDriver

Description

SourceFiles
    autoRefineDriver.C

\*---------------------------------------------------------------------------*/

#ifndef autoRefineDriver_H
#define autoRefineDriver_H

#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class featureEdgeMesh;
class refinementParameters;
class meshRefinement;
class decompositionMethod;
class fvMeshDistribute;

/*---------------------------------------------------------------------------*\
                           Class autoRefineDriver Declaration
\*---------------------------------------------------------------------------*/

class autoRefineDriver
{
    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- Reference to decomposition method
        decompositionMethod& decomposer_;

        //- Reference to mesh distribution engine
        fvMeshDistribute& distributor_;

        //- From surface region to patch
        const labelList globalToPatch_;


    // Private Member Functions

        //- Read explicit feature edges
        label readFeatureEdges
        (
            const PtrList<dictionary>& featDicts,
            PtrList<featureEdgeMesh>& featureMeshes,
            labelList& featureLevel
        ) const;

        //- Refine all cells pierced by explicit feature edges
        label featureEdgeRefine
        (
            const refinementParameters& refineParams,
            const PtrList<dictionary>& featDicts,
            const label maxIter,
            const label minRefine
        );

        //- Refine all cells interacting with the surface
        label surfaceOnlyRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Remove all cells within intersected region
        void removeInsideCells
        (
            const refinementParameters& refineParams,
            const label nBufferLayers
        );

        //- Remove all cells inside/outside shell
        label shellRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Add baffles and remove unreachable cells
        void baffleAndSplitMesh
        (
            const refinementParameters& refineParams,
            const bool handleSnapProblems,
            const dictionary& motionDict
        );

        //- Add zones
        void zonify(const refinementParameters& refineParams);

        void splitAndMergeBaffles
        (
            const refinementParameters& refineParams,
            const bool handleSnapProblems,
            const dictionary& motionDict
        );

        //- Merge refined boundary faces (from exposing coarser cell)
        void mergePatchFaces
        (
            const refinementParameters& refineParams
        );


        //- Disallow default bitwise copy construct
        autoRefineDriver(const autoRefineDriver&);

        //- Disallow default bitwise assignment
        void operator=(const autoRefineDriver&);


public:

    //- Runtime type information
    ClassName("autoRefineDriver");


    // Constructors

        //- Construct from components
        autoRefineDriver
        (
            meshRefinement& meshRefiner,
            decompositionMethod& decomposer,
            fvMeshDistribute& distributor,
            const labelList& globalToPatch
        );


    // Member Functions

        //- Do all the refinement
        void doRefine
        (
            const dictionary& refineDict,
            const refinementParameters& refineParams,
            const bool prepareForSnapping,
            const dictionary& motionDict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
