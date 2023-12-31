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
    Foam::autoSnapDriver

Description
    All to do with snapping to surface

SourceFiles
    autoSnapDriver.C

\*---------------------------------------------------------------------------*/

#ifndef autoSnapDriver_H
#define autoSnapDriver_H

#include "PackedBoolList.H"
#include "meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class motionSmoother;
class snapParameters;

/*---------------------------------------------------------------------------*\
                           Class autoSnapDriver Declaration
\*---------------------------------------------------------------------------*/

class autoSnapDriver
{
    // Private classes

        //- Combine operator class for equalizing displacements.
        class minMagEqOp
        {
        public:

            void operator()(vector& x, const vector& y) const
            {
                if (magSqr(y) < magSqr(x))
                {
                    x = y;
                }
            }
        };


    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- From surface region to patch
        const labelList globalToPatch_;


    // Private Member Functions


        // Snapping

            //- Get faces to repatch. Returns map from face to patch.
            Map<label> getZoneBafflePatches(const bool allowBoundary) const;

            //- Calculates (geometric) shared points
            static label getCollocatedPoints
            (
                const scalar tol,
                const pointField&,
                PackedBoolList&
            );

            //- Calculate displacement per patch point to smooth out patch.
            //  Quite complicated in determining which points to move where.
            pointField smoothPatchDisplacement
            (
                const motionSmoother&,
                const List<labelPair>&
            ) const;

            //- Check that face zones are synced
            void checkCoupledFaceZones() const;

            //- Per edge distance to patch
            static tmp<scalarField> edgePatchDist
            (
                const pointMesh&,
                const indirectPrimitivePatch&
            );

            //- Write displacement as .obj file.
            static void dumpMove
            (
                const fileName&,
                const pointField&,
                const pointField&
            );

            //- Check displacement is outwards pointing
            static bool outwardsDisplacement
            (
                const indirectPrimitivePatch&,
                const vectorField&
            );


        //- Disallow default bitwise copy construct
        autoSnapDriver(const autoSnapDriver&);

        //- Disallow default bitwise assignment
        void operator=(const autoSnapDriver&);


public:

    //- Runtime type information
    ClassName("autoSnapDriver");


    // Constructors

        //- Construct from components
        autoSnapDriver
        (
            meshRefinement& meshRefiner,
            const labelList& globalToPatch
        );


    // Member Functions

        // Snapping

            //- Create baffles for faces straddling zoned surfaces. Return
            //  baffles.
            autoPtr<mapPolyMesh> createZoneBaffles(List<labelPair>&);

            //- Merge baffles.
            autoPtr<mapPolyMesh> mergeZoneBaffles(const List<labelPair>&);

            //- Calculate edge length per patch point.
            scalarField calcSnapDistance
            (
                const snapParameters& snapParams,
                const indirectPrimitivePatch&
            ) const;

            //- Smooth the mesh (patch and internal) to increase visibility
            //  of surface points (on castellated mesh) w.r.t. surface.
            void preSmoothPatch
            (
                const snapParameters& snapParams,
                const label nInitErrors,
                const List<labelPair>& baffles,
                motionSmoother&
            ) const;

            //- Get points both on patch and facezone.
            labelList getZoneSurfacePoints
            (
                const indirectPrimitivePatch&,
                const word& zoneName
            ) const;

            //- Per patch point calculate point on nearest surface. Set as
            //  boundary conditions of motionSmoother displacement field. Return
            //  displacement of patch points.
            vectorField calcNearestSurface
            (
                const scalarField& snapDist,
                motionSmoother& meshMover
            ) const;

            //- Smooth the displacement field to the internal.
            void smoothDisplacement
            (
                const snapParameters& snapParams,
                motionSmoother&
            ) const;

            //- Do the hard work: move the mesh according to displacement,
            //  locally relax the displacement.
            void scaleMesh
            (
                const snapParameters& snapParams,
                const label nInitErrors,
                const List<labelPair>& baffles,
                motionSmoother&
            );

            //- Repatch faces according to surface nearest the face centre
            autoPtr<mapPolyMesh> repatchToSurface
            (
                const snapParameters& snapParams,
                const labelList& adaptPatchIDs
            );


            void doSnap
            (
                const dictionary& snapDict,
                const dictionary& motionDict,
                const snapParameters& snapParams
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
