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
    Foam::meshRefinement

Description
    Helper class which maintains intersections of (changing) mesh with
    (static) surfaces.

    Maintains
    - per face any intersections of the cc-cc segment with any of the surfaces

SourceFiles
    meshRefinement.C
    meshRefinementBaffles.C
    meshRefinementMerge.C
    meshRefinementProblemCells.C
    meshRefinementRefine.C

\*---------------------------------------------------------------------------*/

#ifndef meshRefinement_H
#define meshRefinement_H

#include "hexRef8.H"
#include "mapPolyMesh.H"
#include "autoPtr.H"
#include "labelPair.H"
#include "indirectPrimitivePatch.H"
#include "pointFieldsFwd.H"
#include "Tuple2.H"
#include "pointIndexHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class fvMesh;
class mapDistributePolyMesh;
class decompositionMethod;
class refinementSurfaces;
class shellSurfaces;
class removeCells;
class featureEdgeMesh;
class fvMeshDistribute;
class searchableSurface;
class regionSplit;
class globalIndex;

/*---------------------------------------------------------------------------*\
                           Class meshRefinement Declaration
\*---------------------------------------------------------------------------*/

class meshRefinement
{

public:

    // Public data types

        //- Enumeration for debug dumping
        enum writeFlag
        {
            MESH = 1,
            SCALARLEVELS = 2,
            OBJINTERSECTIONS = 4
        };


        //- Enumeration for how the userdata is to be mapped upon refinement.
        enum mapType
        {
            MASTERONLY = 1, /*!< maintain master only */
            KEEPALL = 2,    /*!< have slaves (upon refinement) from master */
            REMOVE = 4      /*!< set value to -1 any face that was refined */
        };


private:

    // Private data

        //- Reference to mesh
        fvMesh& mesh_;

        //- tolerance used for sorting coordinates (used in 'less' routine)
        const scalar mergeDistance_;

        //- overwrite the mesh?
        const bool overwrite_;

        //- Instance of mesh upon construction. Used when in overwrite_ mode.
        const word oldInstance_;

        //- All surface-intersection interaction
        const refinementSurfaces& surfaces_;

        //- All shell-refinement interaction
        const shellSurfaces& shells_;

        //- refinement engine
        hexRef8 meshCutter_;

        //- per cc-cc vector the index of the surface hit
        labelIOList surfaceIndex_;

        //- user supplied face based data.
        List<Tuple2<mapType, labelList> > userFaceData_;

        //- Meshed patches - are treated differently. Stored as wordList since
        //  order changes.
        wordList meshedPatches_;


    // Private Member Functions

        //- Reorder list according to map.
        template<class T>
        static void updateList
        (
            const labelList& newToOld,
            const T& nullValue,
            List<T>& elems
        );

        //- Add patchfield of given type to all fields on mesh
        template<class GeoField>
        static void addPatchFields(fvMesh&, const word& patchFieldType);

        //- Reorder patchfields of all fields on mesh
        template<class GeoField>
        static void reorderPatchFields(fvMesh&, const labelList& oldToNew);

        //- Find out which faces have changed given cells (old mesh labels)
        //  that were marked for refinement.
        static labelList getChangedFaces
        (
            const mapPolyMesh&,
            const labelList& oldCellsToRefine
        );

        //- Calculate coupled boundary end vector and refinement level
        void calcNeighbourData
        (
            labelList& neiLevel,
            pointField& neiCc
        ) const;

        //- Find any intersection of surface. Store in surfaceIndex_.
        void updateIntersections(const labelList& changedFaces);

        //- Remove cells. Put exposedFaces into exposedPatchIDs.
        autoPtr<mapPolyMesh> doRemoveCells
        (
            const labelList& cellsToRemove,
            const labelList& exposedFaces,
            const labelList& exposedPatchIDs,
            removeCells& cellRemover
        );

        // Get cells which are inside any closed surface. Note that
        // all closed surfaces
        // will have already been oriented to have keepPoint outside.
        labelList getInsideCells(const word&) const;

        // Do all to remove inside cells
        autoPtr<mapPolyMesh> removeInsideCells
        (
            const string& msg,
            const label exposedPatchI
        );

        // For decomposeCombineRegions

            //- Used in decomposeCombineRegions. Given global region per cell
            //  determines master processor/cell for regions straddling
            //  procboundaries.
            void getCoupledRegionMaster
            (
                const globalIndex& globalCells,
                const boolList& blockedFace,
                const regionSplit& globalRegion,
                Map<label>& regionToMaster
            ) const;

            //- Determine regions that are local to me or coupled ones that
            //  are owned by me. Determine representative location.
            void calcLocalRegions
            (
                const globalIndex& globalCells,
                const labelList& globalRegion,
                const Map<label>& coupledRegionToMaster,
                const scalarField& cellWeights,

                Map<label>& globalToLocalRegion,
                pointField& localPoints,
                scalarField& localWeights
            ) const;

            //- Convert region into global index.
            static label getShiftedRegion
            (
                const globalIndex& indexer,
                const Map<label>& globalToLocalRegion,
                const Map<label>& coupledRegionToShifted,
                const label globalRegion
            );

            //- helper: add element if not in list. Linear search.
            static void addUnique(const label, labelList&);

            //- Calculate region connectivity. Major communication.
            void calcRegionRegions
            (
                const labelList& globalRegion,
                const Map<label>& globalToLocalRegion,
                const Map<label>& coupledRegionToMaster,
                labelListList& regionRegions
            ) const;


        // Refinement candidate selection

            //- Mark cell for refinement (if not already marked). Return false
            // if refinelimit hit. Keeps running count (in nRefine) of cells
            // marked for refinement
            static bool markForRefine
            (
                const label markValue,
                const label nAllowRefine,
                label& cellValue,
                label& nRefine
            );

            //- Calculate list of cells to refine based on intersection of
            //  features.
            label markFeatureRefinement
            (
                const point& keepPoint,
                const PtrList<featureEdgeMesh>& featureMeshes,
                const labelList& featureLevels,
                const label nAllowRefine,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cells for refinement-shells based refinement.
            label markInternalRefinement
            (
                const label nAllowRefine,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Collect faces that are intersected and whose neighbours aren't
            //  yet marked  for refinement.
            labelList getRefineCandidateFaces
            (
                const labelList& refineCell
            ) const;

            //- Mark cells for surface intersection based refinement.
            label markSurfaceRefinement
            (
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cell if local curvature > curvature or
            //  markDifferingRegions = true and intersections with different
            //  regions.
            bool checkCurvature
            (
                const scalar curvature,
                const label nAllowRefine,
                const label surfaceLevel,
                const vector& surfaceNormal,
                const label cellI,
                label& cellMaxLevel,
                vector& cellMaxNormal,
                labelList& refineCell,
                label& nRefine
            ) const;


            //- Mark cells for surface curvature based refinement. Marks if
            //  local curvature > curvature or if on different regions
            //  (markDifferingRegions)
            label markSurfaceCurvatureRefinement
            (
                const scalar curvature,
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,
                labelList& refineCell,
                label& nRefine
            ) const;

        // Baffle handling

            //- Determine patches for baffles
            void getBafflePatches
            (
                const labelList& globalToPatch,
                const labelList& neiLevel,
                const pointField& neiCc,
                labelList& ownPatch,
                labelList& neiPatch
            ) const;

            //- Determine patch for baffle using some heuristic (and not
            //  surface)
            label getBafflePatch
            (
                const labelList& facePatch,
                const label faceI
            ) const;

            //- Repatches external face or creates baffle for internal face
            //  with user specified patches (might be different for both sides).
            //  Returns label of added face.
            label createBaffle
            (
                const label faceI,
                const label ownPatch,
                const label neiPatch,
                directTopoChange& meshMod
            ) const;

        // Problem cell handling

            //- Helper function to mark face as being on 'boundary'. Used by
            //  markFacesOnProblemCells
            void markBoundaryFace
            (
                const label faceI,
                boolList& isBoundaryFace,
                boolList& isBoundaryEdge,
                boolList& isBoundaryPoint
            ) const;

            void findNearest
            (
                const labelList& meshFaces,
                List<pointIndexHit>& nearestInfo,
                labelList& nearestSurface,
                labelList& nearestRegion,
                vectorField& nearestNormal
            ) const;

            Map<label> findEdgeConnectedProblemCells
            (
                const scalarField& perpendicularAngle,
                const labelList&
            ) const;

            bool isCollapsedFace
            (
                const pointField&,
                const pointField& neiCc,
                const scalar minFaceArea,
                const scalar maxNonOrtho,
                const label faceI
            ) const;

            bool isCollapsedCell
            (
                const pointField&,
                const scalar volFraction,
                const label cellI
            ) const;

            //- Returns list with for every internal face -1 or the patch
            //  they should be baffled into. If removeEdgeConnectedCells is set
            //  removes cells based on perpendicularAngle.
            labelList markFacesOnProblemCells
            (
                const dictionary& motionDict,
                const bool removeEdgeConnectedCells,
                const scalarField& perpendicularAngle,
                const labelList& globalToPatch
            ) const;

            ////- Initial test of marking faces using geometric information.
            //labelList markFacesOnProblemCellsGeometric
            //(
            //    const dictionary& motionDict
            //) const;


        // Baffle merging

            //- Extract those baffles (duplicate) faces that are on the edge
            //  of a baffle region. These are candidates for merging.
            List<labelPair> filterDuplicateFaces(const List<labelPair>&) const;


        // Zone handling

            //- Finds zone per cell for cells inside closed named surfaces.
            //  (uses geometric test for insideness)
            //  Adapts namedSurfaceIndex so all faces on boundary of cellZone
            //  have corresponding faceZone.
            void findCellZoneGeometric
            (
                const labelList& closedNamedSurfaces,
                labelList& namedSurfaceIndex,
                const labelList& surfaceToCellZone,
                labelList& cellToZone
            ) const;

            //- Determines cell zone from cell region information.
            bool calcRegionToZone
            (
                const label surfZoneI,
                const label ownRegion,
                const label neiRegion,

                labelList& regionToCellZone
            ) const;

            //- Finds zone per cell. Uses topological walk with all faces
            //  marked in namedSurfaceIndex regarded as blocked.
            void findCellZoneTopo
            (
                const point& keepPoint,
                const labelList& namedSurfaceIndex,
                const labelList& surfaceToCellZone,
                labelList& cellToZone
            ) const;

            void makeConsistentFaceIndex
            (
                const labelList& cellToZone,
                labelList& namedSurfaceIndex
            ) const;


        //- Disallow default bitwise copy construct
        meshRefinement(const meshRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const meshRefinement&);

public:

    //- Runtime type information
    ClassName("meshRefinement");


    // Constructors

        //- Construct from components
        meshRefinement
        (
            fvMesh& mesh,
            const scalar mergeDistance,
            const bool overwrite,
            const refinementSurfaces&,
            const shellSurfaces&
        );


    // Member Functions

        // Access

            //- reference to mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }
            fvMesh& mesh()
            {
                return mesh_;
            }

            scalar mergeDistance() const
            {
                return mergeDistance_;
            }

            //- Overwrite the mesh?
            bool overwrite() const
            {
                return overwrite_;
            }

            //- (points)instance of mesh upon construction
            const word& oldInstance() const
            {
                return oldInstance_;
            }

            //- reference to surface search engines
            const refinementSurfaces& surfaces() const
            {
                return surfaces_;
            }

            //- reference to refinement shells (regions)
            const shellSurfaces& shells() const
            {
                return shells_;
            }

            //- reference to meshcutting engine
            const hexRef8& meshCutter() const
            {
                return meshCutter_;
            }

            //- per start-end edge the index of the surface hit
            const labelList& surfaceIndex() const
            {
                return surfaceIndex_;
            }

            labelList& surfaceIndex()
            {
                return surfaceIndex_;
            }

            //- Additional face data that is maintained across
            //  topo changes. Every entry is a list over all faces.
            //  Bit of a hack. Additional flag to say whether to maintain master
            //  only (false) or increase set to account for face-from-face.
            const List<Tuple2<mapType, labelList> >& userFaceData() const
            {
                return userFaceData_;
            }

            List<Tuple2<mapType, labelList> >& userFaceData()
            {
                return userFaceData_;
            }


        // Other

            //- Count number of intersections (local)
            label countHits() const;

            //- Helper function to get decomposition such that all connected
            //  regions get moved onto one processor. Used to prevent baffles
            //  straddling processor boundaries. explicitConnections is to
            //  keep pairs of non-coupled boundary faces together
            //  (e.g. to keep baffles together)
            labelList decomposeCombineRegions
            (
                const scalarField& cellWeights,
                const boolList& blockedFace,
                const List<labelPair>& explicitConnections,
                decompositionMethod&
            ) const;

            //- Redecompose according to cell count
            //  keepZoneFaces : find all faceZones from zoned surfaces and keep
            //                  owner and neighbour together
            //  keepBaffles   : find all baffles and keep them together
            autoPtr<mapDistributePolyMesh> balance
            (
                const bool keepZoneFaces,
                const bool keepBaffles,
                const scalarField& cellWeights,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor
            );

            //- Get faces with intersection.
            labelList intersectedFaces() const;

            //- Get points on surfaces with intersection and boundary faces.
            labelList intersectedPoints() const;

            //- Create patch from set of patches
            static autoPtr<indirectPrimitivePatch> makePatch
            (
                const polyMesh&,
                const labelList&
            );

            //- Helper function to make a pointVectorField with correct
            //  bcs for mesh movement:
            //  - adaptPatchIDs         : fixedValue
            //  - processor             : calculated (so free to move)
            //  - cyclic/wedge/symmetry : slip
            //  - other                 : slip
            static tmp<pointVectorField> makeDisplacementField
            (
                const pointMesh& pMesh,
                const labelList& adaptPatchIDs
            );

            //- Helper function: check that face zones are synced
            static void checkCoupledFaceZones(const polyMesh&);


        // Refinement

            //- Calculate list of cells to refine.
            labelList refineCandidates
            (
                const point& keepPoint,
                const scalar curvature,

                const PtrList<featureEdgeMesh>& featureMeshes,
                const labelList& featureLevels,

                const bool featureRefinement,
                const bool internalRefinement,
                const bool surfaceRefinement,
                const bool curvatureRefinement,
                const label maxGlobalCells,
                const label maxLocalCells
            ) const;

            //- Refine some cells
            autoPtr<mapPolyMesh> refine(const labelList& cellsToRefine);

            //- Refine some cells and rebalance
            autoPtr<mapDistributePolyMesh> refineAndBalance
            (
                const string& msg,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor,
                const labelList& cellsToRefine,
                const scalar maxLoadUnbalance
            );

            //- Balance before refining some cells
            autoPtr<mapDistributePolyMesh> balanceAndRefine
            (
                const string& msg,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor,
                const labelList& cellsToRefine,
                const scalar maxLoadUnbalance
            );


        // Baffle handling

            //- Split off unreachable areas of mesh.
            void baffleAndSplitMesh
            (
                const bool handleSnapProblems,
                const bool removeEdgeConnectedCells,
                const scalarField& perpendicularAngle,
                const bool mergeFreeStanding,
                const dictionary& motionDict,
                Time& runTime,
                const labelList& globalToPatch,
                const point& keepPoint
            );

            //- Split off (with optional buffer layers) unreachable areas
            //  of mesh. Does not introduce baffles.
            autoPtr<mapPolyMesh> splitMesh
            (
                const label nBufferLayers,
                const labelList& globalToPatch,
                const point& keepPoint
            );

            //- Find boundary points that connect to more than one cell
            //  region and split them.
            autoPtr<mapPolyMesh> dupNonManifoldPoints();

            //- Create baffle for every internal face where ownPatch != -1.
            //  External faces get repatched according to ownPatch (neiPatch
            //  should be -1 for these)
            autoPtr<mapPolyMesh> createBaffles
            (
                const labelList& ownPatch,
                const labelList& neiPatch
            );

            //- Return a list of coupled face pairs, i.e. faces that
            //  use the same vertices.
            List<labelPair> getDuplicateFaces(const labelList& testFaces) const;

            //- Merge baffles. Gets pairs of faces.
            autoPtr<mapPolyMesh> mergeBaffles(const List<labelPair>&);

            //- Put faces/cells into zones according to surface specification.
            //  Returns null if no zone surfaces present. Region containing
            //  the keepPoint will not be put into a cellZone.
            autoPtr<mapPolyMesh> zonify
            (
                const point& keepPoint,
                const bool allowFreeStandingZoneFaces
            );


        // Other topo changes

            //- Helper:add patch to mesh. Update all registered fields.
            //  Use addMeshedPatch to add patches originating from surfaces.
            static label addPatch(fvMesh&, const word& name, const word& type);

            //- Add patch originating from meshing. Update meshedPatches_.
            label addMeshedPatch(const word& name, const word& type);

            //- Get patchIDs for patches added in addMeshedPatch.
            labelList meshedPatches() const;

            //- Split mesh. Keep part containing point.
            autoPtr<mapPolyMesh> splitMeshRegions(const point& keepPoint);

            //- Update local numbering for mesh redistribution
            void distribute(const mapDistributePolyMesh&);

            //- Update for external change to mesh. changedFaces are in new mesh
            //  face labels.
            void updateMesh
            (
                const mapPolyMesh&,
                const labelList& changedFaces
            );


            // Restoring : is where other processes delete and reinsert data.

                //- Signal points/face/cells for which to store data
                void storeData
                (
                    const labelList& pointsToStore,
                    const labelList& facesToStore,
                    const labelList& cellsToStore
                );

                //- Update local numbering + undo
                //  Data to restore given as new pointlabel + stored pointlabel
                //  (i.e. what was in pointsToStore)
                void updateMesh
                (
                    const mapPolyMesh&,
                    const labelList& changedFaces,
                    const Map<label>& pointsToRestore,
                    const Map<label>& facesToRestore,
                    const Map<label>& cellsToRestore
                );


            //- Merge faces on the same patch (usually from exposing refinement)
            //  Returns global number of faces merged.
            label mergePatchFaces
            (
                const scalar minCos,
                const scalar concaveCos,
                const labelList& patchIDs
            );

            //- Remove points not used by any face or points used
            // by only two faces where the edges are in line
            autoPtr<mapPolyMesh> mergeEdges(const scalar minCos);


        // Debug/IO

            //- Debugging: check that all faces still obey start()>end()
            void checkData();

            //- Compare two lists over all boundary faces
            template<class T>
            void testSyncBoundaryFaceList
            (
                const scalar mergeDistance,
                const string&,
                const UList<T>&,
                const UList<T>&
            ) const;

            //- Print some mesh stats.
            void printMeshInfo(const bool, const string&) const;

            //- Replacement for Time::timeName() : return oldInstance (if
            //  overwrite_)
            word timeName() const;

            //- Set instance of all local IOobjects
            void setInstance(const fileName&);

            //- Write mesh and all data
            bool write() const;

            //- Write refinement level as volScalarFields for postprocessing
            void dumpRefinementLevel() const;

            //- Debug: Write intersection information to OBJ format
            void dumpIntersections(const fileName& prefix) const;

            //- Do any one of above IO functions. flag is combination of
            //  writeFlag values.
            void write(const label flag, const fileName&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "meshRefinementTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
