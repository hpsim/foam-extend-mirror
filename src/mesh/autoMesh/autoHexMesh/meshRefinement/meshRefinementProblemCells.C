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

\*---------------------------------------------------------------------------*/

#include "meshRefinement.H"
#include "fvMesh.H"
#include "syncTools.H"
#include "foamTime.H"
#include "refinementSurfaces.H"
#include "pointSet.H"
#include "faceSet.H"
#include "indirectPrimitivePatch.H"
#include "OFstream.H"
#include "cellSet.H"
#include "searchableSurfaces.H"
#include "polyMeshGeometry.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshRefinement::markBoundaryFace
(
    const label faceI,
    boolList& isBoundaryFace,
    boolList& isBoundaryEdge,
    boolList& isBoundaryPoint
) const
{
    isBoundaryFace[faceI] = true;

    const labelList& fEdges = mesh_.faceEdges(faceI);

    forAll(fEdges, fp)
    {
        isBoundaryEdge[fEdges[fp]] = true;
    }

    const face& f = mesh_.faces()[faceI];

    forAll(f, fp)
    {
        isBoundaryPoint[f[fp]] = true;
    }
}


void Foam::meshRefinement::findNearest
(
    const labelList& meshFaces,
    List<pointIndexHit>& nearestInfo,
    labelList& nearestSurface,
    labelList& nearestRegion,
    vectorField& nearestNormal
) const
{
    pointField fc(meshFaces.size());
    forAll(meshFaces, i)
    {
        fc[i] = mesh_.faceCentres()[meshFaces[i]];
    }

    const labelList allSurfaces(identity(surfaces_.surfaces().size()));

    surfaces_.findNearest
    (
        allSurfaces,
        fc,
        scalarField(fc.size(), sqr(GREAT)),    // sqr of attraction
        nearestSurface,
        nearestInfo
    );

    // Do normal testing per surface.
    nearestNormal.setSize(nearestInfo.size());
    nearestRegion.setSize(nearestInfo.size());

    forAll(allSurfaces, surfI)
    {
        DynamicList<pointIndexHit> localHits;

        forAll(nearestSurface, i)
        {
            if (nearestSurface[i] == surfI)
            {
                localHits.append(nearestInfo[i]);
            }
        }

        label geomI = surfaces_.surfaces()[surfI];

        pointField localNormals;
        surfaces_.geometry()[geomI].getNormal(localHits, localNormals);

        labelList localRegion;
        surfaces_.geometry()[geomI].getRegion(localHits, localRegion);

        label localI = 0;
        forAll(nearestSurface, i)
        {
            if (nearestSurface[i] == surfI)
            {
                nearestNormal[i] = localNormals[localI];
                nearestRegion[i] = localRegion[localI];
                localI++;
            }
        }
    }
}


Foam::Map<Foam::label> Foam::meshRefinement::findEdgeConnectedProblemCells
(
    const scalarField& perpendicularAngle,
    const labelList& globalToPatch
) const
{
    // Construct addressing engine from all patches added for meshing.
    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh_,
            meshedPatches()
        )
    );
    const indirectPrimitivePatch& pp = ppPtr();


    // 1. Collect faces to test
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    dynamicLabelList candidateFaces(pp.size()/20);

    const labelListList& edgeFaces = pp.edgeFaces();

    const labelList& cellLevel = meshCutter_.cellLevel();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() == 2)
        {
            label face0 = pp.addressing()[eFaces[0]];
            label face1 = pp.addressing()[eFaces[1]];

            label cell0 = mesh_.faceOwner()[face0];
            label cell1 = mesh_.faceOwner()[face1];

            if (cellLevel[cell0] > cellLevel[cell1])
            {
                // cell0 smaller.
                const vector& n0 = pp.faceNormals()[eFaces[0]];
                const vector& n1 = pp.faceNormals()[eFaces[1]];

                if (mag(n0 & n1) < 0.1)
                {
                    candidateFaces.append(face0);
                }
            }
            else if (cellLevel[cell1] > cellLevel[cell0])
            {
                // cell1 smaller.
                const vector& n0 = pp.faceNormals()[eFaces[0]];
                const vector& n1 = pp.faceNormals()[eFaces[1]];

                if (mag(n0 & n1) < 0.1)
                {
                    candidateFaces.append(face1);
                }
            }
        }
    }
    candidateFaces.shrink();

    Info<< "Testing " << returnReduce(candidateFaces.size(), sumOp<label>())
        << " faces on edge-connected cells of differing level."
        << endl;

    if (debug)
    {
        faceSet fSet
        (
            mesh_,
            "edgeConnectedFaces",
            labelHashSet(candidateFaces)
        );
        Pout<< "Writing " << fSet.size()
            << " with problematic topology to faceSet "
            << fSet.objectPath() << endl;
        fSet.write();
    }


    // 2. Find nearest surface on candidate faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    List<pointIndexHit> nearestInfo;
    labelList nearestSurface;
    labelList nearestRegion;
    vectorField nearestNormal;
    findNearest
    (
        candidateFaces,
        nearestInfo,
        nearestSurface,
        nearestRegion,
        nearestNormal
    );


    // 3. Test angle to surface
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    Map<label> candidateCells(candidateFaces.size());

    faceSet perpFaces(mesh_, "perpendicularFaces", pp.size()/100);

    forAll(candidateFaces, i)
    {
        label faceI = candidateFaces[i];

        vector n = mesh_.faceAreas()[faceI];
        n /= mag(n);

        label region = surfaces_.globalRegion
        (
            nearestSurface[i],
            nearestRegion[i]
        );

        scalar angle =
            perpendicularAngle[region]
          / 180.0
          * mathematicalConstant::pi;

        if (angle >= 0)
        {
            if (mag(n & nearestNormal[i]) < Foam::sin(angle))
            {
                perpFaces.insert(faceI);
                candidateCells.insert
                (
                    mesh_.faceOwner()[faceI],
                    globalToPatch[region]
                );
            }
        }
    }

    if (debug)
    {
        Pout<< "Writing " << perpFaces.size()
            << " faces that are perpendicular to the surface to set "
            << perpFaces.objectPath() << endl;
        perpFaces.write();
    }
    return candidateCells;
}


// Check if moving face to new points causes a 'collapsed' face.
// Uses new point position only for the face, not the neighbouring
// cell centres
bool Foam::meshRefinement::isCollapsedFace
(
    const pointField& points,
    const pointField& neiCc,
    const scalar minFaceArea,
    const scalar maxNonOrtho,
    const label faceI
) const
{
    vector s = mesh_.faces()[faceI].normal(points);
    scalar magS = mag(s);

    // Check face area
    if (magS < minFaceArea)
    {
        return true;
    }

    // Check orthogonality
    const point& ownCc = mesh_.cellCentres()[mesh_.faceOwner()[faceI]];

    if (mesh_.isInternalFace(faceI))
    {
        label nei = mesh_.faceNeighbour()[faceI];
        vector d = ownCc - mesh_.cellCentres()[nei];

        scalar dDotS = (d & s)/(mag(d)*magS + VSMALL);

        if (dDotS < maxNonOrtho)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        label patchI = mesh_.boundaryMesh().whichPatch(faceI);

        if (mesh_.boundaryMesh()[patchI].coupled())
        {
            vector d = ownCc - neiCc[faceI-mesh_.nInternalFaces()];

            scalar dDotS = (d & s)/(mag(d)*magS + VSMALL);

            if (dDotS < maxNonOrtho)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            // Collapsing normal boundary face does not cause problems with
            // non-orthogonality
            return true;
        }
    }
}


// Check if moving cell to new points causes it to collapse.
bool Foam::meshRefinement::isCollapsedCell
(
    const pointField& points,
    const scalar volFraction,
    const label cellI
) const
{
    scalar vol = mesh_.cells()[cellI].mag(points, mesh_.faces());

    if (vol/mesh_.cellVolumes()[cellI] < volFraction)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Returns list with for every internal face -1 or the patch they should
// be baffled into. Gets run after createBaffles so all the unzoned surface
// intersections have already been turned into baffles. (Note: zoned surfaces
// are not baffled at this stage)
// Used to remove cells by baffling all their faces and have the
// splitMeshRegions chuck away non used regions.
Foam::labelList Foam::meshRefinement::markFacesOnProblemCells
(
    const dictionary& motionDict,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const labelList& globalToPatch
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Per internal face (boundary faces not used) the patch that the
    // baffle should get (or -1)
    labelList facePatch(mesh_.nFaces(), -1);

    // Mark all points and edges on baffle patches (so not on any inlets,
    // outlets etc.)
    boolList isBoundaryPoint(mesh_.nPoints(), false);
    boolList isBoundaryEdge(mesh_.nEdges(), false);
    boolList isBoundaryFace(mesh_.nFaces(), false);

    // Fill boundary data. All elements on meshed patches get marked.
    // Get the labels of added patches.
    labelList adaptPatchIDs(meshedPatches());

    forAll(adaptPatchIDs, i)
    {
        label patchI = adaptPatchIDs[i];

        const polyPatch& pp = patches[patchI];

        label faceI = pp.start();

        forAll(pp, j)
        {
            markBoundaryFace
            (
                faceI,
                isBoundaryFace,
                isBoundaryEdge,
                isBoundaryPoint
            );

            faceI++;
        }
    }

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);


    // Count of faces marked for baffling
    label nBaffleFaces = 0;

    // Count of faces not baffled since would not cause a collapse
    label nPrevented = 0;

    if (removeEdgeConnectedCells && max(perpendicularAngle) >= 0)
    {
        Info<< "markFacesOnProblemCells :"
            << " Checking for edge-connected cells of highly differing sizes."
            << endl;

        // Pick up the cells that need to be removed and (a guess for)
        // the patch they should be patched with.
        Map<label> problemCells
        (
            findEdgeConnectedProblemCells
            (
                perpendicularAngle,
                globalToPatch
            )
        );

        // Baffle all faces of cells that need to be removed
        forAllConstIter(Map<label>, problemCells, iter)
        {
            const cell& cFaces = mesh_.cells()[iter.key()];

            forAll(cFaces, i)
            {
                label faceI = cFaces[i];

                if (facePatch[faceI] == -1 && mesh_.isInternalFace(faceI))
                {
                    facePatch[faceI] = getBafflePatch(facePatch, faceI);
                    nBaffleFaces++;

                    // Mark face as a 'boundary'
                    markBoundaryFace
                    (
                        faceI,
                        isBoundaryFace,
                        isBoundaryEdge,
                        isBoundaryPoint
                    );
                }
            }
        }
        Info<< "markFacesOnProblemCells : Marked "
            << returnReduce(nBaffleFaces, sumOp<label>())
            << " additional internal faces to be converted into baffles"
            << " due to "
            << returnReduce(problemCells.size(), sumOp<label>())
            << " cells edge-connected to lower level cells." << endl;

        if (debug)
        {
            cellSet problemCellSet
            (
                mesh_,
                "problemCells",
                labelHashSet(problemCells.toc())
            );
            Pout<< "Writing " << problemCellSet.size()
                << " cells that are edge connected to coarser cell to set "
                << problemCellSet.objectPath() << endl;
            problemCellSet.write();
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        isBoundaryPoint,
        orEqOp<bool>(),
        false,              // null value
        false               // no separation
    );

    syncTools::syncEdgeList
    (
        mesh_,
        isBoundaryEdge,
        orEqOp<bool>(),
        false,              // null value
        false               // no separation
    );

    syncTools::syncFaceList
    (
        mesh_,
        isBoundaryFace,
        orEqOp<bool>(),
        false               // no separation
    );


    // See if checking for collapse
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Collapse checking parameters
    scalar volFraction = -1;
    if (motionDict.found("minVolCollapseRatio"))
    {
        volFraction = readScalar(motionDict.lookup("minVolCollapseRatio"));
    }
    const bool checkCollapse = (volFraction > 0);
    scalar minArea = -1;
    scalar maxNonOrtho = -1;


    // Find nearest (non-baffle) surface
    pointField newPoints;

    if (checkCollapse)
    {
        minArea = readScalar(motionDict.lookup("minArea"));
        maxNonOrtho = readScalar(motionDict.lookup("maxNonOrtho"));

        Info<< "markFacesOnProblemCells :"
            << " Deleting all-anchor surface cells only if"
            << "snapping them violates mesh quality constraints:" << nl
            << "    snapped/original cell volume < " << volFraction << nl
            << "    face area                    < " << minArea << nl
            << "    non-orthogonality            > " << maxNonOrtho << nl
            << endl;

        // Construct addressing engine.
        autoPtr<indirectPrimitivePatch> ppPtr
        (
            meshRefinement::makePatch
            (
                mesh_,
                adaptPatchIDs
            )
        );
        const indirectPrimitivePatch& pp = ppPtr();
        const pointField& localPoints = pp.localPoints();
        const labelList& meshPoints = pp.meshPoints();

        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        surfaces_.findNearest
        (
            surfaces_.getUnnamedSurfaces(),
            localPoints,
            scalarField(localPoints.size(), sqr(GREAT)),    // sqr of attraction
            hitSurface,
            hitInfo
        );

        // Start of from current points
        newPoints = mesh_.points();

        forAll(hitInfo, i)
        {
            if (hitInfo[i].hit())
            {
                newPoints[meshPoints[i]] = hitInfo[i].hitPoint();
            }
        }
    }



    // For each cell count the number of anchor points that are on
    // the boundary:
    // 8 : check the number of (baffle) boundary faces. If 3 or more block
    //     off the cell since the cell would get squeezed down to a diamond
    //     (probably; if the 3 or more faces are unrefined (only use the
    //      anchor points))
    // 7 : store. Used to check later on whether there are points with
    //     3 or more of these cells. (note that on a flat surface a boundary
    //     point will only have 4 cells connected to it)


    // Does cell have exactly 7 of its 8 anchor points on the boundary?
    PackedBoolList hasSevenBoundaryAnchorPoints(mesh_.nCells());
    // If so what is the remaining non-boundary anchor point?
    labelHashSet nonBoundaryAnchors(mesh_.nCells()/10000);

    // On-the-fly addressing storage.
    dynamicLabelList dynFEdges;
    dynamicLabelList dynCPoints;

    forAll(cellLevel, cellI)
    {
        const labelList& cPoints = mesh_.cellPoints(cellI, dynCPoints);

        // Get number of anchor points (pointLevel <= cellLevel)

        label nBoundaryAnchors = 0;
        label nNonAnchorBoundary = 0;
        label nonBoundaryAnchor = -1;

        forAll(cPoints, i)
        {
            label pointI = cPoints[i];

            if (pointLevel[pointI] <= cellLevel[cellI])
            {
                // Anchor point
                if (isBoundaryPoint[pointI])
                {
                    nBoundaryAnchors++;
                }
                else
                {
                    // Anchor point which is not on the surface
                    nonBoundaryAnchor = pointI;
                }
            }
            else if (isBoundaryPoint[pointI])
            {
                nNonAnchorBoundary++;
            }
        }

        if (nBoundaryAnchors == 8)
        {
            const cell& cFaces = mesh_.cells()[cellI];

            // Count boundary faces.
            label nBfaces = 0;

            forAll(cFaces, cFaceI)
            {
                if (isBoundaryFace[cFaces[cFaceI]])
                {
                    nBfaces++;
                }
            }

            // If nBfaces > 1 make all non-boundary non-baffle faces baffles.
            // We assume that this situation is where there is a single
            // cell sticking out which would get flattened.

            // Eugene: delete cell no matter what.
            //if (nBfaces > 1)
            {
                if
                (
                    checkCollapse
                && !isCollapsedCell(newPoints, volFraction, cellI)
                )
                {
                    nPrevented++;
                    //Pout<< "Preventing collapse of 8 anchor point cell "
                    //    << cellI << " at " << mesh_.cellCentres()[cellI]
                    //    << " since new volume "
                    //    << mesh_.cells()[cellI].mag(newPoints, mesh_.faces())
                    //    << " old volume " << mesh_.cellVolumes()[cellI]
                    //    << endl;
                }
                else
                {
                    // Block all faces of cell
                    forAll(cFaces, cf)
                    {
                        label faceI = cFaces[cf];

                        if
                        (
                            facePatch[faceI] == -1
                         && mesh_.isInternalFace(faceI)
                        )
                        {
                            facePatch[faceI] = getBafflePatch(facePatch, faceI);
                            nBaffleFaces++;

                            // Mark face as a 'boundary'
                            markBoundaryFace
                            (
                                faceI,
                                isBoundaryFace,
                                isBoundaryEdge,
                                isBoundaryPoint
                            );
                        }
                    }
                }
            }
        }
        else if (nBoundaryAnchors == 7)
        {
            // Mark the cell. Store the (single!) non-boundary anchor point.
            hasSevenBoundaryAnchorPoints.set(cellI, 1u);
            nonBoundaryAnchors.insert(nonBoundaryAnchor);
        }
    }


    // Loop over all points. If a point is connected to 4 or more cells
    // with 7 anchor points on the boundary set those cell's non-boundary faces
    // to baffles

    dynamicLabelList dynPCells;

    forAllConstIter(labelHashSet, nonBoundaryAnchors, iter)
    {
        label pointI = iter.key();

        const labelList& pCells = mesh_.pointCells(pointI, dynPCells);

        // Count number of 'hasSevenBoundaryAnchorPoints' cells.
        label n = 0;

        forAll(pCells, i)
        {
            if (hasSevenBoundaryAnchorPoints.get(pCells[i]) == 1u)
            {
                n++;
            }
        }

        if (n > 3)
        {
            // Point in danger of being what? Remove all 7-cells.
            forAll(pCells, i)
            {
                label cellI = pCells[i];

                if (hasSevenBoundaryAnchorPoints.get(cellI) == 1u)
                {
                    if
                    (
                        checkCollapse
                    && !isCollapsedCell(newPoints, volFraction, cellI)
                    )
                    {
                        nPrevented++;
                        //Pout<< "Preventing collapse of 7 anchor cell "
                        //    << cellI
                        //    << " at " << mesh_.cellCentres()[cellI]
                        //    << " since new volume "
                        //    << mesh_.cells()[cellI].mag
                        //        (newPoints, mesh_.faces())
                        //    << " old volume " << mesh_.cellVolumes()[cellI]
                        //    << endl;
                    }
                    else
                    {
                        const cell& cFaces = mesh_.cells()[cellI];

                        forAll(cFaces, cf)
                        {
                            label faceI = cFaces[cf];

                            if
                            (
                                facePatch[faceI] == -1
                             && mesh_.isInternalFace(faceI)
                            )
                            {
                                facePatch[faceI] = getBafflePatch
                                (
                                    facePatch,
                                    faceI
                                );
                                nBaffleFaces++;

                                // Mark face as a 'boundary'
                                markBoundaryFace
                                (
                                    faceI,
                                    isBoundaryFace,
                                    isBoundaryEdge,
                                    isBoundaryPoint
                                );
                            }
                        }
                    }
                }
            }
        }
    }


    // Sync all. (note that pointdata and facedata not used anymore but sync
    // anyway)

    syncTools::syncPointList
    (
        mesh_,
        isBoundaryPoint,
        orEqOp<bool>(),
        false,              // null value
        false               // no separation
    );

    syncTools::syncEdgeList
    (
        mesh_,
        isBoundaryEdge,
        orEqOp<bool>(),
        false,              // null value
        false               // no separation
    );

    syncTools::syncFaceList
    (
        mesh_,
        isBoundaryFace,
        orEqOp<bool>(),
        false               // no separation
    );


    // Find faces with all edges on the boundary and make them baffles
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (facePatch[faceI] == -1)
        {
            const labelList& fEdges = mesh_.faceEdges(faceI, dynFEdges);
            label nFaceBoundaryEdges = 0;

            forAll(fEdges, fe)
            {
                if (isBoundaryEdge[fEdges[fe]])
                {
                    nFaceBoundaryEdges++;
                }
            }

            if (nFaceBoundaryEdges == fEdges.size())
            {
                if
                (
                    checkCollapse
                && !isCollapsedFace
                    (
                        newPoints,
                        neiCc,
                        minArea,
                        maxNonOrtho,
                        faceI
                    )
                )
                {
                    nPrevented++;
                    //Pout<< "Preventing collapse of face " << faceI
                    //    << " with all boundary edges "
                    //    << " at " << mesh_.faceCentres()[faceI]
                    //    << endl;
                }
                else
                {
                    facePatch[faceI] = getBafflePatch(facePatch, faceI);
                    nBaffleFaces++;

                    // Do NOT update boundary data since this would grow blocked
                    // faces across gaps.
                }
            }
        }
    }

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                if (facePatch[faceI] == -1)
                {
                    const labelList& fEdges = mesh_.faceEdges(faceI, dynFEdges);
                    label nFaceBoundaryEdges = 0;

                    forAll(fEdges, fe)
                    {
                        if (isBoundaryEdge[fEdges[fe]])
                        {
                            nFaceBoundaryEdges++;
                        }
                    }

                    if (nFaceBoundaryEdges == fEdges.size())
                    {
                        if
                        (
                            checkCollapse
                        && !isCollapsedFace
                            (
                                newPoints,
                                neiCc,
                                minArea,
                                maxNonOrtho,
                                faceI
                            )
                        )
                        {
                            nPrevented++;
                            //Pout<< "Preventing collapse of coupled face "
                            //    << faceI
                            //    << " with all boundary edges "
                            //    << " at " << mesh_.faceCentres()[faceI]
                            //    << endl;
                        }
                        else
                        {
                            facePatch[faceI] = getBafflePatch(facePatch, faceI);
                            nBaffleFaces++;

                            // Do NOT update boundary data since this would grow
                            // blocked faces across gaps.
                        }
                    }
                }

                faceI++;
            }
        }
    }

    Info<< "markFacesOnProblemCells : marked "
        << returnReduce(nBaffleFaces, sumOp<label>())
        << " additional internal faces to be converted into baffles."
        << endl;

    if (checkCollapse)
    {
        Info<< "markFacesOnProblemCells : prevented "
            << returnReduce(nPrevented, sumOp<label>())
            << " internal faces fom getting converted into baffles."
            << endl;
    }

    return facePatch;
}


//// Mark faces to be baffled to prevent snapping problems. Does
//// test to find nearest surface and checks which faces would get squashed.
//Foam::labelList Foam::meshRefinement::markFacesOnProblemCellsGeometric
//(
//    const dictionary& motionDict
//) const
//{
//    // Construct addressing engine.
//    autoPtr<indirectPrimitivePatch> ppPtr
//    (
//        meshRefinement::makePatch
//        (
//            mesh_,
//            meshedPatches()
//        )
//    );
//    const indirectPrimitivePatch& pp = ppPtr();
//    const pointField& localPoints = pp.localPoints();
//    const labelList& meshPoints = pp.meshPoints();
//
//    // Find nearest (non-baffle) surface
//    pointField newPoints(mesh_.points());
//    {
//        List<pointIndexHit> hitInfo;
//        labelList hitSurface;
//        surfaces_.findNearest
//        (
//            surfaces_.getUnnamedSurfaces(),
//            localPoints,
//            scalarField(localPoints.size(), sqr(GREAT)),// sqr of attraction
//            hitSurface,
//            hitInfo
//        );
//
//        forAll(hitInfo, i)
//        {
//            if (hitInfo[i].hit())
//            {
//                //label pointI = meshPoints[i];
//                //Pout<< "   " << pointI << " moved from "
//                //    << mesh_.points()[pointI] << " by "
//                //    << mag(hitInfo[i].hitPoint()-mesh_.points()[pointI])
//                //    << endl;
//                newPoints[meshPoints[i]] = hitInfo[i].hitPoint();
//            }
//        }
//    }
//
//    // Per face (internal or coupled!) the patch that the
//    // baffle should get (or -1).
//    labelList facePatch(mesh_.nFaces(), -1);
//    // Count of baffled faces
//    label nBaffleFaces = 0;
//
//    {
//        pointField oldPoints(mesh_.points());
//        mesh_.movePoints(newPoints);
//        faceSet wrongFaces(mesh_, "wrongFaces", 100);
//        {
//            //motionSmoother::checkMesh(false, mesh_, motionDict, wrongFaces);
//
//            // Just check the errors from squashing
//            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//            const labelList allFaces(identity(mesh_.nFaces()));
//            label nWrongFaces = 0;
//
//            scalar minArea(readScalar(motionDict.lookup("minArea")));
//            if (minArea > -SMALL)
//            {
//                polyMeshGeometry::checkFaceArea
//                (
//                    false,
//                    minArea,
//                    mesh_,
//                    mesh_.faceAreas(),
//                    allFaces,
//                    &wrongFaces
//                );
//
//                label nNewWrongFaces = returnReduce
//                (
//                    wrongFaces.size(),
//                    sumOp<label>()
//                );
//
//                Info<< "    faces with area < "
//                    << setw(5) << minArea
//                    << " m^2                            : "
//                    << nNewWrongFaces-nWrongFaces << endl;
//
//                nWrongFaces = nNewWrongFaces;
//            }
//
////            scalar minDet(readScalar(motionDict.lookup("minDeterminant")));
//            scalar minDet = 0.01;
//            if (minDet > -1)
//            {
//                polyMeshGeometry::checkCellDeterminant
//                (
//                    false,
//                    minDet,
//                    mesh_,
//                    mesh_.faceAreas(),
//                    allFaces,
//                    polyMeshGeometry::affectedCells(mesh_, allFaces),
//                    &wrongFaces
//                );
//
//                label nNewWrongFaces = returnReduce
//                (
//                    wrongFaces.size(),
//                    sumOp<label>()
//                );
//
//                Info<< "    faces on cells with determinant < "
//                    << setw(5) << minDet << "                : "
//                    << nNewWrongFaces-nWrongFaces << endl;
//
//                nWrongFaces = nNewWrongFaces;
//            }
//        }
//
//
//        forAllConstIter(faceSet, wrongFaces, iter)
//        {
//            label patchI = mesh_.boundaryMesh().whichPatch(iter.key());
//
//            if (patchI == -1 || mesh_.boundaryMesh()[patchI].coupled())
//            {
//                facePatch[iter.key()] = getBafflePatch(facePatch, iter.key());
//                nBaffleFaces++;
//
//                //Pout<< "    " << iter.key()
//                //    //<< " on patch " << mesh_.boundaryMesh()[patchI].name()
//                //    << " is destined for patch " << facePatch[iter.key()]
//                //    << endl;
//            }
//        }
//        // Restore points.
//        mesh_.movePoints(oldPoints);
//    }
//
//
//    Info<< "markFacesOnProblemCellsGeometric : marked "
//        << returnReduce(nBaffleFaces, sumOp<label>())
//        << " additional internal and coupled faces"
//        << " to be converted into baffles." << endl;
//
//    syncTools::syncFaceList
//    (
//        mesh_,
//        facePatch,
//        maxEqOp<label>(),
//        false               // no separation
//    );
//
//    return facePatch;
//}


// ************************************************************************* //
