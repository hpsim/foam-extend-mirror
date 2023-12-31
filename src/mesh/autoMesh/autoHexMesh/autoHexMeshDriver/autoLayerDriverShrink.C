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

Description
    Shrinking mesh (part of adding cell layers)

\*---------------------------------------------------------------------------*/

#include "autoLayerDriver.H"
#include "fvMesh.H"
#include "foamTime.H"
#include "pointFields.H"
#include "motionSmoother.H"
#include "pointData.H"
#include "PointEdgeWave.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculate inverse sum of edge weights (currently always 1.0)
void Foam::autoLayerDriver::sumWeights
(
    const PackedBoolList& isMasterEdge,
    const labelList& meshEdges,
    const labelList& meshPoints,
    const edgeList& edges,
    scalarField& invSumWeight
) const
{
    invSumWeight = 0;

    forAll(edges, edgeI)
    {
        if (isMasterEdge.get(meshEdges[edgeI]) == 1)
        {
            const edge& e = edges[edgeI];
            //scalar eWeight = edgeWeights[edgeI];
            scalar eWeight = 1.0;

            invSumWeight[e[0]] += eWeight;
            invSumWeight[e[1]] += eWeight;
        }
    }

    syncTools::syncPointList
    (
        meshRefiner_.mesh(),
        meshPoints,
        invSumWeight,
        plusEqOp<scalar>(),
        scalar(0.0),        // null value
        false               // no separation
    );

    forAll(invSumWeight, pointI)
    {
        scalar w = invSumWeight[pointI];

        if (w > 0.0)
        {
            invSumWeight[pointI] = 1.0/w;
        }
    }
}


// Smooth field on moving patch
void Foam::autoLayerDriver::smoothField
(
    const motionSmoother& meshMover,
    const PackedBoolList& isMasterEdge,
    const labelList& meshEdges,
    const scalarField& fieldMin,
    const label nSmoothDisp,
    scalarField& field
) const
{
    const indirectPrimitivePatch& pp = meshMover.patch();
    const edgeList& edges = pp.edges();
    const labelList& meshPoints = pp.meshPoints();

    scalarField invSumWeight(pp.nPoints());
    sumWeights
    (
        isMasterEdge,
        meshEdges,
        meshPoints,
        edges,
        invSumWeight
    );

    // Get smoothly varying patch field.
    Info<< "shrinkMeshDistance : Smoothing field ..." << endl;

    for (label iter = 0; iter < nSmoothDisp; iter++)
    {
        scalarField average(pp.nPoints());
        averageNeighbours
        (
            meshMover.mesh(),
            isMasterEdge,
            meshEdges,
            meshPoints,
            pp.edges(),
            invSumWeight,
            field,
            average
        );

        // Transfer to field
        forAll(field, pointI)
        {
            //full smoothing neighbours + point value
            average[pointI] = 0.5*(field[pointI]+average[pointI]);

            // perform monotonic smoothing
            if
            (
                average[pointI] < field[pointI]
             && average[pointI] >= fieldMin[pointI]
            )
            {
                field[pointI] = average[pointI];
            }
        }

        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            Info<< "    Iteration " << iter << "   residual "
                <<  gSum(mag(field-average))
                   /returnReduce(average.size(), sumOp<label>())
                << endl;
        }
    }
}


// Smooth normals on moving patch.
void Foam::autoLayerDriver::smoothPatchNormals
(
    const motionSmoother& meshMover,
    const PackedBoolList& isMasterEdge,
    const labelList& meshEdges,
    const label nSmoothDisp,
    pointField& normals
) const
{
    const indirectPrimitivePatch& pp = meshMover.patch();
    const edgeList& edges = pp.edges();
    const labelList& meshPoints = pp.meshPoints();

    // Calculate inverse sum of weights

    scalarField invSumWeight(pp.nPoints());
    sumWeights
    (
        isMasterEdge,
        meshEdges,
        meshPoints,
        edges,
        invSumWeight
    );

    // Get smoothly varying internal normals field.
    Info<< "shrinkMeshDistance : Smoothing normals ..." << endl;

    for (label iter = 0; iter < nSmoothDisp; iter++)
    {
        vectorField average(pp.nPoints());
        averageNeighbours
        (
            meshMover.mesh(),
            isMasterEdge,
            meshEdges,
            meshPoints,
            pp.edges(),
            invSumWeight,
            normals,
            average
        );

        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            Info<< "    Iteration " << iter << "   residual "
                <<  gSum(mag(normals-average))
                   /returnReduce(average.size(), sumOp<label>())
                << endl;
        }

        // Transfer to normals vector field
        forAll(average, pointI)
        {
            // full smoothing neighbours + point value
            average[pointI] = 0.5*(normals[pointI]+average[pointI]);
            normals[pointI] = average[pointI];
            normals[pointI] /= mag(normals[pointI]) + VSMALL;
        }
    }
}


// Smooth normals in interior.
void Foam::autoLayerDriver::smoothNormals
(
    const label nSmoothDisp,
    const PackedBoolList& isMasterEdge,
    const labelList& fixedPoints,
    pointVectorField& normals
) const
{
    // Get smoothly varying internal normals field.
    Info<< "shrinkMeshDistance : Smoothing normals ..." << endl;

    const fvMesh& mesh = meshRefiner_.mesh();
    const edgeList& edges = mesh.edges();

    // Points that do not change.
    PackedBoolList isFixedPoint(mesh.nPoints());

    // Internal points that are fixed
    forAll(fixedPoints, i)
    {
        label meshPointI = fixedPoints[i];
        isFixedPoint.set(meshPointI, 1);
    }

    // Correspondence between local edges/points and mesh edges/points
    const labelList meshEdges(identity(mesh.nEdges()));
    const labelList meshPoints(identity(mesh.nPoints()));

    // Calculate inverse sum of weights

    scalarField invSumWeight(mesh.nPoints(), 0);
    sumWeights
    (
        isMasterEdge,
        meshEdges,
        meshPoints,
        edges,
        invSumWeight
    );

    Info<< "shrinkMeshDistance : Smoothing normals in interior ..." << endl;

    for (label iter = 0; iter < nSmoothDisp; iter++)
    {
        vectorField average(mesh.nPoints());
        averageNeighbours
        (
            mesh,
            isMasterEdge,
            meshEdges,
            meshPoints,
            edges,
            invSumWeight,
            normals,
            average
        );

        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            Info<< "    Iteration " << iter << "   residual "
                <<  gSum(mag(normals-average))
                   /returnReduce(average.size(), sumOp<label>())
                << endl;
        }


        // Transfer to normals vector field
        forAll(average, pointI)
        {
            if (isFixedPoint.get(pointI) == 0)
            {
                //full smoothing neighbours + point value
                average[pointI] = 0.5*(normals[pointI]+average[pointI]);
                normals[pointI] = average[pointI];
                normals[pointI] /= mag(normals[pointI]) + VSMALL;
            }
        }
    }
}


// Tries and find a medial axis point. Done by comparing vectors to nearest
// wall point for both vertices of edge.
bool Foam::autoLayerDriver::isMaxEdge
(
    const List<pointData>& pointWallDist,
    const label edgeI,
    const scalar minCos
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const pointField& points = mesh.points();

    // Do not mark edges with one side on moving wall.

    const edge& e = mesh.edges()[edgeI];

    vector v0(points[e[0]] - pointWallDist[e[0]].origin());
    scalar magV0(mag(v0));

    if (magV0 < SMALL)
    {
        return false;
    }

    vector v1(points[e[1]] - pointWallDist[e[1]].origin());
    scalar magV1(mag(v1));

    if (magV1 < SMALL)
    {
        return false;
    }

    v0 /= magV0;
    v1 /= magV1;

    // Test angle.
    if ((v0 & v1) < minCos)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Stop layer growth where mesh wraps around edge with a
// large feature angle
void Foam::autoLayerDriver::handleFeatureAngleLayerTerminations
(
    const indirectPrimitivePatch& pp,
    const scalar minCos,
    List<extrudeMode>& extrudeStatus,
    pointField& patchDisp,
    labelList& patchNLayers,
    label& nPointCounter
) const
{
    // Mark faces that have all points extruded
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    boolList extrudedFaces(pp.size(), true);

    forAll(pp.localFaces(), faceI)
    {
        const face& f = pp.localFaces()[faceI];

        forAll(f, fp)
        {
            if (extrudeStatus[f[fp]] == NOEXTRUDE)
            {
                extrudedFaces[faceI] = false;
                break;
            }
        }
    }


    // Detect situation where two featureedge-neighbouring faces are partly or
    // not extruded and the edge itself is extruded. In this case unmark the
    // edge for extrusion.

    forAll(pp.edgeFaces(), edgeI)
    {
        const labelList& eFaces = pp.edgeFaces()[edgeI];

        if (eFaces.size() == 2)
        {
            const edge& e = pp.edges()[edgeI];
            label v0 = e[0];
            label v1 = e[1];

            if
            (
                extrudeStatus[v0] != NOEXTRUDE
             || extrudeStatus[v1] != NOEXTRUDE
            )
            {
                if (!extrudedFaces[eFaces[0]] || !extrudedFaces[eFaces[1]])
                {
                    const vector& n0 = pp.faceNormals()[eFaces[0]];
                    const vector& n1 = pp.faceNormals()[eFaces[1]];

                    if ((n0 & n1) < minCos)
                    {
                        if
                        (
                            unmarkExtrusion
                            (
                                v0,
                                patchDisp,
                                patchNLayers,
                                extrudeStatus
                            )
                        )
                        {
                            nPointCounter++;
                        }
                        if
                        (
                            unmarkExtrusion
                            (
                                v1,
                                patchDisp,
                                patchNLayers,
                                extrudeStatus
                            )
                        )
                        {
                            nPointCounter++;
                        }
                    }
                }
            }
        }
    }
}


// Find isolated islands (points, edges and faces and layer terminations)
// in the layer mesh and stop any layer growth at these points.
void Foam::autoLayerDriver::findIsolatedRegions
(
    const indirectPrimitivePatch& pp,
    const PackedBoolList& isMasterEdge,
    const labelList& meshEdges,
    const scalar minCosLayerTermination,
    scalarField& field,
    List<extrudeMode>& extrudeStatus,
    pointField& patchDisp,
    labelList& patchNLayers
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< "shrinkMeshDistance : Removing isolated regions ..." << endl;

    // Keep count of number of points unextruded
    label nPointCounter = 0;

    while (true)
    {
        // Stop layer growth where mesh wraps around edge with a
        // large feature angle
        handleFeatureAngleLayerTerminations
        (
            pp,
            minCosLayerTermination,

            extrudeStatus,
            patchDisp,
            patchNLayers,
            nPointCounter
        );


        // Do not extrude from point where all neighbouring
        // faces are not grown
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        boolList extrudedFaces(pp.size(), true);
        forAll(pp.localFaces(), faceI)
        {
            const face& f = pp.localFaces()[faceI];
            forAll(f, fp)
            {
                if (extrudeStatus[f[fp]] == NOEXTRUDE)
                {
                    extrudedFaces[faceI] = false;
                    break;
                }
            }
        }

        const labelListList& pointFaces = pp.pointFaces();

        boolList keptPoints(pp.nPoints(), false);
        forAll(keptPoints, patchPointI)
        {
            const labelList& pFaces = pointFaces[patchPointI];

            forAll(pFaces, i)
            {
                label faceI = pFaces[i];
                if (extrudedFaces[faceI])
                {
                    keptPoints[patchPointI] = true;
                    break;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            keptPoints,
            orEqOp<bool>(),
            false,              // null value
            false               // no separation
        );

        label nChanged = 0;

        forAll(keptPoints, patchPointI)
        {
            if (!keptPoints[patchPointI])
            {
                if
                (
                    unmarkExtrusion
                    (
                        patchPointI,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                   nPointCounter++;
                   nChanged++;
                }
           }
       }


       if (returnReduce(nChanged, sumOp<label>()) == 0)
       {
           break;
       }
   }

    const edgeList& edges = pp.edges();


    // Count number of mesh edges using a point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList isolatedPoint(pp.nPoints(),0);

    forAll(edges, edgeI)
    {
        if (isMasterEdge.get(meshEdges[edgeI]) == 1)
        {
            const edge& e = edges[edgeI];

            label v0 = e[0];
            label v1 = e[1];

            if (extrudeStatus[v1] != NOEXTRUDE)
            {
                isolatedPoint[v0] += 1;
            }
            if (extrudeStatus[v0] != NOEXTRUDE)
            {
                isolatedPoint[v1] += 1;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        isolatedPoint,
        plusEqOp<label>(),
        label(0), // null value
        false    // no separation
    );

    // stop layer growth on isolated faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    forAll(pp, faceI)
    {
        const face& f = pp.localFaces()[faceI];
        bool failed = false;
        forAll(f, fp)
        {
            if (isolatedPoint[f[fp]] > 2)
            {
                failed = true;
                break;
            }
        }
        bool allPointsExtruded = true;
        if (!failed)
        {
            forAll(f, fp)
            {
                if (extrudeStatus[f[fp]] == NOEXTRUDE)
                {
                    allPointsExtruded = false;
                    break;
                }
            }

            if (allPointsExtruded)
            {
                forAll(f, fp)
                {
                    if
                    (
                        unmarkExtrusion
                        (
                            f[fp],
                            patchDisp,
                            patchNLayers,
                            extrudeStatus
                        )
                    )
                    {
                        field[f[fp]] = 0.0;
                        nPointCounter++;
                    }
                }
            }
        }
    }

    reduce(nPointCounter, sumOp<label>());
    Info<< "Number isolated points extrusion stopped : "<< nPointCounter<< endl;
}


// Calculates medial axis fields:
// dispVec     : normal of nearest wall. Where this normal changes direction
//               defines the medial axis
// medialDist  : distance to medial axis
// medialRatio : ratio of medial distance to wall distance.
//               (1 at wall, 0 at medial axis)
void Foam::autoLayerDriver::medialAxisSmoothingInfo
(
    const motionSmoother& meshMover,
    const label nSmoothNormals,
    const label nSmoothSurfaceNormals,
    const scalar minMedianAxisAngleCos,

    pointVectorField& dispVec,
    pointScalarField& medialRatio,
    pointScalarField& medialDist
) const
{

    Info<< "medialAxisSmoothingInfo :"
        << " Calculate distance to Medial Axis ..." << endl;

    const polyMesh& mesh = meshMover.mesh();
    const pointField& points = mesh.points();

    const indirectPrimitivePatch& pp = meshMover.patch();
    const vectorField& faceNormals = pp.faceNormals();
    const labelList& meshPoints = pp.meshPoints();

    // Predetermine mesh edges
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Precalulate master edge (only relevant for shared edges)
    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
    // Precalculate meshEdge per pp edge
    labelList meshEdges(pp.nEdges());

    forAll(meshEdges, patchEdgeI)
    {
        const edge& e = pp.edges()[patchEdgeI];

        label v0 = pp.meshPoints()[e[0]];
        label v1 = pp.meshPoints()[e[1]];
        meshEdges[patchEdgeI] = meshTools::findEdge
        (
            mesh.edges(),
            mesh.pointEdges()[v0],
            v0,
            v1
        );
    }


    // Determine pointNormal
    // ~~~~~~~~~~~~~~~~~~~~~

    pointField pointNormals(pp.nPoints(), vector::zero);
    {
        labelList nPointFaces(pp.nPoints(), 0);

        forAll(faceNormals, faceI)
        {
            const face& f = pp.localFaces()[faceI];

            forAll(f, fp)
            {
                pointNormals[f[fp]] += faceNormals[faceI];
                nPointFaces[f[fp]] ++;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            pointNormals,
            plusEqOp<vector>(),
            vector::zero,       // null value
            false               // no separation
        );

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            nPointFaces,
            plusEqOp<label>(),
            label(0),           // null value
            false               // no separation
        );

        forAll(pointNormals, i)
        {
            pointNormals[i] /= nPointFaces[i];
        }
    }

    // Smooth patch normal vectors
    smoothPatchNormals
    (
        meshMover,
        isMasterEdge,
        meshEdges,
        nSmoothSurfaceNormals,
        pointNormals
    );


    // Calculate distance to pp points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Distance to wall
    List<pointData> pointWallDist(mesh.nPoints());


    // 1. Calculate distance to points where displacement is specified.
    {
        // Seed data.
        List<pointData> wallInfo(meshPoints.size());

        forAll(meshPoints, patchPointI)
        {
            label pointI = meshPoints[patchPointI];
            wallInfo[patchPointI] = pointData
            (
                points[pointI],
                0.0,
                pointI,                       // passive scalar
                pointNormals[patchPointI]     // surface normals
            );
        }

        // Do all calculations
        List<pointData> edgeWallDist(mesh.nEdges());
        PointEdgeWave<pointData> wallDistCalc
        (
            mesh,
            meshPoints,
            wallInfo,
            pointWallDist,
            edgeWallDist,
            mesh.globalData().nTotalPoints()    // max iterations
        );
    }

    // 2. Find points with max distance and transport information back to
    //    wall.
    {
        List<pointData> pointMedialDist(mesh.nPoints());
        List<pointData> edgeMedialDist(mesh.nEdges());

        // Seed point data.
        DynamicList<pointData> maxInfo(meshPoints.size());
        dynamicLabelList maxPoints(meshPoints.size());

        // 1. Medial axis points

        const edgeList& edges = mesh.edges();

        forAll(edges, edgeI)
        {
            if (isMaxEdge(pointWallDist, edgeI, minMedianAxisAngleCos))
            {
                // Both end points of edge have very different nearest wall
                // point. Mark both points as medial axis points.
                const edge& e = edges[edgeI];

                forAll(e, ep)
                {
                    label pointI = e[ep];

                    if (!pointMedialDist[pointI].valid())
                    {
                        maxPoints.append(pointI);
                        maxInfo.append
                        (
                            pointData
                            (
                                points[pointI],
                                0.0,
                                pointI,         // passive data
                                vector::zero    // passive data
                            )
                        );
                        pointMedialDist[pointI] = maxInfo[maxInfo.size()-1];
                    }
                }
            }
        }


        // 2. Seed non-adapt patches
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        labelHashSet adaptPatches(meshMover.adaptPatchIDs());

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if
            (
                !pp.coupled()
             && !isA<emptyPolyPatch>(pp)
             && !adaptPatches.found(patchI)
            )
            {
                Info<< "Inserting points on patch " << pp.name() << endl;

                const labelList& meshPoints = pp.meshPoints();

                forAll(meshPoints, i)
                {
                    label pointI = meshPoints[i];

                    if (!pointMedialDist[pointI].valid())
                    {
                        maxPoints.append(pointI);
                        maxInfo.append
                        (
                            pointData
                            (
                                points[pointI],
                                0.0,
                                pointI,         // passive data
                                vector::zero    // passive data
                            )
                        );
                        pointMedialDist[pointI] = maxInfo[maxInfo.size()-1];
                    }
                }
            }
        }

        maxInfo.shrink();
        maxPoints.shrink();

        // Do all calculations
        PointEdgeWave<pointData> medialDistCalc
        (
            mesh,
            maxPoints,
            maxInfo,

            pointMedialDist,
            edgeMedialDist,
            mesh.globalData().nTotalPoints()    // max iterations
        );

        // Extract medial axis distance as pointScalarField
        forAll(pointMedialDist, pointI)
        {
            medialDist[pointI] = Foam::sqrt(pointMedialDist[pointI].distSqr());
        }
    }

    // Extract transported surface normals as pointVectorField
    forAll(dispVec, i)
    {
        dispVec[i] = pointWallDist[i].v();
    }

    // Smooth normal vectors. Do not change normals on pp.meshPoints
    smoothNormals(nSmoothNormals, isMasterEdge, meshPoints, dispVec);

    // Calculate ratio point medial distance to point wall distance
    forAll(medialRatio, pointI)
    {
        scalar wDist2 = pointWallDist[pointI].distSqr();
        scalar mDist = medialDist[pointI];

        if (wDist2 < sqr(SMALL) && mDist < SMALL)
        {
            medialRatio[pointI] = 0.0;
        }
        else
        {
            medialRatio[pointI] = mDist / (Foam::sqrt(wDist2) + mDist);
        }
    }

    if (debug)
    {
        Info<< "medialAxisSmoothingInfo :"
            << " Writing:" << nl
            << "    " << dispVec.name()
            << " : normalised direction of nearest displacement" << nl
            << "    " << medialDist.name()
            << " : distance to medial axis" << nl
            << "    " << medialRatio.name()
            << " : ratio of medial distance to wall distance" << nl
            << endl;
        dispVec.write();
        medialDist.write();
        medialRatio.write();
    }
}


void Foam::autoLayerDriver::shrinkMeshMedialDistance
(
    motionSmoother& meshMover,
    const dictionary& meshQualityDict,
    const label nSmoothThickness,
    const scalar maxThicknessToMedialRatio,
    const label nAllowableErrors,
    const label nSnap,
    const scalar minCosLayerTermination,

    const scalarField& layerThickness,
    const scalarField& minThickness,

    const pointVectorField& dispVec,
    const pointScalarField& medialRatio,
    const pointScalarField& medialDist,

    List<extrudeMode>& extrudeStatus,
    pointField& patchDisp,
    labelList& patchNLayers
) const
{
    Info<< "shrinkMeshMedialDistance : Smoothing using Medial Axis ..." << endl;

    const polyMesh& mesh = meshMover.mesh();

    const indirectPrimitivePatch& pp = meshMover.patch();
    const labelList& meshPoints = pp.meshPoints();

    // Precalulate master edge (only relevant for shared edges)
    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
    // Precalculate meshEdge per pp edge
    labelList meshEdges(pp.nEdges());

    forAll(meshEdges, patchEdgeI)
    {
        const edge& e = pp.edges()[patchEdgeI];

        label v0 = pp.meshPoints()[e[0]];
        label v1 = pp.meshPoints()[e[1]];
        meshEdges[patchEdgeI] = meshTools::findEdge
        (
            mesh.edges(),
            mesh.pointEdges()[v0],
            v0,
            v1
        );
    }


    scalarField thickness(layerThickness.size());

    thickness = mag(patchDisp);

    forAll(thickness, patchPointI)
    {
        if (extrudeStatus[patchPointI] == NOEXTRUDE)
        {
            thickness[patchPointI] = 0.0;
        }
    }

    label numThicknessRatioExclude = 0;

    // reduce thickness where thickness/medial axis distance large
    forAll(meshPoints, patchPointI)
    {
        if (extrudeStatus[patchPointI] != NOEXTRUDE)
        {
            label pointI = meshPoints[patchPointI];

            scalar mDist = medialDist[pointI];

            scalar thicknessRatio = thickness[patchPointI]/(mDist+VSMALL);

            if (thicknessRatio > maxThicknessToMedialRatio)
            {
                // Truncate thickness.
                thickness[patchPointI] =
                    0.5*(minThickness[patchPointI]+thickness[patchPointI]);

                patchDisp[patchPointI] =
                    thickness[patchPointI]
                  * patchDisp[patchPointI]
                  / (mag(patchDisp[patchPointI]) + VSMALL);
                numThicknessRatioExclude++;
            }
        }
    }

    reduce(numThicknessRatioExclude, sumOp<label>());
    Info<< "shrinkMeshMedialDistance : Reduce layer thickness at "
        << numThicknessRatioExclude
        << " nodes where thickness to medial axis distance is large " << endl;

    // find points where layer growth isolated to a lone point, edge or face

    findIsolatedRegions
    (
        pp,
        isMasterEdge,
        meshEdges,
        minCosLayerTermination,
        thickness,
        extrudeStatus,
        patchDisp,
        patchNLayers
    );

    // smooth layer thickness on moving patch
    smoothField
    (
        meshMover,
        isMasterEdge,
        meshEdges,
        minThickness,
        nSmoothThickness,
        thickness
    );

    List<pointData> pointWallDist(mesh.nPoints());

    const pointField& points = mesh.points();
    // 1. Calculate distance to points where displacement is specified.
    // This wave is used to transport layer thickness
    {
        // Distance to wall and medial axis on edges.
        List<pointData> edgeWallDist(mesh.nEdges());
        labelList wallPoints(meshPoints.size());

        // Seed data.
        List<pointData> wallInfo(meshPoints.size());

        forAll(meshPoints, patchPointI)
        {
            label pointI = meshPoints[patchPointI];
            wallPoints[patchPointI] = pointI;
            wallInfo[patchPointI] = pointData
            (
                points[pointI],
                0.0,
                thickness[patchPointI],       // transport layer thickness
                vector::zero                  // passive vector
            );
        }

        // Do all calculations
        PointEdgeWave<pointData> wallDistCalc
        (
            mesh,
            wallPoints,
            wallInfo,
            pointWallDist,
            edgeWallDist,
            mesh.globalData().nTotalPoints()    // max iterations
        );
    }

    // Calculate scaled displacement vector
    pointVectorField& displacement = meshMover.displacement();

    forAll(displacement, pointI)
    {
        // 1) displacement on nearest wall point, scaled by medialRatio
        //    (wall distance / medial distance)
        // 2) pointWallDist[pointI].s() is layer thickness transported
        //    from closest wall point.
        // 3) shrink in opposite direction of addedPoints
        displacement[pointI] =
            -medialRatio[pointI]
            *pointWallDist[pointI].s()
            *dispVec[pointI];
    }

    // Current faces to check. Gets modified in meshMover.scaleMesh
    labelList checkFaces(identity(mesh.nFaces()));

    Info<< "shrinkMeshMedialDistance : Moving mesh ..." << endl;
    scalar oldErrorReduction = -1;

    for (label iter = 0; iter < 2*nSnap ; iter++)
    {
        Info<< "Iteration " << iter << endl;
        if (iter == nSnap)
        {
            Info<< "Displacement scaling for error reduction set to 0."
                << endl;
            oldErrorReduction = meshMover.setErrorReduction(0.0);
        }

        if
        (
            meshMover.scaleMesh
            (
                checkFaces,
                List<labelPair>(0),
                meshMover.paramDict(),
                meshQualityDict,
                true,
                nAllowableErrors
            )
        )
        {
            Info<< "shrinkMeshMedialDistance : Successfully moved mesh"
                << endl;
            break;
        }
    }

    if (oldErrorReduction >= 0)
    {
        meshMover.setErrorReduction(oldErrorReduction);
    }

    Info<< "shrinkMeshMedialDistance : Finished moving mesh ..." << endl;
}


// ************************************************************************* //
