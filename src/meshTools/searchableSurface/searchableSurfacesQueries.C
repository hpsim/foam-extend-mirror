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

\*----------------------------------------------------------------------------*/

#include "searchableSurfacesQueries.H"
#include "ListOps.H"
#include "OFstream.H"
#include "meshTools.H"
#include "DynamicField.H"
#include "pointConstraint.H"
#include "boundBox.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Foam::searchableSurfacesQueries, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchableSurfacesQueries::mergeHits
(
    const point& start,

    const label testI,                    // index of surface
    const List<pointIndexHit>& surfHits,  // hits on surface

    labelList& allSurfaces,
    List<pointIndexHit>& allInfo,
    scalarList& allDistSqr
)
{
    // Given current set of hits (allSurfaces, allInfo) merge in those coming
    // from surface surfI.

    // Precalculate distances
    scalarList surfDistSqr(surfHits.size());
    forAll(surfHits, i)
    {
        surfDistSqr[i] = magSqr(surfHits[i].hitPoint() - start);
    }

    forAll(surfDistSqr, i)
    {
        label index = findLower(allDistSqr, surfDistSqr[i]);

        // Check if equal to lower.
        if (index >= 0)
        {
            // Same. Do not count.
            //Pout<< "point:" << surfHits[i].hitPoint()
            //    << " considered same as:" << allInfo[index].hitPoint()
            //    << " within tol:" << mergeDist
            //    << endl;
        }
        else
        {
            // Check if equal to higher
            label next = index + 1;

            if (next < allDistSqr.size())
            {
                //Pout<< "point:" << surfHits[i].hitPoint()
                //    << " considered same as:" << allInfo[next].hitPoint()
                //    << " within tol:" << mergeDist
                //    << endl;
            }
            else
            {
                // Insert after index
                label sz = allSurfaces.size();
                allSurfaces.setSize(sz+1);
                allInfo.setSize(allSurfaces.size());
                allDistSqr.setSize(allSurfaces.size());
                // Make space.
                for (label j = sz-1; j > index; --j)
                {
                    allSurfaces[j+1] = allSurfaces[j];
                    allInfo[j+1] = allInfo[j];
                    allDistSqr[j+1] = allDistSqr[j];
                }
                // Insert new value
                allSurfaces[index+1] = testI;
                allInfo[index+1] = surfHits[i];
                allDistSqr[index+1] = surfDistSqr[i];
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Find any intersection
void Foam::searchableSurfacesQueries::findAnyIntersection
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const pointField& start,
    const pointField& end,
    labelList& hitSurfaces,
    List<pointIndexHit>& hitInfo
)
{
    hitSurfaces.setSize(start.size());
    hitSurfaces = -1;
    hitInfo.setSize(start.size());

    // Work arrays
    labelList hitMap(identity(start.size()));
    pointField p0(start);
    pointField p1(end);
    List<pointIndexHit> intersectInfo(start.size());

    forAll(surfacesToTest, testI)
    {
        // Do synchronised call to all surfaces.
        allSurfaces[surfacesToTest[testI]].findLineAny(p0, p1, intersectInfo);

        // Copy all hits into arguments, continue with misses
        label newI = 0;
        forAll(intersectInfo, i)
        {
            if (intersectInfo[i].hit())
            {
                hitInfo[hitMap[i]] = intersectInfo[i];
                hitSurfaces[hitMap[i]] = testI;
            }
            else
            {
                if (i != newI)
                {
                    hitMap[newI] = hitMap[i];
                    p0[newI] = p0[i];
                    p1[newI] = p1[i];
                }
                newI++;
            }
        }

        // All done? Note that this decision should be synchronised
        if (newI == 0)
        {
            break;
        }

        // Trim and continue
        hitMap.setSize(newI);
        p0.setSize(newI);
        p1.setSize(newI);
        intersectInfo.setSize(newI);
    }
}


void Foam::searchableSurfacesQueries::findAllIntersections
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const pointField& start,
    const pointField& end,
    labelListList& hitSurfaces,
    List<List<pointIndexHit> >& hitInfo
)
{
    // Note: maybe move the single-surface all intersections test into
    // searchable surface? Some of the tolerance issues might be
    // lessened.

    // 2. Currently calling searchableSurface::findLine with start==end
    //    is expected to find no intersection. Problem if it does.

    hitSurfaces.setSize(start.size());
    hitInfo.setSize(start.size());

    if (surfacesToTest.empty())
    {
        return;
    }

    // Test first surface
    allSurfaces[surfacesToTest[0]].findLineAll(start, end, hitInfo);

    // Set hitSurfaces and distance
    List<scalarList> hitDistSqr(hitInfo.size());
    forAll(hitInfo, pointI)
    {
        const List<pointIndexHit>& pHits = hitInfo[pointI];

        labelList& pSurfaces = hitSurfaces[pointI];
        pSurfaces.setSize(pHits.size());
        pSurfaces = 0;

        scalarList& pDistSqr = hitDistSqr[pointI];
        pDistSqr.setSize(pHits.size());
        forAll(pHits, i)
        {
            pDistSqr[i] = magSqr(pHits[i].hitPoint() - start[pointI]);
        }
    }


    if (surfacesToTest.size() > 1)
    {
        // Test the other surfaces and merge (according to distance from start).
        for (label testI = 1; testI < surfacesToTest.size(); testI++)
        {
            List<List<pointIndexHit> > surfHits;
            allSurfaces[surfacesToTest[testI]].findLineAll
            (
                start,
                end,
                surfHits
            );

            forAll(surfHits, pointI)
            {
                mergeHits
                (
                    start[pointI],          // Current segment

                    testI,                  // Surface and its hits
                    surfHits[pointI],

                    hitSurfaces[pointI],    // Merge into overall hit info
                    hitInfo[pointI],
                    hitDistSqr[pointI]
                );
            }
        }
    }
}


// Find nearest. Return -1 or nearest point
void Foam::searchableSurfacesQueries::findNearest
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearestSurfaces,
    List<pointIndexHit>& nearestInfo
)
{
    // Find nearest. Return -1 or nearest point

    if (samples.size() != nearestDistSqr.size())
    {
        FatalErrorInFunction
            << "Inconsistent sizes. samples:" << samples.size()
            << " search-radius:" << nearestDistSqr.size()
            << abort(FatalError);
    }

    // Initialise
    nearestSurfaces.setSize(samples.size());
    nearestSurfaces = -1;
    nearestInfo.setSize(samples.size());

    // Work arrays
    scalarField minDistSqr(nearestDistSqr);
    List<pointIndexHit> hitInfo(samples.size());

    forAll(surfacesToTest, testI)
    {
        allSurfaces[surfacesToTest[testI]].findNearest
        (
            samples,
            minDistSqr,
            hitInfo
        );

        // Update minDistSqr and arguments
        forAll(hitInfo, pointI)
        {
            if (hitInfo[pointI].hit())
            {
                minDistSqr[pointI] = magSqr
                (
                    hitInfo[pointI].hitPoint()
                  - samples[pointI]
                );
                nearestInfo[pointI] = hitInfo[pointI];
                nearestSurfaces[pointI] = testI;
            }
        }
    }
}


void Foam::searchableSurfacesQueries::findNearest
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const labelListList& regionIndices,

    const pointField& samples,
    const scalarField& nearestDistSqr,

    labelList& nearestSurfaces,
    List<pointIndexHit>& nearestInfo
)
{
    // Find nearest. Return -1 or nearest point

    if (samples.size() != nearestDistSqr.size())
    {
        FatalErrorInFunction
            << "Inconsistent sizes. samples:" << samples.size()
            << " search-radius:" << nearestDistSqr.size()
            << exit(FatalError);
    }


    if (regionIndices.empty())
    {
        findNearest
        (
            allSurfaces,
            surfacesToTest,
            samples,
            nearestDistSqr,
            nearestSurfaces,
            nearestInfo
        );
    }

    // Initialise
    nearestSurfaces.setSize(samples.size());
    nearestSurfaces = -1;
    nearestInfo.setSize(samples.size());

    // Work arrays
    scalarField minDistSqr(nearestDistSqr);
    List<pointIndexHit> hitInfo(samples.size());

    forAll(surfacesToTest, testI)
    {
        allSurfaces[surfacesToTest[testI]].findNearest
        (
            samples,
            minDistSqr,
            hitInfo
        );

        // Update minDistSqr and arguments
        forAll(hitInfo, pointi)
        {
            if (hitInfo[pointi].hit())
            {
                minDistSqr[pointi] = magSqr
                (
                    hitInfo[pointi].hitPoint()
                  - samples[pointi]
                );
                nearestInfo[pointi] = hitInfo[pointi];
                nearestSurfaces[pointi] = testI;
            }
        }
    }
}


void Foam::searchableSurfacesQueries::findNearest
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const pointField& start,
    const scalarField& distSqr,
    pointField& near,
    List<pointConstraint>& constraint,
    const label nIter
)
{
    // Multi-surface findNearest


    if (start.size() != distSqr.size())
    {
        FatalErrorInFunction << "Inconsistent sizes. samples:" << start.size()
            << " search-radius:" << distSqr.size()
            << exit(FatalError);
    }


    vectorField normal;
    List<pointIndexHit> info;

    allSurfaces[surfacesToTest[0]].findNearest(start, distSqr, info);
    allSurfaces[surfacesToTest[0]].getNormal(info, normal);

    // Extract useful info from initial start point
    near = start;
    forAll(info, i)
    {
        if (info[i].hit())
        {
            near[i] = info[i].hitPoint();
        }
    }

    // Store normal as constraint
    constraint.setSize(near.size());
    constraint = pointConstraint();
    forAll(constraint, i)
    {
        if (info[i].hit())
        {
            constraint[i].applyConstraint(normal[i]);
        }
    }

    if (surfacesToTest.size() >= 2)
    {
        // Work space
        //pointField near1;
        vectorField normal1;

        label surfi = 1;
        for (label iter = 0; iter < nIter; iter++)
        {
            // Find nearest on next surface
            const searchableSurface& s = allSurfaces[surfacesToTest[surfi]];

            // Update: info, normal1
            s.findNearest(near, distSqr, info);
            s.getNormal(info, normal1);

            // Move to intersection of
            //    - previous surface(s) : near+normal
            //    - current surface     : info+normal1
            forAll(near, i)
            {
                if (info[i].hit())
                {
                    if (normal[i] != vector::zero)
                    {
                        // Have previous hit. Find intersection
                        if (mag(normal[i]&normal1[i]) < 1.0-1e-6)
                        {
                            plane pl0(near[i], normal[i]);
                            plane pl1(info[i].hitPoint(), normal1[i]);

                            plane::ray r(pl0.planeIntersect(pl1));
                            vector n = r.dir() / mag(r.dir());

                            // Calculate vector to move onto intersection line
                            vector d(r.refPoint()-near[i]);
                            d -= (d&n)*n;

                            // Trim the max distance
                            scalar magD = mag(d);
                            if (magD > SMALL)
                            {
                                scalar maxDist = Foam::sqrt(distSqr[i]);
                                if (magD > maxDist)
                                {
                                    // Clip
                                    d /= magD;
                                    d *= maxDist;
                                }

                                near[i] += d;
                                normal[i] = normal1[i];
                                constraint[i].applyConstraint(normal1[i]);
                            }
                        }
                    }
                    else
                    {
                        // First hit
                        near[i] = info[i].hitPoint();
                        normal[i] = normal1[i];
                        constraint[i].applyConstraint(normal1[i]);
                    }
                }
            }

            // Step to next surface
            surfi = surfacesToTest.fcIndex(surfi);
        }
    }
}


void Foam::searchableSurfacesQueries::signedDistance
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    const volumeType::type illegalHandling,
    labelList& nearestSurfaces,
    scalarField& distance
)
{
    // Initialise
    distance.setSize(samples.size());
    distance = -GREAT;

    // Find nearest
    List<pointIndexHit> nearestInfo;
    findNearest
    (
        allSurfaces,
        surfacesToTest,
        samples,
        nearestDistSqr,
        nearestSurfaces,
        nearestInfo
    );

    // Determine sign of nearest. Sort by surface to do this.
    DynamicField<point> surfPoints(samples.size());
    DynamicList<label> surfIndices(samples.size());

    forAll(surfacesToTest, testI)
    {
        // Extract samples on this surface
        surfPoints.clear();
        surfIndices.clear();
        forAll(nearestSurfaces, i)
        {
            if (nearestSurfaces[i] == testI)
            {
                surfPoints.append(samples[i]);
                surfIndices.append(i);
            }
        }

        // Calculate sideness of these surface points
        List<volumeType> vType;
        allSurfaces[surfacesToTest[testI]].getVolumeType(surfPoints, vType);

        // Push back to original
        forAll(vType, i)
        {
            label pointi = surfIndices[i];
            scalar dist = mag(samples[pointi] - nearestInfo[pointi].hitPoint());

            const volumeType& vT = vType[i];

            if (vT == volumeType::OUTSIDE)
            {
                distance[pointi] = dist;
            }
            else if (vT == volumeType::INSIDE)
            {
                distance[i] = -dist;
            }
            else
            {
                switch (illegalHandling)
                {
                    case volumeType::OUTSIDE:
                    {
                        distance[pointi] = dist;
                        break;
                    }
                    case volumeType::INSIDE:
                    {
                        distance[pointi] = -dist;
                        break;
                    }
                    default:
                    {
                        FatalErrorInFunction
                            << "getVolumeType failure,"
                            << " neither INSIDE or OUTSIDE."
                            << " point:" << surfPoints[i]
                            << " surface:"
                            << allSurfaces[surfacesToTest[testI]].name()
                            << " vType:" << volumeType::names[vT]
                            << exit(FatalError);
                        break;
                    }
                }
            }
        }
    }
}


Foam::boundBox Foam::searchableSurfacesQueries::bounds
(
    const PtrList<searchableSurface>& allSurfaces,
    const labelUList& surfacesToTest
)
{
    pointField bbPoints(2*surfacesToTest.size());

    forAll(surfacesToTest, testI)
    {
        const searchableSurface& surface(allSurfaces[surfacesToTest[testI]]);

        bbPoints[2*testI] = surface.bounds().min();

        bbPoints[2*testI + 1] = surface.bounds().max();
    }

    return boundBox(bbPoints);
}


// ************************************************************************* //
