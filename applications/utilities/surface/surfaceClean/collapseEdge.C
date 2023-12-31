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

#include "collapseEdge.H"

//- Sets point neighbours of face to val
static void markPointNbrs
(
    const triSurface& surf,
    const label faceI,
    const bool val,
    boolList& okToCollapse
)
{
    const labelledTri& f = surf.localFaces()[faceI];

    forAll(f, fp)
    {
        const labelList& pFaces = surf.pointFaces()[f[fp]];

        forAll(pFaces, i)
        {
            okToCollapse[pFaces[i]] = false;
        }
    }
}


static triSurface pack
(
    const triSurface& surf,
    const pointField& localPoints,
    const labelList& pointMap
)
{
    List<labelledTri> newTriangles(surf.size());
    label newTriangleI = 0;

    forAll(surf, faceI)
    {
        const labelledTri& f = surf.localFaces()[faceI];

        label newA = pointMap[f[0]];
        label newB = pointMap[f[1]];
        label newC = pointMap[f[2]];

        if ((newA != newB) && (newA != newC) && (newB != newC))
        {
            newTriangles[newTriangleI++] =
                labelledTri(newA, newB, newC, f.region());
        }
    }
    newTriangles.setSize(newTriangleI);

    return triSurface(newTriangles, surf.patches(), localPoints);
}


// Collapses small edge to point, thus removing triangle.
label collapseEdge(triSurface& surf, const scalar minLen)
{
    label nTotalCollapsed = 0;

    while (true)
    {
        const pointField& localPoints = surf.localPoints();
        const List<labelledTri>& localFaces = surf.localFaces();


        // Mapping from old to new points
        labelList pointMap(surf.nPoints());
        forAll(pointMap, i)
        {
            pointMap[i] = i;
        }

        // Storage for new points.
        pointField newPoints(localPoints);

        // To protect neighbours of collapsed faces.
        boolList okToCollapse(surf.size(), true);
        label nCollapsed = 0;

        forAll(localFaces, faceI)
        {
            if (okToCollapse[faceI])
            {
                // Check edge lengths.
                const labelledTri& f = localFaces[faceI];

                forAll(f, fp)
                {
                    label v = f[fp];
                    label v1 = f[(fp+1) % 3];

                    if (mag(localPoints[v1] - localPoints[v]) < minLen)
                    {
                        // Collapse f[fp1] onto f[fp].
                        pointMap[v1] = v;
                        newPoints[v] = 0.5*(localPoints[v1] + localPoints[v]);

                        Pout<< "Collapsing triange " << faceI << " to edge mid "
                            << newPoints[v] << endl;

                        nCollapsed++;
                        okToCollapse[faceI] = false;

                        // Protect point neighbours from collapsing.
                        markPointNbrs(surf, faceI, false, okToCollapse);

                        break;
                    }
                }
            }
        }

        Pout<< "collapseEdge : collapsing " << nCollapsed << " triangles"
            << endl;

        nTotalCollapsed += nCollapsed;

        if (nCollapsed == 0)
        {
            break;
        }

        // Pack the triangles
        surf = pack(surf, newPoints, pointMap);
    }

    // Remove any unused vertices
    surf = triSurface(surf.localFaces(), surf.patches(), surf.localPoints());

    return nTotalCollapsed;
}


// ************************************************************************* //
