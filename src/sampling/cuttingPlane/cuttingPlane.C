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

#include "objectRegistry.H"
#include "cuttingPlane.H"
#include "primitiveMesh.H"
#include "linePointRef.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// set values for what is close to zero and what is considered to
// be positive (and not just rounding noise)
//! @cond localScope
const Foam::scalar localSmall = Foam::SMALL;
const Foam::scalar localSmallPositive = 1e-3*Foam::SMALL;
//! @endcond localScope

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Find cut cells
void Foam::cuttingPlane::calcCutCells
(
    const primitiveMesh& mesh,
    const scalarField& dotProducts,
    const UList<label>& cellIdLabels
)
{
    const labelListList& cellEdges = mesh.cellEdges();
    const edgeList& edges = mesh.edges();

    label listSize = cellEdges.size();
    if (!cellIdLabels.empty())
    {
        listSize = cellIdLabels.size();
    }

    cutCells_.setSize(listSize);
    label cutcellI(0);

    // Find the cut cells by detecting any cell that uses points with
    // opposing dotProducts.
    for (label listI = 0; listI < listSize; ++listI)
    {
        label cellI = listI;
        if (!cellIdLabels.empty())
        {
            cellI = cellIdLabels[listI];
        }

        const labelList& cEdges = cellEdges[cellI];

        label nCutEdges = 0;

        forAll(cEdges, i)
        {
            const edge& e = edges[cEdges[i]];

            if
            (
                (
                    dotProducts[e[0]] < localSmall
                 && dotProducts[e[1]] > localSmallPositive
                )
             || (
                    dotProducts[e[1]] < localSmall
                 && dotProducts[e[0]] > localSmallPositive
                )
            )
            {
                nCutEdges++;

                if (nCutEdges > 2)
                {
                    cutCells_[cutcellI++] = cellI;

                    break;
                }
            }
        }
    }

    // Set correct list size
    cutCells_.setSize(cutcellI);
}


// Determine for each edge the intersection point. Calculates
// - cutPoints_ : coordinates of all intersection points
// - edgePoint  : per edge -1 or the index into cutPoints
void Foam::cuttingPlane::intersectEdges
(
    const primitiveMesh& mesh,
    const scalarField& dotProducts,
    labelList& edgePoint
)
{
    // Use the dotProducts to find out the cut edges.
    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();

    // Per edge -1 or the label of the intersection point
    edgePoint.setSize(edges.size());

    DynamicList<point> dynCuttingPoints(4*cutCells_.size());

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        if
        (
            (
                dotProducts[e[0]] < localSmall
             && dotProducts[e[1]] > localSmallPositive
            )
         || (
                dotProducts[e[1]] < localSmall
             && dotProducts[e[0]] > localSmallPositive
            )
        )
        {
            // Edge is cut
            edgePoint[edgeI] = dynCuttingPoints.size();

            const point& p0 = points[e[0]];
            const point& p1 = points[e[1]];

            scalar alpha = lineIntersect(linePointRef(p0, p1));

            if (alpha < localSmall)
            {
                dynCuttingPoints.append(p0);
            }
            else if (alpha >= 1.0)
            {
                dynCuttingPoints.append(p1);
            }
            else
            {
                dynCuttingPoints.append((1-alpha)*p0 + alpha*p1);
            }
        }
        else
        {
            edgePoint[edgeI] = -1;
        }
    }

    this->storedPoints().transfer(dynCuttingPoints);
}


// Coming from startEdgeI cross the edge to the other face
// across to the next cut edge.
bool Foam::cuttingPlane::walkCell
(
    const primitiveMesh& mesh,
    const UList<label>& edgePoint,
    const label cellI,
    const label startEdgeI,
    dynamicLabelList& faceVerts
)
{
    label faceI = -1;
    label edgeI = startEdgeI;

    label nIter = 0;

    faceVerts.clear();
    do
    {
        faceVerts.append(edgePoint[edgeI]);

        // Cross edge to other face
        faceI = meshTools::otherFace(mesh, cellI, faceI, edgeI);

        // Find next cut edge on face.
        const labelList& fEdges = mesh.faceEdges()[faceI];

        label nextEdgeI = -1;

        //Note: here is where we should check for whether there are more
        // than 2 intersections with the face (warped/non-convex face).
        // If so should e.g. decompose the cells on both faces and redo
        // the calculation.

        forAll(fEdges, i)
        {
            label edge2I = fEdges[i];

            if (edge2I != edgeI && edgePoint[edge2I] != -1)
            {
                nextEdgeI = edge2I;
                break;
            }
        }

        if (nextEdgeI == -1)
        {
            // Did not find another cut edge on faceI. Do what?
            WarningIn("Foam::cuttingPlane::walkCell")
                << "Did not find closed walk along surface of cell " << cellI
                << " starting from edge " << startEdgeI
                << " in " << nIter << " iterations." << nl
                << "Collected "  << faceVerts.size() << " cutPoints so far"
                << endl;

            return false;
        }

        edgeI = nextEdgeI;

        nIter++;

        // Note: reduce number of points in the cut before error
        if (nIter > 100)
        {
            WarningIn("Foam::cuttingPlane::walkCell")
                << "Did not find closed walk along surface of cell " << cellI
                << " starting from edge " << startEdgeI
                << " in " << nIter << " iterations." << nl
                << "Collected "  << faceVerts.size() << " cutPoints so far"
                << endl;
            return false;
        }

    } while (edgeI != startEdgeI);


    if (faceVerts.size() >= 3)
    {
        return true;
    }
    else
    {
        WarningIn("Foam::cuttingPlane::walkCell")
            << "Did not find closed walk along surface of cell " << cellI
            << " starting from edge " << startEdgeI << nl
            << "Collected "  << faceVerts.size() << " cutPoints so far"
            << endl;

        return false;
    }
}


// For every cut cell determine a walk through all? its cuts.
void Foam::cuttingPlane::walkCellCuts
(
    const primitiveMesh& mesh,
    const UList<label>& edgePoint
)
{
    const pointField& cutPoints = this->points();

    // use dynamic lists to handle triangulation and/or missed cuts
    DynamicList<face>  dynCutFaces(cutCells_.size());
    dynamicLabelList dynCutCells(cutCells_.size());

    // scratch space for calculating the face vertices
    dynamicLabelList faceVerts(10);

    forAll(cutCells_, i)
    {
        label cellI = cutCells_[i];

        // Find the starting edge to walk from.
        const labelList& cEdges = mesh.cellEdges()[cellI];

        label startEdgeI = -1;

        forAll(cEdges, cEdgeI)
        {
            label edgeI = cEdges[cEdgeI];

            if (edgePoint[edgeI] != -1)
            {
                startEdgeI = edgeI;
                break;
            }
        }

        // Check for the unexpected ...
        if (startEdgeI == -1)
        {
            FatalErrorIn("Foam::cuttingPlane::walkCellCuts")
                << "Cannot find cut edge for cut cell " << cellI
                << abort(FatalError);
        }

        // Walk from starting edge around the circumference of the cell.
        bool okCut = walkCell
        (
            mesh,
            edgePoint,
            cellI,
            startEdgeI,
            faceVerts
        );

        if (okCut)
        {
            face f(faceVerts);

            // Orient face to point in the same direction as the plane normal
            if ((f.normal(cutPoints) && normal()) < 0)
            {
                f = f.reverseFace();
            }

            // the cut faces are usually quite ugly, so always triangulate
            label nTri = f.triangles(cutPoints, dynCutFaces);
            while (nTri--)
            {
                dynCutCells.append(cellI);
            }
        }
    }

    this->storedFaces().transfer(dynCutFaces);
    cutCells_.transfer(dynCutCells);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct without cutting
Foam::cuttingPlane::cuttingPlane(const plane& pln)
:
    plane(pln)
{}


// Construct from plane and mesh reference, restricted to a list of cells
Foam::cuttingPlane::cuttingPlane
(
    const plane& pln,
    const primitiveMesh& mesh,
    const UList<label>& cellIdLabels
)
:
    plane(pln)
{
    reCut(mesh, cellIdLabels);
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// recut mesh with existing planeDesc
void Foam::cuttingPlane::reCut
(
    const primitiveMesh& mesh,
    const UList<label>& cellIdLabels
)
{
    MeshStorage::clear();
    cutCells_.clear();

    scalarField dotProducts = (mesh.points() - refPoint()) & normal();

    // Determine cells that are (probably) cut.
    calcCutCells(mesh, dotProducts, cellIdLabels);

    // Determine cutPoints and return list of edge cuts.
    // per edge -1 or the label of the intersection point
    labelList edgePoint;
    intersectEdges(mesh, dotProducts, edgePoint);

    // Do topological walk around cell to find closed loop.
    walkCellCuts(mesh, edgePoint);
}


// remap action on triangulation
void Foam::cuttingPlane::remapFaces
(
    const UList<label>& faceMap
)
{
    // recalculate the cells cut
    if (!faceMap.empty())
    {
        MeshStorage::remapFaces(faceMap);

        labelList newCutCells(faceMap.size());
        forAll(faceMap, faceI)
        {
            newCutCells[faceI] = cutCells_[faceMap[faceI]];
        }
        cutCells_.transfer(newCutCells);
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::cuttingPlane::operator=(const cuttingPlane& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn ("Foam::cuttingPlane::operator=(const cuttingPlane&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    static_cast<MeshStorage&>(*this) = rhs;
    static_cast<plane&>(*this) = rhs;
    cutCells_ = rhs.cutCells();
}


// ************************************************************************* //
