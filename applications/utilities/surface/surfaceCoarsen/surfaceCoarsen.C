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
    Surface coarsening using 'bunnylod':

    Polygon Reduction Demo
    By Stan Melax (c) 1998
    mailto:melax@cs.ualberta.ca
    http://www.cs.ualberta.ca/~melax


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "triSurface.H"
#include "OFstream.H"
#include "triFace.H"
#include "triFaceList.H"

// From bunnylod
#include "progmesh.h"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int mapVertex(::List<int>& collapse_map, int a, int mx)
{
    if (mx <= 0)
    {
        return 0;
    }
    while (a >= mx)
    {
        a = collapse_map[a];
    }
    return a;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("Foam surface file");
    argList::validArgs.append("reduction factor");
    argList::validArgs.append("Foam output file");
    argList args(argc, argv);

    fileName inFileName(args.additionalArgs()[0]);

    scalar reduction(readScalar(IStringStream(args.additionalArgs()[1])()));

    if (reduction <= 0 || reduction > 1)
    {
        FatalErrorIn(args.executable())
            << "Reduction factor " << reduction
            << " should be within 0..1" << endl
            << "(it is the reduction in number of vertices)"
            << exit(FatalError);
    }

    fileName outFileName(args.additionalArgs()[2]);

    Info<< "Input surface   :" << inFileName << endl
        << "Reduction factor:" << reduction << endl
        << "Output surface  :" << outFileName << endl << endl;

    const triSurface surf(inFileName);

    Info<< "Surface:" << endl;
    surf.writeStats(Info);
    Info<< endl;


    ::List< ::Vector> vert;     // global list of vertices
    ::List< ::tridata> tri;     // global list of triangles


    // Convert triSurface to progmesh format. Note: can use global point
    // numbering since surface read in from file.
    const pointField& pts = surf.points();

    forAll(pts, ptI)
    {
        const point& pt = pts[ptI];

        vert.Add( ::Vector(pt.x(), pt.y(), pt.z()));
    }

    forAll(surf, faceI)
    {
        const labelledTri& f = surf[faceI];

        tridata td;
        td.v[0]=f[0];
        td.v[1]=f[1];
        td.v[2]=f[2];
        tri.Add(td);
    }

    ::List<int> collapse_map;   // to which neighbor each vertex collapses
    ::List<int> permutation;

    ::ProgressiveMesh(vert,tri,collapse_map,permutation);

    // rearrange the vertex list
    ::List< ::Vector> temp_list;
    for(int i=0;i<vert.num;i++)
    {
        temp_list.Add(vert[i]);
    }
    for(int i=0;i<vert.num;i++)
    {
        vert[permutation[i]]=temp_list[i];
    }

    // update the changes in the entries in the triangle list
    for(int i=0;i<tri.num;i++)
    {
        for(int j=0;j<3;j++)
        {
            tri[i].v[j] = permutation[tri[i].v[j]];
        }
    }

    // Only get triangles with non-collapsed edges.
    int render_num = int(reduction * surf.nPoints());

    Info<< "Reducing to " << render_num << " vertices" << endl;


    // Storage for new surface.
    Foam::List<labelledTri> newTris(surf.size());

    label newI = 0;

    for (int i=0; i<tri.num; i++)
    {
        int p0 = mapVertex(collapse_map, tri[i].v[0], render_num);
        int p1 = mapVertex(collapse_map, tri[i].v[1], render_num);
        int p2 = mapVertex(collapse_map, tri[i].v[2], render_num);

        // note:  serious optimization opportunity here,
        //  by sorting the triangles the following "continue"
        //  could have been made into a "break" statement.
        if (p0 == p1 || p1 == p2 || p2 == p0)
        {
            continue;
        }

        newTris[newI++] = labelledTri(p0, p1, p2, 0);
    }
    newTris.setSize(newI);

    // Convert vert into pointField.
    pointField newPoints(vert.num);

    for(int i=0; i<vert.num; i++)
    {
        const ::Vector & v = vert[i];

        newPoints[i] = point(v.x, v.y, v.z);
    }

    triSurface surf2(newTris, newPoints);

    triSurface outSurf
    (
        surf2.localFaces(),
        surf2.patches(),
        surf2.localPoints()
    );

    Info<< "Coarsened surface:" << endl;
    surf2.writeStats(Info);
    Info<< endl;

    Info<< "Writing to file " << outFileName << endl << endl;

    surf2.write(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
