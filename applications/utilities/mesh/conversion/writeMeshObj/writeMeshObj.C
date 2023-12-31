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
    For mesh debugging: writes mesh as three separate OBJ files which can
    be viewed with e.g. javaview.

    meshPoints_XXX.obj : all points and edges as lines.
    meshFaceCentres_XXX.obj : all face centres.
    meshCellCentres_XXX.obj : all cell centres.

    patch_YYY_XXX.obj : all face centres of patch YYY

    Optional: - patch faces (as polygons) : patchFaces_YYY_XXX.obj
              - non-manifold edges : patchEdges_YYY_XXX.obj

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "OFstream.H"
#include "meshTools.H"
#include "cellSet.H"
#include "faceSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeOBJ(const point& pt, Ostream& os)
{
    os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
}

// All edges of mesh
void writePoints(const polyMesh& mesh, const fileName& timeName)
{
    label vertI = 0;

    fileName pointFile(mesh.time().path()/"meshPoints_" + timeName + ".obj");

    Info << "Writing mesh points and edges to " << pointFile << endl;

    OFstream pointStream(pointFile);

    forAll(mesh.points(), pointI)
    {
        writeOBJ(mesh.points()[pointI], pointStream);
        vertI++;
    }

    forAll(mesh.edges(), edgeI)
    {
        const edge& e = mesh.edges()[edgeI];

        pointStream << "l " << e.start() + 1 << ' ' << e.end() + 1 << nl;
    }
}


// Edges for subset of cells
void writePoints
(
    const polyMesh& mesh,
    const labelList& cellLabels,
    const fileName& timeName
)
{
    fileName fName(mesh.time().path()/"meshPoints_" + timeName + ".obj");

    Info << "Writing mesh points and edges to " << fName << endl;

    OFstream str(fName);

    // OBJ file vertex
    label vertI = 0;

    // From point to OBJ file vertex
    Map<label> pointToObj(6*cellLabels.size());

    forAll(cellLabels, i)
    {
        label cellI = cellLabels[i];

        const labelList& cEdges = mesh.cellEdges()[cellI];

        forAll(cEdges, cEdgeI)
        {
            const edge& e = mesh.edges()[cEdges[cEdgeI]];

            label v0;

            Map<label>::iterator e0Fnd = pointToObj.find(e[0]);

            if (e0Fnd == pointToObj.end())
            {
                meshTools::writeOBJ(str, mesh.points()[e[0]]);
                v0 = vertI++;
                pointToObj.insert(e[0], v0);
            }
            else
            {
                v0 = e0Fnd();
            }

            label v1;

            Map<label>::iterator e1Fnd = pointToObj.find(e[1]);

            if (e1Fnd == pointToObj.end())
            {
                meshTools::writeOBJ(str, mesh.points()[e[1]]);
                v1 = vertI++;
                pointToObj.insert(e[1], v1);
            }
            else
            {
                v1 = e1Fnd();
            }


            str << "l " << v0+1 << ' ' << v1+1 << nl;
        }
    }
}


// Edges of single cell
void writePoints
(
    const polyMesh& mesh,
    const label cellI,
    const fileName& timeName
)
{
    fileName fName
    (
        mesh.time().path()
      / "meshPoints_" + timeName + '_' + name(cellI) + ".obj"
    );

    Info << "Writing mesh points and edges to " << fName << endl;

    OFstream pointStream(fName);

    const cell& cFaces = mesh.cells()[cellI];

    meshTools::writeOBJ(pointStream, mesh.faces(), mesh.points(), cFaces);
}



// All face centres
void writeFaceCentres(const polyMesh& mesh,const fileName& timeName)
{
    fileName faceFile
    (
        mesh.time().path()
      / "meshFaceCentres_" + timeName + ".obj"
    );

    Info << "Writing mesh face centres to " << faceFile << endl;

    OFstream faceStream(faceFile);

    forAll(mesh.faceCentres(), faceI)
    {
        writeOBJ(mesh.faceCentres()[faceI], faceStream);
    }
}


void writeCellCentres(const polyMesh& mesh, const fileName& timeName)
{
    fileName cellFile
    (
        mesh.time().path()/"meshCellCentres_" + timeName + ".obj"
    );

    Info << "Writing mesh cell centres to " << cellFile << endl;

    OFstream cellStream(cellFile);

    forAll(mesh.cellCentres(), cellI)
    {
        writeOBJ(mesh.cellCentres()[cellI], cellStream);
    }
}


void writePatchCentres
(
    const polyMesh& mesh,
    const fileName& timeName
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        fileName faceFile
        (
            mesh.time().path()/"patch_" + pp.name() + '_' + timeName + ".obj"
        );

        Info << "Writing patch face centres to " << faceFile << endl;

        OFstream patchFaceStream(faceFile);

        forAll(pp.faceCentres(), faceI)
        {
            writeOBJ(pp.faceCentres()[faceI], patchFaceStream);
        }
    }
}


void writePatchFaces
(
    const polyMesh& mesh,
    const fileName& timeName
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        fileName faceFile
        (
            mesh.time().path()
          / "patchFaces_" + pp.name() + '_' + timeName + ".obj"
        );

        Info << "Writing patch faces to " << faceFile << endl;

        OFstream patchFaceStream(faceFile);

        forAll(pp.localPoints(), pointI)
        {
            writeOBJ(pp.localPoints()[pointI], patchFaceStream);
        }

        forAll(pp.localFaces(), faceI)
        {
            const face& f = pp.localFaces()[faceI];

            patchFaceStream<< 'f';

            forAll(f, fp)
            {
                patchFaceStream << ' ' << f[fp]+1;
            }
            patchFaceStream << nl;
        }
    }
}


void writePatchBoundaryEdges
(
    const polyMesh& mesh,
    const fileName& timeName
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        fileName edgeFile
        (
            mesh.time().path()
          / "patchEdges_" + pp.name() + '_' + timeName + ".obj"
        );

        Info << "Writing patch edges to " << edgeFile << endl;

        OFstream patchEdgeStream(edgeFile);

        forAll(pp.localPoints(), pointI)
        {
            writeOBJ(pp.localPoints()[pointI], patchEdgeStream);
        }

        for (label edgeI = pp.nInternalEdges(); edgeI < pp.nEdges(); edgeI++)
        {
            if (pp.edgeFaces()[edgeI].size() == 1)
            {
                const edge& e = pp.edges()[edgeI];

                patchEdgeStream<< "l " << e[0]+1 << ' ' << e[1]+1 << nl;
            }
        }
    }
}


void writePointCells
(
    const polyMesh& mesh,
    const label pointI,
    const fileName& timeName
)
{
    const labelList& pCells = mesh.pointCells()[pointI];

    labelHashSet allEdges(6*pCells.size());

    forAll(pCells, i)
    {
        const labelList& cEdges = mesh.cellEdges()[pCells[i]];

        forAll(cEdges, i)
        {
            allEdges.insert(cEdges[i]);
        }
    }


    fileName pFile
    (
        mesh.time().path()
      / "pointEdges_" + timeName + '_' + name(pointI) + ".obj"
    );

    Info << "Writing pointEdges to " << pFile << endl;

    OFstream pointStream(pFile);

    label vertI = 0;

    for
    (
        labelHashSet::const_iterator iter = allEdges.begin();
        iter != allEdges.end();
        ++iter
    )
    {
        const edge& e = mesh.edges()[iter.key()];

        meshTools::writeOBJ(pointStream, mesh.points()[e[0]]); vertI++;
        meshTools::writeOBJ(pointStream, mesh.points()[e[1]]); vertI++;
        pointStream<< "l " << vertI-1 << ' ' << vertI << nl;
    }
}


// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validOptions.insert("patchFaces", "");
    argList::validOptions.insert("patchEdges", "");
    argList::validOptions.insert("cell", "cellI");
    argList::validOptions.insert("face", "faceI");
    argList::validOptions.insert("point", "pointI");
    argList::validOptions.insert("cellSet", "setName");
    argList::validOptions.insert("faceSet", "setName");
#   include "addRegionOption.H"

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();

    bool patchFaces = args.optionFound("patchFaces");
    bool patchEdges = args.optionFound("patchEdges");
    bool doCell     = args.optionFound("cell");
    bool doPoint    = args.optionFound("point");
    bool doFace     = args.optionFound("face");
    bool doCellSet  = args.optionFound("cellSet");
    bool doFaceSet  = args.optionFound("faceSet");


    Info<< "Writing mesh objects as .obj files such that the object"
        << " numbering" << endl
        << "(for points, faces, cells) is consistent with"
        << " Foam numbering (starting from 0)." << endl << endl;

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createNamedPolyMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        polyMesh::readUpdateState state = mesh.readUpdate();

        if (!timeI || state != polyMesh::UNCHANGED)
        {
            if (patchFaces)
            {
                writePatchFaces(mesh, runTime.timeName());
            }
            if (patchEdges)
            {
                writePatchBoundaryEdges(mesh, runTime.timeName());
            }
            if (doCell)
            {
                label cellI = args.optionRead<label>("cell");

                writePoints(mesh, cellI, runTime.timeName());
            }
            if (doPoint)
            {
                label pointI = args.optionRead<label>("point");

                writePointCells(mesh, pointI, runTime.timeName());
            }
            if (doFace)
            {
                label faceI = args.optionRead<label>("face");

                fileName fName
                (
                    mesh.time().path()
                  / "meshPoints_"
                  + runTime.timeName()
                  + '_'
                  + name(faceI)
                  + ".obj"
                );

                Info<< "Writing mesh points and edges to " << fName << endl;

                OFstream str(fName);

                const face& f = mesh.faces()[faceI];

                meshTools::writeOBJ(str, faceList(1, f), mesh.points());
            }
            if (doCellSet)
            {
                word setName(args.option("cellSet"));

                cellSet cells(mesh, setName);

                Info<< "Read " << cells.size() << " cells from set " << setName
                    << endl;

                writePoints(mesh, cells.toc(), runTime.timeName());

            }
            if (doFaceSet)
            {
                word setName(args.option("faceSet"));

                faceSet faces(mesh, setName);

                Info<< "Read " << faces.size() << " faces from set " << setName
                    << endl;

                fileName fName
                (
                    mesh.time().path()
                  / "meshPoints_"
                  + runTime.timeName()
                  + '_'
                  + setName
                  + ".obj"
                );

                Info << "Writing mesh points and edges to " << fName << endl;

                OFstream str(fName);

                meshTools::writeOBJ
                (
                    str,
                    mesh.faces(),
                    mesh.points(),
                    faces.toc()
                );
            }
            else if
            (
                !patchFaces
             && !patchEdges
             && !doCell
             && !doPoint
             && !doFace
             && !doCellSet
             && !doFaceSet
            )
            {
                // points & edges
                writePoints(mesh, runTime.timeName());

                // face centres
                writeFaceCentres(mesh, runTime.timeName());

                // cell centres
                writeCellCentres(mesh, runTime.timeName());

                // Patch face centres
                writePatchCentres(mesh, runTime.timeName());
            }
        }
        else
        {
            Info << "No mesh." << endl;
        }

        Info << nl << endl;
    }


    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
