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
    Converts .ele and .node and .face files, written by tetgen.

    Make sure to use add boundary attributes to the smesh file
    (5 fifth column in the element section)
    and run tetgen with -f option.

    Sample smesh file:

        # cube.smesh -- A 10x10x10 cube
        8 3
        1    0 0 0
        2    0 10 0
        3    10 10 0
        4    10 0 0
        5    0 0 10
        6    0 10 10
        7    10 10 10
        8    10 0 10
        6 1                 # 1 for boundary info present
        4    1 2 3 4 11  # region number 11
        4    5 6 7 8 21  # region number 21
        4    1 2 6 5 3
        4    4 3 7 8 43
        4    1 5 8 4 5
        4    2 6 7 3 65
        0
        0

NOTE:
- for some reason boundary faces point inwards. I just reverse them
always. Might use some geometric check instead.
- marked faces might not actually be boundary faces of mesh.
This is hopefully handled now by first constructing without boundaries
and then reconstructing with boundary faces.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "cellModeller.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Find label of face.
label findFace(const primitiveMesh& mesh, const face& f)
{
    const labelList& pFaces = mesh.pointFaces()[f[0]];

    forAll(pFaces, i)
    {
        label faceI = pFaces[i];

        if (mesh.faces()[faceI] == f)
        {
            return faceI;
        }
    }

    FatalErrorIn("findFace(const primitiveMesh&, const face&)")
        << "Cannot find face " << f << " in mesh." << abort(FatalError);

    return -1;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("file prefix");
    argList::validOptions.insert("noFaceFile", "");

#   include "setRootCase.H"
#   include "createTime.H"


    bool readFaceFile = !args.optionFound("noFaceFile");

    fileName prefix(args.additionalArgs()[0]);

    fileName nodeFile(prefix + ".node");
    fileName eleFile(prefix + ".ele");
    fileName faceFile(prefix + ".face");

    if (!readFaceFile)
    {
        Info<< "Files:" << endl
            << "    nodes : " << nodeFile << endl
            << "    elems : " << eleFile << endl
            << endl;
    }
    else
    {
        Info<< "Files:" << endl
            << "    nodes : " << nodeFile << endl
            << "    elems : " << eleFile << endl
            << "    faces : " << faceFile << endl
            << endl;

        Info<< "Reading .face file for boundary information" << nl << endl;
    }

    if (!isFile(nodeFile) || !isFile(eleFile))
    {
        FatalErrorIn(args.executable())
            << "Cannot read " << nodeFile << " or " << eleFile
            << exit(FatalError);
    }

    if (readFaceFile && !isFile(faceFile))
    {
        FatalErrorIn(args.executable())
            << "Cannot read " << faceFile << endl
            << "Did you run tetgen with -f option?" << endl
            << "If you don't want to read the .face file and thus not have"
            << " patches please\nrerun with the -noFaceFile option"
            << exit(FatalError);
    }


    IFstream nodeStream(nodeFile);

    //
    // Read nodes.
    //

    // Read header.
    string line;

    do
    {
        nodeStream.getLine(line);
    }
    while (line.size() && line[0] == '#');

    IStringStream nodeLine(line);

    label nNodes, nDims, nNodeAttr;
    bool hasRegion;

    nodeLine >> nNodes >> nDims >> nNodeAttr >> hasRegion;


    Info<< "Read .node header:" << endl
        << "    nodes     : " << nNodes << endl
        << "    nDims     : " << nDims << endl
        << "    nAttr     : " << nNodeAttr << endl
        << "    hasRegion : " << hasRegion << endl
        << endl;

    //
    // read points
    //

    pointField points(nNodes);
    Map<label> nodeToPoint(nNodes);

    {
        labelList pointIndex(nNodes);

        label pointI = 0;

        while (nodeStream.good())
        {
            nodeStream.getLine(line);

            if (line.size() && line[0] != '#')
            {
                IStringStream nodeLine(line);

                label nodeI;
                scalar x, y, z;
                label dummy;

                nodeLine >> nodeI >> x >> y >> z;

                for (label i = 0; i < nNodeAttr; i++)
                {
                    nodeLine >> dummy;
                }

                if (hasRegion)
                {
                    nodeLine >> dummy;
                }

                // Store point and node number.
                points[pointI] = point(x, y, z);
                nodeToPoint.insert(nodeI, pointI);
                pointI++;
            }
        }
        if (pointI != nNodes)
        {
            FatalIOErrorIn(args.executable().c_str(), nodeStream)
                << "Only " << pointI << " nodes present instead of " << nNodes
                << " from header." << exit(FatalIOError);
        }
    }

    //
    // read elements
    //

    IFstream eleStream(eleFile);

    do
    {
        eleStream.getLine(line);
    }
    while (line.size() && line[0] == '#');

    IStringStream eleLine(line);

    label nTets, nPtsPerTet, nElemAttr;

    eleLine >> nTets >> nPtsPerTet >> nElemAttr;


    Info<< "Read .ele header:" << endl
        << "    tets         : " << nTets << endl
        << "    pointsPerTet : " << nPtsPerTet << endl
        << "    nAttr        : " << nElemAttr << endl
        << endl;

    if (nPtsPerTet != 4)
    {
        FatalIOErrorIn(args.executable().c_str(), eleStream)
            << "Cannot handle tets with "
            << nPtsPerTet << " points per tetrahedron in .ele file" << endl
            << "Can only handle tetrahedra with four points"
            << exit(FatalIOError);
    }

    if (nElemAttr != 0)
    {
        WarningIn(args.executable())
            << "Element attributes (third elemenent in .ele header)"
            << " not used" << endl;
    }



    const cellModel& tet = *(cellModeller::lookup("tet"));

    labelList tetPoints(4);

    cellShapeList cells(nTets);
    label cellI = 0;

    while (eleStream.good())
    {
        eleStream.getLine(line);

        if (line.size() && line[0] != '#')
        {
            IStringStream eleLine(line);

            label elemI;
            eleLine >> elemI;

            for (label i = 0; i < 4; i++)
            {
                label nodeI;
                eleLine >> nodeI;
                tetPoints[i] = nodeToPoint[nodeI];
            }

            cells[cellI++] = cellShape(tet, tetPoints);

            // Skip attributes
            for (label i = 0; i < nElemAttr; i++)
            {
                label dummy;

                eleLine >> dummy;
            }
        }
    }


    //
    // Construct mesh with default boundary only
    //

    autoPtr<polyMesh> meshPtr
    (
        new polyMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime.constant(),
                runTime
            ),
            xferCopy(points),
            cells,
            faceListList(0),
            wordList(0),    // boundaryPatchNames
            wordList(0),    // boundaryPatchTypes
            "defaultFaces",
            polyPatch::typeName,
            wordList(0)
        )
    );
    const polyMesh& mesh = meshPtr;



    if (readFaceFile)
    {
        label nPatches = 0;

        // List of Foam vertices per boundary face
        faceList boundaryFaces;

        // For each boundary faces the Foam patchID
        labelList boundaryPatch;

        //
        // read boundary faces
        //

        IFstream faceStream(faceFile);

        do
        {
            faceStream.getLine(line);
        }
        while (line.size() && line[0] == '#');

        IStringStream faceLine(line);

        label nFaces, nFaceAttr;

        faceLine >> nFaces >> nFaceAttr;


        Info<< "Read .face header:" << endl
            << "    faces : " << nFaces << endl
            << "    nAttr : " << nFaceAttr << endl
            << endl;


        if (nFaceAttr != 1)
        {
            FatalIOErrorIn(args.executable().c_str(), faceStream)
                << "Expect boundary markers to be"
                << " present in .face file." << endl
                << "This is the second number in the header which is now:"
                << nFaceAttr << exit(FatalIOError);
        }

        // List of Foam vertices per boundary face
        boundaryFaces.setSize(nFaces);

        // For each boundary faces the Foam patchID
        boundaryPatch.setSize(nFaces);
        boundaryPatch = -1;

        label faceI = 0;

        // Region to patch conversion
        Map<label> regionToPatch;

        face f(3);

        while (faceStream.good())
        {
            faceStream.getLine(line);

            if (line.size() && line[0] != '#')
            {
                IStringStream faceLine(line);

                label tetGenFaceI, dummy, region;

                faceLine >> tetGenFaceI;

                // Read face and reverse orientation (Foam needs outwards
                // pointing)
                for (label i = 0; i < 3; i++)
                {
                    label nodeI;
                    faceLine >> nodeI;
                    f[2-i] = nodeToPoint[nodeI];
                }


                if (findFace(mesh, f) >= mesh.nInternalFaces())
                {
                    boundaryFaces[faceI] = f;

                    if (nFaceAttr > 0)
                    {
                        // First attribute is the region number
                        faceLine >> region;


                        // Get Foam patchID and update region->patch table.
                        label patchI = 0;

                        Map<label>::iterator patchFind =
                            regionToPatch.find(region);

                        if (patchFind == regionToPatch.end())
                        {
                            patchI = nPatches;

                            Info<< "Mapping tetgen region " << region
                                << " to Foam patch "
                                << patchI << endl;

                            regionToPatch.insert(region, nPatches++);
                        }
                        else
                        {
                            patchI = patchFind();
                        }

                        boundaryPatch[faceI] = patchI;

                        // Skip remaining attributes
                        for (label i = 1; i < nFaceAttr; i++)
                        {
                            faceLine >> dummy;
                        }
                    }

                    faceI++;
                }
            }
        }


        // Trim
        boundaryFaces.setSize(faceI);
        boundaryPatch.setSize(faceI);


         // Print region to patch mapping
        Info<< "Regions:" << endl;

        for
        (
            Map<label>::const_iterator iter = regionToPatch.begin();
            iter != regionToPatch.end();
            ++iter
        )
        {
            Info<< "    region:" << iter.key() << '\t' << "patch:"
                << iter() << endl;
        }
        Info<< endl;


        // Storage for boundary faces
        faceListList patchFaces(nPatches);
        wordList patchNames(nPatches);

        forAll(patchNames, patchI)
        {
            patchNames[patchI] = word("patch") + name(patchI);
        }

        wordList patchTypes(nPatches, polyPatch::typeName);
        word defaultFacesName = "defaultFaces";
        word defaultFacesType = polyPatch::typeName;
        wordList patchPhysicalTypes(nPatches, polyPatch::typeName);


        // Sort boundaryFaces by patch using boundaryPatch.
        List<DynamicList<face> > allPatchFaces(nPatches);

        forAll(boundaryPatch, faceI)
        {
            label patchI = boundaryPatch[faceI];

            allPatchFaces[patchI].append(boundaryFaces[faceI]);
        }

        Info<< "Patch sizes:" << endl;

        forAll(allPatchFaces, patchI)
        {
            Info<< "    " << patchNames[patchI] << " : "
                << allPatchFaces[patchI].size() << endl;

            patchFaces[patchI].transfer(allPatchFaces[patchI]);
        }

        Info<< endl;


        meshPtr.reset
        (
            new polyMesh
            (
                IOobject
                (
                    polyMesh::defaultRegion,
                    runTime.constant(),
                    runTime
                ),
                xferMove(points),
                cells,
                patchFaces,
                patchNames,
                patchTypes,
                defaultFacesName,
                defaultFacesType,
                patchPhysicalTypes
            )
        );
    }


    Info<< "Writing mesh to " << runTime.constant() << endl << endl;

    meshPtr().write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
