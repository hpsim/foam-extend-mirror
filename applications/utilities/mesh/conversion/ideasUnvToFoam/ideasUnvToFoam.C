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
    I-Deas unv format mesh conversion.

    Uses either
    - DOF sets (757) or
    - face groups (2452(Cubit), 2467)
    to do patching.
    Works without but then puts all faces in defaultFaces patch.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "foamTime.H"
#include "IFstream.H"
#include "cellModeller.H"
#include "cellSet.H"
#include "faceSet.H"
#include "DynamicList.H"
#include "triSurface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const string SEPARATOR("    -1");

bool isSeparator(const string& line)
{
    return line.substr(0, 6) == SEPARATOR;
}


// Reads past -1 and reads element type
label readTag(IFstream& is)
{
    string tag;
    do
    {
        if (!is.good())
        {
            return -1;
        }

        string line;

        is.getLine(line);

        if (line.size() < 6)
        {
            return -1;
        }

        tag = line.substr(0, 6);

    } while (tag == SEPARATOR);

    return readLabel(IStringStream(tag)());
}


// Reads and prints header
void readHeader(IFstream& is)
{
    string line;

    while (is.good())
    {
        is.getLine(line);

        if (isSeparator(line))
        {
            break;
        }
        else
        {
            Sout<< line << endl;
        }
    }
}


// Skip
void skipSection(IFstream& is)
{
    Sout<< "Skipping section at line " << is.lineNumber() << '.' << endl;

    string line;

    while (is.good())
    {
        is.getLine(line);

        if (isSeparator(line))
        {
            break;
        }
        else
        {
//            Sout<< line << endl;
        }
    }
}


scalar readUnvScalar(const string& unvString)
{
    string s(unvString);

    s.replaceAll("d", "E");
    s.replaceAll("D", "E");

    return readScalar(IStringStream(s)());
}


// Reads unit section
void readUnits
(
    IFstream& is,
    scalar& lengthScale,
    scalar& forceScale,
    scalar& tempScale,
    scalar& tempOffset
)
{
    Sout<< "Starting reading units at line " << is.lineNumber() << '.' << endl;

    string line;
    is.getLine(line);

    label l = readLabel(IStringStream(line.substr(0, 10))());
    Sout<< "l:" << l << endl;

    string units(line.substr(10, 20));
    Sout<< "units:" << units << endl;

    label unitType = readLabel(IStringStream(line.substr(30, 10))());
    Sout<< "unitType:" << unitType << endl;

    // Read lengthscales
    is.getLine(line);

    lengthScale = readUnvScalar(line.substr(0, 25));
    forceScale = readUnvScalar(line.substr(25, 25));
    tempScale = readUnvScalar(line.substr(50, 25));

    is.getLine(line);
    tempOffset = readUnvScalar(line.substr(0, 25));

    Sout<< "Unit factors:" << nl
        << "    Length scale       : " << lengthScale << nl
        << "    Force scale        : " << forceScale << nl
        << "    Temperature scale  : " << tempScale << nl
        << "    Temperature offset : " << tempOffset << nl
        << endl;
}


// Reads points section. Read region as well?
void readPoints
(
    IFstream& is,
    DynamicList<point>& points,     // coordinates
    dynamicLabelList& unvPointID  // unv index
)
{
    Sout<< "Starting reading points at line " << is.lineNumber() << '.' << endl;

    static bool hasWarned = false;

    while (true)
    {
        string line;
        is.getLine(line);

        label pointI = readLabel(IStringStream(line.substr(0, 10))());

        if (pointI == -1)
        {
            break;
        }
        else if (pointI != points.size()+1 && !hasWarned)
        {
            hasWarned = true;

            IOWarningIn
            (
                "readPoints(IFstream&, label&, DynamicList<point>"
                ", dynamicLabelList&)",
                is
            )   << "Points not in order starting at point " << pointI
                //<< " at line " << is.lineNumber()
                //<< abort(FatalError);
                << endl;
        }

        point pt;
        is.getLine(line);
        pt[0] = readUnvScalar(line.substr(0, 25));
        pt[1] = readUnvScalar(line.substr(25, 25));
        pt[2] = readUnvScalar(line.substr(50, 25));

        unvPointID.append(pointI);
        points.append(pt);
    }

    points.shrink();
    unvPointID.shrink();

    Sout<< "Read " << points.size() << " points." << endl;
}


// Reads cells section. Read region as well? Not handled yet but should just
// be a matter of reading corresponding to boundaryFaces the correct property
// and sorting it later on.
void readCells
(
    IFstream& is,
    DynamicList<cellShape>& cellVerts,
    dynamicLabelList& cellMaterial,
    dynamicLabelList& boundaryFaceIndices,
    DynamicList<face>& boundaryFaces
)
{
    Sout<< "Starting reading cells at line " << is.lineNumber() << '.' << endl;

    const cellModel& hex = *(cellModeller::lookup("hex"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& tet = *(cellModeller::lookup("tet"));

    labelHashSet skippedElements;

    labelHashSet foundFeType;

    while (true)
    {
        string line;
        is.getLine(line);

        if (isSeparator(line))
        {
            break;
        }

        IStringStream lineStr(line);
        label cellI, feID, physProp, matProp, colour, nNodes;
        lineStr >> cellI >> feID >> physProp >> matProp >> colour >> nNodes;

        if (foundFeType.insert(feID))
        {
            Info<< "First occurrence of element type " << feID
                << " for cell " << cellI << " at line "
                << is.lineNumber() << endl;
        }

        if (feID == 11)
        {
            // Rod. Skip.
            is.getLine(line);
            is.getLine(line);
        }
        else if (feID == 171)
        {
            // Rod. Skip.
            is.getLine(line);
        }
        else if (feID == 41 || feID == 91)
        {
            // Triangle. Save - used for patching later on.
            is.getLine(line);

            face cVerts(3);
            IStringStream lineStr(line);
            lineStr >> cVerts[0] >> cVerts[1] >> cVerts[2];
            boundaryFaces.append(cVerts);
            boundaryFaceIndices.append(cellI);
        }
        else if (feID == 44 || feID == 94)
        {
            // Quad. Save - used for patching later on.
            is.getLine(line);

            face cVerts(4);
            IStringStream lineStr(line);
            lineStr >> cVerts[0] >> cVerts[1] >> cVerts[2] >> cVerts[3];
            boundaryFaces.append(cVerts);
            boundaryFaceIndices.append(cellI);
        }
        else if (feID == 111)
        {
            // Tet.
            is.getLine(line);

            labelList cVerts(4);
            IStringStream lineStr(line);
            lineStr >> cVerts[0] >> cVerts[1] >> cVerts[2] >> cVerts[3];

            cellVerts.append(cellShape(tet, cVerts, true));
            cellMaterial.append(physProp);

            if (cellVerts[cellVerts.size()-1].size() != cVerts.size())
            {
                Pout<< "Line:" << is.lineNumber()
                    << " element:" << cellI
                    << " type:" << feID
                    << " collapsed from " << cVerts << nl
                    << " to:" << cellVerts[cellVerts.size()-1]
                    << endl;
            }
        }
        else if (feID == 112)
        {
            // Wedge.
            is.getLine(line);

            labelList cVerts(6);
            IStringStream lineStr(line);
            lineStr >> cVerts[0] >> cVerts[1] >> cVerts[2] >> cVerts[3]
                    >> cVerts[4] >> cVerts[5];

            cellVerts.append(cellShape(prism, cVerts, true));
            cellMaterial.append(physProp);

            if (cellVerts[cellVerts.size()-1].size() != cVerts.size())
            {
                Pout<< "Line:" << is.lineNumber()
                    << " element:" << cellI
                    << " type:" << feID
                    << " collapsed from " << cVerts << nl
                    << " to:" << cellVerts[cellVerts.size()-1]
                    << endl;
            }
        }
        else if (feID == 115)
        {
            // Hex.
            is.getLine(line);

            labelList cVerts(8);
            IStringStream lineStr(line);
            lineStr
                >> cVerts[0] >> cVerts[1] >> cVerts[2] >> cVerts[3]
                >> cVerts[4] >> cVerts[5] >> cVerts[6] >> cVerts[7];

            cellVerts.append(cellShape(hex, cVerts, true));
            cellMaterial.append(physProp);

            if (cellVerts[cellVerts.size()-1].size() != cVerts.size())
            {
                Pout<< "Line:" << is.lineNumber()
                    << " element:" << cellI
                    << " type:" << feID
                    << " collapsed from " << cVerts << nl
                    << " to:" << cellVerts[cellVerts.size()-1]
                    << endl;
            }
        }
        else
        {
            if (skippedElements.insert(feID))
            {
                IOWarningIn("readCells(IFstream&, label&)", is)
                    << "Cell type " << feID << " not supported" << endl;
            }
            is.getLine(line);  //Do nothing
        }
    }

    cellVerts.shrink();
    cellMaterial.shrink();
    boundaryFaces.shrink();
    boundaryFaceIndices.shrink();

    Sout<< "Read " << cellVerts.size() << " cells"
        << " and " << boundaryFaces.size() << " boundary faces." << endl;
}


void readPatches
(
    IFstream& is,
    DynamicList<word>& patchNames,
    DynamicList<labelList>& patchFaceIndices
)
{
    Sout<< "Starting reading patches at line " << is.lineNumber() << '.'
        << endl;

    while(true)
    {
        string line;
        is.getLine(line);

        if (isSeparator(line))
        {
            break;
        }

        IStringStream lineStr(line);
        label group, constraintSet, restraintSet, loadSet, dofSet,
            tempSet, contactSet, nFaces;
        lineStr
            >> group >> constraintSet >> restraintSet >> loadSet
            >> dofSet >> tempSet >> contactSet >> nFaces;

        is.getLine(line);
        word groupName = string::validate<word>(line);

        Info<< "For group " << group
            << " named " << groupName
            << " trying to read " << nFaces << " patch face indices."
            << endl;

        labelList groupIndices(nFaces);
        label groupType = -1;
        nFaces = 0;

        while (nFaces < groupIndices.size())
        {
            is.getLine(line);
            IStringStream lineStr(line);

            // Read one (for last face) or two entries from line.
            label nRead = 2;
            if (nFaces == groupIndices.size() - 1)
            {
                nRead = 1;
            }

            for (label i = 0; i < nRead; i++)
            {
                label tag, nodeLeaf, component;

                lineStr >> groupType >> tag >> nodeLeaf >> component;

                groupIndices[nFaces++] = tag;
            }
        }


        // Store
        if (groupType == 8)
        {
            patchNames.append(groupName);
            patchFaceIndices.append(groupIndices);
        }
        else
        {
            IOWarningIn("readPatches(..)", is)
                << "When reading patches expect entity type code 8"
                << nl << "    Skipping group code " << groupType
                << endl;
        }
    }

    patchNames.shrink();
    patchFaceIndices.shrink();
}



// Read dof set (757)
void readDOFS
(
    IFstream& is,
    DynamicList<word>& patchNames,
    DynamicList<labelList>& dofVertices
)
{
    Sout<< "Starting reading contraints at line " << is.lineNumber() << '.'
        << endl;

    string line;
    is.getLine(line);
    label group;
    {
        IStringStream lineStr(line);
        lineStr >> group;
    }

    is.getLine(line);
    {
        IStringStream lineStr(line);
        patchNames.append(lineStr);
    }

    Info<< "For DOF set " << group
        << " named " << patchNames[patchNames.size()-1]
        << " trying to read vertex indices."
        << endl;

    dynamicLabelList vertices;
    while (true)
    {
        string line;
        is.getLine(line);

        if (isSeparator(line))
        {
            break;
        }

        IStringStream lineStr(line);
        label pointI;
        lineStr >> pointI;

        vertices.append(pointI);
    }

    Info<< "For DOF set " << group
        << " named " << patchNames[patchNames.size()-1]
        << " read " << vertices.size() << " vertex indices." << endl;

    dofVertices.append(vertices.shrink());

    patchNames.shrink();
    dofVertices.shrink();
}


// Returns -1 or group that all of the vertices of f are in,
label findPatch(const List<labelHashSet>& dofGroups, const face& f)
{
    forAll(dofGroups, patchI)
    {
        if (dofGroups[patchI].found(f[0]))
        {
            bool allInGroup = true;

            // Check rest of face
            for (label fp = 1; fp < f.size(); fp++)
            {
                if (!dofGroups[patchI].found(f[fp]))
                {
                    allInGroup = false;
                    break;
                }
            }

            if (allInGroup)
            {
                return patchI;
            }
        }
    }
    return -1;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append(".unv file");
    argList::validOptions.insert("dump", "");

#   include "setRootCase.H"
#   include "createTime.H"

    fileName ideasName(args.additionalArgs()[0]);

    IFstream inFile(ideasName.c_str());

    if (!inFile.good())
    {
        FatalErrorIn(args.executable())
            << "Cannot open file " << ideasName
            << exit(FatalError);
    }


    // Unit scale factors
    scalar lengthScale = 1;
    scalar forceScale = 1;
    scalar tempScale = 1;
    scalar tempOffset = 0;

    // Points
    DynamicList<point> points;
    // Original unv point label
    dynamicLabelList unvPointID;

    // Cells
    DynamicList<cellShape> cellVerts;
    dynamicLabelList cellMat;

    // Boundary faces
    dynamicLabelList boundaryFaceIndices;
    DynamicList<face> boundaryFaces;

    // Patch names and patchFace indices.
    DynamicList<word> patchNames;
    DynamicList<labelList> patchFaceIndices;
    DynamicList<labelList> dofVertIndices;


    while (inFile.good())
    {
        label tag = readTag(inFile);

        if (tag == -1)
        {
            break;
        }

        Sout<< "Processing tag:" << tag << endl;

        switch (tag)
        {
            case 151:
                readHeader(inFile);
            break;

            case 164:
                readUnits
                (
                    inFile,
                    lengthScale,
                    forceScale,
                    tempScale,
                    tempOffset
                );
            break;

            case 2411:
                readPoints(inFile, points, unvPointID);
            break;

            case 2412:
                readCells
                (
                    inFile,
                    cellVerts,
                    cellMat,
                    boundaryFaceIndices,
                    boundaryFaces
                );
            break;

            case 2452:
            case 2467:
                readPatches
                (
                    inFile,
                    patchNames,
                    patchFaceIndices
                );
            break;

            case 757:
                readDOFS
                (
                    inFile,
                    patchNames,
                    dofVertIndices
                );
            break;

            default:
                Sout<< "Skipping tag " << tag << " on line "
                    << inFile.lineNumber() << endl;
                skipSection(inFile);
            break;
        }
        Sout<< endl;
    }


    // Invert point numbering.
    label maxUnvPoint = 0;
    forAll(unvPointID, pointI)
    {
        maxUnvPoint = max(maxUnvPoint, unvPointID[pointI]);
    }
    labelList unvToFoam(invert(maxUnvPoint+1, unvPointID));


    // Renumber vertex numbers on cells

    forAll(cellVerts, cellI)
    {
        labelList foamVerts
        (
            renumber
            (
                unvToFoam,
                static_cast<labelList&>(cellVerts[cellI])
            )
        );

        if (findIndex(foamVerts, -1) != -1)
        {
            FatalErrorIn(args.executable())
                << "Cell " << cellI
                << " unv vertices " << cellVerts[cellI]
                << " has some undefined vertices " << foamVerts
                << abort(FatalError);
        }

        // Bit nasty: replace vertex list.
        cellVerts[cellI].transfer(foamVerts);
    }

    // Renumber vertex numbers on boundaryFaces

    forAll(boundaryFaces, bFaceI)
    {
        labelList foamVerts(renumber(unvToFoam, boundaryFaces[bFaceI]));

        if (findIndex(foamVerts, -1) != -1)
        {
            FatalErrorIn(args.executable())
                << "Boundary face " << bFaceI
                << " unv vertices " << boundaryFaces[bFaceI]
                << " has some undefined vertices " << foamVerts
                << abort(FatalError);
        }

        // Bit nasty: replace vertex list.
        boundaryFaces[bFaceI].transfer(foamVerts);
    }



    // Patchify = create patchFaceVerts for use in cellShape construction
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Works in one of two modes. Either has read boundaryFaces which
    // just need to be sorted according to patch. Or has read point constraint
    // sets (dofVertIndices).

    List<faceList> patchFaceVerts;


    if (dofVertIndices.size())
    {
        // Use the vertex constraints to patch. Is of course bit dodgy since
        // face goes on patch if all its vertices are on a constraint.
        // Note: very inefficient since goes through all faces (including
        // internal ones) twice. Should do a construct without patches
        // and then repatchify.

        Info<< "Using " << dofVertIndices.size()
            << " DOF sets to detect boundary faces."<< endl;

        // Renumber vertex numbers on contraints
        forAll(dofVertIndices, patchI)
        {
            inplaceRenumber(unvToFoam, dofVertIndices[patchI]);
        }


        // Build labelHashSet of points per dofGroup/patch
        List<labelHashSet> dofGroups(dofVertIndices.size());

        forAll(dofVertIndices, patchI)
        {
            const labelList& foamVerts = dofVertIndices[patchI];

            forAll(foamVerts, i)
            {
                dofGroups[patchI].insert(foamVerts[i]);
            }
        }

        List<DynamicList<face> > dynPatchFaces(dofVertIndices.size());

        forAll(cellVerts, cellI)
        {
            const cellShape& shape = cellVerts[cellI];

            const faceList shapeFaces(shape.faces());

            forAll(shapeFaces, i)
            {
                label patchI = findPatch(dofGroups, shapeFaces[i]);

                if (patchI != -1)
                {
                    dynPatchFaces[patchI].append(shapeFaces[i]);
                }
            }
        }

        // Transfer
        patchFaceVerts.setSize(dynPatchFaces.size());

        forAll(dynPatchFaces, patchI)
        {
            patchFaceVerts[patchI].transfer(dynPatchFaces[patchI]);
        }
    }
    else
    {
        // Use the boundary faces.

        // Construct the patch faces by sorting the boundaryFaces according to
        // patch.
        patchFaceVerts.setSize(patchFaceIndices.size());

        Info<< "Sorting boundary faces according to group (patch)" << endl;

        // Construct map from boundaryFaceIndices
        Map<label> boundaryFaceToIndex(boundaryFaceIndices.size());

        forAll(boundaryFaceIndices, i)
        {
            boundaryFaceToIndex.insert(boundaryFaceIndices[i], i);
        }

        forAll(patchFaceVerts, patchI)
        {
            faceList& patchFaces = patchFaceVerts[patchI];
            const labelList& faceIndices = patchFaceIndices[patchI];

            patchFaces.setSize(faceIndices.size());

            forAll(patchFaces, i)
            {
                label bFaceI = boundaryFaceToIndex[faceIndices[i]];

                patchFaces[i] = boundaryFaces[bFaceI];
            }
        }
    }

    pointField polyPoints;
    polyPoints.transfer(points);

    // Length scaling factor
    polyPoints /= lengthScale;


    // For debugging: dump boundary faces as triSurface
    if (args.optionFound("dump"))
    {
        DynamicList<labelledTri> triangles(boundaryFaces.size());

        forAll(boundaryFaces, i)
        {
            const face& f = boundaryFaces[i];

            faceList triFaces(f.nTriangles(polyPoints));
            label nTri = 0;
            f.triangles(polyPoints, nTri, triFaces);

            forAll(triFaces, triFaceI)
            {
                const face& f = triFaces[triFaceI];
                triangles.append(labelledTri(f[0], f[1], f[2], 0));
            }
        }

        // Create globally numbered tri surface
        triSurface rawSurface(triangles.shrink(), polyPoints);

        // Create locally numbered tri surface
        triSurface surface
        (
            rawSurface.localFaces(),
            rawSurface.localPoints()
        );

        Info<< "Writing boundary faces to STL file boundaryFaces.stl"
            << nl << endl;

        surface.write(runTime.path()/"boundaryFaces.stl");
    }


    Info<< "Constructing mesh with non-default patches of size:" << nl;
    forAll(patchNames, patchI)
    {
        Info<< "    " << patchNames[patchI] << '\t'
            << patchFaceVerts[patchI].size() << nl;
    }
    Info<< endl;



    // Construct mesh
    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        xferMove(polyPoints),
        cellVerts,
        patchFaceVerts,             // boundaryFaces,
        patchNames,                 // boundaryPatchNames,
        wordList(patchNames.size(), polyPatch::typeName), // boundaryPatchTypes,
        "defaultFaces",             // defaultFacesName
        polyPatch::typeName,        // defaultFacesType,
        wordList(0)                 // boundaryPatchPhysicalTypes
    );

    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
