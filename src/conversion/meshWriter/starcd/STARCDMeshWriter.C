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

#include "STARCDMeshWriter.H"

#include "foamTime.H"
#include "SortableList.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::meshWriters::STARCD::defaultBoundaryName =
    "Default_Boundary_Region";

const Foam::label Foam::meshWriters::STARCD::foamToStarFaceAddr[4][6] =
{
    { 4, 5, 2, 3, 0, 1 },     // 11 = pro-STAR hex
    { 0, 1, 4, 5, 2, -1 },    // 12 = pro-STAR prism
    { 5, 4, 2, 0, -1, -1 },   // 13 = pro-STAR tetra
    { 0, 4, 3, 5, 2, -1 }     // 14 = pro-STAR pyramid
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::meshWriters::STARCD::findDefaultBoundary() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label id = -1;

    // find Default_Boundary_Region if it exists
    forAll(patches, patchI)
    {
        if (defaultBoundaryName == patches[patchI].name())
        {
            id = patchI;
            break;
        }
    }
    return id;
}


void Foam::meshWriters::STARCD::getCellTable()
{
    // read constant/polyMesh/propertyName
    IOList<label> ioList
    (
        IOobject
        (
            "cellTableId",
            "constant",
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    bool useCellZones = false;
    cellTableId_.setSize(mesh_.nCells(), -1);

    // get information from constant/polyMesh/cellTableId if possible
    if (ioList.headerOk())
    {
        if (ioList.size() == mesh_.nCells())
        {
            cellTableId_.transfer(ioList);

            if (cellTable_.empty())
            {
                Info<< "no cellTable information available" << endl;
            }
        }
        else
        {
            WarningIn("STARCD::getCellTable()")
                << ioList.objectPath() << " has incorrect number of cells "
                << " - use cellZone information"
                << endl;

            ioList.clear();
            useCellZones = true;
        }
    }
    else
    {
        useCellZones = true;
    }


    if (useCellZones)
    {
        if (cellTable_.empty())
        {
            Info<< "created cellTable from cellZones" << endl;
            cellTable_ = mesh_;
        }

        // track if there are unzoned cells
        label nUnzoned = mesh_.nCells();

        // get the cellZone <-> cellTable correspondence
        Info<< "matching cellZones to cellTable" << endl;

        forAll (mesh_.cellZones(), zoneI)
        {
            const cellZone& cZone = mesh_.cellZones()[zoneI];
            if (cZone.size())
            {
                nUnzoned -= cZone.size();

                label tableId = cellTable_.findIndex(cZone.name());
                if (tableId < 0)
                {
                    dictionary dict;

                    dict.add("Label", cZone.name());
                    dict.add("MaterialType", "fluid");
                    tableId = cellTable_.append(dict);
                }

                forAll (cZone, i)
                {
                    cellTableId_[cZone[i]] = tableId;
                }
            }
        }

        if (nUnzoned)
        {
            dictionary dict;

            dict.add("Label", "__unZonedCells__");
            dict.add("MaterialType", "fluid");
            label tableId = cellTable_.append(dict);

            forAll (cellTableId_, i)
            {
                if (cellTableId_[i] < 0)
                {
                    cellTableId_[i] = tableId;
                }
            }
        }
    }
}


void Foam::meshWriters::STARCD::writeHeader(Ostream& os, const char* filetype)
{
    os  << "PROSTAR_" << filetype << nl
        << 4000
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << endl;
}


void Foam::meshWriters::STARCD::writePoints(const fileName& prefix) const
{
    OFstream os(prefix + ".vrt");
    writeHeader(os, "VERTEX");

    // Set the precision of the points data to 10
    os.precision(10);

    // force decimal point for Fortran input
    os.setf(std::ios::showpoint);

    const pointField& points = mesh_.points();

    Info<< "Writing " << os.name() << " : "
        << points.size() << " points" << endl;

    forAll(points, ptI)
    {
        // convert [m] -> [mm]
        os
            << ptI + 1 << " "
            << scaleFactor_ * points[ptI].x() << " "
            << scaleFactor_ * points[ptI].y() << " "
            << scaleFactor_ * points[ptI].z() << nl;
    }
    os.flush();

}


void Foam::meshWriters::STARCD::writeCells(const fileName& prefix) const
{
    OFstream os(prefix + ".cel");
    writeHeader(os, "CELL");

    // this is what we seem to need
    // map foam cellModeller index -> star shape
    Map<label> shapeLookupIndex;
    shapeLookupIndex.insert(hexModel->index(), 11);
    shapeLookupIndex.insert(prismModel->index(), 12);
    shapeLookupIndex.insert(tetModel->index(), 13);
    shapeLookupIndex.insert(pyrModel->index(), 14);

    const cellShapeList& shapes = mesh_.cellShapes();
    const cellList& cells  = mesh_.cells();
    const faceList& faces  = mesh_.faces();
    const labelList& owner = mesh_.faceOwner();

    Info<< "Writing " << os.name() << " : "
        << cells.size() << " cells" << endl;

    forAll(cells, cellId)
    {
        label tableId = cellTableId_[cellId];
        label materialType  = 1;        // 1(fluid)
        if (cellTable_.found(tableId))
        {
            const dictionary& dict = cellTable_[tableId];
            if (dict.found("MaterialType"))
            {
                word matType;
                dict.lookup("MaterialType") >> matType;
                if (matType == "solid")
                {
                    materialType = 2;
                }

            }
        }

        const cellShape& shape = shapes[cellId];
        label mapIndex = shape.model().index();

        // a registered primitive type
        if (shapeLookupIndex.found(mapIndex))
        {
            label shapeId = shapeLookupIndex[mapIndex];
            const labelList& vrtList = shapes[cellId];

            os  << cellId + 1
                << " " << shapeId
                << " " << vrtList.size()
                << " " << tableId
                << " " << materialType;

            // primitives have <= 8 vertices, but prevent overrun anyhow
            // indent following lines for ease of reading
            label count = 0;
            forAll(vrtList, i)
            {
                if ((count % 8) == 0)
                {
                    os  << nl
                        << "  " << cellId + 1;
                }
                os << " " << vrtList[i] + 1;
                count++;
            }
            os << endl;

        }
        else
        {
            label shapeId = 255;        // treat as general polyhedral
            const labelList& cFaces  = cells[cellId];

            // create (beg,end) indices
            labelList indices(cFaces.size() + 1);
            indices[0] = indices.size();

            label count = indices.size();
            // determine the total number of vertices
            forAll(cFaces, faceI)
            {
                count += faces[cFaces[faceI]].size();
                indices[faceI+1] = count;
            }

            os  << cellId + 1
                << " " << shapeId
                << " " << count
                << " " << tableId
                << " " << materialType;

            // write indices - max 8 per line
            // indent following lines for ease of reading
            count = 0;
            forAll(indices, i)
            {
                if ((count % 8) == 0)
                {
                    os  << nl
                        << "  " << cellId + 1;
                }
                os << " " << indices[i];
                count++;
            }

            // write faces - max 8 per line
            forAll(cFaces, faceI)
            {
                label meshFace = cFaces[faceI];
                face f;

                if (owner[meshFace] == cellId)
                {
                    f = faces[meshFace];
                }
                else
                {
                    f = faces[meshFace].reverseFace();
                }

                forAll(f, i)
                {
                    if ((count % 8) == 0)
                    {
                        os  << nl
                            << "  " << cellId + 1;
                    }

                    os << " " << f[i] + 1;
                    count++;
                }
            }

            os << endl;
        }
    }
}


void Foam::meshWriters::STARCD::writeBoundary(const fileName& prefix) const
{
    OFstream os(prefix + ".bnd");
    writeHeader(os, "BOUNDARY");

    const cellShapeList& shapes = mesh_.cellShapes();
    const cellList& cells  = mesh_.cells();
    const faceList& faces  = mesh_.faces();
    const labelList& owner = mesh_.faceOwner();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // this is what we seem to need
    // these MUST correspond to foamToStarFaceAddr
    //
    Map<label> faceLookupIndex;
    faceLookupIndex.insert(hexModel->index(), 0);
    faceLookupIndex.insert(prismModel->index(), 1);
    faceLookupIndex.insert(tetModel->index(), 2);
    faceLookupIndex.insert(pyrModel->index(), 3);

    Info<< "Writing " << os.name() << " : "
        << (mesh_.nFaces() - patches[0].start()) << " boundaries" << endl;


    label defaultId = findDefaultBoundary();

    //
    // write boundary faces - skip Default_Boundary_Region entirely
    //
    label boundId = 0;
    forAll(patches, patchI)
    {
        label regionId = patchI;
        if (regionId == defaultId)
        {
            continue;    // skip - already written
        }
        else if (defaultId == -1 || regionId < defaultId)
        {
            regionId++;
        }

        label patchStart = patches[patchI].start();
        label patchSize  = patches[patchI].size();
        word  bndType = boundaryRegion_.boundaryType(patches[patchI].name());

        for
        (
            label faceI = patchStart;
            faceI < (patchStart + patchSize);
            ++faceI
        )
        {
            label cellId = owner[faceI];
            const labelList& cFaces  = cells[cellId];
            const cellShape& shape = shapes[cellId];
            label cellFaceId = findIndex(cFaces, faceI);

            //      Info<< "cell " << cellId + 1 << " face " << faceI
            //          << " == " << faces[faceI]
            //          << " is index " << cellFaceId << " from " << cFaces;

            // Unfortunately, the order of faces returned by
            //   primitiveMesh::cells() is not necessarily the same
            //   as defined by primitiveMesh::cellShapes()
            // Thus, for registered primitive types, do the lookup ourselves.
            // Finally, the cellModel face number is re-mapped to the
            // STAR-CD local face number

            label mapIndex = shape.model().index();

            // a registered primitive type
            if (faceLookupIndex.found(mapIndex))
            {
                const faceList sFaces = shape.faces();
                forAll(sFaces, sFaceI)
                {
                    if (faces[faceI] == sFaces[sFaceI])
                    {
                        cellFaceId = sFaceI;
                        break;
                    }
                }

                mapIndex = faceLookupIndex[mapIndex];
                cellFaceId = foamToStarFaceAddr[mapIndex][cellFaceId];
            }
            // Info<< endl;

            boundId++;

            os
                << boundId
                << " " << cellId + 1
                << " " << cellFaceId + 1
                << " " << regionId
                << " " << 0
                << " " << bndType.c_str()
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshWriters::STARCD::STARCD
(
    const polyMesh& mesh,
    const scalar scaleFactor
)
:
    meshWriter(mesh, scaleFactor)
{
    boundaryRegion_.readDict(mesh_);
    cellTable_.readDict(mesh_);
    getCellTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshWriters::STARCD::~STARCD()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshWriters::STARCD::rmFiles(const fileName& baseName) const
{
    rm(baseName + ".vrt");
    rm(baseName + ".cel");
    rm(baseName + ".bnd");
    rm(baseName + ".inp");
}


bool Foam::meshWriters::STARCD::write(const fileName& meshName) const
{
    fileName baseName(meshName);

    if (baseName.empty())
    {
        baseName = meshWriter::defaultMeshName;

        if
        (
            mesh_.time().timeName() != "0"
         && mesh_.time().timeName() != "constant"
        )
        {
            baseName += "_" + mesh_.time().timeName();
        }
    }

    rmFiles(baseName);
    writePoints(baseName);
    writeCells(baseName);

    if (writeBoundary_)
    {
        writeBoundary(baseName);
    }

    return true;
}


bool Foam::meshWriters::STARCD::writeSurface
(
    const fileName& meshName,
    const bool& triangulate
) const
{
    fileName baseName(meshName);

    if (baseName.empty())
    {
        baseName = meshWriter::defaultSurfaceName;

        if
        (
            mesh_.time().timeName() != "0"
         && mesh_.time().timeName() != "constant"
        )
        {
            baseName += "_" + mesh_.time().timeName();
        }
    }

    rmFiles(baseName);

    OFstream celFile(baseName + ".cel");
    writeHeader(celFile, "CELL");

    Info << "Writing " << celFile.name() << endl;

    // mesh and patch info
    const pointField& points = mesh_.points();
    const labelList& owner = mesh_.faceOwner();
    const faceList& meshFaces = mesh_.faces();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label shapeId = 3;  // shell/baffle element
    label typeId  = 4;  // 4(shell)

    // remember which points need to be written
    labelHashSet pointHash;

    // write boundary faces as normal STAR-CD mesh
    if (triangulate)
    {
        // cell Id has no particular meaning - just increment
        // use the cellTable id from the patch Number
        label cellId = 0;

        forAll(patches, patchI)
        {
            label patchStart = patches[patchI].start();
            label patchSize  = patches[patchI].size();

            label ctableId = patchI + 1;

            for
            (
                label faceI = patchStart;
                faceI < (patchStart + patchSize);
                ++faceI
            )
            {
                const face& f = meshFaces[faceI];

                label nTri = f.nTriangles(points);
                faceList triFaces;

                // triangulate polygons, but not quads
                if (nTri <= 2)
                {
                    triFaces.setSize(1);
                    triFaces[0] = f;
                }
                else
                {
                    triFaces.setSize(nTri);
                    nTri = 0;
                    f.triangles(points, nTri, triFaces);
                }

                forAll(triFaces, faceI)
                {
                    const labelList& vrtList = triFaces[faceI];

                    celFile
                        << cellId + 1 << " "
                        << shapeId << " "
                        << vrtList.size() << " "
                        << ctableId << " "
                        << typeId;

                    // must be 3 (triangle) but could be quad
                    label count = 0;
                    forAll(vrtList, i)
                    {
                        if ((count % 8) == 0)
                        {
                            celFile
                                << nl
                                << "  " << cellId + 1;
                        }
                        // remember which points we'll need to write
                        pointHash.insert(vrtList[i]);
                        celFile << " " << vrtList[i] + 1;
                        count++;
                    }
                    celFile << endl;

                    cellId++;
                }
            }
        }
    }
    else
    {
        // cell Id is the foam face Id
        // use the cellTable id from the face owner
        // - allows separation of parts
        forAll(patches, patchI)
        {
            label patchStart = patches[patchI].start();
            label patchSize  = patches[patchI].size();

            for
            (
                label faceI = patchStart;
                faceI < (patchStart + patchSize);
                ++faceI
            )
            {
                const labelList& vrtList = meshFaces[faceI];
                label cellId = faceI;

                celFile
                    << cellId + 1 << " "
                    << shapeId << " "
                    << vrtList.size() << " "
                    << cellTableId_[owner[faceI]] << " "
                    << typeId;

                // likely <= 8 vertices, but prevent overrun anyhow
                label count = 0;
                forAll(vrtList, i)
                {
                    if ((count % 8) == 0)
                    {
                        celFile
                            << nl
                            << "  " << cellId + 1;
                    }
                    // remember which points we'll need to write
                    pointHash.insert(vrtList[i]);
                    celFile << " " << vrtList[i] + 1;
                    count++;
                }
                celFile << endl;
            }
        }
    }

    OFstream vrtFile(baseName + ".vrt");
    writeHeader(vrtFile, "VERTEX");

    vrtFile.precision(10);
    vrtFile.setf(std::ios::showpoint);  // force decimal point for Fortran

    Info << "Writing " << vrtFile.name() << endl;

    // build sorted table of contents
    SortableList<label> toc(pointHash.size());
    {
        label i = 0;
        forAllConstIter(labelHashSet, pointHash, iter)
        {
            toc[i++] = iter.key();
        }
    }
    toc.sort();
    toc.shrink();
    pointHash.clear();

    // write points in sorted order
    forAll(toc, i)
    {
        label vrtId = toc[i];
        vrtFile
            << vrtId + 1
            << " " << scaleFactor_ * points[vrtId].x()
            << " " << scaleFactor_ * points[vrtId].y()
            << " " << scaleFactor_ * points[vrtId].z()
            << endl;
    }

    return true;
}


// ************************************************************************* //
