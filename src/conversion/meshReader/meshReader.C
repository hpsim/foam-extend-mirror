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

#include "meshReader.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "faceSet.H"
#include "emptyPolyPatch.H"
#include "cellModeller.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::cellModel* Foam::meshReader::unknownModel = Foam::cellModeller::
lookup
(
    "unknown"
);

const Foam::cellModel* Foam::meshReader::tetModel = Foam::cellModeller::
lookup
(
    "tet"
);

const Foam::cellModel* Foam::meshReader::pyrModel = Foam::cellModeller::
lookup
(
    "pyr"
);

const Foam::cellModel* Foam::meshReader::prismModel = Foam::cellModeller::
lookup
(
    "prism"
);

const Foam::cellModel* Foam::meshReader::hexModel = Foam::cellModeller::
lookup
(
    "hex"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshReader::addCellZones(polyMesh& mesh) const
{
    cellTable_.addCellZones(mesh, cellTableId_);
    warnDuplicates("cellZones", mesh.cellZones().names());
}


void Foam::meshReader::addFaceZones(polyMesh& mesh) const
{
    label nZone = monitoringSets_.size();
    mesh.faceZones().setSize(nZone);

    if (!nZone)
    {
        return;
    }

    nZone = 0;
    for
    (
        HashTable<List<label>, word, string::hash>::const_iterator
        iter = monitoringSets_.begin();
        iter != monitoringSets_.end();
        ++iter
    )
    {
        Info<< "faceZone " << nZone
            << " (size: " << iter().size() << ") name: "
            << iter.key() << endl;

        mesh.faceZones().set
        (
            nZone,
            new faceZone
            (
                iter.key(),
                iter(),
                boolList(iter().size(), false),
                nZone,
                mesh.faceZones()
            )
        );

        nZone++;
    }
    mesh.faceZones().writeOpt() = IOobject::AUTO_WRITE;
    warnDuplicates("faceZones", mesh.faceZones().names());
}


Foam::autoPtr<Foam::polyMesh> Foam::meshReader::mesh
(
    const objectRegistry& registry
)
{
    readGeometry();

    Info<< "Creating a polyMesh" << endl;
    createPolyCells();

    Info<< "Number of internal faces: " << nInternalFaces_ << endl;

    createPolyBoundary();
    clearExtraStorage();

    autoPtr<polyMesh> mesh
    (
        new polyMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                "constant",
                registry
            ),
            xferMove(points_),
            xferMove(meshFaces_),
            xferMove(cellPolys_)
        )
    );

    // adding patches also checks the mesh
    mesh().addPatches(polyBoundaryPatches(mesh));

    warnDuplicates("boundaries", mesh().boundaryMesh().names());

    addCellZones(mesh());
    addFaceZones(mesh());

    return mesh;
}


void Foam::meshReader::writeMesh
(
    const polyMesh& mesh,
    IOstream::streamFormat fmt
) const
{
    mesh.removeFiles();

    Info<< "Writing polyMesh" << endl;
    mesh.writeObject
    (
        fmt,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED
    );
    writeAux(mesh);
}


void Foam::meshReader::clearExtraStorage()
{
    cellFaces_.clear();
    baffleFaces_.clear();
    boundaryIds_.clear();
    baffleIds_.clear();

    deleteDemandDrivenData(pointCellsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshReader::meshReader
(
    const fileName& fileOrPrefix,
    const scalar scaleFactor
)
    :
    pointCellsPtr_(nullptr),
    nInternalFaces_(0),
    patchStarts_(0),
    patchSizes_(0),
    interfaces_(0),
    baffleIds_(0),
    meshFaces_(0),
    cellPolys_(0),
    geometryFile_(fileOrPrefix),
    scaleFactor_(scaleFactor),
    points_(0),
    origCellId_(0),
    boundaryIds_(0),
    patchTypes_(0),
    patchNames_(0),
    patchPhysicalTypes_(0),
    cellFaces_(0),
    baffleFaces_(0),
    cellTableId_(0),
    cellTable_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshReader::~meshReader()
{
    deleteDemandDrivenData(pointCellsPtr_);
}


// ************************************************************************* //
