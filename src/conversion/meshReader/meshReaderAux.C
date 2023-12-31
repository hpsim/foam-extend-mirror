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

#include "meshReader.H"
#include "IOMap.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

void Foam::meshReader::warnDuplicates
(
    const word& context,
    const wordList& list
)
{
    HashTable<label> hashed(list.size());
    bool duplicates = false;

    forAll(list, listI)
    {
        // check duplicate name
        HashTable<label>::iterator iter = hashed.find(list[listI]);
        if (iter != hashed.end())
        {
            (*iter)++;
            duplicates = true;
        }
        else
        {
            hashed.insert(list[listI], 1);
        }
    }

    // warn about duplicate names
    if (duplicates)
    {
        Info << nl << "WARNING: " << context << " with identical names:";
        forAllConstIter(HashTable<label>, hashed, iter)
        {
            if (*iter > 1)
            {
                Info << "  " << iter.key();
            }
        }
        Info << nl << endl;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshReader::writeInterfaces(const objectRegistry& registry) const
{
    // write constant/polyMesh/interface
    IOList<labelList> ioObj
    (
        IOobject
        (
            "interfaces",
            "constant",
            polyMesh::meshSubDir,
            registry,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    ioObj.note() = "as yet unsupported interfaces (baffles)";

    Info<< "Writing " << ioObj.name() << " to " << ioObj.objectPath() << endl;

    OFstream os(ioObj.objectPath());
    ioObj.writeHeader(os);

    os << interfaces_;
    ioObj.writeEndDivider(os);
}


void Foam::meshReader::writeMeshLabelList
(
    const objectRegistry& registry,
    const word& propertyName,
    const labelList& list,
    IOstream::streamFormat fmt
) const
{
    // write constant/polyMesh/propertyName
    IOList<label> ioObj
    (
        IOobject
        (
            propertyName,
            "constant",
            polyMesh::meshSubDir,
            registry,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        list
    );


    ioObj.note() = "persistent data for star-cd <-> foam translation";
    Info<< "Writing " << ioObj.name() << " to " << ioObj.objectPath() << endl;

    // NOTE:
    // the cellTableId is an integer and almost always < 1000, thus ASCII
    // will be compacter than binary and makes external scripting easier
    //
    ioObj.writeObject
    (
        fmt,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshReader::writeAux(const objectRegistry& registry) const
{
    cellTable_.writeDict(registry);
    writeInterfaces(registry);

    // write origCellId as labelList
    writeMeshLabelList
    (
        registry,
        "origCellId",
        origCellId_,
        IOstream::BINARY
    );

    // write cellTableId as labelList
    // this is crucial for later conversion back to ccm/starcd
    writeMeshLabelList
    (
        registry,
        "cellTableId",
        cellTableId_,
        IOstream::ASCII
    );
}


// ************************************************************************* //
