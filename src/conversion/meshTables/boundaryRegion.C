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

#include "boundaryRegion.H"
#include "IOMap.H"
#include "OFstream.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryRegion::boundaryRegion()
:
    Map<dictionary>()
{}


Foam::boundaryRegion::boundaryRegion
(
    const objectRegistry& registry,
    const word& name,
    const fileName& instance
)
:
    Map<dictionary>()
{
    readDict(registry, name, instance);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryRegion::~boundaryRegion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::boundaryRegion::append(const dictionary& dict)
{
    label maxId = -1;
    forAllConstIter(Map<dictionary>, *this, iter)
    {
        if (maxId < iter.key())
        {
            maxId = iter.key();
        }
    }

    insert(++maxId, dict);
    return maxId;
}


Foam::Map<Foam::word> Foam::boundaryRegion::names() const
{
    Map<word> lookup;

    forAllConstIter(Map<dictionary>, *this, iter)
    {
        lookup.insert
        (
            iter.key(),
            iter().lookupOrDefault<word>
            (
                "Label",
                "boundaryRegion_" + Foam::name(iter.key())
            )
        );
    }

    return lookup;
}


Foam::Map<Foam::word> Foam::boundaryRegion::names
(
    const List<wordRe>& patterns
) const
{
    Map<word> lookup;

    forAllConstIter(Map<dictionary>, *this, iter)
    {
        word lookupName = iter().lookupOrDefault<word>
        (
            "Label",
            "boundaryRegion_" + Foam::name(iter.key())
        );

        if (findStrings(patterns, lookupName))
        {
            lookup.insert(iter.key(), lookupName);
        }
    }

    return lookup;
}


Foam::Map<Foam::word> Foam::boundaryRegion::boundaryTypes() const
{
    Map<word> lookup;

    forAllConstIter(Map<dictionary>, *this, iter)
    {
        lookup.insert
        (
            iter.key(),
            iter().lookupOrDefault<word>("BoundaryType", "patch")
        );
    }

    return lookup;
}


Foam::label Foam::boundaryRegion::findIndex(const word& name) const
{
    if (name.empty())
    {
        return -1;
    }

    forAllConstIter(Map<dictionary>, *this, iter)
    {
        if (iter().lookupOrDefault<word>("Label", word::null) == name)
        {
            return iter.key();
        }
    }

    return -1;
}


Foam::word Foam::boundaryRegion::boundaryType(const word& name) const
{
    word bndType("patch");

    label id = this->findIndex(name);
    if (id >= 0)
    {
        operator[](id).readIfPresent<word>("BoundaryType", bndType);
    }

    return bndType;
}


void Foam::boundaryRegion::readDict
(
    const objectRegistry& registry,
    const word& name,
    const fileName& instance
)
{
    clear();

    // read constant/dictName
    IOMap<dictionary> ioObj
    (
        IOobject
        (
            name,
            instance,
            registry,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    if (ioObj.headerOk())
    {
        *this = ioObj;
    }
    else
    {
        Info<< "no constant/boundaryRegion information available" << endl;
    }
}


void Foam::boundaryRegion::writeDict
(
    const objectRegistry& registry,
    const word& name,
    const fileName& instance
) const
{
    // write constant/dictName
    IOMap<dictionary> ioObj
    (
        IOobject
        (
            name,
            instance,
            registry,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    ioObj.note() =
        "persistent data for thirdParty mesh <-> foam translation";

    Info<< "Writing " << ioObj.name() << " to " << ioObj.objectPath() << endl;

    OFstream os(ioObj.objectPath());
    ioObj.writeHeader(os);
    os << *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::boundaryRegion::operator=(const boundaryRegion& rhs)
{
    Map<dictionary>::operator=(rhs);
}


void Foam::boundaryRegion::operator=(const Map<dictionary>& rhs)
{
    Map<dictionary>::operator=(rhs);
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

void Foam::boundaryRegion::rename(const dictionary& mapDict)
{
    if (mapDict.empty())
    {
        return;
    }

    // Use 1st pass to collect all the regions to be changed
    // and 2nd pass to relabel regions.
    // This avoid re-matching any renamed regions

    Map<word> mapping;
    forAllConstIter(dictionary, mapDict, iter)
    {
        word oldName(iter().stream());

        label id = this->findIndex(oldName);
        if (id >= 0)
        {
            mapping.insert(id, iter().keyword());
        }
    }

    forAllConstIter(Map<word>, mapping, iter)
    {
        dictionary& dict = operator[](iter.key());

        Info<< "rename patch: " << iter()
            << " <- " << word(dict.lookup("Label")) << nl;

        dict.set("Label", iter());
    }
}


// ************************************************************************* //
