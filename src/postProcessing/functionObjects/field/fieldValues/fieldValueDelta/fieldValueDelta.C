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

#include "fieldValueDelta.H"
#include "ListOps.H"
#include "foamTime.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fieldValues
    {
        defineTypeNameAndDebug(fieldValueDelta, 0);
    }

    template<>
    const char*
    NamedEnum<fieldValues::fieldValueDelta::operationType, 5>::names[] =
    {
        "add",
        "subtract",
        "min",
        "max",
        "average"
    };

    const NamedEnum<fieldValues::fieldValueDelta::operationType, 5>
        fieldValues::fieldValueDelta::operationTypeNames_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldValues::fieldValueDelta::fieldValueDelta
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, typeName),
    name_(name),
    obr_(obr),
    loadFromFiles_(loadFromFiles),
    log_(true),
    operation_(opSubtract),
    source1Ptr_(nullptr),
    source2Ptr_(nullptr)
{
    read(dict);
}


void Foam::fieldValues::fieldValueDelta::writeFileHeader(const label i)
{
    const wordList& fields1 = source1Ptr_->fields();
    const wordList& fields2 = source2Ptr_->fields();

    DynamicList<word> commonFields(fields1.size());
    forAll(fields1, i)
    {
        label index = findIndex(fields2, fields1[i]);
        if (index != -1)
        {
            commonFields.append(fields1[i]);
        }
    }

    Ostream& os = file();

    writeHeaderValue(os, "Source1", source1Ptr_->name());
    writeHeaderValue(os, "Source2", source2Ptr_->name());
    writeHeaderValue(os, "Operation", operationTypeNames_[operation_]);
    writeCommented(os, "Time");

    forAll(commonFields, i)
    {
        os  << tab << commonFields[i];
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldValues::fieldValueDelta::~fieldValueDelta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldValues::fieldValueDelta::read(const dictionary& dict)
{
    log_ = dict.lookupOrDefault<Switch>("log", true);
    source1Ptr_.reset
    (
        fieldValue::New
        (
            name_ + ".source1",
            obr_,
            dict.subDict("source1"),
            loadFromFiles_,
            false
        ).ptr()
    );
    source2Ptr_.reset
    (
        fieldValue::New
        (
            name_ + ".source2",
            obr_,
            dict.subDict("source2"),
            loadFromFiles_,
            false
        ).ptr()
    );

    operation_ = operationTypeNames_.read(dict.lookup("operation"));
}


void Foam::fieldValues::fieldValueDelta::write()
{
    functionObjectFile::write();

    source1Ptr_->write();
    source2Ptr_->write();

    if (Pstream::master())
    {
        file()<< obr_.time().value();
    }

    if (log_) Info<< type() << " " << name_ << " output:" << endl;

    bool found = false;
    processFields<scalar>(found);
    processFields<vector>(found);
    processFields<sphericalTensor>(found);
    processFields<symmTensor>(found);
    processFields<tensor>(found);

    if (Pstream::master())
    {
        file()<< endl;
    }

    if (log_)
    {
        if (!found)
        {
            Info<< "    none" << endl;
        }
        else
        {
            Info<< endl;
        }
    }
}


void Foam::fieldValues::fieldValueDelta::execute()
{
    // Do nothing
}


void Foam::fieldValues::fieldValueDelta::end()
{
    // Do nothing
}


void Foam::fieldValues::fieldValueDelta::timeSet()
{
    // Do nothing
}


void Foam::fieldValues::fieldValueDelta::updateMesh(const mapPolyMesh&)
{
    // Do nothing
}


void Foam::fieldValues::fieldValueDelta::movePoints(const pointField&)
{
    // Do nothing
}


// ************************************************************************* //
