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
    Output for ensightPart

\*---------------------------------------------------------------------------*/

#include "ensightPart.H"
#include "dictionary.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightPart::writeHeader
(
    ensightFile& os,
    bool withDescription
) const
{
    os.write("part");
    os.newline();

    os.write(number() + 1);   // Ensight starts with 1
    os.newline();

    if (withDescription)
    {
        os.write(name());
        os.newline();
    }
}


void Foam::ensightPart::writeFieldList
(
    ensightFile& os,
    const scalarList& field,
    const labelList& idList
) const
{
    forAll(idList, i)
    {
#       ifndef Intel
        if (idList[i] >= field.size() || std::isnan(field[idList[i]]))
#       else
        if (idList[i] >= field.size())
#       endif
        {
            os.writeUndef();
        }
        else
        {
            os.write(field[idList[i]]);
        }

        os.newline();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ensightPart::writeSummary(Ostream& os) const
{
    os  << indent << type() << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;

    // Ensight starts with 1
    os.writeKeyword("id") << (number() + 1) << token::END_STATEMENT << nl;
    os.writeKeyword("name") << name() << token::END_STATEMENT << nl;
    os.writeKeyword("offset") << offset() << token::END_STATEMENT << nl;
    os.writeKeyword("size") << size() << token::END_STATEMENT << nl;

    os   << decrIndent << indent << token::END_BLOCK << nl << endl;

    return true;
}


bool Foam::ensightPart::writeData(Ostream& os) const
{
    os  << indent << type() << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;

    os.writeKeyword("id") << number() << token::END_STATEMENT << nl;
    os.writeKeyword("name") << name() << token::END_STATEMENT << nl;
    os.writeKeyword("offset") << offset() << token::END_STATEMENT << nl;

    forAll(elementTypes(), typeI)
    {
        word key(elementTypes()[typeI]);
        if (elemLists_[typeI].size())
        {
            elemLists_[typeI].writeEntry(key, os);
        }
    }

    os   << decrIndent << indent << token::END_BLOCK << nl << endl;

    return true;
}


void Foam::ensightPart::writeGeometry(ensightGeoFile& os) const
{
    if (size() && meshPtr_)
    {
        const polyMesh& mesh = *meshPtr_;
        const pointField& meshPoints = mesh.points();

        localPoints ptList = calcLocalPoints();
        labelList& pointMap = ptList.list;

        writeHeader(os, true);

        // write points
        os.writeKeyword("coordinates");
        os.write(ptList.nPoints);
        os.newline();

        for (direction cmpt=0; cmpt < vector::nComponents; cmpt++)
        {
            forAll(pointMap, ptI)
            {
                if (pointMap[ptI] > -1)
                {
                    os.write( meshPoints[ptI].component(cmpt) );
                    os.newline();
                }
            }
        }

        // write parts
        forAll(elementTypes(), elemI)
        {
            if (elemLists_[elemI].size())
            {
                writeConnectivity
                (
                    os,
                    elementTypes()[elemI],
                    elemLists_[elemI],
                    pointMap
                );
            }
        }
    }
}


void Foam::ensightPart::writeScalarField
(
    ensightFile& os,
    const scalarList& field
) const
{
    if (size() && field.size() && (os.allowUndef() || isFieldDefined(field)))
    {
        writeHeader(os);

        forAll(elementTypes(), elemI)
        {
            const labelList& idList = elemLists_[elemI];

            if (idList.size())
            {
                os.writeKeyword( elementTypes()[elemI] );
                writeFieldList(os, field, idList);
            }
        }
    }
}


void Foam::ensightPart::writeVectorField
(
    ensightFile& os,
    const scalarList& field0,
    const scalarList& field1,
    const scalarList& field2
) const
{
    if (size() && field0.size() && (os.allowUndef() || isFieldDefined(field0)))
    {
        writeHeader(os);

        forAll(elementTypes(), elemI)
        {
            const labelList& idList = elemLists_[elemI];

            if (idList.size())
            {
                os.writeKeyword( elementTypes()[elemI] );
                writeFieldList(os, field0, idList);
                writeFieldList(os, field1, idList);
                writeFieldList(os, field2, idList);
            }
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ensightPart& part
)
{
    part.writeData(os);
    return os;
}


Foam::ensightGeoFile& Foam::operator<<
(
    ensightGeoFile& os,
    const ensightPart& part
)
{
    part.writeGeometry(os);
    return os;
}


// ************************************************************************* //
