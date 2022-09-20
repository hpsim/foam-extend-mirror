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

#include "volumeType.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::volumeType,
        4
    >::names[] =
    {
        "unknown",
        "mixed",
        "inside",
        "outside"
    };
}

const Foam::NamedEnum<Foam::volumeType, 4> Foam::volumeType::names;


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, volumeType& vt)
{
    // Read beginning of volumeType
    is.readBegin("volumeType");

    int type;
    is  >> type;

    vt.t_ = static_cast<volumeType::type>(type);

    // Read end of volumeType
    is.readEnd("volumeType");

    // Check state of Istream
    is.check("operator>>(Istream&, volumeType&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const volumeType& vt)
{
    os  << static_cast<int>(vt.t_);

    return os;
}


// ************************************************************************* //
