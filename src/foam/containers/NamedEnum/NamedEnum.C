/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "stringList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Enum, unsigned int nEnum>
Foam::NamedEnum<Enum, nEnum>::NamedEnum()
:
    HashTable<unsigned int>(2*nEnum)
{
    for (unsigned int enumI = 0; enumI < nEnum; ++enumI)
    {
        if (!names[enumI] || names[enumI][0] == '\0')
        {
            stringList goodNames(enumI);

            for (unsigned int i = 0; i < enumI; ++i)
            {
                goodNames[i] = names[i];
            }

            FatalErrorInFunction
                << "Illegal enumeration name at position " << enumI << endl
                << "after entries " << goodNames << ".\n"
                << "Possibly your NamedEnum<Enum, nEnum>::names array"
                << " is not of size " << nEnum << endl
                << abort(FatalError);
        }
        insert(names[enumI], enumI);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Enum, unsigned int nEnum>
Enum Foam::NamedEnum<Enum, nEnum>::read(Istream& is) const
{
    const word name(is);

    HashTable<unsigned int>::const_iterator iter = find(name);

    if (iter == HashTable<unsigned int>::end())
    {
        FatalIOErrorInFunction(is)
            << name << " is not in enumeration: "
            << sortedToc() << exit(FatalIOError);
    }

    return Enum(iter());
}


template<class Enum, unsigned int nEnum>
void Foam::NamedEnum<Enum, nEnum>::write(const Enum e, Ostream& os) const
{
    os  << operator[](e);
}


template<class Enum, unsigned int nEnum>
const char* Foam::NamedEnum<Enum, nEnum>::operator[](const Enum e) const
{
    if (e < nEnum)
    {
        return names[e];
    }
    else
    {
        FatalErrorInFunction
            << "names array index " << e << " out of range 0-"
            << nEnum - 1
            << exit(FatalError);

        return names[0];
    }
}


template<class Enum, unsigned int nEnum>
const char* Foam::NamedEnum<Enum, nEnum>::operator[](const int e) const
{
    unsigned int ue = unsigned(e);

    return this->operator[](ue);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Enum, unsigned int nEnum>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const NamedEnum<Enum, nEnum>& n
)
{
    for (unsigned int e = 0; e < nEnum; e++)
    {
        os << e << " " << n.names[e] << nl;
    }

    return os;
}


// ************************************************************************* //
