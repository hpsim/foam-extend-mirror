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

Class
    Foam::Table

Description
    Templated table container data entry. Items are stored in a list of
    Tuple2's. First column is always stored as scalar entries. Data is read
    in the form, e.g. for an entry <entryName> that is (scalar, vector):

    @verbatim
        <entryName>   table
        (
            0.0 (1 2 3)
            1.0 (4 5 6)
        );
    @endverbatim

SourceFiles
    Table.C

\*---------------------------------------------------------------------------*/

#ifndef Table_H
#define Table_H

#include "DataEntry.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class Table;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const Table<Type>&
);

/*---------------------------------------------------------------------------*\
                           Class Table Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Table
:
    public DataEntry<Type>
{
    // Private data

        //- Table data
        List<Tuple2<scalar, Type> > table_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const Table<Type>&);


public:

    //- Runtime type information
    TypeName("table");


    // Constructors

        //- Construct from entry name and Istream
        Table(const word& entryName, Istream& is);

        //- Copy constructor
        Table(const Table<Type>& tbl);


    //- Destructor
    virtual ~Table();


    // Member Functions

        //- Return Table value
        Type value(const scalar x) const;

        //- Integrate between two (scalar) values
        Type integrate(const scalar x1, const scalar x2) const;


    //- Ostream Operator
    friend Ostream& operator<< <Type>
    (
        Ostream&,
        const Constant<Type>&
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Table.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
