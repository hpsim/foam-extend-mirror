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
    Foam::jplotSetWriter

SourceFiles
    jplotSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef jplotSetWriter_H
#define jplotSetWriter_H

#include "writer.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class jplotSetWriter Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class jplotSetWriter
:
    public writer<Type>
{

    // Private Member Functions

        //- Write header
        Ostream& writeHeader(Ostream&) const;

public:

    //- Runtime type information
    TypeName("jplot");


    // Constructors

        //- Construct null
        jplotSetWriter();


    //- Destructor
    virtual ~jplotSetWriter();


    // Member Functions

        virtual fileName getFileName
        (
            const coordSet&,
            const wordList&
        ) const;

        virtual void write
        (
            const coordSet&,
            const wordList&,
            const List<const Field<Type>*>&,
            Ostream&
        ) const;

        virtual void write
        (
            const bool writeTracks,
            const PtrList<coordSet>&,
            const wordList& valueSetNames,
            const List<List<Field<Type> > >&,
            Ostream&
        ) const
        {
            notImplemented
            (
                "jplotSetWriter<Type>::write\n"
                "(\n"
                "    const bool,\n"
                "    const PtrList<coordSet>&,\n"
                "    const wordList&,\n"
                "    const List<List<Field<Type> > >&,\n"
                "    Ostream&\n"
                ") const"
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "jplotSetWriter.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
