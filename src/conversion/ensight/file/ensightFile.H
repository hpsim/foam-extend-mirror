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
    Foam::ensightFile

Description
    Ensight output with specialized write() for strings, integers and floats.
    Correctly handles binary write as well.

\*---------------------------------------------------------------------------*/

#ifndef ensightFile_H
#define ensightFile_H

#include "OFstream.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ensightFile Declaration
\*---------------------------------------------------------------------------*/

class ensightFile
:
    public OFstream
{
    // Private data

        //- allow undef in results
        static bool allowUndef_;

        //- value to represent undef in results
        static scalar undefValue_;


    // Private member functions

        //- Disallow default bitwise assignment
        void operator=(const ensightFile&);

        //- Disallow default copy constructor
        ensightFile(const ensightFile&);


public:

    // Constructors

        //- Construct from pathname
        ensightFile
        (
            const fileName& pathname,
            IOstream::streamFormat format=IOstream::BINARY
        );


    // Destructor
    ~ensightFile();


    // Access

        //- Return setting for whether 'undef' values are allowed in results
        static bool allowUndef();

        //- '*' mask appropriate for subDir
        static string mask();

        //- consistent zero-padded numbers for subdirectories
        static string subDir(const label);


    // Edit

        static bool allowUndef(bool);

        //- Assign the value to represent undef in the results
        //  Returns the previous value
        //  NB: do not use values larger than floatScalarVGREAT
        static scalar undefValue(const scalar&);


    // Output

        //- binary write
        virtual Ostream& write(const char* buf, std::streamsize count);

        //- write element keyword with trailing newline, optionally with undef
        virtual Ostream& writeKeyword(const string& key);

        //- write "C Binary" for binary files (eg, geometry/measured)
        Ostream& writeBinaryHeader();

        //- write undef value
        Ostream& writeUndef();

        //- write string as "%80s" or as binary
        Ostream& write(const string& value);

        //- write integer as "%10d" or as binary
        Ostream& write(const label& value);

        //- write integer with specified width or as binary
        Ostream& write(const label& value, const label fieldWidth);

        //- write float as "%12.5e" or as binary
        Ostream& write(const scalar& value);

        //- Add carriage return to ascii stream
        void newline();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
