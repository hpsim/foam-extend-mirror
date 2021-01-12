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

Class
    Foam::volumeType

Description
    An enumeration wrapper for classification of a location as being
    inside/outside of a volume.

SourceFiles
    volumeType.C

\*---------------------------------------------------------------------------*/

#ifndef volumeType_H
#define volumeType_H

#include "contiguous.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class volumeType;
Istream& operator>>(Istream& is, volumeType&);
Ostream& operator<<(Ostream& os, const volumeType& C);


/*---------------------------------------------------------------------------*\
                         Class volumeType Declaration
\*---------------------------------------------------------------------------*/

class volumeType
{
public:

    //- Volume classification types
    enum type
    {
        UNKNOWN = 0,    //!< Unknown state
        INSIDE  = 1,    //!< A location inside the volume
        OUTSIDE = 2,    //!< A location outside the volume
        MIXED   = 3     //!< A location that is partly inside and outside
    };

private:

    // Private data

        //- Volume type
        type t_;


public:

    // Static data

        static const NamedEnum<volumeType, 4> names;


    // Constructors

        //- Construct null as UNKNOWN state
        volumeType()
        :
            t_(UNKNOWN)
        {}

        //- Construct from enumeration
        volumeType(type t)
        :
            t_(t)
        {}


    // Member Functions

        //- Return the enumeration
        operator type() const
        {
            return t_;
        }


    // IOstream operators

        friend Istream& operator>>(Istream& is, volumeType& vt);
        friend Ostream& operator<<(Ostream& os, const volumeType& vt);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for volumeType
template<>
inline bool contiguous<volumeType>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
