/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
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
    partTet

Description
    A partTet is a primitive needed during mesh smoothing

SourceFiles
    partTetI.H

\*---------------------------------------------------------------------------*/

#ifndef partTet_H
#define partTet_H

#include "pointField.H"
#include "edge.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                           class partTet Declaration
\*---------------------------------------------------------------------------*/

class partTet
{
protected:

    // Protected data
        label data_[4];

public:

    // Constructors
        //- Null construct
        inline partTet();

        //- Construct from point labels
        inline partTet
        (
            const label a,
            const label b,
            const label c,
            const label d
        );

    // Destructor
        ~partTet();


    // Member Functions

        // Access

            //- Return vertices
            inline label a() const;

            inline label b() const;

            inline label c() const;

            inline label d() const;

            inline label size() const;

        // Searches
            //- find position of the node in the partTet
            inline label whichPosition(const label) const;

        // Properties

            //- Return face normal
            template<class PointField>
            inline vector Sa(const PointField&) const;

            template<class PointField>
            inline vector Sb(const PointField&) const;

            template<class PointField>
            inline vector Sc(const PointField&) const;

            template<class PointField>
            inline vector Sd(const PointField&) const;


            //- Return volume
            template<class PointField>
            inline scalar mag(const PointField&) const;

            //- Return circum-centre
            template<class PointField>
            inline point crcmCentre(const PointField&) const;

            template<class PointField>
            inline scalar crcmRadius(const PointField&) const;

            //- Return centroid of the tetrahedron
            template<class PointField>
            inline point centroid(const PointField&) const;

            //- Return edges
            inline FixedList<edge, 6> edges() const;

            //- Return edge
            inline edge getEdge(const label) const;

    // Member operators

            inline label operator[](const label) const;

            inline void operator=(const partTet&);

            inline bool operator==(const partTet&) const;
            inline bool operator!=(const partTet&) const;

    // IOstream operators

            friend Ostream& operator<<(Ostream&, const partTet&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "partTetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
