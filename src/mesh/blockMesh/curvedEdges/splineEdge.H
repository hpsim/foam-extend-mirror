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
    Foam::splineEdge

Description
    A curvedEdge interface for Catmull-Rom splines.

SourceFiles
    splineEdge.C

\*---------------------------------------------------------------------------*/

#ifndef splineEdge_H
#define splineEdge_H

#include "curvedEdge.H"
#include "CatmullRomSpline.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class splineEdge Declaration
\*---------------------------------------------------------------------------*/

class splineEdge
:
    public curvedEdge,
    public CatmullRomSpline
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        splineEdge(const splineEdge&);

        //- Disallow default bitwise assignment
        void operator=(const splineEdge&);


public:

    //- Runtime type information
    TypeName("spline");


    // Constructors

        //- Construct from components
        splineEdge
        (
            const pointField&,
            const label start,
            const label end,
            const pointField& internalPoints
        );

        //- Construct from Istream, setting pointsList
        splineEdge(const pointField&, Istream&);


    //- Destructor
    virtual ~splineEdge();


    // Member Functions

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        virtual point position(const scalar) const;

        //- Return the length of the spline curve (not implemented)
        virtual scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
