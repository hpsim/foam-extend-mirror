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
    Foam::lineDivide

Description
    Divides a line into segments using grading descriptors

SourceFiles
    lineDivide.C

\*---------------------------------------------------------------------------*/

#ifndef lineDivide_H
#define lineDivide_H

#include "pointField.H"
#include "scalarList.H"
#include "gradingDescriptors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class curvedEdge;

/*---------------------------------------------------------------------------*\
                         Class lineDivide Declaration
\*---------------------------------------------------------------------------*/

class lineDivide
{
    // Private data

        //- Point list
        pointField points_;

        //- Divisions
        scalarList lambdaDivisions_;


    // Private member functions

        //- Calculate the geometric expension factor from the expansion ratio
        inline scalar calcGexp(const scalar expRatio, const label nDiv)
        {
            return nDiv > 1 ? pow(expRatio, 1.0/(nDiv - 1)) : 0.0;
        }


public:

    // Constructors

        //- Construct from components with discretization and expansion ratio
        lineDivide
        (
            const curvedEdge&,
            const label ndiv,
            const gradingDescriptors& gd = gradingDescriptors()
        );


    // Member Functions

        //- Return the points
        inline const pointField& points() const
        {
            return points_;
        }

        //- Return the list of lambda values
        const scalarList& lambdaDivisions() const
        {
            return lambdaDivisions_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
