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

#include "normal.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(normal, 0);
    addToRunTimeSelectionTable(pdf, normal, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::normal::normal(const dictionary& dict, Random& rndGen)
:
    pdf(dict, rndGen),
    pdfDict_(dict.subDict(typeName + "PDF")),
    minValue_(readScalar(pdfDict_.lookup("minValue"))),
    maxValue_(readScalar(pdfDict_.lookup("maxValue"))),
    expectation_(pdfDict_.lookup("expectation")),
    variance_(pdfDict_.lookup("variance")),
    strength_(pdfDict_.lookup("strength")),
    range_(maxValue_-minValue_)
{
    scalar sMax = 0;
    label n = strength_.size();
    for (label i=0; i<n; i++)
    {
        scalar x = expectation_[i];
        scalar s = strength_[i];
        for (label j=0; j<n; j++)
        {
            if (i!=j)
            {
                scalar x2 = (x-expectation_[j])/variance_[j];
                scalar y = exp(-0.5*x2*x2);
                s += strength_[j]*y;
            }
        }

        sMax = max(sMax, s);
    }

    for (label i=0; i<n; i++)
    {
        strength_[i] /= sMax;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::normal::~normal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::normal::sample() const
{
    scalar y = 0;
    scalar x = 0;
    label n = expectation_.size();
    bool success = false;

    while (!success)
    {
        x = minValue_ + range_*rndGen_.scalar01();
        y = rndGen_.scalar01();
        scalar p = 0.0;

        for (label i=0; i<n; i++)
        {
            scalar nu = expectation_[i];
            scalar sigma = variance_[i];
            scalar s = strength_[i];
            scalar v = (x-nu)/sigma;
            p += s*exp(-0.5*v*v);
        }

        if (y<p)
        {
            success = true;
        }
    }

    return x;
}


Foam::scalar Foam::normal::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::normal::maxValue() const
{
    return maxValue_;
}


// ************************************************************************* //
