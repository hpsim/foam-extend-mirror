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

#include "general.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(general, 0);
    addToRunTimeSelectionTable(pdf, general, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::general::general(const dictionary& dict, Random& rndGen)
:
    pdf(dict, rndGen),
    pdfDict_(dict.subDict(typeName + "PDF")),
    xy_(pdfDict_.lookup("distribution")),
    nEntries_(xy_.size()),
    minValue_(xy_[0][0]),
    maxValue_(xy_[nEntries_-1][0]),
    range_(maxValue_-minValue_)
{
    // normalize the pdf
    scalar yMax = 0;

    for (label i=0; i<nEntries_; i++)
    {
        yMax = max(yMax, xy_[i][1]);
    }

    for (label i=0; i<nEntries_; i++)
    {
        xy_[i][1] /= yMax;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::general::~general()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::general::sample() const
{
    scalar y = 0;
    scalar x = 0;

    bool success = false;

    while (!success)
    {
        x = minValue_ + range_*rndGen_.scalar01();
        y = rndGen_.scalar01();

        bool intervalFound = false;
        label i = -1;
        while (!intervalFound)
        {
            i++;
            if ( (x>xy_[i][0]) && (x<xy_[i+1][0]) )
            {
                intervalFound = true;
            }
        }

        scalar p =
            xy_[i][1]
          + (x-xy_[i][0])
           *(xy_[i+1][1]-xy_[i][1])
           /(xy_[i+1][0]-xy_[i][0]);

        if (y<p)
        {
            success = true;
        }
    }

    return x;
}


Foam::scalar Foam::general::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::general::maxValue() const
{
    return maxValue_;
}


// ************************************************************************* //
