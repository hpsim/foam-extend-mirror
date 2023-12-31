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
    RiddersRoot

Description
    Ridder's method of root finiding given a function, bracketed root
    and accuracy.  Based on Numerical Recipes in C++, Section 9.2,
    page 362.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*----------------------------------------------------------------------------*/

#include "RiddersRoot.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Func>
const Foam::label Foam::RiddersRoot<Func>::maxIter = 60;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Func>
Foam::RiddersRoot<Func>::RiddersRoot(const Func& f, const scalar eps)
:
    f_(f),
    eps_(eps)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Func>
Foam::scalar Foam::RiddersRoot<Func>::root
(
    const scalar x0,
    const scalar x1
) const
{
    scalar fl = f_(x0);
    scalar fh = f_(x1);

    // Check bracketing of the root
    if ((fl > 0 && fh < 0) || (fl < 0 && fh > 0))
    {
        scalar xl = x0;
        scalar xh = x1;

        // Bad guess for answer, to simplify logic below
        scalar ans = -1.11e-30;

        scalar xm, fm, s, xNew, fNew;

        for (label nIter = 0; nIter < maxIter; nIter++)
        {
            // Create a mean value and evaluate function
            xm = 0.5*(xl + xh);
            fm = f_(xm);

            // Solve quadratic equation if well-posed
            s = sqrt(sqr(fm) - fl*fh);

            if (s < SMALL)
            {
                return ans;
            }

            // Updating formula
            if (fl >= fh)
            {
                xNew = xm + (xm - xl)*fm/s;
            }
            else
            {
                xNew = xm - (xm - xl)*fm/s;
            }

            // Check answer
            if (mag(xNew - ans) <= eps_)
            {
                return ans;
            }

            ans = xNew;

            fNew = f_(ans);

            if (mag(fNew) < SMALL)
            {
                return ans;
            }

            if (checkSign(fm, fNew))
            {
                xl = xm;
                fl = fm;
                xh = ans;
                fh = fNew;
            }
            else if (checkSign(fl, fNew))
            {
                xh = ans;
                fh = fNew;
            }
            else if (checkSign(fh, fNew))
            {
                xl = ans;
                fl = fNew;
            }
            else
            {
                // Never get here
                FatalErrorIn
                (
                    "Foam::scalar Foam::RiddersRoot<Func>::root\n"
                    "(\n"
                    "    const scalar x0,\n"
                    "    const scalar x1\n"
                    ") const"
                )   << "Error in search logic" << abort(FatalError);
            }

            if (mag(xh - xl) <= eps_)
            {
                return ans;
            }
        }

        FatalErrorIn
        (
            "Foam::scalar Foam::RiddersRoot<Func>::root\n"
            "(\n"
            "    const scalar x0,\n"
            "    const scalar x1\n"
            ") const"
        )   << "Maximum number of iterations exceeded" << abort(FatalError);
    }
    else if (mag(fl) < SMALL)
    {
        return x0;
    }
    else if (mag(fh) < SMALL)
    {
        return x1;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::scalar Foam::RiddersRoot<Func>::root\n"
            "(\n"
            "    const scalar x0,\n"
            "    const scalar x1\n"
            ") const"
        )   << "Root is not bracketed.  f(x0) = " << fl << " f(x1) = " << fh
            << abort(FatalError);

    }

    // Dummy return to keep compiler happy
    return x0;
}


// ************************************************************************* //
