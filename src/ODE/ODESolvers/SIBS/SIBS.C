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

Description
    Semi-Implicit Burlisch-Stoer scheme for stiff equation sets
    Numerical Recipes in C, Secn 16.6 page 742-747
    Alternative reference in Numerical Recipes in C++

\*---------------------------------------------------------------------------*/

#include "SIBS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::SIBS, 0);

namespace Foam
{
    addToRunTimeSelectionTable(ODESolver, SIBS, ODE);

    const label SIBS::nSeq_[iMaxX_] = {2, 6, 10, 14, 22, 34, 50, 70};

    const scalar
        SIBS::safe1 = 0.25,
        SIBS::safe2 = 0.7,
        SIBS::redMax = 1.0e-5,
        SIBS::redMin = 0.7,
        SIBS::scaleMX = 0.1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SIBS::SIBS(ODE& ode)
:
    ODESolver(ode),
    a_(iMaxX_),
    alpha_(kMaxX_, scalar(0)),
    d_p_(ode_.nEqns(), kMaxX_),
    x_p_(kMaxX_),
    err_(kMaxX_),

    yTemp_(ode_.nEqns()),
    ySeq_(ode_.nEqns()),
    yErr_(ode_.nEqns()),
    dfdx_(ode_.nEqns()),
    dfdy_(ode_.nEqns(), scalar(0)),
    first_(1),
    epsOld_(-1.0),
    xNew_(0.0)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::SIBS::solve
(
    scalar& x,
    scalarField& y,
    scalarField& dydx,
    const scalar eps,
    const scalarField& yScale,
    const scalar hTry,
    scalar& hDid,
    scalar& hNext
) const
{
    bool exitflag = false;
    const label nEqns = ode_.nEqns();

    if (eps != epsOld_)
    {
        hNext = xNew_ = -GREAT;
        scalar eps1 = safe1*eps;
        a_[0] = nSeq_[0] + 1;

        for (label k=0; k<kMaxX_; k++)
        {
            a_[k + 1] = a_[k] + nSeq_[k + 1];
        }

        for (label iq = 1; iq<kMaxX_; iq++)
        {
            for (label k=0; k<iq; k++)
            {
                alpha_[k][iq] =
                    pow(eps1, (a_[k + 1] - a_[iq + 1])
                   /((a_[iq + 1] - a_[0] + 1.0)*(2*k + 3)));
            }
        }

        epsOld_ = eps;
        a_[0] += nEqns;

        for (label k=0; k<kMaxX_; k++)
        {
            a_[k + 1] = a_[k] + nSeq_[k + 1];
        }

        for (kOpt_ = 1; kOpt_<kMaxX_ - 1; kOpt_++)
        {
            if (a_[kOpt_ + 1] > a_[kOpt_]*alpha_[kOpt_ - 1][kOpt_])
            {
                break;
            }
        }

        kMax_ = kOpt_;
    }

    label k = 0;
    scalar h = hTry;
    yTemp_ = y;

    ode_.jacobian(x, y, dfdx_, dfdy_);

    if (x != xNew_ || h != hNext)
    {
        first_ = 1;
        kOpt_ = kMax_;
    }

    label km=0;
    label reduct=0;
    scalar maxErr = SMALL;

    for (;;)
    {
        scalar red = 1.0;

        for (k=0; k <= kMax_; k++)
        {
            xNew_ = x + h;

            if (xNew_ == x)
            {
                FatalErrorIn("ODES::SIBS")
                    << "step size underflow"
                    << exit(FatalError);
            }

            SIMPR(x, yTemp_, dydx, dfdx_, dfdy_, h, nSeq_[k], ySeq_);
            scalar xest = sqr(h/nSeq_[k]);

            polyExtrapolate(k, xest, ySeq_, y, yErr_, x_p_, d_p_);

            if (k != 0)
            {
                maxErr = SMALL;
                for (label i=0; i<nEqns; i++)
                {
                    maxErr = max(maxErr, mag(yErr_[i]/yScale[i]));
                }
                maxErr /= eps;
                km = k - 1;
                err_[km] = pow(maxErr/safe1, 1.0/(2*km + 3));
            }

            if (k != 0 && (k >= kOpt_ - 1 || first_))
            {
                if (maxErr < 1.0)
                {
                    exitflag = true;
                    break;
                }

                if (k == kMax_ || k == kOpt_ + 1)
                {
                    red = safe2/err_[km];
                    break;
                }
                else if (k == kOpt_ && alpha_[kOpt_ - 1][kOpt_] < err_[km])
                {
                    red = 1.0/err_[km];
                    break;
                }
                else if (kOpt_ == kMax_ && alpha_[km][kMax_ - 1] < err_[km])
                {
                    red = alpha_[km][kMax_ - 1]*safe2/err_[km];
                    break;
                }
                else if (alpha_[km][kOpt_] < err_[km])
                {
                    red = alpha_[km][kOpt_ - 1]/err_[km];
                    break;
                }
            }
        }

        if (exitflag)
        {
            break;
        }

        red = min(red, redMin);
        red = max(red, redMax);
        h *= red;
        reduct = 1;
    }

    x = xNew_;
    hDid = h;
    first_=0;
    scalar wrkmin = GREAT;
    scalar scale = 1.0;

    for (label kk=0; kk<=km; kk++)
    {
        scalar fact = max(err_[kk], scaleMX);
        scalar work = fact*a_[kk + 1];
        if (work < wrkmin)
        {
            scale = fact;
            wrkmin = work;
            kOpt_ = kk + 1;
        }
    }

    hNext = h/scale;

    if (kOpt_ >= k && kOpt_ != kMax_ && !reduct)
    {
        scalar fact = max(scale/alpha_[kOpt_ - 1][kOpt_], scaleMX);
        if (a_[kOpt_ + 1]*fact <= wrkmin)
        {
            hNext = h/fact;
            kOpt_++;
        }
    }
}


// ************************************************************************* //
