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
    Fourth-order Kaps-Rentrop scheme with adjustive time-step size
    Numerical Recipes in C, Secn 16.6 page 739-742
    Alternative reference in Numerical Recipes in C++

\*---------------------------------------------------------------------------*/

#include "KRR4.H"
#include "simpleMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::KRR4, 0);

namespace Foam
{
    addToRunTimeSelectionTable(ODESolver, KRR4, ODE);

const scalar
    KRR4::safety = 0.9, KRR4::grow = 1.5, KRR4::pgrow = -0.25,
    KRR4::shrink = 0.5, KRR4::pshrink = (-1.0/3.0), KRR4::errcon = 0.1296;

const scalar
    KRR4::gamma = 1.0/2.0,
    KRR4::a21 = 2.0, KRR4::a31 = 48.0/25.0, KRR4::a32 = 6.0/25.0,
    KRR4::c21 = -8.0, KRR4::c31 = 372.0/25.0, KRR4::c32 = 12.0/5.0,
    KRR4::c41 = -112.0/125.0, KRR4::c42 = -54.0/125.0, KRR4::c43 = -2.0/5.0,
    KRR4::b1 = 19.0/9.0, KRR4::b2 = 1.0/2.0, KRR4::b3 = 25.0/108.0,
    KRR4::b4 = 125.0/108.0,
    KRR4::e1 = 17.0/54.0, KRR4::e2 = 7.0/36.0, KRR4::e3 = 0.0,
    KRR4::e4 = 125.0/108.0,
    KRR4::c1X = 1.0/2.0, KRR4::c2X = -3.0/2.0, KRR4::c3X = 121.0/50.0,
    KRR4::c4X = 29.0/250.0,
    KRR4::a2X = 1.0, KRR4::a3X = 3.0/5.0;
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::KRR4::KRR4(ODE& ode)
:
    ODESolver(ode),
    yTemp_(ode_.nEqns()),
    dydxTemp_(ode_.nEqns()),
    g1_(ode_.nEqns()),
    g2_(ode_.nEqns()),
    g3_(ode_.nEqns()),
    g4_(ode_.nEqns()),
    yErr_(ode_.nEqns()),
    dfdx_(ode_.nEqns()),
    dfdy_(ode_.nEqns(), scalar(0)),
    a_(ode_.nEqns(), scalar(0)),
    pivotIndices_(ode_.nEqns())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::KRR4::solve
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
    scalar xTemp = x;
    yTemp_ = y;
    dydxTemp_ = dydx;
    const label nEqns = ode_.nEqns();

    ode_.jacobian(xTemp, yTemp_, dfdx_, dfdy_);

    scalar h = hTry;

    for (label jtry=0; jtry<maxtry; jtry++)
    {
        for (label i=0; i<nEqns; i++)
        {
            for (label j=0; j<nEqns; j++)
            {
                a_[i][j] = -dfdy_[i][j];
            }

            a_[i][i] += 1.0/(gamma*h);
        }

        simpleMatrix<scalar>::LUDecompose(a_, pivotIndices_);

        for (label i=0; i<nEqns; i++)
        {
            g1_[i] = dydxTemp_[i] + h*c1X*dfdx_[i];
        }

        simpleMatrix<scalar>::LUBacksubstitute(a_, pivotIndices_, g1_);

        for (label i=0; i<nEqns; i++)
        {
            y[i] = yTemp_[i] + a21*g1_[i];
        }

        x = xTemp + a2X*h;
        ode_.derivatives(x, y, dydx_);

        for (label i=0; i<nEqns; i++)
        {
            g2_[i] = dydx_[i] + h*c2X*dfdx_[i] + c21*g1_[i]/h;
        }

        simpleMatrix<scalar>::LUBacksubstitute(a_, pivotIndices_, g2_);

        for (label i=0; i<nEqns; i++)
        {
            y[i] = yTemp_[i] + a31*g1_[i] + a32*g2_[i];
        }

        x = xTemp + a3X*h;
        ode_.derivatives(x, y, dydx_);

        for (label i=0; i<nEqns; i++)
        {
            g3_[i] = dydx[i] + h*c3X*dfdx_[i] + (c31*g1_[i] + c32*g2_[i])/h;
        }

        simpleMatrix<scalar>::LUBacksubstitute(a_, pivotIndices_, g3_);

        for (label i=0; i<nEqns; i++)
        {
            g4_[i] = dydx_[i] + h*c4X*dfdx_[i]
                + (c41*g1_[i] + c42*g2_[i] + c43*g3_[i])/h;
        }

        simpleMatrix<scalar>::LUBacksubstitute(a_, pivotIndices_, g4_);

        for (label i=0; i<nEqns; i++)
        {
            y[i] = yTemp_[i] + b1*g1_[i] + b2*g2_[i] + b3*g3_[i] + b4*g4_[i];
            yErr_[i] = e1*g1_[i] + e2*g2_[i] + e3*g3_[i] + e4*g4_[i];
        }

        x = xTemp + h;

        if (x == xTemp)
        {
            FatalErrorIn("ODES::KRR4")
                << "stepsize not significant"
                << exit(FatalError);
        }

        scalar maxErr = 0.0;
        for (label i=0; i<nEqns; i++)
        {
            maxErr = max(maxErr, mag(yErr_[i]/yScale[i]));
        }
        maxErr /= eps;

        if (maxErr <= 1.0)
        {
            hDid = h;
            hNext = (maxErr > errcon ? safety*h*pow(maxErr, pgrow) : grow*h);
            return;
        }
        else
        {
            hNext = safety*h*pow(maxErr, pshrink);
            h = (h >= 0.0 ? max(hNext, shrink*h) : min(hNext, shrink*h));
        }
    }

    FatalErrorIn("ODES::KRR4")
        << "exceeded maxtry"
        << exit(FatalError);
}


// ************************************************************************* //
