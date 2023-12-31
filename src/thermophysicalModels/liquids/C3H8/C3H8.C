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

#include "C3H8.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C3H8, 0);
    addToRunTimeSelectionTable(liquid, C3H8,);
    addToRunTimeSelectionTable(liquid, C3H8, Istream);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C3H8::C3H8()
:
    liquid
    (
        44.096,
        369.83,
        4.248e+6,
        0.2, 0.276,
        85.47,
        1.685e-4,
        231.11,
        0.0,
        0.1523,
        1.31e+4
    ),
    rho_(60.6628672, 0.27453, 369.83, 0.29359),
    pv_(59.078, -3492.6, -6.0669, 1.0919e-05, 2.0),
    hl_(369.83, 662395.682148041, 0.78237, -0.77319, 0.39246, 0.0),
    cp_
    (
        369.83,
        9.48470319647089,
        2576.87772133527,
        95.3560311677331,
       -131.535634282099
    ),
    h_(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    cpg_(1177.43105950653, 4364.34143686502, 1626.5, 2648.76632801161, 723.6),
    B_
    (
        0.00255578737300435,
       -2.24963715529753,
       -102276.850507983,
        7.00743831640058e+15,
       -1.59878447024673e+18
    ),
    mu_(-6.9281, 420.76, -0.63276, -1.713e-26, 10.0),
    mug_(2.4993e-07, 0.68612, 179.34, -8254.6),
    K_(0.26755, -0.00066457, 2.774e-07, 0.0, 0.0, 0.0),
    Kg_(-1.12, 0.10972, -9834.6, -7535800),
    sigma_(369.83, 0.05092, 1.2197, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 44.096, 28) // note: Same as nHeptane
{}


Foam::C3H8::C3H8
(
    const liquid& l,
    const NSRDSfunc5& density,
    const NSRDSfunc1& vapourPressure,
    const NSRDSfunc6& heatOfVapourisation,
    const NSRDSfunc14& heatCapacity,
    const NSRDSfunc0& enthalpy,
    const NSRDSfunc7& idealGasHeatCapacity,
    const NSRDSfunc4& secondVirialCoeff,
    const NSRDSfunc1& dynamicViscosity,
    const NSRDSfunc2& vapourDynamicViscosity,
    const NSRDSfunc0& thermalConductivity,
    const NSRDSfunc2& vapourThermalConductivity,
    const NSRDSfunc6& surfaceTension,
    const APIdiffCoefFunc& vapourDiffussivity
)
:
    liquid(l),
    rho_(density),
    pv_(vapourPressure),
    hl_(heatOfVapourisation),
    cp_(heatCapacity),
    h_(enthalpy),
    cpg_(idealGasHeatCapacity),
    B_(secondVirialCoeff),
    mu_(dynamicViscosity),
    mug_(vapourDynamicViscosity),
    K_(thermalConductivity),
    Kg_(vapourThermalConductivity),
    sigma_(surfaceTension),
    D_(vapourDiffussivity)
{}


Foam::C3H8::C3H8(Istream& is)
:
    liquid(is),
    rho_(is),
    pv_(is),
    hl_(is),
    cp_(is),
    h_(is),
    cpg_(is),
    B_(is),
    mu_(is),
    mug_(is),
    K_(is),
    Kg_(is),
    sigma_(is),
    D_(is)
{}


// ************************************************************************* //
