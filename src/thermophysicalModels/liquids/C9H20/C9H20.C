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

#include "C9H20.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C9H20, 0);
    addToRunTimeSelectionTable(liquid, C9H20,);
    addToRunTimeSelectionTable(liquid, C9H20, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C9H20::C9H20()
:
    liquid
    (
        128.258,
        594.60,
        2.29e+6,
        0.544,
        0.252,
        219.66,
        4.3058e-1,
        423.97,
        0.0,
        0.4435,
        1.56e+4
    ),
    rho_(62.06019846, 0.26147, 594.6, 0.28281),
    pv_(109.35, -90304.0, -12.882, 7.8544e-06, 2.0),
    hl_(594.60, 470691.886665939, 0.38522, 0.0, 0.0, 0.0),
    cp_
    (
        2986.79224687739,
       -8.88677509395125,
        0.0211300659607978,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
       -2825628.50868792,
        2986.79224687739,
       -4.44338754697563,
        0.00704335532026592,
        0.0,
        0.0
    ),
    cpg_(1183.16206396482, 3832.11963386299, 1644.8, 2705.48425829188, 749.6),
    B_
    (
        0.00304542406711472,
       -3.65357326638494,
       -520825.211682702,
       -6.15400208953827e+18,
        1.41901479829718e+21
    ),
    mu_(-21.149, 1658, 1.454, 0.0, 0.0),
    mug_(1.0344e-07, 0.77301, 220.47, 0.0),
    K_(0.209, -0.000264, 0.0, 0.0, 0.0, 0.0),
    Kg_(-0.065771, 0.27198, -3482.3, -1580300.0),
    sigma_(594.60, 0.054975, 1.2897, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 128.258, 28.0) // note: Same as nHeptane
{}


Foam::C9H20::C9H20
(
    const liquid& l,
    const NSRDSfunc5& density,
    const NSRDSfunc1& vapourPressure,
    const NSRDSfunc6& heatOfVapourisation,
    const NSRDSfunc0& heatCapacity,
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


Foam::C9H20::C9H20(Istream& is)
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
