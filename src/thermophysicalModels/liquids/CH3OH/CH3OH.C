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

#include "CH3OH.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CH3OH, 0);
    addToRunTimeSelectionTable(liquid, CH3OH,);
    addToRunTimeSelectionTable(liquid, CH3OH, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CH3OH::CH3OH()
:
    liquid
    (
        32.042,
        512.58,
        8.0959e+6,
        0.1178,
        0.224,
        175.47,
        1.054e-1,
        337.85,
        5.6706e-30,
        0.5656,
        2.9523e+4
    ),
    rho_(73.952936, 0.27192, 512.58, 0.2331),
    pv_(109.93, -7471.3, -13.988, 0.015281, 1.0),
    hl_(512.58, 1644716.30984333, 0.3766, 0.0, 0.0, 0.0),
    cp_
    (
        3358.09250358904,
       -11.8781599151114,
        0.0305536483365583,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
       -8190474.32066862,
        3358.09250358904,
       -5.93907995755571,
        0.0101845494455194,
        0.0,
        0.0
    ),
    cpg_(1226.9521253355, 2772.92303851195, 1963, 1733.66206853505, 909.6),
    B_
    (
       -0.0199737844079645,
        19.3496036452157,
       -3342487.98452032,
        2.40808938268523e+19,
       -6.85787404032208e+21
    ),
    mu_(-7.288, 1065.3, -0.6657, 0.0, 0.0),
    mug_(3.0663e-07, 0.69655, 205.0, 0.0),
    K_(0.2837, -0.000281, 0.0, 0.0, 0.0, 0.0),
    Kg_(-7.763, 1.0279, -74360000.0, 6770000000.0),
    sigma_(512.58, 0.056, -0.00014583, 1.08e-07, 0.0, 0.0),
    D_(147.18, 20.1, 32.042, 28.0) // note: Same as nHeptane
{}


Foam::CH3OH::CH3OH
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


Foam::CH3OH::CH3OH(Istream& is)
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
