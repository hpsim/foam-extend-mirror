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

#include "C4H10O.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C4H10O, 0);
    addToRunTimeSelectionTable(liquid, C4H10O,);
    addToRunTimeSelectionTable(liquid, C4H10O, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C4H10O::C4H10O()
:
    liquid
    (
        74.123,
        466.70,
        3.6376e+6,
        0.28,
        0.262,
        156.85,
        4.0709e-1,
        307.58,
        3.836e-30,
        0.2846,
        1.5532e+4
    ),
    rho_(75.2793188, 0.27608, 466.7, 0.29358),
    pv_(101.03, -6311.5, -12.27, 1.377e-05, 2),
    hl_(466.70, 566355.921913576, 0.40717, 0, 0, 0),
    cp_
    (
        599.004357621791,
        17.5519069654493,
       -0.0742009902459426,
        0.00011822241409549,
        0.0,
        0.0
    ),
    h_
    (
       -4312350.92187216,
        599.004357621791,
        8.77595348272466,
       -0.0247336634153142,
        2.95556035238725e-05,
        0.0
    ),
    cpg_(1163.06679438231, 3441.57683849817, 1541.3, 1938.66950878944, -688.9),
    B_
    (
        0.00215992337061371,
       -1.810504162001,
       -276972.0599544,
       -2.12349742994752e+17,
        3.1016013922804e+19
    ),
    mu_(10.197, -63.8, -3.226, 0.0, 0.0),
    mug_(1.948e-06, 0.41, 495.8, 0.0),
    K_(0.249, -0.0004005, 0.0, 0.0, 0.0, 0.0),
    Kg_(-0.0044894, 0.6155, -3266.3, 0.0),
    sigma_(466.70, 0.057356, 1.288, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 74.123, 28) // note: Same as nHeptane
{}


Foam::C4H10O::C4H10O
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


Foam::C4H10O::C4H10O(Istream& is)
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
