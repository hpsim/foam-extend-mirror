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

#include "C14H30.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C14H30, 0);
    addToRunTimeSelectionTable(liquid, C14H30,);
    addToRunTimeSelectionTable(liquid, C14H30, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C14H30::C14H30()
:
    liquid
    (
        198.392,
        692.40,
        1.6212e+6,
        0.8428,
        0.237,
        279.01,
        1.8849e-1,
        526.73,
        0.0,
        0.6617,
        1.6173e+4
    ),
    rho_(60.92023144, 0.2582, 692.4, 0.26628),
    pv_(249.21, -16915, -35.195, 0.028451, 1.0),
    hl_(692.40, 455764.345336506, 0.428, 0.0, 0.0, 0.0),
    cp_
    (
        2565.72845679261,
       -4.78114036856325,
        0.0120362716238558,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
       -2690601.01887934,
        2565.72845679261,
       -2.39057018428162,
        0.00401209054128527,
        0.0,
        0.0
    ),
    cpg_(1134.11831122223, 3629.17859591113, -1440.3, 2275.29335860317, -682),
    B_
    (
        0.00247837614419936,
       -2.62692044034034,
       -1427174.48284205,
       -1.68288035807895e+19,
        3.48854792531957e+21
    ),
    mu_(-18.964, 2010.9, 1.0648, 0.0, 0.0),
    mug_(4.4565e-08, 0.8684, 228.16, -4347.2),
    K_(0.1957, -0.0001993, 0.0, 0.0, 0.0, 0.0),
    Kg_(-0.000628, 0.944, -5490, 0.0),
    sigma_(692.40, 0.056436, 1.3658, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 198.392, 28.0) // note: Same as nHeptane
{}


Foam::C14H30::C14H30
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


Foam::C14H30::C14H30(Istream& is)
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
