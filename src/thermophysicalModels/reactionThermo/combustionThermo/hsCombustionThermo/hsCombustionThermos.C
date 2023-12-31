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

#include "makeHsCombustionThermo.H"

#include "hsCombustionThermo.H"
#include "hsPsiMixtureThermo.H"

#include "perfectGas.H"

#include "hConstThermo.H"
#include "janafThermo.H"
#include "specieThermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "dieselMixture.H"
#include "homogeneousMixture.H"
#include "inhomogeneousMixture.H"
#include "veryInhomogeneousMixture.H"

#include "reactingMixture.H"
#include "multiComponentMixture.H"

#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeHsCombustionThermo
(
    hsCombustionThermo,
    hsPsiMixtureThermo,
    homogeneousMixture,
    constTransport,
    hConstThermo,
    perfectGas
);

makeHsCombustionThermo
(
    hsCombustionThermo,
    hsPsiMixtureThermo,
    inhomogeneousMixture,
    constTransport,
    hConstThermo,
    perfectGas
);

makeHsCombustionThermo
(
    hsCombustionThermo,
    hsPsiMixtureThermo,
    veryInhomogeneousMixture,
    constTransport,
    hConstThermo,
    perfectGas
);

makeHsCombustionThermo
(
    hsCombustionThermo,
    hsPsiMixtureThermo,
    homogeneousMixture,
    sutherlandTransport,
    janafThermo,
    perfectGas
);

makeHsCombustionThermo
(
    hsCombustionThermo,
    hsPsiMixtureThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    janafThermo,
    perfectGas
);

makeHsCombustionThermo
(
    hsCombustionThermo,
    hsPsiMixtureThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    janafThermo,
    perfectGas
);

makeHsCombustionThermo
(
    hsCombustionThermo,
    hsPsiMixtureThermo,
    dieselMixture,
    sutherlandTransport,
    janafThermo,
    perfectGas
);

// Multi-component thermo

makeHsCombustionMixtureThermo
(
    hsCombustionThermo,
    hsPsiMixtureThermo,
    multiComponentMixture,
    gasThermoPhysics
);


// Multi-component reaction thermo

makeHsCombustionMixtureThermo
(
    hsCombustionThermo,
    hsPsiMixtureThermo,
    reactingMixture,
    constGasThermoPhysics
);

makeHsCombustionMixtureThermo
(
    hsCombustionThermo,
    hsPsiMixtureThermo,
    reactingMixture,
    gasThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
