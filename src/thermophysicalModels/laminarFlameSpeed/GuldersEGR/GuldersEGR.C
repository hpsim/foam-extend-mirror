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

#include "GuldersEGR.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarFlameSpeedModels
{
    defineTypeNameAndDebug(GuldersEGR, 0);

    addToRunTimeSelectionTable
    (
        laminarFlameSpeed,
        GuldersEGR,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::GuldersEGR::GuldersEGR
(
    const dictionary& dict,
    const hhuCombustionThermo& ct
)
:
    laminarFlameSpeed(dict, ct),

    coeffsDict_(dict.subDict(typeName + "Coeffs").subDict(fuel_)),
    W_(readScalar(coeffsDict_.lookup("W"))),
    eta_(readScalar(coeffsDict_.lookup("eta"))),
    xi_(readScalar(coeffsDict_.lookup("xi"))),
    f_(readScalar(coeffsDict_.lookup("f"))),
    alpha_(readScalar(coeffsDict_.lookup("alpha"))),
    beta_(readScalar(coeffsDict_.lookup("beta")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::GuldersEGR::~GuldersEGR()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar Foam::laminarFlameSpeedModels::GuldersEGR::SuRef
(
    scalar phi
) const
{
    if (phi > SMALL)
    {
        return W_*pow(phi, eta_)*exp(-xi_*sqr(phi - 1.075));
    }
    else
    {
        return 0.0;
    }
}


inline Foam::scalar Foam::laminarFlameSpeedModels::GuldersEGR::Su0pTphi
(
    scalar p,
    scalar Tu,
    scalar phi,
    scalar Yres
) const
{
    static const scalar Tref = 300.0;
    static const scalar pRef = 1.013e5;

    return SuRef(phi)*pow((Tu/Tref), alpha_)*pow((p/pRef), beta_)*(1 - f_*Yres);
}


Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::GuldersEGR::Su0pTphi
(
    const volScalarField& p,
    const volScalarField& Tu,
    scalar phi
) const
{
    tmp<volScalarField> tSu0
    (
        new volScalarField
        (
            IOobject
            (
                "Su0",
                p.time().timeName(),
                p.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p.mesh(),
            dimensionedScalar("Su0", dimVelocity, 0.0)
        )
    );

    volScalarField& Su0 = tSu0();

    forAll(Su0, celli)
    {
        Su0[celli] = Su0pTphi(p[celli], Tu[celli], phi, 0.0);
    }

    forAll(Su0.boundaryField(), patchi)
    {
        forAll(Su0.boundaryField()[patchi], facei)
        {
            Su0.boundaryField()[patchi][facei] =
                Su0pTphi
                (
                    p.boundaryField()[patchi][facei],
                    Tu.boundaryField()[patchi][facei],
                    phi,
                    0.0
                );
        }
    }

    return tSu0;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::GuldersEGR::Su0pTphi
(
    const volScalarField& p,
    const volScalarField& Tu,
    const volScalarField& phi,
    const volScalarField& egr
) const
{
    tmp<volScalarField> tSu0
    (
        new volScalarField
        (
            IOobject
            (
                "Su0",
                p.time().timeName(),
                p.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p.mesh(),
            dimensionedScalar("Su0", dimVelocity, 0.0)
        )
    );

    volScalarField& Su0 = tSu0();

    forAll(Su0, celli)
    {
        Su0[celli] = Su0pTphi(p[celli], Tu[celli], phi[celli], egr[celli]);
    }

    forAll(Su0.boundaryField(), patchi)
    {
        forAll(Su0.boundaryField()[patchi], facei)
        {
            Su0.boundaryField()[patchi][facei] =
                Su0pTphi
                (
                    p.boundaryField()[patchi][facei],
                    Tu.boundaryField()[patchi][facei],
                    phi.boundaryField()[patchi][facei],
                    egr.boundaryField()[patchi][facei]
                );
        }
    }

    return tSu0;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::GuldersEGR::operator()() const
{
    if
    (
        hhuCombustionThermo_.composition().contains("ft")
     && hhuCombustionThermo_.composition().contains("egr")
    )
    {
        return Su0pTphi
        (
            hhuCombustionThermo_.p(),
            hhuCombustionThermo_.Tu(),
            dimensionedScalar
            (
                hhuCombustionThermo_.lookup("stoichiometricAirFuelMassRatio")
            )/
            (
                scalar(1)/hhuCombustionThermo_.composition().Y("ft")
              - scalar(1)
            ),
            hhuCombustionThermo_.composition().Y("egr")
        );
    }
    else
    {
        return Su0pTphi
        (
            hhuCombustionThermo_.p(),
            hhuCombustionThermo_.Tu(),
            equivalenceRatio_
        );
    }
}


// ************************************************************************* //
