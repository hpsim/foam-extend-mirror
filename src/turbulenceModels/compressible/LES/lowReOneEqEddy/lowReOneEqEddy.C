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

#include "lowReOneEqEddy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(lowReOneEqEddy, 0);
addToRunTimeSelectionTable(LESModel, lowReOneEqEddy, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void lowReOneEqEddy::updateSubGridScaleFields()
{
    // High Re eddy viscosity
    muSgs_ = ck_*rho()*sqrt(k_)*delta();

    // low Re no corrected eddy viscosity
    muSgs_ -= (mu()/beta_)*(scalar(1) - exp(-beta_*muSgs_/mu()));
    muSgs_.correctBoundaryConditions();

    alphaSgs_ = muSgs_/Prt_;
    alphaSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lowReOneEqEddy::lowReOneEqEddy
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermophysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, rho, U, phi, thermophysicalModel, turbulenceModelName),
    GenEddyVisc(rho, U, phi, thermophysicalModel),

    ck_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.07
        )
    ),
    beta_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "beta",
            coeffDict_,
            0.01
        )
    )
{
    updateSubGridScaleFields();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void lowReOneEqEddy::correct(const tmp<volTensorField>& tgradU)
{
    const volTensorField& gradU = tgradU();

    GenEddyVisc::correct(gradU);

    volScalarField divU = fvc::div(phi()/fvc::interpolate(rho()));
    volScalarField G = 2*muSgs_*(gradU && dev(symm(gradU)));

    solve
    (
        fvm::ddt(rho(), k_)
      + fvm::div(phi(), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::SuSp(2.0/3.0*rho()*divU, k_)
      - fvm::Sp(ce_*rho()*sqrt(k_)/delta(), k_)
    );

    bound(k_, k0());

    updateSubGridScaleFields();
}


bool lowReOneEqEddy::read()
{
    if (GenEddyVisc::read())
    {
        ck_.readIfPresent(coeffDict());
        beta_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
