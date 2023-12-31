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

#include "DeardorffDiffStress.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(DeardorffDiffStress, 0);
addToRunTimeSelectionTable(LESModel, DeardorffDiffStress, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void DeardorffDiffStress::updateSubGridScaleFields(const volScalarField& K)
{
    nuSgs_ = ck_*sqrt(K)*delta();
    nuSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DeardorffDiffStress::DeardorffDiffStress
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),
    GenSGSStress(U, phi, transport),

    ck_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.094
        )
    ),
    cm_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cm",
            coeffDict_,
            4.13
        )
    )
{
    updateSubGridScaleFields(0.5*tr(B_));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DeardorffDiffStress::correct(const tmp<volTensorField>& tgradU)
{
    const volTensorField& gradU = tgradU();

    GenSGSStress::correct(gradU);

    volSymmTensorField D = symm(gradU);

    volSymmTensorField P = -twoSymm(B_ & gradU);

    volScalarField K = 0.5*tr(B_);
    volScalarField Epsilon = 2*nuEff()*magSqr(D);

    fvSymmTensorMatrix BEqn
    (
        fvm::ddt(B_)
      + fvm::div(phi(), B_)
      - fvm::laplacian(DBEff(), B_)
      + fvm::Sp(cm_*sqrt(K)/delta(), B_)
     ==
        P
      + 0.8*K*D
      - (2*ce_ - 0.667*cm_)*I*Epsilon
    );

    BEqn.relax();
    BEqn.solve();

    // Bounding the component kinetic energies

    forAll(B_, celli)
    {
        B_[celli].component(symmTensor::XX) =
            max(B_[celli].component(symmTensor::XX), k0().value());
        B_[celli].component(symmTensor::YY) =
            max(B_[celli].component(symmTensor::YY), k0().value());
        B_[celli].component(symmTensor::ZZ) =
            max(B_[celli].component(symmTensor::ZZ), k0().value());
    }

    K = 0.5*tr(B_);
    bound(K, k0());

    updateSubGridScaleFields(K);
}


bool DeardorffDiffStress::read()
{
    if (GenSGSStress::read())
    {
        ck_.readIfPresent(coeffDict());
        cm_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
