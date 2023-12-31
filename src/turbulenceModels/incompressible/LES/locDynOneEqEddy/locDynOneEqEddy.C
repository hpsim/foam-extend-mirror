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

#include "locDynOneEqEddy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(locDynOneEqEddy, 0);
addToRunTimeSelectionTable(LESModel, locDynOneEqEddy, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void locDynOneEqEddy::updateSubGridScaleFields
(
    const volSymmTensorField& D,
    const volScalarField& KK
)
{
    nuSgs_ = ck(D, KK)*sqrt(k_)*delta();
    nuSgs_.correctBoundaryConditions();
}


volScalarField locDynOneEqEddy::ck
(
    const volSymmTensorField& D,
    const volScalarField& KK
) const
{
    volSymmTensorField LL =
        simpleFilter_(dev(filter_(sqr(U())) - (sqr(filter_(U())))));

    volSymmTensorField MM = simpleFilter_(-2.0*delta()*pow(KK, 0.5)*filter_(D));

    volScalarField ck =
        simpleFilter_(0.5*(LL && MM))
       /(
            simpleFilter_(magSqr(MM))
          + dimensionedScalar("small", sqr(MM.dimensions()), VSMALL)
        );

    return 0.5*(mag(ck) + ck);
}


volScalarField locDynOneEqEddy::ce
(
    const volSymmTensorField& D,
    const volScalarField& KK
) const
{
    volScalarField ce =
        simpleFilter_(nuEff()*(filter_(magSqr(D)) - magSqr(filter_(D))))
       /simpleFilter_(pow(KK, 1.5)/(2.0*delta()));

    return 0.5*(mag(ce) + ce);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

locDynOneEqEddy::locDynOneEqEddy
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),
    GenEddyVisc(U, phi, transport),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    simpleFilter_(U.mesh()),
    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    volScalarField KK = 0.5*(filter_(magSqr(U)) - magSqr(filter_(U)));
    updateSubGridScaleFields(symm(fvc::grad(U)), KK);

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void locDynOneEqEddy::correct(const tmp<volTensorField>& gradU)
{
    LESModel::correct(gradU);

    volSymmTensorField D = symm(gradU);

    volScalarField KK = 0.5*(filter_(magSqr(U())) - magSqr(filter_(U())));
    KK.max(dimensionedScalar("small", KK.dimensions(), SMALL));

    volScalarField P = 2.0*nuSgs_*magSqr(D);

    fvScalarMatrix kEqn
    (
       fvm::ddt(k_)
     + fvm::div(phi(), k_)
     - fvm::laplacian(DkEff(), k_)
    ==
       P
     - fvm::Sp(ce(D, KK)*sqrt(k_)/delta(), k_)
    );

    kEqn.relax();
    kEqn.solve();

    bound(k_, k0());

    updateSubGridScaleFields(D, KK);
}


bool locDynOneEqEddy::read()
{
    if (GenEddyVisc::read())
    {
        filter_.read(coeffDict());

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
