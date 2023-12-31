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

#include "turbulenceModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbulenceModel, 0);
defineRunTimeSelectionTable(turbulenceModel, turbulenceModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulenceModel::turbulenceModel
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermophysicalModel,
    const word& turbulenceModelName
)
:
    regIOobject
    (
        IOobject
        (
            turbulenceModelName,
            U.time().constant(),
            U.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(U.time()),
    mesh_(U.mesh()),

    rho_(rho),
    U_(U),
    phi_(phi),
    thermophysicalModel_(thermophysicalModel),

    y_(mesh_)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<turbulenceModel> turbulenceModel::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermophysicalModel,
    const word& turbulenceModelName
)
{
    word modelName;

    // Enclose the creation of the dictionary to ensure it is deleted
    // before the turbulenceModel is created otherwise the dictionary is
    // entered in the database twice
    {
        IOdictionary dict
        (
            IOobject
            (
                "turbulenceProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

        dict.lookup("simulationType") >> modelName;
    }

    Info<< "Selecting turbulence model type " << modelName << endl;

    turbulenceModelConstructorTable::iterator cstrIter =
        turbulenceModelConstructorTablePtr_->find(modelName);

    if (cstrIter == turbulenceModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "turbulenceModel::New(const volScalarField&, "
            "const volVectorField&, const surfaceScalarField&, "
            "basicThermo&)"
        )   << "Unknown turbulenceModel type " << modelName
            << endl << endl
            << "Valid turbulenceModel types are :" << endl
            << turbulenceModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<turbulenceModel>
    (
        cstrIter()(rho, U, phi, thermophysicalModel, turbulenceModelName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar turbulenceModel::yPlusLam(const scalar kappa, const scalar E) const
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(E*ypl)/kappa;
    }

    return ypl;
}


tmp<volScalarField> turbulenceModel::rhoEpsilonEff() const
{
    tmp<volTensorField> tgradU = fvc::grad(U_);
    return mu()*(tgradU() && dev(twoSymm(tgradU()))) + rho_*epsilon();
}


void turbulenceModel::correct()
{
    if (mesh_.changing())
    {
        y_.correct();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
