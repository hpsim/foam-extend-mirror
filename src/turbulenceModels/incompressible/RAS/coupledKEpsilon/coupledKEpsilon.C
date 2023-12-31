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

#include "coupledKEpsilon.H"
#include "fvBlockMatrix.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledKEpsilon, 0);
addToRunTimeSelectionTable(RASModel, coupledKEpsilon, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledKEpsilon::coupledKEpsilon
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    Cmu_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),
    sigmaEps_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateEpsilon("epsilon", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    ),
    kEpsilon_
    (
        IOobject
        (
            "kEpsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector2("zero", dimless, vector2::zero)
    )
{
    nut_ = Cmu_*sqr(k_)/(epsilon_ + epsilonSmall_);
    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> coupledKEpsilon::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> coupledKEpsilon::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> coupledKEpsilon::divDevReff() const
{
    const volScalarField nuEffective = nuEff();

    return
    (
      - fvm::laplacian(nuEffective, U_)
      - (fvc::grad(U_) & fvc::grad(nuEffective))
    );
}


bool coupledKEpsilon::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void coupledKEpsilon::correct()
{
    // Bound in case of topological change
    // HJ, 22/Aug/2007
    if (mesh_.changing())
    {
        bound(k_, k0_);
        bound(epsilon_, epsilon0_);
    }

    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    // Make coupled matrix
    fvBlockMatrix<vector2> keEqn(kEpsilon_);

    volScalarField G("RASModel::G", nut_*2*magSqr(symm(fvc::grad(U_))));

    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    // Dissipation equation
    {
        fvScalarMatrix epsEqn
        (
            fvm::ddt(epsilon_)
          + fvm::div(phi_, epsilon_)
          + fvm::SuSp(-fvc::div(phi_), epsilon_)
          - fvm::laplacian(DepsilonEff(), epsilon_)
         ==
            C1_*G*epsilon_/k_
          - fvm::Sp(2*C2_*epsilon_/k_, epsilon_)
          + C2_*sqr(epsilon_)/k_
        );

        epsEqn.relax();
        epsEqn.completeAssembly();

        keEqn.insertEquation(1, epsEqn);

        // Coupling term
        volScalarField coupling
        (
            "coupling",
            -C2_*sqr(epsilon_/k_)
        );
        scalarField& couplingIn = coupling.internalField();

        // Eliminate coupling in wall function cells
        labelList eliminated = epsEqn.eliminatedEqns().toc();

        forAll (eliminated, cellI)
        {
            couplingIn[eliminated[cellI]] *= 0;
        }

        // Insert coupling
        keEqn.insertEquationCoupling(1, 0, coupling);
    }

    // Turbulent kinetic energy equation
    {
        dimensionedScalar nutSmall("nutSmall", nut_.dimensions(), SMALL);

        fvScalarMatrix kEqn
        (
            fvm::ddt(k_)
          + fvm::div(phi_, k_)
          + fvm::SuSp(-fvc::div(phi_), k_)
          - fvm::laplacian(DkEff(), k_)
         ==
            G
          + Cmu_*sqr(k_)/(nut_+nutSmall)
          - fvm::Sp(2*Cmu_*k_/(nut_+nutSmall), k_)
        );

        kEqn.relax();

        keEqn.insertEquation(0, kEqn);
    }

    // Update source coupling: coupling terms eliminated from source
    keEqn.updateSourceCoupling();

    // Solve the block matrix
    keEqn.solve();

    // Retrieve solution
    keEqn.retrieveSolution(0, k_.internalField());
    keEqn.retrieveSolution(1, epsilon_.internalField());

    bound(epsilon_, epsilon0_);
    bound(k_, k0_);

    k_.correctBoundaryConditions();
    epsilon_.correctBoundaryConditions();

    // Re-calculate viscosity
    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
