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

#include "LaunderSharmaKE.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LaunderSharmaKE, 0);
addToRunTimeSelectionTable(RASModel, LaunderSharmaKE, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> LaunderSharmaKE::fMu() const
{
    return exp(-3.4/sqr(scalar(1) + rho_*sqr(k_)/(mu()*epsilon_)/50.0));
}


tmp<volScalarField> LaunderSharmaKE::f2() const
{
    return
        scalar(1)
      - 0.3*exp(-min(sqr(rho_*sqr(k_)/(mu()*epsilon_)), scalar(50.0)));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LaunderSharmaKE::LaunderSharmaKE
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermophysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, rho, U, phi, thermophysicalModel, turbulenceModelName),

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
    C3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C3",
            coeffDict_,
            -0.33
        )
    ),
    sigmak_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmak",
            coeffDict_,
            1.0
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
    Prt_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Prt",
            coeffDict_,
            1.0
        )
    ),

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

    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    mut_
    (
        IOobject
        (
            "mut",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateLowReMut("mut", mesh_)
    ),

    alphat_
    (
        IOobject
        (
            "alphat",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateAlphat("alphat", mesh_)
    )
{
    mut_ = rho_*Cmu_*fMu()*sqr(k_)/(epsilon_ + epsilonSmall_);
    mut_ = min(mut_, muRatio()*mu());
    mut_.correctBoundaryConditions();

    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> LaunderSharmaKE::R() const
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
            ((2.0/3.0)*I)*k_ - (mut_/rho_)*dev(twoSymm(fvc::grad(U_))),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> LaunderSharmaKE::devRhoReff() const
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
           -muEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> LaunderSharmaKE::divDevRhoReff() const
{
    return
    (
      - fvm::laplacian(muEff(), U_)
      - fvc::div(muEff()*dev2(T(fvc::grad(U_))))
    );
}


bool LaunderSharmaKE::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        C3_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        Prt_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void LaunderSharmaKE::correct()
{
    // Bound in case of topological change
    // HJ, 22/Aug/2007
    if (mesh_.changing())
    {
        bound(k_, k0_);
        bound(epsilon_, epsilon0_);
    }

    if (!turbulence_)
    {
        // Re-calculate viscosity
        mut_ == rho_*Cmu_*fMu()*sqr(k_)/(epsilon_ + epsilonSmall_);
        mut_ == min(mut_, muRatio()*mu());

        // Re-calculate thermal diffusivity
        alphat_ = mut_/Prt_;
        alphat_.correctBoundaryConditions();

        return;
    }

    RASModel::correct();

    // Calculate parameters and coefficients for Launder-Sharma low-Reynolds
    // number model

    volScalarField E = 2.0*mu()*mut_*fvc::magSqrGradGrad(U_)/rho_;
    volScalarField D = 2.0*mu()*magSqr(fvc::grad(sqrt(k_)))/rho_;

    volScalarField divU = fvc::div(phi_/fvc::interpolate(rho_));

    if (mesh_.moving())
    {
        divU += fvc::div(mesh_.phi());
    }

    // Change return type due to gradient cacheing.  HJ, 22/Apr/2016
    const tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();

    volScalarField G("RASModel::G", mut_*(gradU && dev(twoSymm(gradU))));

    // Dissipation equation

    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(rho_, epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        C1_*G*epsilon_/k_ + fvm::SuSp((C3_ - 2.0/3.0*C1_)*rho_*divU, epsilon_)
      - fvm::Sp(C2_*f2()*rho_*epsilon_/k_, epsilon_)
    //+ 0.75*1.5*flameKproduction*epsilon_/k_
      + E
    );

    epsEqn().relax();
    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(rho_, k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G - fvm::SuSp(2.0/3.0*rho_*divU, k_)
      - fvm::Sp(rho_*(epsilon_ + D)/k_, k_)
    //+ flameKproduction
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
    mut_ == Cmu_*fMu()*rho_*sqr(k_)/(epsilon_ + epsilonSmall_);
    mut_ == min(mut_, muRatio()*mu());


    // Re-calculate thermal diffusivity
    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
