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

#include "LRR.H"
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

defineTypeNameAndDebug(LRR, 0);
addToRunTimeSelectionTable(RASModel, LRR, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LRR::LRR
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
    Clrr1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Clrr1",
            coeffDict_,
            1.8
        )
    ),
    Clrr2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Clrr2",
            coeffDict_,
            0.6
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
    Cs_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cs",
            coeffDict_,
            0.25
        )
    ),
    Ceps_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Ceps",
            coeffDict_,
            0.15
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
    couplingFactor_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "couplingFactor",
            coeffDict_,
            0.0
        )
    ),

    R_
    (
        IOobject
        (
            "R",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateR("R", mesh_)
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
    )
{
    nut_ = Cmu_*sqr(k_)/(epsilon_ + epsilonSmall_);
    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();

    if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
    {
        FatalErrorIn
        (
            "LRR::LRR"
            "(const volVectorField& U, const surfaceScalarField& phi,"
            "transportModel& lamTransportModel)"
        )   << "couplingFactor = " << couplingFactor_
            << " is not in range 0 - 1" << nl
            << exit(FatalError);
    }

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> LRR::devReff() const
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
            R_ - nu()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> LRR::divDevReff() const
{
    if (couplingFactor_.value() > 0.0)
    {
        return
        (
            fvc::div(R_ + couplingFactor_*nut_*fvc::grad(U_), "div(R)")
          + fvc::laplacian
            (
                 (1.0 - couplingFactor_)*nut_,
                 U_,
                 "laplacian(nuEff,U)"
            )
          - fvm::laplacian(nuEff(), U_)
        );
    }
    else
    {
        return
        (
            fvc::div(R_)
          + fvc::laplacian(nut_, U_, "laplacian(nuEff,U)")
          - fvm::laplacian(nuEff(), U_)
        );
    }
}


bool LRR::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        Clrr1_.readIfPresent(coeffDict());
        Clrr2_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        Cs_.readIfPresent(coeffDict());
        Ceps_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());

        couplingFactor_.readIfPresent(coeffDict());

        if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
        {
            FatalErrorIn("LRR::read()")
                << "couplingFactor = " << couplingFactor_
                << " is not in range 0 - 1"
                << exit(FatalError);
        }

        return true;
    }
    else
    {
        return false;
    }
}


void LRR::correct()
{
    // Bound in case of topological change
    // HJ, 22/Aug/2007
    if (mesh_.changing())
    {
        R_.correctBoundaryConditions();
        bound(epsilon_, epsilon0_);
    }

    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    volSymmTensorField P = -twoSymm(R_ & fvc::grad(U_));
    volScalarField G("RASModel::G", 0.5*mag(tr(P)));

    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      + fvm::SuSp(-fvc::div(phi_), epsilon_)
    //- fvm::laplacian(Ceps*(K/epsilon_)*R, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        C1_*G*epsilon_/k_
      - fvm::Sp(C2_*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

    // No longer needed: matrix completes at the point of solution
    // HJ, 17/Apr/2012
//     epsEqn().completeAssembly();

    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Reynolds stress equation

    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (curPatch.isWall())
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                P[faceCelli] *=
                    min
                    (
                        G[faceCelli]/(0.5*mag(tr(P[faceCelli])) + SMALL),
                        1.0
                    );
            }
        }
    }


    tmp<fvSymmTensorMatrix> REqn
    (
        fvm::ddt(R_)
      + fvm::div(phi_, R_)
      + fvm::SuSp(-fvc::div(phi_), R_)
    //- fvm::laplacian(Cs*(k_/epsilon_)*R_, R_)
      - fvm::laplacian(DREff(), R_)
      + fvm::Sp(Clrr1_*epsilon_/k_, R_)
      ==
        P
      - (2.0/3.0*(1 - Clrr1_)*I)*epsilon_
      - Clrr2_*dev(P)
    );

    REqn().relax();
    solve(REqn);

    R_.max
    (
        dimensionedSymmTensor
        (
            "zero",
            R_.dimensions(),
            symmTensor
            (
                k0_.value(), -GREAT, -GREAT,
                k0_.value(), -GREAT,
                k0_.value()
            )
        )
    );

    k_ = 0.5*tr(R_);
    bound(k_, k0_);


    // Re-calculate viscosity
    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();


    // Correct wall shear stresses

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (curPatch.isWall())
        {
            symmTensorField& Rw = R_.boundaryField()[patchi];

            const scalarField& nutw = nut_.boundaryField()[patchi];

            vectorField snGradU = U_.boundaryField()[patchi].snGrad();

            const vectorField& faceAreas
                = mesh_.Sf().boundaryField()[patchi];

            const scalarField& magFaceAreas
                = mesh_.magSf().boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                // Calculate near-wall velocity gradient
                tensor gradUw
                    = (faceAreas[facei]/magFaceAreas[facei])*snGradU[facei];

                // Calculate near-wall shear-stress tensor
                tensor tauw = -nutw[facei]*2*symm(gradUw);

                // Reset the shear components of the stress tensor
                Rw[facei].xy() = tauw.xy();
                Rw[facei].xz() = tauw.xz();
                Rw[facei].yz() = tauw.yz();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
