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

#include "kOmegaSSTSAS.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kOmegaSSTSAS, 0);
addToRunTimeSelectionTable(LESModel, kOmegaSSTSAS, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void kOmegaSSTSAS::updateSubGridScaleFields(const volScalarField& S2)
{
    nuSgs_ == a1_*k_/max(a1_*omega_, F2()*sqrt(S2));
    nuSgs_.correctBoundaryConditions();
}


tmp<volScalarField> kOmegaSSTSAS::F1(const volScalarField& CDkOmega) const
{
    volScalarField CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    volScalarField arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*nu()/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


tmp<volScalarField> kOmegaSSTSAS::F2() const
{
    volScalarField arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*nu()/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


tmp<volScalarField> kOmegaSSTSAS::Lvk2
(
    const volScalarField& S2
) const
{
    return max
    (
        kappa_*sqrt(S2)
       /(
            mag(fvc::laplacian(U()))
          + dimensionedScalar
            (
                "ROOTVSMALL",
                dimensionSet(0, -1 , -1, 0, 0, 0, 0),
                ROOTVSMALL
            )
        ),
        Cs_*delta()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOmegaSSTSAS::kOmegaSSTSAS
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport),

    alphaK1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaK1",
            coeffDict_,
            0.85034
        )
    ),
    alphaK2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaK2",
            coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaOmega1",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaOmega2",
            coeffDict_,
            0.85616
        )
    ),
    gamma1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "gamma1",
            coeffDict_,
            0.5532
        )
    ),
    gamma2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "gamma2",
            coeffDict_,
            0.4403
        )
    ),
    beta1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),
    c1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
        )
    ),
    Cs_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cs",
            coeffDict_,
            0.262
        )
    ),
    alphaPhi_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaPhi",
            coeffDict_,
            0.666667
        )
    ),
    zetaTilda2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "zetaTilda2",
            coeffDict_,
            1.755
        )
    ),
    FSAS_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "FSAS",
            coeffDict_,
            1.25
        )
    ),

    omega0_("omega0", dimless/dimTime, SMALL),
    omegaSmall_("omegaSmall", dimless/dimTime, SMALL),
    y_(mesh_),
    Cmu_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
         )
    ),
    kappa_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "kappa",
            *this,
            0.41
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

    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    nuSgs_
    (
        IOobject
        (
            "nuSgs",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    updateSubGridScaleFields(magSqr(symm(fvc::grad(U))));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kOmegaSSTSAS::correct(const tmp<volTensorField>& gradU)
{
    LESModel::correct(gradU);

    if (mesh_.changing())
    {
        y_.correct();
    }

    volScalarField S2 = magSqr(symm(gradU()));
    gradU.clear();

    volVectorField gradK = fvc::grad(k_);
    volVectorField gradOmega = fvc::grad(omega_);
    volScalarField L = sqrt(k_)/(pow(Cmu_, 0.25)*(omega_ + omegaSmall_));
    volScalarField CDkOmega =
        (2.0*alphaOmega2_)*(gradK & gradOmega)/(omega_ + omegaSmall_);
    volScalarField F1 = this->F1(CDkOmega);
    volScalarField G = nuSgs_*2.0*S2;

    // Turbulent kinetic energy equation
    {
        fvScalarMatrix kEqn
        (
            fvm::ddt(k_)
          + fvm::div(phi(), k_)
          - fvm::Sp(fvc::div(phi()), k_)
          - fvm::laplacian(DkEff(F1), k_)
        ==
            min(G, c1_*betaStar_*k_*omega_)
          - fvm::Sp(betaStar_*omega_, k_)
        );

        kEqn.relax();
        kEqn.solve();
    }
    bound(k_, k0());

    volScalarField grad_omega_k = max
    (
        magSqr(gradOmega)/
        sqr(omega_ + omegaSmall_),
        magSqr(gradK)/
        sqr(k_ + k0())
    );

    // Turbulent frequency equation
    {
        fvScalarMatrix omegaEqn
        (
            fvm::ddt(omega_)
          + fvm::div(phi(), omega_)
          - fvm::Sp(fvc::div(phi()), omega_)
          - fvm::laplacian(DomegaEff(F1), omega_)
        ==
            gamma(F1)*2.0*S2
          - fvm::Sp(beta(F1)*omega_, omega_)
          - fvm::SuSp       // cross diffusion term
            (
                (F1 - scalar(1))*CDkOmega/omega_,
                omega_
            )
          + FSAS_
           *max
            (
                dimensionedScalar("zero",dimensionSet(0, 0 , -2, 0, 0),0. ),
                zetaTilda2_*kappa_*S2*(L/Lvk2(S2))
              - 2.0/alphaPhi_*k_*grad_omega_k
            )
        );

        omegaEqn.relax();
        omegaEqn.solve();
    }
    bound(omega_, omega0_);

    updateSubGridScaleFields(S2);
}


tmp<volScalarField> kOmegaSSTSAS::epsilon() const
{
    return 2.0*nuEff()*magSqr(symm(fvc::grad(U())));
}


tmp<volSymmTensorField> kOmegaSSTSAS::B() const
{
    return ((2.0/3.0)*I)*k() - nuSgs()*twoSymm(fvc::grad(U()));
}


tmp<volSymmTensorField> kOmegaSSTSAS::devBeff() const
{
    return -nuEff()*dev(twoSymm(fvc::grad(U())));
}


tmp<fvVectorMatrix> kOmegaSSTSAS::divDevBeff() const
{
    return
    (
      - fvm::laplacian(nuEff(), U_) - fvc::div(nuEff()*dev(T(fvc::grad(U_))))
    );
}


bool kOmegaSSTSAS::read()
{
    if (LESModel::read())
    {
        alphaK1_.readIfPresent(coeffDict());
        alphaK2_.readIfPresent(coeffDict());
        alphaOmega1_.readIfPresent(coeffDict());
        alphaOmega2_.readIfPresent(coeffDict());
        gamma1_.readIfPresent(coeffDict());
        gamma2_.readIfPresent(coeffDict());
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        betaStar_.readIfPresent(coeffDict());
        a1_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());
        Cs_.readIfPresent(coeffDict());
        alphaPhi_.readIfPresent(coeffDict());
        zetaTilda2_.readIfPresent(coeffDict());
        FSAS_.readIfPresent(coeffDict());

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
