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

Application
    rhoCentralFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedRhoFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "readThermophysicalProperties.H"
#   include "createTimeControls.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero",dimVolume/dimTime, 0.0);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // --- upwind interpolation of primitive fields on faces

        surfaceScalarField rho_pos =
            fvc::interpolate(rho, pos, "reconstruct(rho)");
        surfaceScalarField rho_neg =
            fvc::interpolate(rho, neg, "reconstruct(rho)");

        surfaceVectorField rhoU_pos =
            fvc::interpolate(rhoU, pos, "reconstruct(U)");
        surfaceVectorField rhoU_neg =
            fvc::interpolate(rhoU, neg, "reconstruct(U)");

        volScalarField rPsi = 1.0/psi;
        surfaceScalarField rPsi_pos =
            fvc::interpolate(rPsi, pos, "reconstruct(T)");
        surfaceScalarField rPsi_neg =
            fvc::interpolate(rPsi, neg, "reconstruct(T)");

        surfaceScalarField e_pos =
            fvc::interpolate(e, pos, "reconstruct(T)");
        surfaceScalarField e_neg =
            fvc::interpolate(e, neg, "reconstruct(T)");

        surfaceVectorField U_pos = rhoU_pos/rho_pos;
        surfaceVectorField U_neg = rhoU_neg/rho_neg;

        surfaceScalarField p_pos = rho_pos*rPsi_pos;
        surfaceScalarField p_neg = rho_neg*rPsi_neg;

        surfaceScalarField phiv_pos = U_pos & mesh.Sf();
        surfaceScalarField phiv_neg = U_neg & mesh.Sf();

        volScalarField c = sqrt(thermo.Cp()/thermo.Cv()*rPsi);
        surfaceScalarField cSf_pos = fvc::interpolate(c, pos, "reconstruct(T)")*mesh.magSf();
        surfaceScalarField cSf_neg = fvc::interpolate(c, neg, "reconstruct(T)")*mesh.magSf();

        surfaceScalarField ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);
        surfaceScalarField am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);

        surfaceScalarField a_pos = ap/(ap - am);

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf = am*a_pos;

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg = (1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos = phiv_pos - aSf;
        surfaceScalarField aphiv_neg = phiv_neg + aSf;

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        #include "compressibleCourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        surfaceScalarField phi("phi", aphiv_pos*rho_pos + aphiv_neg*rho_neg);

        surfaceVectorField phiUp =
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
          + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf();

        surfaceScalarField phiEp =
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg;

        volTensorField tauMC("tauMC", mu*dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

        U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField();
        U.correctBoundaryConditions();
        rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

        volScalarField rhoBydt(rho/runTime.deltaT());

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(mu, U)
              - fvc::div(tauMC)
            );
            rhoU = rho*U;
        }

        // --- Solve energy
        surfaceScalarField sigmaDotU =
        (
            (
                fvc::interpolate(mu)*mesh.magSf()*fvc::snGrad(U)
              + (mesh.Sf() & fvc::interpolate(tauMC))
            )
            & (a_pos*U_pos + a_neg*U_neg)
        );

        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
        );

        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        thermo.correct();
        rhoE.boundaryField() =
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );

        if (!inviscid)
        {
            volScalarField k("k", thermo.Cp()*mu/Pr);
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e)
              - fvm::laplacian(thermo.alpha(), e)
              + fvc::laplacian(thermo.alpha(), e)
              - fvc::laplacian(k, T)
            );
            thermo.correct();
            rhoE = rho*(e + 0.5*magSqr(U));
        }

        p.dimensionedInternalField() =
            rho.dimensionedInternalField()
           /psi.dimensionedInternalField();
        p.correctBoundaryConditions();
        rho.boundaryField() = psi.boundaryField()*p.boundaryField();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
