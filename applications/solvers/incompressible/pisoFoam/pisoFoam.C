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
    pisoFoam

Description
    Transient solver for incompressible, turbulent flow.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    Consistent formulation without time-step and relaxation dependence by Jasak
    and Tukovic.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    pisoControl piso(mesh);

#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            // Momentum predictor

            // Time-derivative matrix
            fvVectorMatrix ddtUEqn(fvm::ddt(U));

            // Convection-diffusion matrix
            fvVectorMatrix HUEqn
            (
                fvm::div(phi, U)
              + turbulence->divDevReff()
            );

            if (piso.momentumPredictor())
            {
                solve(relax(ddtUEqn + HUEqn) == -fvc::grad(p));
            }

            // --- PISO loop

            while (piso.correct())
            {
                // Prepare clean 1/a_p without time derivative and
                // under-relaxation contribution
                volScalarField rAU = 1.0/HUEqn.A();

                // Calculate U from convection-diffusion matrix
                U = rAU*HUEqn.H();

                // Consistently calculate flux
                piso.calcTransientConsistentFlux(phi, U, rAU, ddtUEqn);

                // Global flux balance
                adjustPhi(phi, U, p);

                // Non-orthogonal pressure corrector loop
                while (piso.correctNonOrthogonal())
                {
                    // Pressure corrector

                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian
                        (
                            fvc::interpolate(rAU)/piso.aCoeff(U.name()),
                            p,
                            "laplacian(rAU," + p.name() + ')'
                        )
                     ==
                        fvc::div(phi)
                    );

                    pEqn.setReference(pRefCell, pRefValue);
                    pEqn.solve
                    (
                        mesh.solutionDict().solver
                        (
                            p.select(piso.finalInnerIter())
                        )
                    );

                    if (piso.finalNonOrthogonalIter())
                    {
                        phi -= pEqn.flux();
                    }
                }

#               include "continuityErrs.H"

                // Consistently reconstruct velocity after pressure
                // equation. Note: flux is made relative inside the function
                piso.reconstructTransientVelocity(U, phi, ddtUEqn, rAU, p);
            }
        }

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
