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
    shallowWaterFoam

Description
    Transient solver for inviscid shallow-water equations with rotation.

    If the geometry is 3D then it is assumed to be one layers of cells and the
    component of the velocity normal to gravity is removed.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "readGravitationalAcceleration.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "\n Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            surfaceScalarField phiv("phiv", phi/fvc::interpolate(h));

            fvVectorMatrix hUEqn
            (
                fvm::ddt(hU)
              + fvm::div(phiv, hU)
            );

            hUEqn.relax();

            if (pimple.momentumPredictor())
            {
                if (rotating)
                {
                    solve(hUEqn + (F ^ hU) == -magg*h*fvc::grad(h + h0));
                }
                else
                {
                    solve(hUEqn == -magg*h*fvc::grad(h + h0));
                }

                // Constrain the momentum to be in the geometry if 3D geometry
                if (mesh.nGeometricD() == 3)
                {
                    hU -= (gHat & hU)*gHat;
                    hU.correctBoundaryConditions();
                }
            }

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                surfaceScalarField hf = fvc::interpolate(h);
                volScalarField rUA = 1.0/hUEqn.A();
                surfaceScalarField ghrUAf = magg*fvc::interpolate(h*rUA);

                surfaceScalarField phih0 = ghrUAf*mesh.magSf()*fvc::snGrad(h0);

                if (rotating)
                {
                    hU = rUA*(hUEqn.H() - (F ^ hU));
                }
                else
                {
                    hU = rUA*hUEqn.H();
                }

                phi = (fvc::interpolate(hU) & mesh.Sf())
                    + fvc::ddtPhiCorr(rUA, h, hU, phi)
                    - phih0;

                while (pimple.correctNonOrthogonal())
                {
                    fvScalarMatrix hEqn
                    (
                        fvm::ddt(h)
                      + fvc::div(phi)
                      - fvm::laplacian(ghrUAf, h)
                    );

                    hEqn.solve
                    (
                        mesh.solutionDict().solver
                        (
                            h.select(pimple.finalInnerIter())
                        )
                    );

                    if (pimple.finalNonOrthogonalIter())
                    {
                        phi += hEqn.flux();
                    }
                }

                hU -= rUA*h*magg*fvc::grad(h + h0);

                // Constrain the momentum to be in the geometry if 3D geometry
                if (mesh.nGeometricD() == 3)
                {
                    hU -= (gHat & hU)*gHat;
                }

                hU.correctBoundaryConditions();
            }
        }

        U == hU/h;
        hTotal == h + h0;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
