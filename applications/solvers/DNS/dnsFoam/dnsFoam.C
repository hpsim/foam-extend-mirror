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
    dnsFoam

Description
    Direct numerical simulation solver for boxes of isotropic turbulence

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Kmesh.H"
#include "UOprocess.H"
#include "fft.H"
#include "calcEk.H"
#include "graph.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMeshNoClear.H"

    pisoControl piso(mesh);

    #include "readTransportProperties.H"
    #include "createFields.H"
    #include "readTurbulenceProperties.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        force.internalField() = ReImSum
        (
            fft::reverseTransform
            (
                K/(mag(K) + 1.0e-6) ^ forceGen.newField(), K.nn()
            )
        );

        #include "globalProperties.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
         ==
            force
        );

        solve(UEqn == -fvc::grad(p));


        // --- PISO loop

        while (piso.correct())
        {
            volScalarField rUA = 1.0/UEqn.A();

            U = rUA*UEqn.H();
            phi = (fvc::interpolate(U) & mesh.Sf())
                + fvc::ddtPhiCorr(rUA, U, phi);

            fvScalarMatrix pEqn
            (
                fvm::laplacian(rUA, p) == fvc::div(phi)
            );

            pEqn.solve();

            phi -= pEqn.flux();

#           include "continuityErrs.H"

            U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        runTime.write();

        if (runTime.outputTime())
        {
            calcEk(U, K).write(runTime.timePath()/"Ek", runTime.graphFormat());
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
