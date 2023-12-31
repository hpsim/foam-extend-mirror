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
    conjugateHeatFoam

Description
    Transient solver for buoyancy-driven turbulent flow of incompressible
    Newtonian fluids with conjugate heat transfer, complex heat conduction
    and radiation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "coupledFvMatrices.H"
#include "regionCouplePolyPatch.H"
#include "radiationModel.H"
#include "fluidThermalModel.H"
#include "thermalModel.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "pisoControl.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createSolidMesh.H"

    pisoControl piso(mesh);
    simpleControl simpleSolid(solidMesh);

#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "createSolidFields.H"
#   include "initContinuityErrs.H"
#   include "createTimeControls.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readTimeControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        // Detach coupled CHT patches
#       include "detachPatches.H"

#       include "UEqn.H"

        p_rgh.storePrevIter();

        while (piso.correct())
        {
#           include "pEqn.H"
        }

        // Update turbulent quantities
        turbulence->correct();

        radiation->correct();

        // Update thermal conductivity in the fluid
        kappaEff = rho*Cp*(turbulence->nu()/Pr + turbulence->nut()/Prt);

        // Update thermal conductivity in the solid
        solidThermalModel.correct();
        kSolid = solidThermalModel.k();

        rhoCpsolid.oldTime();
        rhoCpsolid = solidThermalModel.rho()*solidThermalModel.C();

        // Attached coupled CHT patches
#       include "attachPatches.H"

        kappaEff.correctBoundaryConditions();
        kSolid.correctBoundaryConditions();

        // Interpolate to the faces and add thermal resistance
        surfaceScalarField kSolidf = fvc::interpolate(kSolid);
        solidThermalModel.modifyResistance(kSolidf);

#       include "solveEnergy.H"

        // Update density according to Boussinesq approximation
        rhok = 1.0 - beta*(T - TRef);
        rhok.correctBoundaryConditions();

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
