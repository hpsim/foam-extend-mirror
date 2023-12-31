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
    twoPhaseEulerFoam

Description
    Solver for a system of 2 incompressible fluid phases with one phase
    dispersed, e.g. gas bubbles in a liquid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "nearWallDist.H"
#include "wallFvPatch.H"
#include "Switch.H"

#include "IFstream.H"
#include "OFstream.H"

#include "dragModel.H"
#include "phaseModel.H"
#include "kineticTheoryModel.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    pisoControl piso(mesh);

#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "readPPProperties.H"
#   include "initContinuityErrs.H"
#   include "createTimeControls.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTwoPhaseEulerFoamControls.H"
#       include "CourantNos.H"
#       include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "alphaEqn.H"

#       include "liftDragCoeffs.H"

#       include "UEqns.H"

        // --- PISO loop
        while (piso.correct())
        {
#           include "pEqn.H"

            if (correctAlpha && !piso.finalInnerIter())
            {
#               include "alphaEqn.H"
            }
        }

#       include "DDtU.H"

#       include "kEpsilon.H"

        if (kineticTheory.on())
        {
            kineticTheory.solve();
            nuEffa += kineticTheory.mua()/rhoa;
        }
#       include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
