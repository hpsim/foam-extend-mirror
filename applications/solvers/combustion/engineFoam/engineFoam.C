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
    engineFoam

Description
    Solver for internal combustion engines.

    Combusting RANS code using the b-Xi two-equation model.
    Xi may be obtained by either the solution of the Xi transport
    equation or from an algebraic exression.  Both approaches are
    based on Gulder's flame speed correlation which has been shown
    to be appropriate by comparison with the results from the
    spectral model.

    Strain effects are encorporated directly into the Xi equation
    but not in the algebraic approximation.  Further work need to be
    done on this issue, particularly regarding the enhanced removal rate
    caused by flame compression.  Analysis using results of the spectral
    model will be required.

    For cases involving very lean Propane flames or other flames which are
    very strain-sensitive, a transport equation for the laminar flame
    speed is present.  This equation is derived using heuristic arguments
    involving the strain time scale and the strain-rate at extinction.
    the transport velocity is the same as that for the Xi equation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "engineTime.H"
#include "engineMesh.H"
#include "hhuCombustionThermo.H"
#include "turbulenceModel.H"
#include "laminarFlameSpeed.H"
#include "ignition.H"
#include "Switch.H"
#include "OFstream.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createEngineTime.H"
#   include "createEngineMesh.H"

    pisoControl piso(mesh);

#   include "readCombustionProperties.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "createTimeControls.H"
#   include "readEngineTimeControls.H"
#   include "compressibleCourantNo.H"
#   include "setInitialDeltaT.H"
#   include "startSummary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readEngineTimeControls.H"
#       include "compressibleCourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Crank angle = " << runTime.theta() << " CA-deg" << endl;

        mesh.move();

#       include "rhoEqn.H"

#       include "UEqn.H"

        // --- PISO loop
        while (piso.correct())
        {
#           include "ftEqn.H"
#           include "bEqn.H"
#           include "huEqn.H"
#           include "hEqn.H"

            if (!ign.ignited())
            {
                hu == h;
            }

#           include "pEqn.H"
        }

        turbulence->correct();

#       include "logSummary.H"

        rho = thermo.rho();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
