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
    porousExplicitSourceReactingParcelFoam

Description
    Transient PISO solver for compressible, laminar or turbulent flow with
    reacting multiphase Lagrangian parcels for porous media, including explicit
    sources for mass, momentum and energy

    The solver includes:
    - reacting multiphase parcel cloud
    - porous media
    - mass, momentum and energy sources
    - polynomial based, incompressible thermodynamics (f(T))

    Note: ddtPhiCorr not used here when porous zones are active
    - not well defined for porous calculations

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hReactionThermo.H"
#include "turbulenceModel.H"
#include "BasicReactingMultiphaseCloud.H"
#include "rhoChemistryModel.H"
#include "chemistrySolver.H"
#include "radiationModel.H"
#include "porousZones.H"
#include "timeActivatedExplicitSource.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "readChemistryProperties.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createRadiationModel.H"
    #include "createClouds.H"
    #include "createExplicitSources.H"
    #include "createPorousZones.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readAdditionalSolutionControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();

        #include "chemistry.H"
        #include "rhoEqn.H"
        #include "UEqn.H"
        #include "YEqn.H"
        #include "hsEqn.H"

        // --- PISO loop
        while (piso.correct())
        {
            #include "pEqn.H"
        }

        turbulence->correct();

        rho = thermo.rho();

        if (runTime.write())
        {
            chemistry.dQ()().write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
