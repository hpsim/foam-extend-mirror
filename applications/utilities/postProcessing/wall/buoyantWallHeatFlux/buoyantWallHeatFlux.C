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
    buoyantWallHeatFlux

Description
    Calculates and writes the heat flux in incompressible flow with Boussinesq's
    buoyancy assumption (e.g. buoyantBoussineqSimpleFoam) for all patches as the
    boundary field of a volScalarField and also prints the integrated flux for
    all wall patches.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "RASModel.H"
#include "LESModel.H"
#include "singlePhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        #include "createFields.H"

        // Calculate effective kappa
        const volScalarField kappaEff
        (
            "kappaEff",
            turbulence->nu()/Pr + turbulence->nut()/Prt
        );

        // Calculate the heat flux
        const surfaceScalarField heatFlux =
            rhoRef*c*fvc::interpolate(kappaEff)*fvc::snGrad(T);

        // Get the heat flux at the patch
        const surfaceScalarField::GeometricBoundaryField& patchHeatFlux =
            heatFlux.boundaryField();

        Info<< "\nWall heat fluxes [W]" << endl;
        forAll(patchHeatFlux, patchi)
        {
            if (mesh.boundary()[patchi].isWall())
            {
                Info<< mesh.boundary()[patchi].name()
                    << " "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux[patchi]
                       )
                    << endl;
            }
        }
        Info<< endl;

        // Create the volScalarField which will have heat fluxes at the boundary
        volScalarField buoyantWallHeatFlux
        (
            IOobject
            (
                "buoyantWallHeatFlux",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("buoyantWallHeatFlux", heatFlux.dimensions(), 0.0)
        );

        forAll(buoyantWallHeatFlux.boundaryField(), patchi)
        {
            buoyantWallHeatFlux.boundaryField()[patchi] = patchHeatFlux[patchi];
        }

        buoyantWallHeatFlux.write();
    }

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
