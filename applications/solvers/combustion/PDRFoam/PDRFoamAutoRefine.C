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
    PDRFoam

Description
    Compressible premixed/partially-premixed combustion solver with turbulence
    modelling.

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

    For large flames e.g. explosions additional modelling for the flame
    wrinkling due to surface instabilities may be applied.

    PDR (porosity/distributed resistance) modelling is included to handle
    regions containing blockages which cannot be resolved by the mesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "hhuCombustionThermo.H"
#include "RASModel.H"
#include "laminarFlameSpeed.H"
#include "XiModel.H"
#include "PDRDragModel.H"
#include "ignition.H"
#include "Switch.H"
#include "bound.H"
#include "dynamicRefineFvMesh.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createDynamicFvMesh.H"

    pisoControl piso(mesh);

#   include "readCombustionProperties.H"
#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "createTimeControls.H"
#   include "setInitialDeltaT.H"

    scalar StCoNum = 0.0;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTimeControls.H"
#       include "CourantNo.H"

#       include "setDeltaT.H"

        runTime++;

        Info<< "\n\nTime = " << runTime.timeName() << endl;

        // Indicators for refinement. Note: before runTime++
        // only for postprocessing reasons.
        tmp<volScalarField> tmagGradP = mag(fvc::grad(p));
        volScalarField normalisedGradP
        (
            "normalisedGradP",
            tmagGradP()/max(tmagGradP())
        );
        normalisedGradP.writeOpt() = IOobject::AUTO_WRITE;
        tmagGradP.clear();

        bool meshChanged = false;
        {
            // Make the fluxes absolute
            fvc::makeAbsolute(phi, rho, U);

            // Test : disable refinement for some cells
            PackedBoolList& protectedCell =
                refCast<dynamicRefineFvMesh>(mesh).protectedCell();

            if (protectedCell.empty())
            {
                protectedCell.setSize(mesh.nCells());
                protectedCell = 0;
            }

            forAll(betav, cellI)
            {
                if (betav[cellI] < 0.99)
                {
                    protectedCell[cellI] = 1;
                }
            }

            //volScalarField pIndicator("pIndicator",
            //    p*(fvc::laplacian(p))
            //  / (
            //        magSqr(fvc::grad(p))
            //      + dimensionedScalar
            //        (
            //            "smallish",
            //            sqr(p.dimensions()/dimLength),
            //            1E-6
            //        )
            //    ));
            //pIndicator.writeOpt() = IOobject::AUTO_WRITE;

            // Flux estimate for introduced faces.
            volVectorField rhoU("rhoU", rho*U);

            // Do any mesh changes
            meshChanged = mesh.update();

//        if (mesh.moving() || meshChanged)
//        {
//#           include "correctPhi.H"
//        }

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, rho, U);
        }


#       include "rhoEqn.H"
#       include "UEqn.H"

        // --- PISO loop
        while (piso.correct())
        {
#           include "bEqn.H"
#           include "ftEqn.H"
#           include "huEqn.H"
#           include "hEqn.H"

            if (!ign.ignited())
            {
                hu == h;
            }

#           include "pEqn.H"
        }

        turbulence->correct();

        runTime.write();

        Info<< "\nExecutionTime = "
             << runTime.elapsedCpuTime()
             << " s\n" << endl;
    }

    Info<< "\n end\n";

    return 0;
}


// ************************************************************************* //
