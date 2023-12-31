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
    interFoamPressure

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Description
    Calculate static pressure from interFoam results

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"

    pisoControl piso(mesh);

#   include "readGravitationalAcceleration.H"

    label pRefCell = 0;
    scalar pRefValue = 0.0;

    for (label i = startTime; i < endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject pdHeader
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject gammaHeader
        (
            "gamma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject phiHeader
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check all fields exists
        if
        (
            pdHeader.headerOk()
         && gammaHeader.headerOk()
         && Uheader.headerOk()
         && phiHeader.headerOk()
        )
        {
            mesh.readUpdate();

            Info<< "    Reading pd" << endl;
            volScalarField pd(pdHeader, mesh);

            Info<< "    Reading gamma" << endl;
            volScalarField gamma(gammaHeader, mesh);

            Info<< "    Reading U" << endl;
            volVectorField U(Uheader, mesh);

            Info<< "    Reading phi" << endl;
            surfaceScalarField phi(phiHeader, mesh);

            Info<< "Reading transportProperties\n" << endl;
            twoPhaseMixture twoPhaseProperties(U, phi, "gamma");

            twoPhaseProperties.correct();

            // Construct interface from gamma distribution
            interfaceProperties interface(gamma, U, twoPhaseProperties);

            // Create momentum matrix

            const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
            const dimensionedScalar& rho2 = twoPhaseProperties.rho2();

            volScalarField rho
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT
                ),
                gamma*rho1 + (scalar(1) - gamma)*rho2,
                gamma.boundaryField().types()
            );

            surfaceScalarField rhoPhi
            (
                IOobject
                (
                    "rho*phi",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::interpolate(rho)*phi
            );

            surfaceScalarField muf = twoPhaseProperties.muf();

            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U)
              + fvm::div(rhoPhi, U)
              - fvm::laplacian(muf, U)
              - (fvc::grad(U) & fvc::grad(muf))
             ==
                interface.sigmaK()*fvc::grad(gamma)
              + rho*g
            );

            // Solve for static pressure
            volScalarField p
            (
                IOobject
                (
                    "p",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                pd
            );

            setRefCell(p, piso.dict(), pRefCell, pRefValue);
            mesh.schemesDict().setFluxRequired(p.name());

            volScalarField rUA = 1.0/UEqn.A();
            surfaceScalarField rUAf = fvc::interpolate(rUA);

            U = rUA*UEqn.H();

            phi = fvc::interpolate(U) & mesh.Sf();

            while (piso.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rUAf, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();
            }

            Info << "Writing p" << endl;
            p.write();
        }
        else
        {
            Info << "Not all fields are present.  " << endl;

            if (!pdHeader.headerOk())
            {
                Info << "pd ";
            }

            if (!gammaHeader.headerOk())
            {
                Info << "gamma ";
            }

            if (!Uheader.headerOk())
            {
                Info << "U ";
            }

            if (!phiHeader.headerOk())
            {
                Info << "phi ";
            }

            Info << "missing." << endl;
        }
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
