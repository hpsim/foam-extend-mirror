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
    estimateScalarError

Description
    Estimates the error in the solution for a scalar transport equation in the
    standard form

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "errorEstimate.H"
#include "resError.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    Info<< "\nEstimating error in scalar transport equation\n"
        << "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    Info<< "Reading diffusivity DT\n" << endl;

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        IOobject THeader
        (
            "T",
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

        if (THeader.headerOk() && Uheader.headerOk())
        {
            Info << "Reading T" << endl;
            volScalarField T(THeader, mesh);

            Info << "Reading U" << endl;
            volVectorField U(Uheader, mesh);

#           include "createPhi.H"

            errorEstimate<scalar> ee
            (
                resError::div(phi, T)
              - resError::laplacian(DT, T)
            );

            ee.residual()().write();
            volScalarField e = ee.error();
            e.write();
            mag(e)().write();
        }
        else
        {
            Info<< "    No T or U" << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
