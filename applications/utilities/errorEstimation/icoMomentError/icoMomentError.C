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
    icoMomentError

Description
    Estimates error for the incompressible laminar CFD application icoFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "linear.H"
#include "gaussConvectionScheme.H"
#include "gaussLaplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    Info<< "\nEstimating error in the incompressible momentum equation\n"
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

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        IOobject pHeader
        (
            "p",
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

        if (pHeader.headerOk() && Uheader.headerOk())
        {
            Info << "Reading p" << endl;
            volScalarField p(pHeader, mesh);

            Info << "Reading U" << endl;
            volVectorField U(Uheader, mesh);

#           include "createPhi.H"

            volScalarField ek = 0.5*magSqr(U);
            volTensorField gradU = fvc::grad(U);

            // Divergence of the error in U squared

            volScalarField L
            (
                IOobject
                (
                    "L",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("one", dimLength, 1.0)
            );

            L.internalField() =
                mesh.V()/fvc::surfaceSum(mesh.magSf())().internalField();

            // Warning: 4th row of this equation specially modified
            // for the momentum equation. The "real" formulation would
            // have diffusivity*(gradV && gradV)
            volScalarField momError
            (
                IOobject
                (
                    "momErrorL" + U.name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                sqrt
                (
                    2.0*mag
                    (
                        (
                            fv::gaussConvectionScheme<scalar>
                            (
                                mesh,
                                phi,
                                tmp<surfaceInterpolationScheme<scalar> >
                                (
                                    new linear<scalar>(mesh)
                                )
                            ).fvcDiv(phi, ek)

                          - nu*
                            fv::gaussLaplacianScheme<scalar, scalar>(mesh)
                           .fvcLaplacian
                            (
                                ek
                            )
                          - (U & fvc::grad(p))
//                        + nu*(gradU && gradU)
                          + 0.5*nu*
                            (
                                gradU && (gradU + gradU.T())
                            )
                        )*L/(mag(U) + nu/L)
                    )
                )
            );

            momError.boundaryField() = 0.0;
            momError.write();
        }
        else
        {
            Info<< "    No p or U" << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
