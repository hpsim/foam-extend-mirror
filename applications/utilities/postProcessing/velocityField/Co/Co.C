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
    Co

Description
    Calculates and writes the Co number as a surfaceScalarField obtained
    from field phi.

    The -noWrite option just outputs the max values without writing the
    field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    tmp<volScalarField> Co(const surfaceScalarField& Cof)
    {
        const fvMesh& mesh = Cof.mesh();

        tmp<volScalarField> tCo
        (
            new volScalarField
            (
                IOobject
                (
                    "Co",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("0", Cof.dimensions(), 0)
            )
        );

        volScalarField& Co = tCo();

        // Set local references to mesh data
        const unallocLabelList& owner = mesh.owner();
        const unallocLabelList& neighbour = mesh.neighbour();

        forAll(owner, facei)
        {
            label own = owner[facei];
            label nei = neighbour[facei];

            Co[own] = max(Co[own], Cof[facei]);
            Co[nei] = max(Co[nei], Cof[facei]);
        }

        forAll(Co.boundaryField(), patchi)
        {
            Co.boundaryField()[patchi] = Cof.boundaryField()[patchi];
        }

        return tCo;
    }
}


void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject phiHeader
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (phiHeader.headerOk())
    {
        autoPtr<surfaceScalarField> CoPtr;

        Info<< "    Reading phi" << endl;
        surfaceScalarField phi(phiHeader, mesh);
        Info<< "    Calculating Co" << endl;

        if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
        {
            // compressible
            volScalarField rho
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ
                ),
                mesh
            );

            CoPtr.set
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "Cof",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    (
                        mesh.surfaceInterpolation::deltaCoeffs()
                      * (mag(phi)/(fvc::interpolate(rho)*mesh.magSf()))
                      * runTime.deltaT()
                    )
                )
            );
        }
        else if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
        {
            // incompressible
            CoPtr.set
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "Cof",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    (
                        mesh.surfaceInterpolation::deltaCoeffs()
                      * (mag(phi)/mesh.magSf())
                      * runTime.deltaT()
                    )
                )
            );
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Incorrect dimensions of phi: " << phi.dimensions()
                    << abort(FatalError);
        }

        Info<< "Co max : " << max(CoPtr()).value() << endl;

        if (writeResults)
        {
            CoPtr().write();
            Co(CoPtr())().write();
        }
    }
    else
    {
        Info<< "    No phi" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
