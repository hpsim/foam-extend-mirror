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
    Pe

Description
    Calculates and writes the Pe number as a surfaceScalarField obtained from
    field phi.

    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "incompressible/LES/LESModel/LESModel.H"
#include "basicPsiThermo.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "compressible/LES/LESModel/LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
        autoPtr<surfaceScalarField> PePtr;

        Info<< "    Reading phi" << endl;
        surfaceScalarField phi(phiHeader, mesh);

        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        IOobject RASPropertiesHeader
        (
            "RASProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        IOobject LESPropertiesHeader
        (
            "LESProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        Info<< "    Calculating Pe" << endl;

        if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
        {
            if (RASPropertiesHeader.headerOk())
            {
                IOdictionary RASProperties(RASPropertiesHeader);

                singlePhaseTransportModel laminarTransport(U, phi);

                autoPtr<incompressible::RASModel> RASModel
                (
                    incompressible::RASModel::New
                    (
                        U,
                        phi,
                        laminarTransport
                    )
                );

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mag(phi)
                       /(
                            mesh.magSf()
                          * mesh.surfaceInterpolation::deltaCoeffs()
                          * fvc::interpolate(RASModel->nuEff())
                         )
                    )
                );
            }
            else if (LESPropertiesHeader.headerOk())
            {
                IOdictionary LESProperties(LESPropertiesHeader);

                singlePhaseTransportModel laminarTransport(U, phi);

                autoPtr<incompressible::LESModel> sgsModel
                (
                    incompressible::LESModel::New(U, phi, laminarTransport)
                );

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mag(phi)
                       /(
                            mesh.magSf()
                          * mesh.surfaceInterpolation::deltaCoeffs()
                          * fvc::interpolate(sgsModel->nuEff())
                        )
                    )
                );
            }
            else
            {
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

                dimensionedScalar nu(transportProperties.lookup("nu"));

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mesh.surfaceInterpolation::deltaCoeffs()
                      * (mag(phi)/mesh.magSf())*(runTime.deltaT()/nu)
                    )
                );
            }
        }
        else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
        {
            if (RASPropertiesHeader.headerOk())
            {
                IOdictionary RASProperties(RASPropertiesHeader);

                autoPtr<basicPsiThermo> thermo(basicPsiThermo::New(mesh));

                volScalarField rho
                (
                    IOobject
                    (
                        "rho",
                        runTime.timeName(),
                        mesh
                    ),
                    thermo->rho()
                );

                autoPtr<compressible::RASModel> RASModel
                (
                    compressible::RASModel::New
                    (
                        rho,
                        U,
                        phi,
                        thermo()
                    )
                );

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mag(phi)
                       /(
                            mesh.magSf()
                          * mesh.surfaceInterpolation::deltaCoeffs()
                          * fvc::interpolate(RASModel->muEff())
                        )
                    )
                );
            }
            else if (LESPropertiesHeader.headerOk())
            {
                IOdictionary LESProperties(LESPropertiesHeader);

                autoPtr<basicPsiThermo> thermo(basicPsiThermo::New(mesh));

                volScalarField rho
                (
                    IOobject
                    (
                        "rho",
                        runTime.timeName(),
                        mesh
                    ),
                    thermo->rho()
                );

                autoPtr<compressible::LESModel> sgsModel
                (
                    compressible::LESModel::New(rho, U, phi, thermo())
                );

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mag(phi)
                       /(
                            mesh.magSf()
                          * mesh.surfaceInterpolation::deltaCoeffs()
                          * fvc::interpolate(sgsModel->muEff())
                        )
                    )
                );
            }
            else
            {
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

                dimensionedScalar mu(transportProperties.lookup("mu"));

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mesh.surfaceInterpolation::deltaCoeffs()
                      * (mag(phi)/(mesh.magSf()))*(runTime.deltaT()/mu)
                    )
                );
            }
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Incorrect dimensions of phi: " << phi.dimensions()
                    << abort(FatalError);
        }


        // can also check how many cells exceed a particular Pe limit
        /*
        {
            label count = 0;
            label PeLimit = 200;
            forAll(PePtr(), i)
            {
                if (PePtr()[i] > PeLimit)
                {
                    count++;
                }

            }

            Info<< "Fraction > " << PeLimit << " = "
                << scalar(count)/Pe.size() << endl;
        }
        */

        Info << "Pe max : " << max(PePtr()).value() << endl;

        if (writeResults)
        {
            PePtr().write();
        }
    }
    else
    {
        Info<< "    No phi" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
