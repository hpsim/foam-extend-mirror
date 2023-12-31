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

\*---------------------------------------------------------------------------*/

#include "heatFlux.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "chtRcTemperatureFvPatchScalarField.H"
#include "processorFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::heatFlux, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatFlux::heatFlux
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    Kfluid_(dict.lookup("K")),
    obr_(obr),
    active_(true)
{
    // Only active if a fvMesh is available
    if (isA<fvMesh>(obr_))
    {
        read(dict);
    }
    else
    {
        active_ = false;
        WarningIn
        (
            "heatFlux::heatFlux\n"
            "(\n"
                "const word&,\n"
                "const objectRegistry&,\n"
                "const dictionary&,\n"
                "const bool\n"
            ")"
        )   << "No fvMesh available, deactivating."
            << nl << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heatFlux::calcAndPrint()
{
    const volScalarField& T =
        obr_.lookupObject<volScalarField>("T");

    const volScalarField& kappaEff =
        obr_.lookupObject<volScalarField>(Kfluid_);

    //const surfaceScalarField& kappaEff =
    //    obr_.lookupObject<surfaceScalarField>(Kfluid_);

    scalar rho = 1.2;
    scalar Cp = 1000;

    const fvMesh& mesh = T.mesh();

    surfaceScalarField heatFluxD =
        -fvc::interpolate(kappaEff)*fvc::snGrad(T);
        // -kappaEff*fvc::snGrad(T);

    const surfaceScalarField::GeometricBoundaryField& patchHeatFluxD =
        heatFluxD.boundaryField();
    //const surfaceScalarField::GeometricBoundaryField& patchHeatFluxD =
    //    heatFluxC.boundaryField();

    //surfaceScalarField::GeometricBoundaryField patchHeatFluxC =

    scalar sumConduction = 0.0;
    scalar sumConvection = 0.0;
    scalar sumRadiation = 0.0;

    Info<< "\nWall heat fluxes [W]" << endl;
    forAll(patchHeatFluxD, patchI)
    {
        if (isA<processorFvPatchScalarField>(T.boundaryField()[patchI]))
        {
            continue;
        }

        scalar conduction = gSum
        (
            mesh.magSf().boundaryField()[patchI]
           *heatFluxD.boundaryField()[patchI]
        );

        // Account for heat sources at region couple BCs
        if (isA<chtRcTemperatureFvPatchScalarField>(T.boundaryField()[patchI]))
        {
            const chtRcTemperatureFvPatchScalarField& pT =
                dynamic_cast<const chtRcTemperatureFvPatchScalarField&>
                (
                    T.boundaryField()[patchI]
                );

            conduction -= gSum
            (
                pT.source()*mesh.magSf().boundaryField()[patchI]
            );
        }

        scalar convection = 0.0;
        scalar radiation = 0.0;

        if (obr_.foundObject<surfaceScalarField>("phi"))
        {
            const surfaceScalarField& phi =
                obr_.lookupObject<surfaceScalarField>("phi");

            convection = gSum
            (
                 rho*Cp*T.boundaryField()[patchI]*
                 phi.boundaryField()[patchI]
            );
        }

        if (obr_.foundObject<volScalarField>("Qr"))
        {
            const volScalarField& Qr =
                obr_.lookupObject<volScalarField>("Qr");

            radiation = gSum
            (
                 Qr.boundaryField()[patchI]*
                 mesh.magSf().boundaryField()[patchI]
            );
        }

        Info<< mesh.boundary()[patchI].name()
            << " conduction = "
            << conduction
            << " convection = "
            << convection
            << " radiation = "
            << radiation
            << " total = "
            << conduction + convection + radiation
            << endl;

        sumConduction += conduction;
        sumConvection += convection;
        sumRadiation += radiation;
    }

    Info<< "sum conduction ="
        << sumConduction
        << " sum convection = "
        << sumConvection
        << " sum radiation = "
        << sumRadiation
        << " sum total = "
        << sumConduction + sumConvection + sumRadiation
        << nl << endl;
}

void Foam::heatFlux::read(const dictionary& dict)
{
    if (active_)
    {
    }
}


void Foam::heatFlux::execute()
{
    if (active_)
    {
        calcAndPrint();
    }
}


void Foam::heatFlux::end()
{
    if (active_)
    {
        calcAndPrint();
    }
}


void Foam::heatFlux::timeSet()
{
    // Do nothing - only valid on write
}


void Foam::heatFlux::write()
{
    if (active_)
    {
    }
}


void Foam::heatFlux::updateMesh(const mapPolyMesh&)
{
    // Do nothing
}


void Foam::heatFlux::movePoints(const pointField&)
{
    // Do nothing
}


// ************************************************************************* //
