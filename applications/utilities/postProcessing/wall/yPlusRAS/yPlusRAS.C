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
    yPlusRAS

Description
    Calculates and reports yPlus for all wall patches, for the specified times
    when using RAS turbulence models.

    Default behaviour assumes operating in incompressible mode. To apply to
    compressible RAS cases, use the -compressible option.

    Extended version for being able to handle two phase flows using the
    -twoPhase option.

    Frank Albina, 16/Nov/2009

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/incompressibleTwoPhaseMixture/twoPhaseMixture.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "nutWallFunction/nutWallFunctionFvPatchScalarField.H"

#include "basicPsiThermo.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "mutWallFunction/mutWallFunctionFvPatchScalarField.H"

#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void calcIncompressibleYPlus
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& yPlus
)
{
    typedef incompressible::RASModels::nutWallFunctionFvPatchScalarField
        wallFunctionPatchField;

#   include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> RASModel
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    const volScalarField::GeometricBoundaryField nutPatches =
        RASModel->nut()().boundaryField();

    bool foundNutPatch = false;
    forAll(nutPatches, patchi)
    {
        if (isA<wallFunctionPatchField>(nutPatches[patchi]))
        {
            foundNutPatch = true;

            const wallFunctionPatchField& nutPw =
                dynamic_cast<const wallFunctionPatchField&>
                    (nutPatches[patchi]);

            yPlus.boundaryField()[patchi] = nutPw.yPlus();
            const scalarField& Yp = yPlus.boundaryField()[patchi];

            Info<< "Patch " << patchi
                << " named " << nutPw.patch().name()
                << " y+ : min: " << gMin(Yp) << " max: " << gMax(Yp)
                << " average: " << gAverage(Yp) << nl << endl;
        }
    }

    if (!foundNutPatch)
    {
        Info<< "    no " << wallFunctionPatchField::typeName << " patches"
            << endl;
    }
}


void calcCompressibleYPlus
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& yPlus
)
{
    typedef compressible::RASModels::mutWallFunctionFvPatchScalarField
        wallFunctionPatchField;

    IOobject rhoHeader
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!rhoHeader.headerOk())
    {
        Info<< "    no rho field" << endl;
        return;
    }

    Info << "Reading field rho\n" << endl;
    volScalarField rho(rhoHeader, mesh);

#   include "compressibleCreatePhi.H"

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

    autoPtr<compressible::RASModel> RASModel
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    const volScalarField::GeometricBoundaryField mutPatches =
        RASModel->mut()().boundaryField();

    bool foundMutPatch = false;
    forAll(mutPatches, patchi)
    {
        if (isA<wallFunctionPatchField>(mutPatches[patchi]))
        {
            foundMutPatch = true;

            const wallFunctionPatchField& mutPw =
                dynamic_cast<const wallFunctionPatchField&>
                    (mutPatches[patchi]);

            yPlus.boundaryField()[patchi] = mutPw.yPlus();
            const scalarField& Yp = yPlus.boundaryField()[patchi];

            Info<< "Patch " << patchi
                << " named " << mutPw.patch().name()
                << " y+ : min: " << gMin(Yp) << " max: " << gMax(Yp)
                << " average: " << gAverage(Yp) << nl << endl;
        }
    }

    if (!foundMutPatch)
    {
        Info<< "    no " << wallFunctionPatchField::typeName << " patches"
            << endl;
    }
}


// Calculate two phase Y+
void calcTwoPhaseYPlus
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& yPlus
)
{
    typedef incompressible::RASModels::nutWallFunctionFvPatchScalarField
        wallFunctionPatchField;

#   include "createPhi.H"

    IOobject alphaHeader
    (
        "alpha1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!alphaHeader.headerOk())
    {
        Info<< "    no alpha1 field" << endl;
        return;
    }

    Info << "Reading field alpha1\n" << endl;
    volScalarField alpha(alphaHeader, mesh);

    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi, "alpha1");

    autoPtr<incompressible::RASModel> RASModel
    (
        incompressible::RASModel::New(U, phi, twoPhaseProperties)
    );

    const volScalarField::GeometricBoundaryField nutPatches =
        RASModel->nut()().boundaryField();

    bool foundNutPatch = false;
    forAll(nutPatches, patchi)
    {
        if (isA<wallFunctionPatchField>(nutPatches[patchi]))
        {
            foundNutPatch = true;

            const wallFunctionPatchField& nutPw =
                dynamic_cast<const wallFunctionPatchField&>
                    (nutPatches[patchi]);

            yPlus.boundaryField()[patchi] = nutPw.yPlus();
            const scalarField& Yp = yPlus.boundaryField()[patchi];

            Info<< "Patch " << patchi
                << " named " << nutPw.patch().name()
                << " y+ : min: " << gMin(Yp) << " max: " << gMax(Yp)
                << " average: " << gAverage(Yp) << nl << endl;
        }
    }

    if (!foundNutPatch)
    {
        Info<< "    no " << wallFunctionPatchField::typeName << " patches"
            << endl;
    }
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    argList::validOptions.insert("compressible","");
    argList::validOptions.insert("twoPhase","");

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    bool compressible = args.optionFound("compressible");

    // Check if two phase model was selected
    bool twoPhase = args.optionFound("twoPhase");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        fvMesh::readUpdateState state = mesh.readUpdate();

        // Wall distance
        if (timeI == 0 || state != fvMesh::UNCHANGED)
        {
            Info<< "Calculating wall distance\n" << endl;
            wallDist y(mesh, true);
            Info<< "Writing wall distance to field "
                << y.name() << nl << endl;
            y.write();
        }

        volScalarField yPlus
        (
            IOobject
            (
                "yPlus",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlus", dimless, 0.0)
        );

        Info << "Reading field U\n" << endl;
        IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (UHeader.headerOk())
        {
            Info << "Reading field U\n" << endl;
            volVectorField U(UHeader, mesh);

            if (compressible)
            {
                calcCompressibleYPlus(mesh, runTime, U, yPlus);
            }
            else if (twoPhase)
            {
                calcTwoPhaseYPlus(mesh, runTime, U, yPlus);
            }
            else
            {
                calcIncompressibleYPlus(mesh, runTime, U, yPlus);
            }
        }
        else
        {
            Info<< "    no U field" << endl;
        }

        Info<< "Writing yPlus to field " << yPlus.name() << nl << endl;

        yPlus.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
