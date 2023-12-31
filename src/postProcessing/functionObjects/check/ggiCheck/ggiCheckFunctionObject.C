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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "ggiCheckFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "ggiFvsPatchFields.H"
#include "cyclicGgiFvsPatchFields.H"
#include "overlapGgiFvsPatchFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ggiCheckFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        ggiCheckFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ggiCheckFunctionObject::ggiCheckFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    phiName_(dict.lookup("phi"))
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info << "Creating ggi check function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ggiCheckFunctionObject::start()
{
    return true;
}


bool Foam::ggiCheckFunctionObject::execute(const bool forceWrite)
{
    const objectRegistry& mesh =
        time_.lookupObject<objectRegistry>(regionName_);

    if (mesh.foundObject<surfaceScalarField>(phiName_))
    {
        Info<< "Checking flux " << phiName_ << " GGI balance." << endl;
        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>(phiName_);

        boolList visited(phi.boundaryField().size(), false);

        forAll (phi.boundaryField(), patchI)
        {
            if (visited[patchI])
            {
                continue;
            }
            else if (isA<ggiFvsPatchScalarField>(phi.boundaryField()[patchI]))
            {
                const ggiPolyPatch& ggiPatch =
                    refCast<const ggiPolyPatch>
                    (
                        phi.boundaryField()[patchI].patch().patch()
                    );

                const label shadowPatchI = ggiPatch.shadowIndex();

                visited[patchI] = true;
                visited[shadowPatchI] = true;

                // Calculate local flux and area
                const scalar localArea =
                    gSum(phi.mesh().magSf().boundaryField()[patchI]);

                const scalar localFlux = gSum(phi.boundaryField()[patchI]);
                const scalar localFluxMag = mag(localFlux);

                // Calculate shadow flux and area
                const scalar shadowArea =
                    gSum(phi.mesh().magSf().boundaryField()[shadowPatchI]);

                const scalar shadowFlux =
                    gSum(phi.boundaryField()[shadowPatchI]);
                const scalar shadowFluxMag = mag(shadowFlux);

                // Report
                Info<< "GGI pair (" << ggiPatch.name() << ", "
                    << ggiPatch.shadow().name() << ")" << nl
                    << "Area: " << localArea << " " << shadowArea
                    << " Diff = " << localArea - shadowArea << " or "
                    << mag(localArea - shadowArea)/
                       (Foam::max(localArea, shadowArea) + SMALL)*100
                    << " % " << nl
                    << "Flux: " << localFluxMag << " " << shadowFluxMag
                    << " Diff = " << localFlux + shadowFlux << " or "
                    << mag(localFlux + shadowFlux)/
                       (Foam::max(localFluxMag, shadowFluxMag) + SMALL)*100
                    << " %" << endl;
            }
            else if
            (
                isA<cyclicGgiFvsPatchScalarField>(phi.boundaryField()[patchI])
            )
            {
                const cyclicGgiPolyPatch& ggiPatch =
                    refCast<const cyclicGgiPolyPatch>
                    (
                        phi.boundaryField()[patchI].patch().patch()
                    );

                const label shadowPatchI = ggiPatch.shadowIndex();

                visited[patchI] = true;
                visited[shadowPatchI] = true;

                // Calculate local flux and area
                const scalar localArea =
                    gSum(phi.mesh().magSf().boundaryField()[patchI]);

                const scalar localFlux = gSum(phi.boundaryField()[patchI]);
                const scalar localFluxMag = mag(localFlux);

                // Calculate shadow flux and area
                const scalar shadowArea =
                    gSum(phi.mesh().magSf().boundaryField()[shadowPatchI]);

                const scalar shadowFlux =
                    gSum(phi.boundaryField()[shadowPatchI]);
                const scalar shadowFluxMag = mag(shadowFlux);

                // Report
                Info<< "Cyclic GGI pair (" << ggiPatch.name() << ", "
                    << ggiPatch.shadow().name() << ")" << nl
                    << "Area: " << localArea << " " << shadowArea
                    << " Diff = " << localArea - shadowArea << " or "
                    << mag(localArea - shadowArea)/
                       (Foam::max(localArea, shadowArea) + SMALL)*100
                    << " % " << nl
                    << "Flux: " << localFluxMag << " " << shadowFluxMag
                    << " Diff = " << localFlux + shadowFlux << " or "
                    << mag(localFlux + shadowFlux)/
                       (Foam::max(localFluxMag, shadowFluxMag) + SMALL)*100
                    << " %" << endl;
            }
            else if
            (
                isA<overlapGgiFvsPatchScalarField>(phi.boundaryField()[patchI])
            )
            {
                const overlapGgiPolyPatch& ggiPatch =
                    refCast<const overlapGgiPolyPatch>
                    (
                        phi.boundaryField()[patchI].patch().patch()
                    );

                const label shadowPatchI = ggiPatch.shadowIndex();

                const overlapGgiPolyPatch& ggiShadowPatch =
                    refCast<const overlapGgiPolyPatch>
                    (
                        phi.boundaryField()[shadowPatchI].patch().patch()
                    );

                visited[patchI] = true;
                visited[shadowPatchI] = true;

                // Calculate local flux and area
                const scalar localArea =
                    gSum(phi.mesh().magSf().boundaryField()[patchI])
                   *ggiPatch.nCopies();

                const scalar localFlux =
                    gSum(phi.boundaryField()[patchI])
                   *ggiPatch.nCopies();
                const scalar localFluxMag = mag(localFlux);

                // Calculate shadow flux and area
                const scalar shadowArea =
                    gSum(phi.mesh().magSf().boundaryField()[shadowPatchI])
                   *ggiShadowPatch.nCopies();

                const scalar shadowFlux =
                    gSum(phi.boundaryField()[shadowPatchI])
                   *ggiShadowPatch.nCopies();
                const scalar shadowFluxMag = mag(shadowFlux);

                Info<< "Overlap GGI pair (" << ggiPatch.name() << ", "
                    << ggiPatch.shadow().name() << ")" << nl
                    << "Area: " << localArea << " " << shadowArea
                    << " Diff = " << localArea - shadowArea << " or "
                    << mag(localArea - shadowArea)/
                       (Foam::max(localArea, shadowArea) + SMALL)*100
                    << " % " << nl
                    << "Flux: " << localFlux << " " << shadowFlux
                    << " Diff = " << localFlux + shadowFlux << " or "
                    << mag(localFlux + shadowFlux)/
                       (localFluxMag + SMALL)*100
                    << " %" << endl;
            }
        }

        return true;
    }
    else
    {
        InfoIn("bool ggiCheckFunctionObject::execute(const bool forceWrite)")
            << "Cannot find flux field " << phiName_
            << endl;

        return false;
    }
}


bool Foam::ggiCheckFunctionObject::timeSet()
{
    return true;
}


bool Foam::ggiCheckFunctionObject::read(const dictionary& dict)
{
    return false;
}

// ************************************************************************* //
