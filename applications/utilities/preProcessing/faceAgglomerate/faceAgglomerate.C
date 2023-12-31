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
    faceAgglomerate

Description
    Agglomerate boundary faces using the pairPatchAgglomeration algorithm.
    It writes a map from the fine to coarse grid.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "volFields.H"
#include "unitConversion.H"
#include "pairPatchAgglomeration.H"
#include "labelIOList.H"
#include "syncTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
    argList::validOptions.insert
    (
        "dict",
        "name of dictionary to provide patch agglomeration controls"
    );

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createNamedMesh.H"

    word agglomDictName("faceAgglomerateDict");
    args.optionReadIfPresent("dict", agglomDictName);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    labelListIOList finalAgglom
    (
        IOobject
        (
            "finalAgglom",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        patches.size()
    );


    // Read control dictionary
    IOdictionary agglomDict
    (
        IOobject
        (
            agglomDictName,
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    bool writeAgglom = readBool(agglomDict.lookup("writeFacesAgglomeration"));

    const polyBoundaryMesh& boundary = mesh.boundaryMesh();

    forAll(boundary, patchId)
    {
        const polyPatch& pp = boundary[patchId];

        label patchI = pp.index();
        finalAgglom[patchI].setSize(pp.size(), 0);

        if (!pp.coupled())
        {
            if (agglomDict.found(pp.name()))
            {
                Info << "\nAgglomerating patch : " << pp.name() << endl;
                pairPatchAgglomeration agglomObject
                (
                    pp,
                    agglomDict.subDict(pp.name())
                );
                agglomObject.agglomerate();
                finalAgglom[patchI] =
                    agglomObject.restrictTopBottomAddressing();
            }
            else
            {
                FatalErrorIn(args.executable())
                    << "Patch " << pp.name() << " not found in dictionary: "
                    << agglomDict.name() << exit(FatalError);
            }
        }
    }

    // Sync agglomeration across coupled patches
    labelList nbrAgglom(mesh.nFaces() - mesh.nInternalFaces(), -1);

    forAll(boundary, patchId)
    {
        const polyPatch& pp = boundary[patchId];
        if (pp.coupled())
        {
            finalAgglom[patchId] = identity(pp.size());
            forAll(pp, i)
            {
                nbrAgglom[pp.start() - mesh.nInternalFaces() + i] =
                    finalAgglom[patchId][i];
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh, nbrAgglom, false);
    forAll(boundary, patchId)
    {
        const polyPatch& pp = boundary[patchId];
        if (pp.coupled() && !refCast<const coupledPolyPatch>(pp).master())
        {
            forAll(pp, i)
            {
                finalAgglom[patchId][i] =
                    nbrAgglom[pp.start() - mesh.nInternalFaces() + i];
            }
        }
    }

    finalAgglom.write();

    if (writeAgglom)
    {
        volScalarField facesAgglomeration
        (
            IOobject
            (
                "facesAgglomeration",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("facesAgglomeration", dimless, 0)
        );

        forAll(boundary, patchId)
        {

            fvPatchScalarField& bFacesAgglomeration =
                facesAgglomeration.boundaryField()[patchId];

            forAll(bFacesAgglomeration, j)
            {
                bFacesAgglomeration[j] = finalAgglom[patchId][j];
            }
        }

        Info << "\nWriting  facesAgglomeration" << endl;
        facesAgglomeration.write();
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
