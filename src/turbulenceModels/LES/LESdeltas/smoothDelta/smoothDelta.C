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

#include "smoothDelta.H"
#include "addToRunTimeSelectionTable.H"
#include "FaceCellWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(smoothDelta, 0);
addToRunTimeSelectionTable(LESdelta, smoothDelta, dictionary);

scalar smoothDelta::deltaData::maxDeltaRatio = 1.2;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Fill changedFaces (with face labels) and changedFacesInfo (with delta)
// This is the initial set of faces from which to start the waves.
// Since there might be lots of places with delta jumps we can follow various
// strategies for this initial 'seed'.
// - start from single cell/face and let FaceCellWave pick up all others
//   from there. might be quite a few waves before everything settles.
// - start from all faces. Lots of initial transfers.
// We do something inbetween:
// - start from all faces where there is a jump. Since we cannot easily
//   determine this across coupled patches (cyclic, processor) introduce
//   all faces of these and let FaceCellWave sort it out.
void smoothDelta::setChangedFaces
(
    const polyMesh& mesh,
    const volScalarField& delta,
    dynamicLabelList& changedFaces,
    DynamicList<deltaData>& changedFacesInfo
)
{
    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        scalar ownDelta = delta[mesh.faceOwner()[faceI]];

        scalar neiDelta = delta[mesh.faceNeighbour()[faceI]];

        // Check if owner delta much larger than neighbour delta or vice versa

        if (ownDelta > deltaData::maxDeltaRatio * neiDelta)
        {
            changedFaces.append(faceI);
            changedFacesInfo.append(deltaData(ownDelta));
        }
        else if (neiDelta > deltaData::maxDeltaRatio * ownDelta)
        {
            changedFaces.append(faceI);
            changedFacesInfo.append(deltaData(neiDelta));
        }
    }

    // Insert all faces of coupled patches no matter what. Let FaceCellWave
    // sort it out.
    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        if (patch.coupled())
        {
            forAll(patch, patchFaceI)
            {
                label meshFaceI = patch.start() + patchFaceI;

                scalar ownDelta = delta[mesh.faceOwner()[meshFaceI]];

                changedFaces.append(meshFaceI);
                changedFacesInfo.append(deltaData(ownDelta));
            }
        }
    }

    changedFaces.shrink();
    changedFacesInfo.shrink();
}


void smoothDelta::calcDelta()
{
    deltaData::maxDeltaRatio = maxDeltaRatio_;
    const volScalarField& geometricDelta = geometricDelta_();

    // Fill changed faces with info
    dynamicLabelList changedFaces(mesh_.nFaces()/100 + 100);
    DynamicList<deltaData> changedFacesInfo(changedFaces.size());

    setChangedFaces(mesh_, geometricDelta, changedFaces, changedFacesInfo);

    // Set initial field on cells.
    List<deltaData> cellDeltaData(mesh_.nCells());

    forAll(geometricDelta, cellI)
    {
        cellDeltaData[cellI] = geometricDelta[cellI];
    }

    // Set initial field on faces.
    List<deltaData> faceDeltaData(mesh_.nFaces());


    // Propagate information over whole domain.
    FaceCellWave<deltaData> deltaCalc
    (
        mesh_,
        changedFaces,
        changedFacesInfo,
        faceDeltaData,
        cellDeltaData,
        mesh_.nCells()  // max iterations
    );

    forAll(delta_, cellI)
    {
        delta_[cellI] = cellDeltaData[cellI].delta();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

smoothDelta::smoothDelta
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dd
)
:
    LESdelta(name, mesh),
    geometricDelta_
    (
        LESdelta::New("geometricDelta", mesh, dd.subDict(type() + "Coeffs"))
    ),
    maxDeltaRatio_
    (
        readScalar(dd.subDict(type() + "Coeffs").lookup("maxDeltaRatio"))
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void smoothDelta::read(const dictionary& d)
{
    const dictionary& dd(d.subDict(type() + "Coeffs"));

    geometricDelta_().read(dd);
    dd.lookup("maxDeltaRatio") >> maxDeltaRatio_;
    calcDelta();
}


void smoothDelta::correct()
{
    geometricDelta_().correct();

    if (mesh_.changing())
    {
        calcDelta();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
