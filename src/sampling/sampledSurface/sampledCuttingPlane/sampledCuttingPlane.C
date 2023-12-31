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

#include "sampledCuttingPlane.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "isoSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledCuttingPlane, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledCuttingPlane,
        word,
        cuttingPlane
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledCuttingPlane::createGeometry()
{
    if (debug)
    {
        Pout<< "sampledCuttingPlane::createGeometry :updating geometry."
            << endl;
    }

    // Clear any stored topologies
    facesPtr_.clear();
    isoSurfPtr_.ptr();
    pointDistance_.clear();
    cellDistancePtr_.clear();

    // Clear derived data
    clearGeom();

    // Get any subMesh
    if (zoneID_.index() != -1 && !subMeshPtr_.valid())
    {
        const polyBoundaryMesh& patches = mesh().boundaryMesh();

        // Patch to put exposed internal faces into
        label exposedPatchI = patches.findPatchID(exposedPatchName_);

        if (debug)
        {
            Info<< "Allocating subset of size "
                << mesh().cellZones()[zoneID_.index()].size()
                << " with exposed faces into patch "
                << patches[exposedPatchI].name() << endl;
        }

        const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

        subMeshPtr_.reset
        (
            new fvMeshSubset
            (
                IOobject
                (
                    "set",
                    fvm.time().timeName(),
                    fvm,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvm
            )
        );
        subMeshPtr_().setLargeCellSubset
        (
            labelHashSet(mesh().cellZones()[zoneID_.index()]),
            exposedPatchI
        );
    }


    // Select either the submesh or the underlying mesh
    const fvMesh& fvm =
    (
        subMeshPtr_.valid()
      ? subMeshPtr_().subMesh()
      : static_cast<const fvMesh&>(mesh())
    );


    // Distance to cell centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    cellDistancePtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "cellDistance",
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fvm,
            dimensionedScalar("zero", dimLength, 0)
        )
    );
    volScalarField& cellDistance = cellDistancePtr_();

    // Internal field
    {
        const pointField& cc = fvm.cellCentres();
        scalarField& fld = cellDistance.internalField();

        forAll(cc, i)
        {
            // Signed distance
            fld[i] = (cc[i] - plane_.refPoint()) & plane_.normal();
        }
    }

    // Patch fields
    {
        forAll(cellDistance.boundaryField(), patchI)
        {
            if
            (
                isA<emptyFvPatchScalarField>
                (
                    cellDistance.boundaryField()[patchI]
                )
            )
            {
                cellDistance.boundaryField().set
                (
                    patchI,
                    new calculatedFvPatchScalarField
                    (
                        fvm.boundary()[patchI],
                        cellDistance
                    )
                );

                const polyPatch& pp = fvm.boundary()[patchI].patch();
                pointField::subField cc = pp.patchSlice(fvm.faceCentres());

                fvPatchScalarField& fld = cellDistance.boundaryField()[patchI];
                fld.setSize(pp.size());
                forAll(fld, i)
                {
                    fld[i] = (cc[i] - plane_.refPoint()) & plane_.normal();
                }
            }
            else
            {
                const pointField& cc = fvm.C().boundaryField()[patchI];
                fvPatchScalarField& fld = cellDistance.boundaryField()[patchI];

                forAll(fld, i)
                {
                    fld[i] =  (cc[i] - plane_.refPoint()) & plane_.normal();
                }
            }
        }
    }


    // On processor patches the mesh.C() will already be the cell centre
    // on the opposite side so no need to swap cellDistance.


    // Distance to points
    pointDistance_.setSize(fvm.nPoints());
    {
        const pointField& pts = fvm.points();

        forAll(pointDistance_, i)
        {
            pointDistance_[i] = (pts[i] - plane_.refPoint()) & plane_.normal();
        }
    }


    if (debug)
    {
        Pout<< "Writing cell distance:" << cellDistance.objectPath() << endl;
        cellDistance.write();
        pointScalarField pDist
        (
            IOobject
            (
                "pointDistance",
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(fvm),
            dimensionedScalar("zero", dimLength, 0)
        );
        pDist.internalField() = pointDistance_;

        Pout<< "Writing point distance:" << pDist.objectPath() << endl;
        pDist.write();
    }


    //- Direct from cell field and point field.
    isoSurfPtr_.reset
    (
        new isoSurface
        (
            cellDistance,
            pointDistance_,
            0.0,
            regularise_
        )
    );

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledCuttingPlane::sampledCuttingPlane
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    plane_(dict),
    mergeTol_(dict.lookupOrDefault("mergeTol", 1E-6)),
    regularise_(dict.lookupOrDefault("regularise", true)),
    average_(dict.lookupOrDefault("average", false)),
    zoneID_(dict.lookupOrDefault("zone", word::null), mesh.cellZones()),
    exposedPatchName_(word::null),
    needsUpdate_(true),
    subMeshPtr_(nullptr),
    cellDistancePtr_(nullptr),
    isoSurfPtr_(nullptr),
    facesPtr_(nullptr)
{
    if (zoneID_.index() != -1)
    {
        dict.lookup("exposedPatchName") >> exposedPatchName_;

        if (mesh.boundaryMesh().findPatchID(exposedPatchName_) == -1)
        {
            FatalErrorIn
            (
                "sampledCuttingPlane::sampledCuttingPlane"
                "(const word&, const polyMesh&, const dictionary&)"
            )   << "Cannot find patch " << exposedPatchName_
                << " in which to put exposed faces." << endl
                << "Valid patches are " << mesh.boundaryMesh().names()
                << exit(FatalError);
        }

        if (debug && zoneID_.index() != -1)
        {
            Info<< "Restricting to cellZone " << zoneID_.name()
                << " with exposed internal faces into patch "
                << exposedPatchName_ << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledCuttingPlane::~sampledCuttingPlane()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledCuttingPlane::needsUpdate() const
{
    // Update for changing mesh
    if (mesh().changing())
    {
        needsUpdate_ = true;
    }

    return needsUpdate_;
}


bool Foam::sampledCuttingPlane::expire()
{
    if (debug)
    {
        Pout<< "sampledCuttingPlane::expire :"
            << " have-facesPtr_:" << facesPtr_.valid()
            << " needsUpdate_:" << needsUpdate_ << endl;
    }

    // Clear any stored topologies
    facesPtr_.clear();

    // Clear derived data
    clearGeom();

    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledCuttingPlane::update()
{
    if (debug)
    {
        Pout<< "sampledCuttingPlane::update :"
            << " have-facesPtr_:" << facesPtr_.valid()
            << " needsUpdate_:" << needsUpdate_ << endl;
    }

    if (!needsUpdate_)
    {
        return false;
    }

    createGeometry();

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField>
Foam::sampledCuttingPlane::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::sampledCuttingPlane::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledCuttingPlane::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledCuttingPlane::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::sampledCuttingPlane::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledCuttingPlane::print(Ostream& os) const
{
    os  << "sampledCuttingPlane: " << name() << " :"
        << "  plane:" << plane_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
