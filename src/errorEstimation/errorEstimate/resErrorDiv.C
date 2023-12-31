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

Description
    Residual error estimate for the fv convection operators.

\*---------------------------------------------------------------------------*/

#include "resErrorDiv.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace resError
{

template<class Type>
tmp<errorEstimate<Type> >
div
(
    const surfaceScalarField& flux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    const scalarField& vols = mesh.V();
    const surfaceVectorField& faceCentres = mesh.Cf();
    const volVectorField& cellCentres = mesh.C();
    const fvPatchList& patches = mesh.boundary();
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    Field<Type> res(vols.size(), pTraits<Type>::zero);
    scalarField aNorm(vols.size(), 0.0);

    // Get sign of flux
    const surfaceScalarField signF = pos(flux);

    // Calculate gradient of the solution
    // Change of return type due to gradient cacheing.  HJ, 22/Apr/2016
    const tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvPatchField,
            volMesh
        >
    > tgradVf = fvc::grad(vf);

    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& gradVf = tgradVf();

    // Internal faces
    forAll (owner, faceI)
    {
        // Calculate the centre of the face
        const vector& curFaceCentre = faceCentres[faceI];

        // Owner
        vector ownD = curFaceCentre - cellCentres[owner[faceI]];

        // Subtract convection
        res[owner[faceI]] -=
        (
            vf[owner[faceI]]
          + (ownD & gradVf[owner[faceI]])
        )*flux[faceI];

        aNorm[owner[faceI]] += signF[faceI]*flux[faceI];

        // Neighbour
        vector neiD = curFaceCentre - cellCentres[neighbour[faceI]];

        // Subtract convection
        res[neighbour[faceI]] +=
        (
            vf[neighbour[faceI]]
          + (neiD & gradVf[neighbour[faceI]])
        )*flux[faceI];

        aNorm[neighbour[faceI]] -= (1.0 - signF[faceI])*flux[faceI];
    }

    forAll (patches, patchI)
    {
        const vectorField& patchFaceCentres =
            faceCentres.boundaryField()[patchI];

        const scalarField& patchFlux = flux.boundaryField()[patchI];
        const scalarField& patchSignFlux = signF.boundaryField()[patchI];

        const labelList& fCells = patches[patchI].faceCells();

        forAll (fCells, faceI)
        {
            vector d =
                patchFaceCentres[faceI] - cellCentres[fCells[faceI]];

            // Subtract convection
            res[fCells[faceI]] -=
            (
                vf[fCells[faceI]]
              + (d & gradVf[fCells[faceI]])
            )*patchFlux[faceI];

            aNorm[fCells[faceI]] += patchSignFlux[faceI]*patchFlux[faceI];
        }
    }

    res /= vols;
    aNorm /= vols;

    return tmp<errorEstimate<Type> >
    (
        new errorEstimate<Type>
        (
            vf,
            flux.dimensions()*vf.dimensions(),
            res,
            aNorm
        )
    );
}


template<class Type>
tmp<errorEstimate<Type> >
div
(
    const tmp<surfaceScalarField>& tflux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<errorEstimate<Type> > Div(resError::div(tflux(), vf));
    tflux.clear();
    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace resError

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

