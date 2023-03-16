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

#include "gaussGrad.H"
#include "extrapolatedCalculatedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
gaussGrad<Type>::gradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tgGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "0",
                ssf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();

    Field<GradType>& igGrad = gGrad;
    const Field<Type>& issf = ssf;

    forAll(owner, faceI)
    {
        GradType Sfssf = Sf[faceI]*issf[faceI];

        igGrad[owner[faceI]] += Sfssf;
        igGrad[neighbour[faceI]] -= Sfssf;
    }

    forAll(mesh.boundary(), patchI)
    {
        const unallocLabelList& pFaceCells =
            mesh.boundary()[patchI].faceCells();

        const vectorField& pSf = mesh.Sf().boundaryField()[patchI];

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchI];

        forAll(mesh.boundary()[patchI], faceI)
        {
            igGrad[pFaceCells[faceI]] += pSf[faceI]*pssf[faceI];
        }
    }

    igGrad /= mesh.V();

    gGrad.correctBoundaryConditions();

    return tgGrad;
}


// template<class Type>
// tmp
// <
//     GeometricField
//     <
//         typename outerProduct<vector, Type>::type, fvPatchField, volMesh
//     >
// >
// gaussGrad<Type>::calcGrad
// (
//     const GeometricField<Type, fvPatchField, volMesh>& vsf,
//     const word& name
// ) const
// {
//     typedef typename outerProduct<vector, Type>::type GradType;

//     tmp<GeometricField<GradType, fvPatchField, volMesh> > tgGrad
//     (
//         gradf(tinterpScheme_().interpolate(vsf), name)
//     );
//     GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad();

//     gGrad.rename("grad(" + vsf.name() + ')');
//     this->correctBoundaryConditions(vsf, gGrad);

//     return tgGrad;
// }


// NEW FORMULATION: deltas.  See CJ Marooney OFW17.  HJ, 8/Dec/2022
template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
gaussGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tgGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                vsf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "zero",
                vsf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad();

    // Get weights
    surfaceScalarField w = this->tinterpScheme_().weights(vsf);

    const surfaceVectorField& Sf = mesh.Sf();
    
    // updateCoupledPatchFields for patchNeighbourField update
    // HJ, 10/Sep/2021
    vsf.boundaryField().updateCoupledPatchFields();

    // Get access to internal fields

    const Field<Type>& vsfIn = vsf.internalField();

    Field<GradType>& gGradIn = gGrad.internalField();

    const scalarField& wIn = w.internalField();
    const vectorField& SfIn = Sf.internalField();

    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();

    label ownFaceI, neiFaceI;

    forAll (own, faceI)
    {
        ownFaceI = own[faceI];
        neiFaceI = nei[faceI];

        Type deltaVsf = vsfIn[neiFaceI] - vsfIn[ownFaceI];

        // Both Sf and own-nei swap values: sign remains the same
        // HJ, 8/Dec/2022
        gGradIn[ownFaceI] += (1 - wIn[faceI])*SfIn[faceI]*deltaVsf;
        gGradIn[neiFaceI] += wIn[faceI]*SfIn[faceI]*deltaVsf;
    }

    // Boundary faces
    forAll (vsf.boundaryField(), patchI)
    {
        const scalarField& patchW = w.boundaryField()[patchI];

        const vectorField& patchSf = Sf.boundaryField()[patchI];

        const unallocLabelList& faceCells =
            gGrad.boundaryField()[patchI].patch().faceCells();

        if (vsf.boundaryField()[patchI].coupled())
        {
            Field<Type> neiVsf =
                vsf.boundaryField()[patchI].patchNeighbourField();

            forAll (neiVsf, patchFaceI)
            {
                gGrad[faceCells[patchFaceI]] +=
                    (1 - patchW[patchFaceI])*patchSf[patchFaceI]*
                    (neiVsf[patchFaceI] - vsfIn[faceCells[patchFaceI]]);
            }
        }
        else
        {
            const fvPatchField<Type>& patchVsf = vsf.boundaryField()[patchI];

            forAll (patchVsf, patchFaceI)
            {
                gGrad[faceCells[patchFaceI]] +=
                    (1 - patchW[patchFaceI])*patchSf[patchFaceI]*
                    (patchVsf[patchFaceI] - vsfIn[faceCells[patchFaceI]]);
            }
        }
    }

    gGradIn /= mesh.V();

    gGrad.correctBoundaryConditions();
    this->correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
> gaussGrad<Type>::fvmGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    FatalErrorInFunction
        << "Implicit gradient operator defined only for scalar."
        << abort(FatalError);

    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<BlockLduSystem<vector, GradType> > tbs
    (
        new BlockLduSystem<vector, GradType>(vf.mesh())
    );

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
