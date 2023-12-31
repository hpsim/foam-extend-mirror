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

#include "fvMesh.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::GeometricField
<
    typename Foam::outerProduct<Foam::vector, Type>::type,
    Foam::fvPatchField,
    Foam::volMesh
>&
Foam::calcFvcGrad::gradField(const word& gradName, const dimensionSet& dims)
{
    Info<< "gradField" << endl;

    typedef typename outerProduct<vector, Type>::type gradType;
    typedef GeometricField<gradType, fvPatchField, volMesh> vfGradType;

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    if (!mesh.foundObject<vfGradType>(gradName))
    {
        vfGradType* gradFieldPtr
        (
            new vfGradType
            (
                IOobject
                (
                    gradName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<gradType>
                (
                    "zero",
                    dims/dimLength,
                    pTraits<gradType>::zero
                )
            )
        );

        mesh.objectRegistry::store(gradFieldPtr);
    }

    const vfGradType& field = mesh.lookupObject<vfGradType>(gradName);

    return const_cast<vfGradType&>(field);
}


template<class Type>
void Foam::calcFvcGrad::calcGrad
(
    const word& fieldName,
    const word& resultName,
    bool& processed
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> vfType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sfType;

    typedef typename outerProduct<vector, Type>::type gradType;
    typedef GeometricField<gradType, fvPatchField, volMesh> vfGradType;

    const fvMesh& mesh = refCast<const fvMesh>(obr_);


    if (mesh.foundObject<vfType>(fieldName))
    {
        const vfType& vf = mesh.lookupObject<vfType>(fieldName);

        vfGradType& field = gradField<Type>(resultName, vf.dimensions());

        // De-reference the tmp to avoid a clash with the cached grad field
        field = fvc::grad(vf)();

        processed = true;
    }
    else if (mesh.foundObject<sfType>(fieldName))
    {
        const sfType& sf = mesh.lookupObject<sfType>(fieldName);

        vfGradType& field = gradField<Type>(resultName, sf.dimensions());

        // De-reference the tmp to avoid a clash with the cached grad field
        field = fvc::grad(sf)();

        processed = true;
    }
}


// ************************************************************************* //
