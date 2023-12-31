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

#include "zeroGradientIbFvPatchField.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void zeroGradientIbFvPatchField<Type>::updateIbValues()
{
    // // Interpolate the values from tri surface using nearest triangle
    // const labelList& nt = this->ibPatch().ibPolyPatch().nearestTri();

    Field<Type>::operator=(this->patchInternalField());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
zeroGradientIbFvPatchField<Type>::zeroGradientIbFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    zeroGradientFvPatchField<Type>(p, iF),
    immersedBoundaryFieldBase<Type>(p, false, pTraits<Type>::zero)
{}


template<class Type>
zeroGradientIbFvPatchField<Type>::zeroGradientIbFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchField<Type>(p, iF),
    immersedBoundaryFieldBase<Type>
    (
        p,
        Switch(dict.lookup("setDeadValue")),
        pTraits<Type>(dict.lookup("deadValue"))
    )
{
    // Since patch does not read a dictionary, the patch type needs to be read
    // manually.  HJ, 6/Sep/2018
    this->readPatchType(dict);

    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    zeroGradientFvPatchField<Type>::evaluate();
}


template<class Type>
zeroGradientIbFvPatchField<Type>::zeroGradientIbFvPatchField
(
    const zeroGradientIbFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchField<Type>(p, iF),
    immersedBoundaryFieldBase<Type>
    (
        p,
        ptf.setDeadValue(),
        ptf.deadValue()
    )
{
    // Note: NO MAPPING.  Fields are created on the immersed boundary
    // HJ, 12/Apr/2012
    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    this->setPatchType(ptf);

    // On creation of the mapped field, the internal field is dummy and
    // cannot be used.  Initialise the value to avoid errors
    // HJ, 1/Dec/2017
    Field<Type>::operator=(pTraits<Type>::zero);
}


template<class Type>
zeroGradientIbFvPatchField<Type>::zeroGradientIbFvPatchField
(
    const zeroGradientIbFvPatchField<Type>& ptf
)
:
    zeroGradientFvPatchField<Type>(ptf),
    immersedBoundaryFieldBase<Type>
    (
        ptf.ibPatch(),
        ptf.setDeadValue(),
        ptf.deadValue()
    )
{
    this->setPatchType(ptf);
}


template<class Type>
zeroGradientIbFvPatchField<Type>::zeroGradientIbFvPatchField
(
    const zeroGradientIbFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    zeroGradientFvPatchField<Type>(ptf, iF),
    immersedBoundaryFieldBase<Type>
    (
        ptf.ibPatch(),
        ptf.setDeadValue(),
        ptf.deadValue()
    )
{
    this->setPatchType(ptf);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void zeroGradientIbFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    // Base fields do not map: re-interpolate them from tri data
    this->updateIbValues();
}


template<class Type>
void zeroGradientIbFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList&
)
{
    // Base fields do not rmap: re-interpolate them from tri data
    this->updateIbValues();
}


template<class Type>
void zeroGradientIbFvPatchField<Type>::updateOnMotion()
{
    if (this->size() != this->ibPatch().size())
    {
        this->updateIbValues();
    }
}


template<class Type>
void zeroGradientIbFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    this->updateIbValues();

    // Set dead value
    this->setDeadValues(*this);

    // Evaluate fixed value condition
    zeroGradientFvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::zeroGradientIbFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix
)
{
    this->setDeadValues(matrix);
}


template<class Type>
void zeroGradientIbFvPatchField<Type>::write(Ostream& os) const
{
    // Resolve post-processing issues.  HJ, 1/Dec/2017
    zeroGradientFvPatchField<Type>::write(os);
    // triValue_.writeEntry("triValue", os);
    immersedBoundaryFieldBase<Type>::writeDeadData(os);

    // The value entry needs to be written with zero size
    Field<Type>::null().writeEntry("value", os);
    // this->writeEntry("value", os);

    this->writeField(*this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
