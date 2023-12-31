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

#include "mixedFixedValueSlipFvPatchField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
mixedFixedValueSlipFvPatchField<Type>::mixedFixedValueSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(p, iF),
    refValue_(p.size(), pTraits<Type>::zero),
    valueFraction_(p.size(), 1.0)
{}


template<class Type>
mixedFixedValueSlipFvPatchField<Type>::mixedFixedValueSlipFvPatchField
(
    const mixedFixedValueSlipFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<Type>(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    valueFraction_(ptf.valueFraction_, mapper)
{}


template<class Type>
mixedFixedValueSlipFvPatchField<Type>::mixedFixedValueSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<Type>(p, iF),
    refValue_("refValue", dict, p.size()),
    valueFraction_("valueFraction", dict, p.size())
{}


template<class Type>
mixedFixedValueSlipFvPatchField<Type>::mixedFixedValueSlipFvPatchField
(
    const mixedFixedValueSlipFvPatchField<Type>& ptf
)
:
    transformFvPatchField<Type>(ptf),
    refValue_(ptf.refValue_),
    valueFraction_(ptf.valueFraction_)
{}

template<class Type>
mixedFixedValueSlipFvPatchField<Type>::mixedFixedValueSlipFvPatchField
(
    const mixedFixedValueSlipFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(ptf, iF),
    refValue_(ptf.refValue_),
    valueFraction_(ptf.valueFraction_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
template<class Type>
void mixedFixedValueSlipFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
    refValue_.autoMap(m);
    valueFraction_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
template<class Type>
void mixedFixedValueSlipFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    transformFvPatchField<Type>::rmap(ptf, addr);

    const mixedFixedValueSlipFvPatchField<Type>& dmptf =
        refCast<const mixedFixedValueSlipFvPatchField<Type> >(ptf);

    refValue_.rmap(dmptf.refValue_, addr);
    valueFraction_.rmap(dmptf.valueFraction_, addr);
}


// Return gradient at boundary
template<class Type>
tmp<Field<Type> > mixedFixedValueSlipFvPatchField<Type>::snGrad() const
{
    vectorField nHat = this->patch().nf();
    Field<Type> pif = this->patchInternalField();

    return
    (
        valueFraction_*refValue_
      + (1.0 - valueFraction_)*transform(I - sqr(nHat), pif) - pif
    )*this->patch().deltaCoeffs();
}


// Evaluate the field on the patch
template<class Type>
void mixedFixedValueSlipFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField nHat = this->patch().nf();

    Field<Type>::operator=
    (
        valueFraction_*refValue_
      +
        (1.0 - valueFraction_)
       *transform(I - nHat*nHat, this->patchInternalField())
    );

    transformFvPatchField<Type>::evaluate();
}


// Return defining fields
template<class Type>
tmp<Field<Type> > mixedFixedValueSlipFvPatchField<Type>::snGradTransformDiag() const
{
    vectorField nHat = this->patch().nf();
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return
        valueFraction_*Type(pTraits<Type>::one)
      + (1.0 - valueFraction_)*transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


// Write
template<class Type>
void mixedFixedValueSlipFvPatchField<Type>::write(Ostream& os) const
{
    transformFvPatchField<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    valueFraction_.writeEntry("valueFraction", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
