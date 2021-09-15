/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "processorFaPatchField.H"
#include "processorFaPatch.H"
#include "IPstream.H"
#include "OPstream.H"
#include "demandDrivenData.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFaPatch>(p)),
    pnf_()
{}


template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& f
)
:
    coupledFaPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFaPatch>(p)),
    pnf_()
{}


// Construct by mapping given processorFaPatchField<Type>
template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const processorFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    coupledFaPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFaPatch>(p)),
    pnf_()
{
    if (!isType<processorFaPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    coupledFaPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorFaPatch>(p)),
    pnf_()
{
    if (!isType<processorFaPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const processorFaPatchField<Type>& ptf
)
:
    processorLduInterfaceField(),
    coupledFaPatchField<Type>(ptf),
    procPatch_(refCast<const processorFaPatch>(ptf.patch())),
    pnf_()
{}


template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const processorFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFaPatch>(ptf.patch())),
    pnf_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::processorFaPatchField<Type>::patchNeighbourField() const
{
    if (Pstream::parRun())
    {
        if (pnf_.empty())
        {
            FatalErrorInFunction
                << "Processor patchNeighbourField not ready"
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Attempting to access processor patchNeighbourField "
            << "in a serial run"
            << abort(FatalError);
    }

    return pnf_;
}


template<class Type>
void Foam::processorFaPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        procPatch_.send(commsType, this->patchInternalField()());
    }
}


template<class Type>
void Foam::processorFaPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        pnf_.setSize(this->size());

        procPatch_.receive<Type>(commsType, pnf_);

        if (doTransform())
        {
            transform(*this, procPatch_.forwardT(), pnf_);
        }
    }

    // Evaluate patch field to carry interpolated values
    Field<Type>::operator=
    (
        this->patch().weights()*this->patchInternalField()
      + (1 - this->patch().weights())*pnf_
    );

    // Signal completion
    faPatchField<Type>::evaluate(commsType);
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::processorFaPatchField<Type>::snGrad() const
{
    return this->patch().deltaCoeffs()*(*this - this->patchInternalField());
}


template<class Type>
void Foam::processorFaPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    procPatch_.send
    (
        commsType,
        this->patch().patchInternalField(psiInternal)()
    );
}


template<class Type>
void Foam::processorFaPatchField<Type>::updateInterfaceMatrix
(
    const scalarField&,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    scalarField pnf
    (
        procPatch_.receive<scalar>(commsType, this->size())()
    );

    // Transform according to the transformation tensor
    transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result

    const unallocLabelList& edgeFaces = this->patch().edgeFaces();

    if (switchToLhs)
    {
        forAll(edgeFaces, elemI)
        {
            result[edgeFaces[elemI]] += coeffs[elemI]*pnf[elemI];
        }
    }
    else
    {
        forAll(edgeFaces, elemI)
        {
            result[edgeFaces[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::processorFaPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    pnf_.clear();
    faPatchField<Type>::operator=(ul);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator=
(
    const faPatchField<Type>& ptf
)
{
    pnf_.clear();
    faPatchField<Type>::operator=(ptf);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator+=
(
    const faPatchField<Type>& ptf
)
{
    pnf_.clear();
    faPatchField<Type>::operator+=(ptf);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator-=
(
    const faPatchField<Type>& ptf
)
{
    pnf_.clear();
    faPatchField<Type>::operator-=(ptf);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator*=
(
    const faPatchField<scalar>& ptf
)
{
    pnf_.clear();
    faPatchField<Type>::operator*=(ptf);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator/=
(
    const faPatchField<scalar>& ptf
)
{
    pnf_.clear();
    faPatchField<Type>::operator/=(ptf);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    pnf_.clear();
    faPatchField<Type>::operator+=(tf);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    pnf_.clear();
    faPatchField<Type>::operator-=(tf);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    pnf_.clear();
    faPatchField<Type>::operator*=(tf);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    pnf_.clear();
    faPatchField<Type>::operator/=(tf);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator=
(
    const Type& t
)
{
    pnf_.clear();
    faPatchField<Type>::operator=(t);
}


// Note: it is possible to manipulate pnf_ in operations with primitive
// types, but I am currently avoiding this because this assumes synchronised
// field functions across multiple processors
// HJ, 8/Sep/2021


template<class Type>
void Foam::processorFaPatchField<Type>::operator+=
(
    const Type& t
)
{
    pnf_.clear();
    faPatchField<Type>::operator+=(t);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator-=
(
    const Type& t
)
{
    pnf_.clear();
    faPatchField<Type>::operator-=(t);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator*=
(
    const scalar s
)
{
    pnf_.clear();
    faPatchField<Type>::operator*=(s);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator/=
(
    const scalar s
)
{
    pnf_.clear();
    faPatchField<Type>::operator/=(s);
}


// Force an assignment, overriding fixedValue status
template<class Type>
void Foam::processorFaPatchField<Type>::operator==
(
    const faPatchField<Type>& ptf
)
{
    pnf_.clear();
    faPatchField<Type>::operator=(ptf);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    pnf_.clear();
    faPatchField<Type>::operator=(tf);
}


template<class Type>
void Foam::processorFaPatchField<Type>::operator==
(
    const Type& t
)
{
    pnf_.clear();
    faPatchField<Type>::operator=(t);
}


// ************************************************************************* //
