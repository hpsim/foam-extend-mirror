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

#include "cellSource.H"
#include "volFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::fieldValues::cellSource::validField(const word& fieldName) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> vf;

    if (obr_.foundObject<vf>(fieldName))
    {
        return true;
    }

    return false;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fieldValues::cellSource::setFieldValues
(
    const word& fieldName,
    const bool mustGet
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> vf;

    if (obr_.foundObject<vf>(fieldName))
    {
        return filterField(obr_.lookupObject<vf>(fieldName));
    }

    if (mustGet)
    {
        FatalErrorIn
        (
            "Foam::tmp<Foam::Field<Type> > "
            "Foam::fieldValues::cellSource::setFieldValues"
            "("
                "const word&, "
                "const bool"
            ") const"
        )   << "Field " << fieldName << " not found in database"
            << abort(FatalError);
    }

    return tmp<Field<Type> >(new Field<Type>(0.0));
}


template<class Type>
Type Foam::fieldValues::cellSource::processValues
(
    const Field<Type>& values,
    const scalarField& V,
    const scalarField& weightField
) const
{
    Type result = pTraits<Type>::zero;
    switch (operation_)
    {
        case opSum:
        {
            result = sum(values);
            break;
        }
        case opSumMag:
        {
            result = sum(cmptMag(values));
            break;
        }
        case opAverage:
        {
            result = sum(values)/values.size();
            break;
        }
        case opWeightedAverage:
        {
            result = sum(weightField*values)/sum(weightField);
            break;
        }
        case opVolAverage:
        {
            result = sum(V*values)/sum(V);
            break;
        }
        case opWeightedVolAverage:
        {
            result = sum(weightField*V*values)/sum(weightField*V);
            break;
        }
        case opVolIntegrate:
        {
            result = sum(V*values);
            break;
        }
        case opMin:
        {
            result = min(values);
            break;
        }
        case opMax:
        {
            result = max(values);
            break;
        }
        case opCoV:
        {
            Type meanValue = sum(values*V)/sum(V);

            const label nComp = pTraits<Type>::nComponents;

            for (direction d=0; d<nComp; ++d)
            {
                scalarField vals(values.component(d));
                scalar mean = component(meanValue, d);
                scalar& res = setComponent(result, d);

                res = sqrt(sum(V*sqr(vals - mean))/sum(V))/mean;
            }

            break;
        }
        default:
        {
            // Do nothing
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fieldValues::cellSource::writeValues(const word& fieldName)
{
    const bool ok = validField<Type>(fieldName);

    if (ok)
    {
        Field<Type> values(setFieldValues<Type>(fieldName));
        scalarField V(filterField(mesh().V()));
        scalarField weightField(values.size(), 1.0);

        if (weightFieldName_ != "none")
        {
            weightField = setFieldValues<scalar>(weightFieldName_, true);
        }

        // Combine onto master
        combineFields(values);
        combineFields(V);
        combineFields(weightField);

        if (Pstream::master())
        {
            Type result = processValues(values, V, weightField);

            // Add to result dictionary, over-writing any previous entry
            resultDict_.add(fieldName, result, true);

            if (valueOutput_)
            {
                IOField<Type>
                (
                    IOobject
                    (
                        fieldName + "_" + sourceTypeNames_[source_] + "-"
                            + sourceName_,
                        obr_.time().timeName(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    weightField*values
                ).write();
            }


            file()<< tab << result;

            if (log_) Info<< "    " << operationTypeNames_[operation_]
                << "(" << sourceName_ << ") of " << fieldName
                <<  " = " << result << endl;
        }
    }

    return ok;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fieldValues::cellSource::filterField
(
    const Field<Type>& field
) const
{
    return tmp<Field<Type> >(new Field<Type>(field, cellId_));
}


// ************************************************************************* //
