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

#include "heatFluxItem.H"
#include "volFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::heatFlux::addMeanField
(
    const label fieldI,
    wordList& meanFieldList
) const
{
    if (faItems_[fieldI].mean())
    {
        typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

        const word& fieldName = faItems_[fieldI].fieldName();

        const word meanFieldName = fieldName + EXT_MEAN;

        Info<< "Reading/calculating field " << meanFieldName << nl << endl;

        if (obr_.foundObject<fieldType>(meanFieldName))
        {
            meanFieldList[fieldI] = meanFieldName;
        }
        else if (obr_.found(meanFieldName))
        {
            Info<< "Cannot allocate average field " << meanFieldName
                << " since an object with that name already exists."
                << " Disabling averaging." << nl << endl;
            meanFieldList[fieldI] = word::null;
        }
        else
        {
            const fieldType& baseField =
                obr_.lookupObject<fieldType>(fieldName);

            // Store on registry
            obr_.store
            (
                new fieldType
                (
                    IOobject
                    (
                        meanFieldName,
                        obr_.time().timeName(),
                        obr_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    1*baseField
                )
            );

            meanFieldList[fieldI] = meanFieldName;
        }
    }
}


template<class Type1, class Type2>
void Foam::heatFlux::addPrime2MeanField
(
    const label fieldI,
    const wordList& meanFieldList,
    wordList& prime2MeanFieldList
) const
{
    if (faItems_[fieldI].mean() && meanFieldList[fieldI].size())
    {
        typedef GeometricField<Type1, fvPatchField, volMesh> fieldType1;
        typedef GeometricField<Type2, fvPatchField, volMesh> fieldType2;

        const word& fieldName = faItems_[fieldI].fieldName();

        const word meanFieldName = fieldName + EXT_PRIME2MEAN;
        Info<< "Reading/calculating field " << meanFieldName << nl << endl;

        if (obr_.foundObject<fieldType2>(meanFieldName))
        {
            prime2MeanFieldList[fieldI] = meanFieldName;
        }
        else if (obr_.found(meanFieldName))
        {
            Info<< "Cannot allocate average field " << meanFieldName
                << " since an object with that name already exists."
                << " Disabling averaging." << nl << endl;
            prime2MeanFieldList[fieldI] = word::null;
        }
        else
        {
            const fieldType1& baseField =
                obr_.lookupObject<fieldType1>(fieldName);
            const fieldType1& meanField =
                obr_.lookupObject<fieldType1>(meanFieldList[fieldI]);

            obr_.store
            (
                new fieldType2
                (
                    IOobject
                    (
                        meanFieldName,
                        obr_.time().timeName(),
                        obr_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    sqr(baseField) - sqr(meanField)
                )
            );

            prime2MeanFieldList[fieldI] = meanFieldName;
        }
    }
}


template<class Type>
void Foam::heatFlux::calculateMeanFields(const wordList& meanFieldList)
const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const scalar dt = obr_.time().deltaTValue();

    forAll(faItems_, i)
    {
        if (faItems_[i].mean() && meanFieldList[i].size())
        {
            const word& fieldName = faItems_[i].fieldName();
            const fieldType& baseField =
                obr_.lookupObject<fieldType>(fieldName);
            fieldType& meanField = const_cast<fieldType&>
            (
                obr_.lookupObject<fieldType>(meanFieldList[i])
            );

            scalar alpha = 0.0;
            scalar beta = 0.0;
            if (faItems_[i].timeBase())
            {
                 alpha = (totalTime_[i] - dt)/totalTime_[i];
                 beta = dt/totalTime_[i];
            }
            else
            {
                alpha = scalar(totalIter_[i] - 1)/scalar(totalIter_[i]);
                beta = 1.0/scalar(totalIter_[i]);
            }

            meanField = alpha*meanField + beta*baseField;
        }
    }
}


template<class Type1, class Type2>
void Foam::heatFlux::calculatePrime2MeanFields
(
    const wordList& meanFieldList,
    const wordList& prime2MeanFieldList
) const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> fieldType1;
    typedef GeometricField<Type2, fvPatchField, volMesh> fieldType2;

    const scalar dt = obr_.time().deltaTValue();

    forAll(faItems_, i)
    {
        if
        (
            faItems_[i].prime2Mean()
         && meanFieldList[i].size()
         && prime2MeanFieldList[i].size()
        )
        {
            const word& fieldName = faItems_[i].fieldName();
            const fieldType1& baseField =
                obr_.lookupObject<fieldType1>(fieldName);
            const fieldType1& meanField =
                obr_.lookupObject<fieldType1>(meanFieldList[i]);
            fieldType2& prime2MeanField = const_cast<fieldType2&>
            (
                obr_.lookupObject<fieldType2>(prime2MeanFieldList[i])
            );

            scalar alpha = 0.0;
            scalar beta = 0.0;
            if (faItems_[i].timeBase())
            {
                alpha = (totalTime_[i] - dt)/totalTime_[i];
                beta = dt/totalTime_[i];
            }
            else
            {
                alpha = scalar(totalIter_[i] - 1)/scalar(totalIter_[i]);
                beta = 1.0/scalar(totalIter_[i]);
            }

            prime2MeanField =
                alpha*prime2MeanField
              + beta*sqr(baseField)
              - sqr(meanField);
        }
    }
}


template<class Type1, class Type2>
void Foam::heatFlux::addMeanSqrToPrime2Mean
(
    const wordList& meanFieldList,
    const wordList& prime2MeanFieldList
) const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> fieldType1;
    typedef GeometricField<Type2, fvPatchField, volMesh> fieldType2;

    forAll(faItems_, i)
    {
        if
        (
            faItems_[i].prime2Mean()
         && meanFieldList[i].size()
         && prime2MeanFieldList[i].size()
        )
        {
            const fieldType1& meanField =
                obr_.lookupObject<fieldType1>(meanFieldList[i]);
            fieldType2& prime2MeanField = const_cast<fieldType2&>
            (
                obr_.lookupObject<fieldType2>(prime2MeanFieldList[i])
            );

            prime2MeanField += sqr(meanField);
        }
    }
}


template<class Type>
void Foam::heatFlux::writeFieldList(const wordList& fieldList) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    forAll(fieldList, i)
    {
        if (fieldList[i].size())
        {
            const fieldType& f = obr_.lookupObject<fieldType>(fieldList[i]);
            f.write();
        }
    }
}


// ************************************************************************* //
