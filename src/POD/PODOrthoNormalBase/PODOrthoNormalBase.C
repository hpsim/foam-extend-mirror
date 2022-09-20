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

Class
    PODOrthoNormalBase

Description
    Establish POD ortho-normal base and interpolation coefficients give a list
    of fields. Size of ortho-normal base is calculated from the desired
    accuracy, e.g. 0.99-0.99999 (in energy terms)

\*---------------------------------------------------------------------------*/

#include "PODOrthoNormalBase.H"
#include "POD.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class GeoField>
void Foam::PODOrthoNormalBase<Type>::calcOrthoBase
(
    const FieldField<Field, scalar> eigenVectors,
    const PtrList<GeoField>& snapshots,
    PtrList<GeoField>& orthoFields
)
{
    // Check if there are less requested orthoFields than eigen vectors
    if
    (
        eigenVectors.empty()
     || snapshots.empty()
     || (orthoFields.size() > snapshots.size())
     || (eigenVectors[0].size() != snapshots.size())
    )
    {
        FatalErrorInFunction
            << "Incompatible snapshots, eigenVectors and orthoFields: "
            << snapshots.size() << ", " << orthoFields.size()
            << " and " << eigenVectors[0].size()
            << abort(FatalError);
    }

    forAll (orthoFields, baseI)
    {
        // Scale the eigenvector before establishing the ortho-normal base
        const scalarField& eigenVector = eigenVectors[baseI];

        // Use calculated boundary conditions on the eigenbase
        GeoField* onBasePtr
        (
            new GeoField
            (
                IOobject
                (
                    snapshots[0].name() + "POD" + name(baseI),
                    snapshots[0].time().timeName(),
                    snapshots[0].mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                snapshots[0].mesh(),
                dimensioned<typename GeoField::PrimitiveType>
                (
                    "zero",
                    snapshots[0].dimensions(),
                    pTraits<typename GeoField::PrimitiveType>::zero
                )
            )
        );
        GeoField& onBase = *onBasePtr;

        forAll (eigenVector, eigenI)
        {
            onBase += eigenVector[eigenI]*snapshots[eigenI];
        }

        // Note: onBase field is not normalised
        // This is required when the fields are rebuild using consistent
        // scaling factors, eg. flux fields from velocity
        // HJ, 20/Sep/2021
        orthoFields.set(baseI, onBasePtr);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PODOrthoNormalBase<Type>::PODOrthoNormalBase
(
    const PtrList<GeometricField<Type, fvPatchField, volMesh> >& snapshots,
    const scalar accuracy
)
:
    eigenBase_(snapshots),
    scaledEigenVectors_(eigenBase_.eigenVectors()), // To be scaled later
    orthoFields_(),
    interpolationCoeffsPtr_(nullptr)
{
    label baseSize = 0;

    const scalarField& cumEigenValues = eigenBase_.cumulativeEigenValues();

    forAll (cumEigenValues, i)
    {
        baseSize++;

        if (cumEigenValues[i] > accuracy)
        {
            break;
        }
    }

    Info<< setprecision(14)
        << "eigen-values: " << eigenBase_.eigenValues() << nl
        << "Cumulative eigen-values: " << cumEigenValues << nl
        << "Base size: " << baseSize << endl;

    // Establish orthonormal base
    orthoFields_.setSize(baseSize);

    // Calculate orthogonal base
    calcOrthoBase
    (
        scaledEigenVectors_,
        snapshots,
        orthoFields_
    );

    // Scale the eigen-vectors and orthoFields
    forAll (orthoFields_, baseI)
    {
        GeoTypeField& onBase = orthoFields_[baseI];

        scalarField& eigenVector = scaledEigenVectors_[baseI];

        // Re-normalise the ortho-normal vector and eigenVector
        scalar magSumSquare = Foam::sqrt(sumSqr(onBase));

        if (magSumSquare > SMALL)
        {
            // Rescale ortho base
            onBase /= magSumSquare;
            onBase.correctBoundaryConditions();

            // Also rescale scaledEigenVectors_ for future use
            eigenVector /= magSumSquare;
        }
    }

    // Calculate interpolation coefficients
    interpolationCoeffsPtr_ =
        new RectangularMatrix<scalar>(snapshots.size(), orthoFields_.size());
    RectangularMatrix<scalar>& coeffs = *interpolationCoeffsPtr_;

    forAll (snapshots, snapshotI)
    {
        forAll (orthoFields_, baseI)
        {
            coeffs[snapshotI][baseI] =
                POD::projection
                (
                    snapshots[snapshotI],
                    orthoFields_[baseI]
                );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::PODOrthoNormalBase<Type>::~PODOrthoNormalBase()
{
    deleteDemandDrivenData(interpolationCoeffsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::PODOrthoNormalBase<Type>::checkBase() const
{
    // Check orthogonality and magnitude of snapshots
    Info<< "Check orthogonal base: dot-products: " << baseSize() << endl;
    for (label i = 0; i < baseSize(); i++)
    {
        for (label j = 0; j < baseSize(); j++)
        {
            Info<< "(" << i << ", " << j << ") = "
                << POD::projection
                   (
                       orthoField(i),
                       orthoField(j)
                   )
                << endl;
        }
    }
}


template<class Type>
template<class GeoField>
void Foam::PODOrthoNormalBase<Type>::getOrthoBase
(
    const PtrList<GeoField>& snapshots,
    PtrList<GeoField>& orthoFields
) const
{
    calcOrthoBase
    (
        this->scaledEigenVectors(),
        snapshots,
        orthoFields
    );
}


// ************************************************************************* //
