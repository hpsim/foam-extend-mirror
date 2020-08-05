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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PODOrthoNormalBase<Type>::PODOrthoNormalBase
(
    const PtrList<GeometricField<Type, fvPatchField, volMesh> >& snapshots,
    const scalar accuracy
)
:
    eigenBase_(snapshots),
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

    Info<< "Cumulative eigen-values: "
        << setprecision(14) << cumEigenValues << nl
        << "Base size: " << baseSize << endl;

    // Establish orthonormal base
    orthoFields_.setSize(baseSize);

    calcOrthoBase(snapshots, orthoFields_);

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
void Foam::PODOrthoNormalBase<Type>::calcOrthoBase
(
    const PtrList<GeoField>& snapshots,
    PtrList<GeoField>& orthoFields
)
{
    // Check if there are less requested orthoFields than eigen vectors
    if (orthoFields.size() > eigenBase_.eigenValues().size())
    {
        FatalErrorInFunction
            << "Requested " << orthoFields.size()
            << " orthogonal fields in snapshot size of "
            << eigenBase_.eigenValues().size() << "(" << snapshots.size() << ")"
            << abort(FatalError);
    }

    forAll (orthoFields, baseI)
    {
        const scalarField& eigenVector = eigenBase_.eigenVectors()[baseI];

        // Reconsider boundary conditions on ortho-normal base.
        // Construct from zeroth snapshot to carry patch types?
        // HJ, 5/Aug/2020
        
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

        // Re-normalise ortho-normal vector
        scalar magSumSquare = Foam::sqrt(sumSqr(onBase));
        if (magSumSquare > SMALL)
        {
            onBase /= magSumSquare;
            onBase.correctBoundaryConditions();
        }

        orthoFields.set(baseI, onBasePtr);
    }
}


// ************************************************************************* //
