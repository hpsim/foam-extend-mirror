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
    PODEigenBase

\*---------------------------------------------------------------------------*/

#include "PODEigenBase.H"
#include "volFields.H"
#include "POD.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::PODEigenBase<Type>::calcEigenBase
(
    const scalarSquareMatrix& orthMatrix
)
{
    // Calculate eigen-values

    EigenSolver<scalar> eigenSolver(orthMatrix);

    // Sort and assemble

    SortableList<scalar> sortedList(orthMatrix.n());

    forAll (sortedList, i)
    {
        // Note:
        // For stability under normalisation, make sure that very small
        // eigen-values remain positive
        // HJ, 19/Feb/2021
        sortedList[i] = mag(eigenSolver.eigenValue(i));
    }

    // Sort will sort the values in descending order and insert values
    sortedList.sort();

    // Add normalisation to eigenvectors
    // HJ, 19/Feb/2021
    label n = 0;
    forAllReverse (sortedList, i)
    {
        eigenValues_[n] = sortedList[i];
        eigenVectors_.set
        (
            n,
            new scalarField
            (
                eigenSolver.eigenVector(sortedList.indices()[i])*
                sqrt(eigenVectors_.size()*mag(eigenValues_[n]))
            )
        );

        n++;
    }

    // Info<< "Check products of base vectors" << endl;
    // forAll (eigenVectors_, i)
    // {
    //     Info<< " i " << i << " eigenValue " << eigenValues_[i]
    //         << " scale " << eigenVectors_.size()*mag(eigenValues_[i])
    //         << nl;

    //         forAll (eigenVectors_, j)
    //         {
    //             Info<< "(" << i << ", " << j << ") = "
    //                 <<
    //                 POD::projection
    //                 (
    //                     eigenVectors_[i],
    //                     eigenVectors_[j]
    //                 )
    //                 << nl;
    //         }
    //         Info<< nl << endl;
    // }

    // Assemble cumulative relative eigen-values
    cumEigenValues_[0] = eigenValues_[0];

    // Assemble accumulated normalised eigenvalues
    for (label i = 1; i < cumEigenValues_.size(); i++)
    {
        cumEigenValues_[i] = cumEigenValues_[i - 1] + eigenValues_[i];
    }

    // Renormalise
    cumEigenValues_ /= sum(eigenValues_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given a list of fields
template<class Type>
Foam::PODEigenBase<Type>::PODEigenBase
(
    const PtrList<GeometricField<Type, fvPatchField, volMesh> >& snapshots
)
:
    eigenValues_(snapshots.size()),
    cumEigenValues_(snapshots.size()),
    eigenVectors_(snapshots.size())
{
    // Calculate the snapshot products of the field with all available fields
    label nSnapshots = snapshots.size();

    scalarSquareMatrix orthMatrix(nSnapshots);

    for (label snapI = 0; snapI < nSnapshots; snapI++)
    {
        for (label snapJ = 0; snapJ <= snapI; snapJ++)
        {
            // Calculate the inner product and insert it into the matrix
            // Added normalisation with the number of snapshots
            // Note projection with cell volumes
            orthMatrix[snapI][snapJ] =
                1.0/nSnapshots*
                POD::projection
                (
                    snapshots[snapI].internalField(),
                    snapshots[snapJ].internalField()
                );

            // Fill in the symmetric part of the matrix
            if (snapI != snapJ)
            {
                orthMatrix[snapJ][snapI] = orthMatrix[snapI][snapJ];
            }
        }
    }

    // Calculate eigenbase: this fills in the eigenvalues and eigenvectors
    this->calcEigenBase(orthMatrix);
}


// ************************************************************************* //
