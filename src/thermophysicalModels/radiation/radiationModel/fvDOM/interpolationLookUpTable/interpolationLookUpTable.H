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
    Foam::interpolationLookUpTable

Description
    A list of lists. Interpolates based on the first dimension.
    The values must be positive and monotonically increasing in each dimension

Note
    - Accessing an empty list results in an error.
    - Accessing a list with a single element always returns the same value.

SourceFiles
    interpolationLookUpTable.C

\*---------------------------------------------------------------------------*/

#ifndef interpolationLookUpTable_H
#define interpolationLookUpTable_H

#include "List.H"
#include "ListOps.H"
#include "scalarField.H"
#include "HashTable.H"
#include "IOdictionary.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class interpolationLookUpTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolationLookUpTable
:
    public List<scalarField>
{
private:

    // Privsate data

        //- File name
        fileName fileName_;

        //- Table dimensions
        labelList dim_;

        //- Min on each dimension
        scalarList min_;

        //- Deltas on each dimension
        scalarList delta_;

        //- Maximum on each dimension
        scalarList max_;

        //- Dictionary entries
        List<dictionary> entries_;

        //- Output dictionaries
        List<dictionary> output_;

        //- Input indices from the look up table
        labelList entryIndices_;

        //- Output Indeces from the Look Up Table
        labelList outputIndices_;

        //- Field names and indices
        HashTable<label> fieldIndices_;

        //- Output list containing input and interpolation values of outputs
        scalarList interpOutput_;


    // Private Member Functions

        //- Read the table of data from file
        void readTable(const word& instance, const fvMesh& mesh);

        //- Dimension table from dictionaries input and output
        void dimensionTable();

        //- Find table index by scalarList
        label index(const scalarList&, const bool lastDim=true) const;

        //- Find table index by scalar
        label index(const scalar) const;

        //- Check range of lookup value
        bool checkRange(const scalar, const label) const;

        //- Interpolate function return an scalar
        scalar interpolate
        (
            const label lo,
            const label hi,
            const scalar lookUpValue,
            const label ofield,
            const label interfield
        ) const;

        // Check list is monotonically increasing
        void check() const;

        // find hi index, interpolate and populate interpOutput_
        void findHi(const label lo, const scalar retvals);


public:

    // Constructors

        //- Construct null
        interpolationLookUpTable();

        //- Construct given the name of the file containing the table of data
        interpolationLookUpTable
        (
            const fileName& fn,
            const word& instance,
            const fvMesh& mesh
        );

         //- Construct from dictionary
        interpolationLookUpTable(const dictionary& dict);

        //- Construct copy
        interpolationLookUpTable(const interpolationLookUpTable& interpTable);


    // Member Functions

        //- Return true if the filed exists in the table
        bool found(const word& fieldName) const;

        //- Return the output list given a single input scalar
        const scalarList& lookUp(const scalar);

        //- Write Look Up Table to filename.
        void write
        (
            Ostream& os,
            const fileName& fn,
            const word& instance,
            const fvMesh& mesh
        ) const;


    // Access

        //- Return the index of a field by name
        inline label findFieldIndex(const word& fieldName) const;

        //- Return const access to the output dictionaries
        inline const List<dictionary>& output() const;

        //- Return const access tp the dictionary entries
        inline const List<dictionary>& entries() const;

        //- Return const access to the list of min dimensions
        inline const scalarList& min() const;

        //- Return const access to the list of dimensions
        inline const labelList& dim() const;

        //- Return const access to the deltas in each dimension
        inline const scalarList& delta() const;

        //- Return const access to the list of max dimensions
        inline const scalarList& max() const;

        //- Return const access to the table name
        inline word tableName() const;


     // Member Operators

        //- Return an element of constant List<scalar, Type>
        const scalarField& operator[](const label) const;

        //- Return an element of List<scalar, Type>
        scalarField& operator[](const label);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interpolationLookUpTableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
#   include "interpolationLookUpTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
