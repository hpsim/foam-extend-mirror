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

#include "error.H"

#include "dispersionModel.H"
#include "noDispersion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<dispersionModel> dispersionModel::New
(
    const dictionary& dict,
    spray& sm
)
{
    word dispersionModelType
    (
        dict.lookup("dispersionModel")
    );

    Info<< "Selecting dispersionModel "
         << dispersionModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dispersionModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "dispersionModel::New(const dictionary&, const spray&) : "
            << endl
            << "    unknown dispersionModelType type "
            << dispersionModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid dispersionModel types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc() << abort(FatalError);
    }

    return autoPtr<dispersionModel>(cstrIter()(dict, sm));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
