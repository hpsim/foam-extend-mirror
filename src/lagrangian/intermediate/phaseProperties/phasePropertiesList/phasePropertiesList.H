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
    Foam::phasePropertiesList

Description
    Simple container for a list of phase properties

SourceFiles
    phasePropertiesList.C

\*---------------------------------------------------------------------------*/

#ifndef phasePropertiesList_H
#define phasePropertiesList_H

#include "Istream.H"
#include "volFields.H"
#include "wordList.H"
#include "phaseProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class phasePropertiesList Declaration
\*---------------------------------------------------------------------------*/

class phasePropertiesList
{
    // Private data

        //- List of phase properties
        List<phaseProperties> props_;

        //- List of word representation of phase types
        wordList phaseTypeNames_;

        //- List of state labels
        wordList stateLabels_;


public:

    //- Constructor
    phasePropertiesList
    (
        Istream& is,
        const wordList& gasNames,
        const wordList& liquidNames,
        const wordList& solidNames
    );

    //- Destructor
    ~phasePropertiesList();


    // Public member functions

        //- Return the list of phase properties
        const List<phaseProperties>& props() const;

        //- Return the list of word representation of phase types
        const wordList& phaseTypes() const;

        //- Return the list of state labels
        const wordList& stateLabels() const;

        //- Return the size (number of phases)
        label size() const;


    // Member operators

        const phaseProperties& operator[](const label) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
