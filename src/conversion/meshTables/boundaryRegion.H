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
    Foam::boundaryRegion

Description
    The boundaryRegion persistent data saved as a Map<dictionary>.

    The meshReader supports boundaryRegion information.

    The <tt>constant/boundaryRegion</tt> file is an @c IOMap<dictionary>
    that is used to save the information persistently.
    It contains the boundaryRegion information of the following form:

    @verbatim
        (
            INT
            {
                BoundaryType    WORD;
                Label           WORD;
            }
            ...
        )
    @endverbatim

SourceFiles
    boundaryRegion.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryRegion_H
#define boundaryRegion_H

#include "polyMesh.H"
#include "Map.H"
#include "dictionary.H"
#include "labelList.H"
#include "wordList.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class boundaryRegion Declaration
\*---------------------------------------------------------------------------*/

class boundaryRegion
:
    public Map<dictionary>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        boundaryRegion(const boundaryRegion&);


public:

    // Constructors

        //- Construct null
        boundaryRegion();

        //- Construct read from registry, name. instance
        boundaryRegion
        (
            const objectRegistry&,
            const word& name = "boundaryRegion",
            const fileName& instance = "constant"
        );


    //- Destructor
    ~boundaryRegion();


    // Member Functions

        //- Append to the end, return index
        label append(const dictionary&);

        //- Return index corresponding to patch 'name'
        //  returns -1 if not found
        label findIndex(const word& name) const;

        //- Return a Map of (id => name)
        Map<word> names() const;

        //- Return a Map of (id => names) selected by patterns
        Map<word> names(const List<wordRe>& patterns) const;

        //- Return a Map of (id => type)
        Map<word> boundaryTypes() const;

        //- Return BoundaryType corresponding to patch 'name'
        word boundaryType(const word& name) const;

        //- Read constant/boundaryRegion
        void readDict
        (
            const objectRegistry&,
            const word& name = "boundaryRegion",
            const fileName& instance = "constant"
        );

        //- Write constant/boundaryRegion for later reuse
        void writeDict
        (
            const objectRegistry&,
            const word& name = "boundaryRegion",
            const fileName& instance = "constant"
        ) const;


    // Member Operators

        //- Assignment
        void operator=(const boundaryRegion&);

        //- Assign from Map<dictionary>
        void operator=(const Map<dictionary>&);


    // Friend Functions

        //- Rename regions
        //  each dictionary entry is a single word:
        //  @verbatim
        //      newPatchName    originalName;
        //  @endverbatim
        void rename(const dictionary&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
