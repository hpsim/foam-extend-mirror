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
    thermalSourceModel

Description
    Thermal source model for fluids and solids.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    thermalSourceModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermalSourceModel_H
#define thermalSourceModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "tmp.H"
#include "thermalLaw.H"
#include "thermalSource.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class thermalSourceModel Declaration
\*---------------------------------------------------------------------------*/

class thermalSourceModel
{
    // Private data

        //- Dictionary
        const dictionary& dict_;

        //- Reference to temperature field
        const volScalarField& T_;

        //- Thermal source
        PtrList<thermalSource> sources_;


    // Private Member Functions

        //- Disallow copy construct
        thermalSourceModel(const thermalSourceModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const thermalSourceModel&) = delete;


public:

    //- Runtime type information
    TypeName("thermalSourceModel");

    // Constructors

        //- Construct from dictionary
        thermalSourceModel
        (
            const dictionary& dict,
            const volScalarField& T
        );


    //- Destructor
    virtual ~thermalSourceModel() = default;


    // Member Functions

        //- Is the thermal model active?
        inline bool active() const
        {
            return !sources_.empty();
        }

        //- Return the aggregated source
        tmp<volScalarField> S() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
