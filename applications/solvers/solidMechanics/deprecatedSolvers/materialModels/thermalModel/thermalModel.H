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
    thermalModel

Description
    Thermal material properties for solids.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    thermalModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermalModel_H
#define thermalModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "thermalLaw.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalModel Declaration
\*---------------------------------------------------------------------------*/

class thermalModel
:
    public IOdictionary
{
    // Private data

        //- Reference to temperature field
        const volScalarField& T_;

        //- Thermal law
        autoPtr<thermalLaw> lawPtr_;


    // Private Member Functions

        //- Disallow copy construct
        thermalModel(const thermalModel&);

        //- Disallow default bitwise assignment
        void operator=(const thermalModel&);


public:

    //- Runtime type information
    TypeName("thermalModel");

    // Constructors

        //- Construct from dictionary
        thermalModel(const volScalarField& T);


    // Destructor

        virtual ~thermalModel()
        {}


    // Member Functions

        //- Return thermal law
        const thermalLaw& law() const
        {
            return lawPtr_();
        }

        //- Return specific heat capacity
        tmp<volScalarField> C() const
        {
            return lawPtr_->C();
        }

        //- Return thermal conductivity
        tmp<volScalarField> k() const
        {
            return lawPtr_->k();
        }

        //- Return thermal expansion coefficient
        tmp<volScalarField> alpha() const
        {
            return lawPtr_->alpha();
        }

        //- Return reference temperature
        tmp<volScalarField> T0() const
        {
            return lawPtr_->T0();
        }

        //- Correct the thermal model
        void correct()
        {
            lawPtr_->correct();
        }

        //- Read thermalProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
