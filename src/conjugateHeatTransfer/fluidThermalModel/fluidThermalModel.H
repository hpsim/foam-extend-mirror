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
    fluidThermalModel

Description
    Thermal model for fluids.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    fluidThermalModel.C

\*---------------------------------------------------------------------------*/

#ifndef fluidThermalModel_H
#define fluidThermalModel_H

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
                      Class fluidThermalModel Declaration
\*---------------------------------------------------------------------------*/

class fluidThermalModel
{
    // Private data

        //- Dictionary
        const dictionary& dict_;

        //- Reference to temperature field
        const volScalarField& T_;

        //- Thermal law
        autoPtr<thermalLaw> lawPtr_;

        //- Thermal source
        PtrList<thermalSource> sources_;


    // Private Member Functions

        //- Disallow copy construct
        fluidThermalModel(const fluidThermalModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const fluidThermalModel&) = delete;


public:

    //- Runtime type information
    TypeName("fluidThermalModel");

    // Constructors

        //- Construct from dictionary
        fluidThermalModel
        (
            const dictionary& dict,
            const volScalarField& T
        );


    //- Destructor
    virtual ~fluidThermalModel() = default;


    // Member Functions

        //- Is the thermal model active?
        inline bool active() const
        {
            return lawPtr_.valid();
        }

        //- Return thermal law
        const thermalLaw& law() const
        {
            return lawPtr_();
        }

        //- Return density
        tmp<volScalarField> rho() const
        {
            return lawPtr_->rho();
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

        //- Return the aggregated source
        tmp<volScalarField> S() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
