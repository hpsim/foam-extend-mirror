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
    Foam::granularPressureModel

SourceFiles
    granularPressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef granularPressureModel_H
#define granularPressureModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class granularPressureModel Declaration
\*---------------------------------------------------------------------------*/

class granularPressureModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        granularPressureModel(const granularPressureModel&);

        //- Disallow default bitwise assignment
        void operator=(const granularPressureModel&);


protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("granularPressureModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        granularPressureModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        granularPressureModel(const dictionary& dict);


    // Selectors

        static autoPtr<granularPressureModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~granularPressureModel();


    // Member Functions

        //- Granular pressure coefficient
        virtual tmp<volScalarField> granularPressureCoeff
        (
            const volScalarField& alpha,
            const volScalarField& g0,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& e
        ) const = 0;

        //- Derivative of the granular pressure coefficient
        virtual tmp<volScalarField> granularPressureCoeffPrime
        (
            const volScalarField& alpha,
            const volScalarField& g0,
            const volScalarField& g0prime,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& e
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
