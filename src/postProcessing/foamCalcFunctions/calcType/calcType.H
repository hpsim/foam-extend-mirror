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

Namespace
    Foam::calcTypes

Description
    Namespace for post-processing calculation functions


Class
    Foam::calcType

Description
    Base class for post-processing calculation functions

SourceFiles
    calcType.C

\*---------------------------------------------------------------------------*/

#ifndef calcType_H
#define calcType_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class calcType Declaration
\*---------------------------------------------------------------------------*/

class calcType
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        calcType(const calcType&);

        //- Disallow default bitwise assignment
        void operator=(const calcType&);


protected:

    // Protected member functions

        // Calculation routines

            //- Initialise - typically setting static variables,
            //  e.g. command line arguments
            virtual void init();

            //- Pre-time loop calculations
            virtual void preCalc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );

            //- Time loop calculations
            virtual void calc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );

            //- Post-time loop calculations
            virtual void postCalc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );


public:

    //- Runtime type information
    TypeName("calcType");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            calcType,
            dictionary,
            (),
            ()
        );


    // Constructors

        //- Construct null
        calcType();


    // Selectors

        static autoPtr<calcType> New(const word& calcTypeName);


    // Destructor

        virtual ~calcType();


    // Member Functions

        // Calculation routines - wrapped by exception handling loop

            //- Initialise - typically setting static variables,
            //  e.g. command line arguments
            void tryInit();

            //- Pre-time loop calculations
            void tryPreCalc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );

            //- Time loop calculations
            void tryCalc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );

            //- Post-time loop calculations
            void tryPostCalc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
