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
    Foam::potentialControl

Description
    potential control class to supply convergence information/checks for
    the potential flow solver.

\*---------------------------------------------------------------------------*/

#ifndef potentialControl_H
#define potentialControl_H

#include "solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class potentialControl Declaration
\*---------------------------------------------------------------------------*/

class potentialControl
:
    public solutionControl
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        potentialControl(const potentialControl&);

        //- Disallow default bitwise assignment
        void operator=(const potentialControl&);


    // Protected Member Functions

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();


public:


    // Static Data Members

        //- Run-time type information
        TypeName("potentialControl");


    // Constructors

        //- Construct from mesh
        potentialControl(fvMesh& mesh);


    //- Destructor
    virtual ~potentialControl();


    // Member Functions

        // Solution control

            //- Loop loop
            virtual bool loop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
