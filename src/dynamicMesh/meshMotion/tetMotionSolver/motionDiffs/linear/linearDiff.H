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
    linearDiff

Description
    Linear distance-based motion diffusion.

SourceFiles
    linearDiff.C

\*---------------------------------------------------------------------------*/

#ifndef linearDiff_H
#define linearDiff_H

#include "motionDiff.H"
#include "elementFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linearDiff Declaration
\*---------------------------------------------------------------------------*/

class linearDiff
:
    public motionDiff
{
private:

    // Private data

        //- Patches selected to base the distance on
        wordList patchNames_;

        //- Diffusion field
        elementScalarField motionGamma_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearDiff(const linearDiff&);

        //- Disallow default bitwise assignment
        void operator=(const linearDiff&);


protected:

        //- Return non constant diffusivity field
        elementScalarField& motionGamma()
        {
            return motionGamma_;
        }

public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from components
        linearDiff
        (
            const tetMotionSolver& mSolver
        );


    // Destructor

        virtual ~linearDiff();


    // Member Functions

        //- Return cell distance field
        tmp<scalarField> L() const;

        //- Return diffusivity field
        virtual const elementScalarField& motionGamma() const
        {
            return motionGamma_;
        }

        //- Correct the motion diffusivity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
