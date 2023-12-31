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
    Foam::radiation::constantAbsorptionEmission

Description
    Constant radiation absorption and emission coefficients for continuous
    phase

SourceFiles
    constantAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef radiationConstantAbsorptionEmission_H
#define radiationConstantAbsorptionEmission_H

#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class constantAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class constantAbsorptionEmission
:
    public absorptionEmissionModel
{
    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Absorption coefficient / [1/m]
        dimensionedScalar a_;

        //- Emission coefficient / [1/m]
        dimensionedScalar e_;

        //- Emission contribution / [kg/(m s^3)]
        dimensionedScalar E_;


public:

    //- Runtime type information
    TypeName("constantAbsorptionEmission");


    // Constructors

        //- Construct from components
        constantAbsorptionEmission
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Destructor
    virtual ~constantAbsorptionEmission();


    // Member Operators

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;


    // Member Functions

        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
