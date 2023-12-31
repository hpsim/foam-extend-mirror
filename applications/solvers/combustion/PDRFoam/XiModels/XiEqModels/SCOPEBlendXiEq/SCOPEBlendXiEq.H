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
    Foam::XiEqModels::SCOPEBlend

Description
    Simple SCOPEBlendXiEq model for XiEq based on SCOPEXiEqs correlation
    with a linear correction function to give a plausible profile for XiEq.
    See @link SCOPELaminarFlameSpeed.H @endlink for details on the SCOPE
    laminar flame speed model.

SourceFiles
    SCOPEBlend.C

\*---------------------------------------------------------------------------*/

#ifndef SCOPEBlend_H
#define SCOPEBlend_H

#include "XiEqModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{

/*---------------------------------------------------------------------------*\
                          Class SCOPEBlend Declaration
\*---------------------------------------------------------------------------*/

class SCOPEBlend
:
    public XiEqModel
{
    // Private data

        //- Low turbulence intensity equilibrium Xi model
        autoPtr<XiEqModel> XiEqModelL_;

        //- High turbulence intensity equilibrium Xi model
        autoPtr<XiEqModel> XiEqModelH_;


    // Private Member Functions

        //- Disallow copy construct
        SCOPEBlend(const SCOPEBlend&);

        //- Disallow default bitwise assignment
        void operator=(const SCOPEBlend&);


public:

    //- Runtime type information
    TypeName("SCOPEBlend");


    // Constructors

        //- Construct from components
        SCOPEBlend
        (
            const dictionary& XiEqProperties,
            const hhuCombustionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    // Destructor

        virtual ~SCOPEBlend();


    // Member Functions

        //- Return the flame-wrinking XiEq
        virtual tmp<volScalarField> XiEq() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiEqProperties)
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiEqModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
