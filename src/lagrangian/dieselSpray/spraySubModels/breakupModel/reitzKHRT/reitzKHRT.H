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
    Foam::reitzKHRT

Description
    secondary breakup model which uses the Kelvin-Helmholtz
    instability theory to predict the 'stripped' droplets... and
    the Raleigh-Taylor instability as well.

\*---------------------------------------------------------------------------*/

#ifndef reitzKHRT_H
#define reitzKHRT_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class reitzKHRT Declaration
\*---------------------------------------------------------------------------*/

class reitzKHRT
:
    public breakupModel
{

private:

    // Private data

        dictionary coeffsDict_;

        // reference to gravity
        const vector& g_;

        // model constants
        scalar b0_;
        scalar b1_;
        scalar cTau_;
        scalar cRT_;
        scalar msLimit_;
        scalar weberLimit_;

public:

    //- Runtime type information
    TypeName("ReitzKHRT");


    // Constructors

        //- Construct from components
        reitzKHRT
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        ~reitzKHRT();


    // Member Operators

        void breakupParcel
        (
            parcel& parcel,
            const scalar deltaT,
            const vector& vel,
            const liquidMixture& fuels
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
