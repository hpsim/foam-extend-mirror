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
    Foam::DualCompetingRateDevolatilisation

Group
    grpLagrangianIntermediateDevolatilisationSubModels

Description
    Double competing rate kinetic rate devolatisation model.
    - Rate given by Arrhenius eqn

        kappa1 = A1.exp(- E1/R.T)
        kappa2 = A1.exp(- E1/R.T)

    - Competing reaction rate

       r = min((a1*k1 + a2*k2)

    where:
        a1, a2    = stoichiometric coefficient
        A1, A2    = activation energy (user input)
        E1, E2    = pre-exponential factor (user input)
        R         = universal gas constant
        T         = temperature

    Usage:

        dualCompetingRateDevolatilisationCoeffs
        {
            volatileData
            (
                (CH4     0.l5    1e3     1e5
             1.00    1e5    1e5)   // (name a1 A1 E1 a2 A2 E2)
                (CO2     0.l5    1e3    1e5
             1.00    1e5    1e5)   // (name a1 A1 E1 a2 A2 E2)
            );

            volatileResidualCoeff 1e-6;
        }

    NOTE:

    R uses kmol insteasd of mol.

\*---------------------------------------------------------------------------*/

#ifndef DualCompetingRateDevolatilisation_H
#define DualCompetingRateDevolatilisation_H

#include "DevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class DualCompetingRateDevolatilisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DualCompetingRateDevolatilisation
:
    public DevolatilisationModel<CloudType>
{
    // Private data

        // Model constants

            //- Stoichiometric coefficient - 1
            const scalar a1_;

            //- Pre-exponential factor - 1
            const scalar A1_;

            //- Activiation Energy - 1
            const scalar E1_;

            //- Stoichiometric coefficient - 2
            const scalar a2_;

            //- Pre-exponential factor - 2
            const scalar A2_;

            //- Activiation Energy - 2
            const scalar E2_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar volatileResidualCoeff_;


public:

    //- Runtime type information
    TypeName("DualCompetingRateDevolatilisation");


    // Constructors

        //- Construct from dictionary
        DualCompetingRateDevolatilisation
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~DualCompetingRateDevolatilisation() = default;


    // Member Functions

        //- Flag to indicate whether model activates devolatilisation model
        virtual bool active() const;

        //- Update model
        virtual scalar calculate
        (
            const scalar dt,
            const scalar mass0,
            const scalar mass,
            const scalar T,
            const scalar YVolatile0,
            const scalar YVolatile,
            bool& canCombust
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DualCompetingRateDevolatilisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
