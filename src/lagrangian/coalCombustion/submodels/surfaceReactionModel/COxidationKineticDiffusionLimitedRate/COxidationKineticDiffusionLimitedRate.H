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
    COxidationKineticDiffusionLimitedRate

Description
    Kinetic/diffusion limited rate surface reaction model for coal parcels.
    Limited to:

        C(s) + Sb*O2 -> CO2

    where Sb is the stoichiometry of the reaction

\*---------------------------------------------------------------------------*/

#ifndef COxidationKineticDiffusionLimitedRate_H
#define COxidationKineticDiffusionLimitedRate_H

#include "SurfaceReactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
template<class CloudType>
class COxidationKineticDiffusionLimitedRate;

/*---------------------------------------------------------------------------*\
            Class COxidationKineticDiffusionLimitedRate Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class COxidationKineticDiffusionLimitedRate
:
    public SurfaceReactionModel<CloudType>
{
    // Private data

        // Model constants

            //- Stoichiometry of reaction
            const scalar Sb_;

            //- Mass diffusion limited rate constant, C1
            const scalar C1_;

            //- Kinetics limited rate pre-exponential constant, C2
            const scalar C2_;

            //- Kinetics limited rate activation energy
            const scalar E_;


        // Addressing

            //- Cs positions in global/local lists
            label CsLocalId_;

            //- O2 position in global list
            label O2GlobalId_;

            //- CO2 positions in global list
            label CO2GlobalId_;


        // Local copies of thermo properties

            //- Molecular weight of C [kg/kmol]
            scalar WC_;

            //- Molecular weight of O2 [kg/kmol]
            scalar WO2_;

            //- Chemical enthalpy of CO2 [J/kg]
            scalar HcCO2_;


public:

    //- Runtime type information
    TypeName("COxidationKineticDiffusionLimitedRate");


    // Constructors

        //- Construct from dictionary
        COxidationKineticDiffusionLimitedRate
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~COxidationKineticDiffusionLimitedRate();


    // Member Functions

        //- Flag to indicate whether model activates surface reaction model
        virtual bool active() const;

        //- Update surface reactions
        virtual scalar calculate
        (
            const scalar dt,
            const label cellI,
            const scalar d,
            const scalar T,
            const scalar Tc,
            const scalar pc,
            const scalar rhoc,
            const scalar mass,
            const scalarField& YGas,
            const scalarField& YLiquid,
            const scalarField& YSolid,
            const scalarField& YMixture,
            const scalar N,
            scalarField& dMassGas,
            scalarField& dMassLiquid,
            scalarField& dMassSolid,
            scalarField& dMassSRCarrier
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "COxidationKineticDiffusionLimitedRate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
