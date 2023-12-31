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
    Foam::SingleMixtureFraction

Description
    Templated parcel multi-phase, multi-component class

SourceFiles
    SingleMixtureFraction.C

\*---------------------------------------------------------------------------*/

#ifndef SingleMixtureFraction_H
#define SingleMixtureFraction_H

#include "CompositionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SingleMixtureFraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SingleMixtureFraction
:
    public CompositionModel<CloudType>
{
    // Private data

        // Indices of the phases

            //- Gas
            label idGas_;

            //- Liquid
            label idLiquid_;

            //- Solid
            label idSolid_;


       // Mixture properties

            //- Phase component total fractions
            scalarField YMixture0_;


    // Private member functions

        //- Construct the indices and check correct specification of
        //  1 gas, 1 liquid and 1 solid
        void constructIds();


public:

    //- Runtime type information
    TypeName("SingleMixtureFraction");


    // Constructors

        //- Construct from dictionary
        SingleMixtureFraction(const dictionary& dict, CloudType& owner);


    //- Destructor
    virtual ~SingleMixtureFraction();


    // Member Functions

        // Access

            // Gas properties

                //- Return the list of gas mass fractions
                const scalarField& YGas0() const;

                //- Return the total gas mass fraction
                scalar YGasTot0() const;


            // Liquid properties

                //- Return the list of liquid mass fractions
                const scalarField& YLiquid0() const;

                //- Return the total liquid mass fraction
                scalar YLiquidTot0() const;


            // Solid properties

                //- Return the list of solid mass fractions
                const scalarField& YSolid0() const;

                //- Return the total solid mass fraction
                scalar YSolidTot0() const;


            // Mixture properties

                //- Return the list of mixture mass fractions
                virtual const scalarField& YMixture0() const;

                // Indices of gas, liquid and solid phases in phase properties
                // list

                    //- Gas id
                    virtual label idGas() const;

                    //- Liquid id
                    virtual label idLiquid() const;

                    //- Solid id
                    virtual label idSolid() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SingleMixtureFraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
