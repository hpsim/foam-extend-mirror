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
    Foam::pairPotential

Description

SourceFiles
    pairPotential.C
    newPairPotential.C

\*---------------------------------------------------------------------------*/

#ifndef pairPotential_H
#define pairPotential_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "scalarList.H"
#include "Pair.H"
#include "Switch.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class energyScalingFunction;

/*---------------------------------------------------------------------------*\
                     Class pairPotential Declaration
\*---------------------------------------------------------------------------*/

class pairPotential
{

protected:

    // Protected data

        word name_;
        dictionary pairPotentialProperties_;

        scalar rCut_;
        scalar rCutSqr_;

        scalar rMin_;
        scalar dr_;

        scalarList forceLookup_;
        scalarList energyLookup_;

        mutable energyScalingFunction* esfPtr_;

        bool writeTables_;


    // Private Member Functions

        void scaleEnergy(scalar& e, const scalar r) const;

        //- Disallow copy construct
        pairPotential(const pairPotential&);

        //- Disallow default bitwise assignment
        void operator=(const pairPotential&);


public:

    //- Runtime type information
        TypeName("pairPotential");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            pairPotential,
            dictionary,
            (
                const word& name,
                const dictionary& pairPotentialProperties
            ),
            (name, pairPotentialProperties)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<pairPotential> New
        (
            const word& name,
            const dictionary& pairPotentialProperties
        );


    // Constructors

        //- Construct from components
        pairPotential
        (
            const word& name,
            const dictionary& pairPotentialProperties
        );


    // Destructor

        virtual ~pairPotential()
        {}


    // Member Functions

        void setLookupTables();

        inline scalar rMin() const;

        inline scalar dr() const;

        inline scalar rCut() const;

        inline scalar rCutSqr() const;

        scalar energy (const scalar r) const;

        scalar force (const scalar r) const;

        List<Pair<scalar> > energyTable() const;

        List<Pair<scalar> > forceTable() const;

        inline bool writeTables() const;

        virtual scalar unscaledEnergy(const scalar r) const = 0;

        scalar scaledEnergy(const scalar r) const;

        scalar energyDerivative
        (
            const scalar r,
            const bool scaledEnergyDerivative = true
        ) const;

        const dictionary& pairPotentialProperties() const
        {
            return pairPotentialProperties_;
        }

        bool writeEnergyAndForceTables(Ostream& os) const;

        //- Read pairPotential dictionary
        virtual bool read(const dictionary& pairPotentialProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pairPotentialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
