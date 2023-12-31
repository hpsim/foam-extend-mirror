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
    Foam::ode

Description
    An ODE solver for chemistry

SourceFiles
    odeChemistrySolver.C

\*---------------------------------------------------------------------------*/

#ifndef odeChemistrySolver_H
#define odeChemistrySolver_H

#include "chemistrySolver.H"
#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CompType, class ThermoType>
class ode;

/*---------------------------------------------------------------------------*\
                            Class ode Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class ode
:
    public chemistrySolver<CompType, ThermoType>
{
    // Private data

        dictionary coeffsDict_;
        const word solverName_;
        autoPtr<ODESolver> odeSolver_;

        // Model constants

            scalar eps_;
            scalar scale_;


public:

    //- Runtime type information
    TypeName("ODE");


    // Constructors

        //- Construct from components
        ode
        (
            ODEChemistryModel<CompType, ThermoType>& model,
            const word& modelName
        );


    //- Destructor
    virtual ~ode();


    // Member Functions

        scalar solve
        (
            scalarField& c,
            const scalar T,
            const scalar p,
            const scalar t0,
            const scalar dt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "odeChemistrySolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
