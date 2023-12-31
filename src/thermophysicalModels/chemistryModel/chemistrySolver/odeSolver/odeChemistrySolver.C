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

\*---------------------------------------------------------------------------*/

#include "odeChemistrySolver.H"
#include "ODEChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::ode<CompType, ThermoType>::ode
(
    ODEChemistryModel<CompType, ThermoType>& model,
    const word& modelName
)
:
    chemistrySolver<CompType, ThermoType>(model, modelName),
    coeffsDict_(model.subDict(modelName + "Coeffs")),
    solverName_(coeffsDict_.lookup("ODESolver")),
    odeSolver_(ODESolver::New(solverName_, model)),
    eps_(readScalar(coeffsDict_.lookup("eps"))),
    scale_(readScalar(coeffsDict_.lookup("scale")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::ode<CompType, ThermoType>::~ode()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::scalar Foam::ode<CompType, ThermoType>::solve
(
    scalarField& c,
    const scalar T,
    const scalar p,
    const scalar t0,
    const scalar dt
) const
{
    label nSpecie = this->model_.nSpecie();
    scalarField& c1 = this->model_.coeffs();

    // copy the concentration, T and P to the total solve-vector
    for (label i=0; i<nSpecie; i++)
    {
        c1[i] = c[i];
    }
    c1[nSpecie] = T;
    c1[nSpecie+1] = p;

    scalar dtEst = dt;

    odeSolver_->solve
    (
        t0,
        t0 + dt,
        eps_,
        dtEst
    );

    for (label i=0; i<c.size(); i++)
    {
        c[i] = max(0.0, c1[i]);
    }

    return dtEst;
}


// ************************************************************************* //
