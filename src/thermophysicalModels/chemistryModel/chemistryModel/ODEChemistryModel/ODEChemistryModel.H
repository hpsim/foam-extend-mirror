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
    Foam::ODEChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    ODEChemistryModelI.H
    ODEChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef ODEChemistryModel_H
#define ODEChemistryModel_H

#include "Reaction.H"
#include "ODE.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class objectRegistry;

template<class CompType, class ThermoType>
class chemistrySolver;

/*---------------------------------------------------------------------------*\
                      Class ODEChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class ODEChemistryModel
:
    public CompType,
    public ODE
{
    // Private data

        //- Reference to the field of specie mass fractions
        PtrList<volScalarField>& Y_;

        //- Reactions
        const PtrList<Reaction<ThermoType> >& reactions_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermo_;

        //- Number of species
        label nSpecie_;

        //- Number of reactions
        label nReaction_;

        //- Chemistry solver
        autoPtr<chemistrySolver<CompType, ThermoType> > solver_;

        //- Chemical source term [kg/m3/s]
        PtrList<scalarField> RR_;

        //- ODE coefficients
        scalarField coeffs_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const ODEChemistryModel&);


protected:

    // Protected Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<scalarField>& RR();


public:

    //- Runtime type information
    TypeName("ODEChemistryModel");


    // Constructors

        //- Construct from components
        ODEChemistryModel
        (
            const fvMesh& mesh,
            const objectRegistry& obj,
            const word& compTypeName,
            const word& thermoTypeName
        );


    //- Destructor
    virtual ~ODEChemistryModel();


    // Member Functions

        //- Access the fractions of the species
        inline const PtrList<volScalarField>& Y() const;

        //- The reactions
        inline const PtrList<Reaction<ThermoType> >& reactions() const;

        //- Thermodynamic data of the species
        inline const PtrList<ThermoType>& specieThermo() const;

        //- The number of species
        inline label nSpecie() const;

        //- The number of reactions
        inline label nReaction() const;

        //- Return the chemisty solver
        inline const chemistrySolver<CompType, ThermoType>& solver() const;

        //- dc/dt = omega, rate of change in concentration, for each species
        virtual scalarField omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p
        ) const;

        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
        (
            const Reaction<ThermoType>& r,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;

        //- Calculates the reaction rates
        virtual void calculate();


        // Chemistry model functions (overriding abstract functions in
        // basicChemistryModel.H)

            //- Return const access to the chemical source terms
            inline tmp<volScalarField> RR(const label i) const;

            //- Solve the reaction system for the given start time and time
            //  step and return the characteristic time
            virtual scalar solve(const scalar t0, const scalar deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return source for enthalpy equation [kg/m/s3]
            virtual tmp<volScalarField> Sh() const;

            //- Return the heat release, i.e. enthalpy/sec [m2/s3]
            virtual tmp<volScalarField> dQ() const;


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            virtual label nEqns() const;

            //- Return access to coefficients
            virtual scalarField& coeffs();

            //- Return reference to coefficients
            virtual const scalarField& coeffs() const;

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt,
                scalarSquareMatrix& dfdc
            ) const;

            //- Update ODE after the solution, advancing by delta
            virtual void update(const scalar delta)
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ODEChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ODEChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
