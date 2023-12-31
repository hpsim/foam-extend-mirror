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

Namespace
    Foam::compressible::LESModels

Description
    Namespace for compressible LES models.


Class
    Foam::compressible::LESModel

Description
    Base class for all compressible flow LES SGS models.

    This class defines the basic interface for a compressible flow SGS
    model, and encapsulates data of value to all possible models.
    In particular this includes references to all the dependent fields
    (rho, U, phi), the physical viscosity mu, and the LESProperties
    dictionary, which contains the model selection and model coefficients.

SourceFiles
    LESModel.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleLESModel_H
#define compressibleLESModel_H

#include "compressible/turbulenceModel/turbulenceModel.H"
#include "LESdelta.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "bound.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
                           Class LESModel Declaration
\*---------------------------------------------------------------------------*/

class LESModel
:
    public turbulenceModel,
    public IOdictionary
{

protected:

    // Protected data

        Switch printCoeffs_;
        dictionary coeffDict_;

        dimensionedScalar k0_;

        autoPtr<LESdelta> delta_;


    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        LESModel(const LESModel&);

        //- Disallow default bitwise assignment
        LESModel& operator=(const LESModel&);


public:

    //- Runtime type information
    TypeName("LESModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LESModel,
            dictionary,
            (
                const volScalarField& rho,
                const volVectorField& U,
                const surfaceScalarField& phi,
                const basicThermo& thermoPhysicalModel,
                const word& turbulenceModelName
            ),
            (rho, U, phi, thermoPhysicalModel, turbulenceModelName)
        );


    // Constructors

        //- Construct from components
        LESModel
        (
            const word& type,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermoPhysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName
        );


    // Selectors

        //- Return a reference to the selected LES model
        static autoPtr<LESModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermoPhysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName
        );


    //- Destructor

        virtual ~LESModel()
        {}


    // Member Functions

        // Access

            //- Const access to the coefficients dictionary,
            //  which provides info. about choice of models,
            //  and all related data (particularly model coefficients).
            inline const dictionary& coeffDict() const
            {
                return coeffDict_;
            }

            //- Return the value of k0 which k is not allowed to be less than
            const dimensionedScalar& k0() const
            {
                return k0_;
            }

            //- Allow k0 to be changed
            dimensionedScalar& k0()
            {
                return k0_;
            }

            //- Access function to filter width
            inline const volScalarField& delta() const
            {
                return delta_();
            }


        //- Return the SGS turbulent kinetic energy.
        virtual tmp<volScalarField> k() const = 0;

        //- Return the SGS turbulent dissipation.
        virtual tmp<volScalarField> epsilon() const = 0;

        //- Return the SGS turbulent viscosity
        virtual tmp<volScalarField> muSgs() const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> muEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("muEff", muSgs() + mu())
            );
        }

        //- Return the SGS turbulent thermal diffusivity
        virtual tmp<volScalarField> alphaSgs() const = 0;

        //- Return the SGS thermal conductivity.
        virtual tmp<volScalarField> alphaEff() const = 0;

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const = 0;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devRhoBeff() const = 0;

        //- Returns div(rho*dev(B)).
        // This is the additional term due to the filtering of the NSE.
        virtual tmp<fvVectorMatrix> divDevRhoBeff() const = 0;


        // RAS compatibility functions for the turbulenceModel base class

            //- Return the turbulence viscosity
            virtual tmp<volScalarField> mut() const
            {
                return muSgs();
            }

            //- Return the turbulence thermal diffusivity
            virtual tmp<volScalarField> alphat() const
            {
                return alphaSgs();
            }

            //- Return the Reynolds stress tensor
            virtual tmp<volSymmTensorField> R() const
            {
                return B();
            }

            //- Return the effective stress tensor including the laminar stress
            virtual tmp<volSymmTensorField> devRhoReff() const
            {
                return devRhoBeff();
            }

            //- Return the source term for the momentum equation
            virtual tmp<fvVectorMatrix> divDevRhoReff() const
            {
                return divDevRhoBeff();
            }


        //- Correct Eddy-Viscosity and related properties.
        //  This calls correct(const tmp<volTensorField>& gradU) by supplying
        //  gradU calculated locally.
        void correct();

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
