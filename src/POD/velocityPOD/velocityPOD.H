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
    velocityPOD

Description
    POD solver for laminar Navier-Stokes equations, based on
    velocity decomposition

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    velocityPOD.C

\*---------------------------------------------------------------------------*/

#ifndef velocityPOD_H
#define velocityPOD_H

#include "PODODE.H"
#include "PODOrthoNormalBases.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class velocityPOD Declaration
\*---------------------------------------------------------------------------*/

class velocityPOD
:
    public PODODE
{
    // Private data

        //- Velocity field name
        word UName_;

        //- Pressure field name
        word pName_;

        //- Flux field name
        word phiName_;

        //- Transport properties
        IOdictionary transportProperties_;

        //- Kinematic viscosity
        dimensionedScalar nu_;

        //- Lagrange multiplier pre-factor
        scalar beta_;

        //- Use zero field as snapshot
        Switch useZeroField_;

        //- Drift correction
        Switch driftCorrection_;

        //- POD coefficients
        scalarField coeffs_;


        //- List of valid times for snapshots
        mutable instantList* validTimesPtr_;

        //- Convection derivative coefficient matrix
        mutable PtrList<scalarSquareMatrix>* convectionDerivativePtr_;

        //- Diffusion and pressure derivative coefficient matrix
        mutable scalarSquareMatrix* derivativePtr_;

        //- Lagrange multiplier derivative (enforcing boundary conditions)
        mutable scalarSquareMatrix* lagrangeDerPtr_;

        //- Lagrange multiplier source (enforcing boundary conditions)
        mutable scalarField* lagrangeSrcPtr_;

        //- Ortho-normal base pointer, carrying velocity fields
        mutable vectorPODOrthoNormalBase* orthoBasePtr_;

        //- Pressure field base
        mutable PtrList<volScalarField>* pBasePtr_;

        //- Flux field base
        mutable PtrList<surfaceScalarField>* phiBasePtr_;

        //- Reconstructed velocity field pointer
        mutable volVectorField* reconUPtr_;

        //- Reconstructed pressure field pointer
        mutable volScalarField* reconPPtr_;

        //- Time index of current field update
        mutable label fieldUpdateTimeIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        velocityPOD(const velocityPOD&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const velocityPOD&) = delete;


        //- Calculated list valid times for snapshots
        void calcValidTimes() const;

        //- Return valid times
        const instantList& validTimes() const;

        //- Calculate ortho-normal base for velocity
        void calcOrthoBase() const;

        //- Calculate coefficients
        void calcDerivativeCoeffs() const;

        //- Update reconstructed fields
        void updateFields() const;


public:

    //- Runtime type information
    TypeName("velocity");


    // Constructors

        //- Construct from components
        velocityPOD
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~velocityPOD();


    // Member Functions

        // ODE parameters

            //- Return number of equations
            virtual label nEqns() const;

            //- Return reference to interpolation coefficients
            virtual scalarField& coeffs();

            //- Return reference to interpolation coefficients
            virtual const scalarField& coeffs() const;

            //- Return derivatives
            virtual void derivatives
            (
                const scalar x,
                const scalarField& y,
                scalarField& dydx
            ) const;

            //- Return Jacobian
            virtual void jacobian
            (
                const scalar x,
                const scalarField& y,
                scalarField& dfdx,
                scalarSquareMatrix& dfdy
            ) const;


        // Orthogonalisation and fields

            //- Return ortho-normal base: velocity fields
            const vectorPODOrthoNormalBase& orthoBase() const;

            //- Return pressure field base
            const PtrList<volScalarField>& pBase() const;

            //- Return flux field base
            const PtrList<surfaceScalarField>& phiBase() const;


            //- Return reconstructed velocity field
            const volVectorField& reconU() const;

            //- Return reconstructed pressure field
            const volScalarField& reconP() const;


        // Update

            //- Update ODE after the solution, advancing by delta
            virtual void update(const scalar delta);


        // Write

            //- Write snapshots
            virtual void writeSnapshots() const;

            //- Write reconstructed fields
            virtual void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
