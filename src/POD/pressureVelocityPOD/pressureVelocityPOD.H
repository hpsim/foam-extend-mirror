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
    pressureVelocityPOD

Description
    POD solver for laminar Navier-Stokes equations, based on
    the coupled p-U solver formulation

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    pressureVelocityPOD.C

\*---------------------------------------------------------------------------*/

#ifndef pressureVelocityPOD_H
#define pressureVelocityPOD_H

#include "PODODE.H"
#include "PODOrthoNormalBases.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class pressureVelocityPOD Declaration
\*---------------------------------------------------------------------------*/

class pressureVelocityPOD
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

        //- POD coefficients
        scalarField coeffs_;


        //- List of valid times for snapshots
        mutable instantList* validTimesPtr_;

        //- Convection derivative coefficient matrix
        mutable PtrList<scalarSquareMatrix>* convectionDerivativePtr_;

        //- Diffusion derivative coefficient matrix
        mutable scalarSquareMatrix* diffusionDerivativePtr_;

        //- Pressure gradient derivative coefficient matrix
        mutable scalarRectangularMatrix* pressureGradDerivativePtr_;

        //- Pressure laplacian derivative coefficient matrix.
        //  Also used for source term
        mutable scalarSquareMatrix* pressureLaplaceDerivativePtr_;

        //- Pressure source derivative coefficient matrix div(U)
        mutable scalarRectangularMatrix* pressureSourceDerivativePtr_;

        //- Lagrange multiplier derivative (enforcing boundary conditions)
        mutable scalarField* lagrangeDerPtr_;

        //- Lagrange multiplier source (enforcing boundary conditions)
        mutable scalarField* lagrangeSrcPtr_;

        //- Velocity field ortho-normal base
        mutable vectorPODOrthoNormalBase* UBasePtr_;

        //- Pressure field ortho-normal base
        mutable scalarPODOrthoNormalBase* pBasePtr_;

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
        pressureVelocityPOD(const pressureVelocityPOD&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const pressureVelocityPOD&) = delete;


        //- Calculated list valid times for snapshots
        void calcValidTimes() const;

        //- Return valid times
        const instantList& validTimes() const;

        //- Calculate ortho-normal base for velocity and pressure
        void calcOrthoBase() const;

        //- Calculate coefficients
        void calcDerivativeCoeffs() const;

        //- Update reconstructed fields
        void updateFields() const;


public:

    //- Runtime type information
    TypeName("pressureVelocity");


    // Constructors

        //- Construct from components
        pressureVelocityPOD
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~pressureVelocityPOD();


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

            //- Update ODE after the solution, advancing by delta
            virtual void update(const scalar delta);


    // Orthogonalisation and fields

        //- Return velocity ortho-normal base
        const vectorPODOrthoNormalBase& UBase() const;

        //- Return pressure ortho-normal base
        const scalarPODOrthoNormalBase& pBase() const;

        //- Return flux field base
        const PtrList<surfaceScalarField>& phiBase() const;


        //- Return reconstructed velocity field
        const volVectorField& reconU() const;

        //- Return reconstructed pressure field
        const volScalarField& reconP() const;


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
