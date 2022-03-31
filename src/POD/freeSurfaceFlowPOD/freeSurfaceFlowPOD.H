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
    freeSurfaceFlowPOD

Description
    POD solver for laminar free surface flows using modified U-p decomposition

    M. Bergmann and C.-H. Bruneau and A. Iollo
    Enablers for robust POD models
    Journal of Computational Physics v228 n2, pp 516-538
    doi: 10.1016/j.jcp.2008.09.024

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    freeSurfaceFlowPOD.C

\*---------------------------------------------------------------------------*/

#ifndef freeSurfaceFlowPOD_H
#define freeSurfaceFlowPOD_H

#include "PODODE.H"
#include "PODOrthoNormalBases.H"
#include "VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class freeSurfaceFlowPOD Declaration
\*---------------------------------------------------------------------------*/

class freeSurfaceFlowPOD
:
    public PODODE
{
    // Private data

        //- Velocity field name
        word UName_;

        //- Pressure field name (dynamic pressure)
        word pdName_;

        //- VOF field name
        word alphaName_;

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

        //- Ortho-normal base pointer
        mutable PODOrthoNormalBase<vector4>* orthoBasePtr_;

        //- Pressure field base
        mutable PtrList<volVectorField>* UBasePtr_;

        //- Pressure field base
        mutable PtrList<volScalarField>* pdBasePtr_;

        //- VOF field base
        mutable PtrList<volScalarField>* alphaBasePtr_;

        //- Flux field base
        mutable PtrList<surfaceScalarField>* phiBasePtr_;

        //- Reconstructed velocity field pointer
        mutable volVectorField* reconUPtr_;

        //- Reconstructed pressure field pointer
        mutable volScalarField* reconPdPtr_;

        //- Reconstructed VOF field pointer
        mutable volScalarField* reconAlphaPtr_;

        //- Time index of current field update
        mutable label fieldUpdateTimeIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        freeSurfaceFlowPOD(const freeSurfaceFlowPOD&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const freeSurfaceFlowPOD&) = delete;


        //- Calculated list valid times for snapshots
        void calcValidTimes() const;

        //- Return valid times
        const instantList& validTimes() const;

        //- Calculate ortho-normal base
        void calcOrthoBase() const;

        //- Calculate coefficients
        void calcDerivativeCoeffs() const;

        //- Update reconstructed fields
        void updateFields() const;


public:

    //- Runtime type information
    TypeName("freeSurfaceFlow");


    // Constructors

        //- Construct from components
        freeSurfaceFlowPOD
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~freeSurfaceFlowPOD();


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

            //- Return ortho-normal base
            const PODOrthoNormalBase<vector4>& orthoBase() const;

            //- Return velocity field base
            const PtrList<volVectorField>& UBase() const;

            //- Return pressure field base
            const PtrList<volScalarField>& pdBase() const;

            //- Return VOF field base
            const PtrList<volScalarField>& alphaBase() const;

            //- Return flux field base
            const PtrList<surfaceScalarField>& phiBase() const;


            //- Return reconstructed velocity field
            const volVectorField& reconU() const;

            //- Return reconstructed pressure field
            const volScalarField& reconPd() const;

            //- Return reconstructed VOF field
            const volScalarField& reconAlpha() const;


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
