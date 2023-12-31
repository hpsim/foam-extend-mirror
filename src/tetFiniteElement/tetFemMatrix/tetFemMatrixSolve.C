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

Description
    Tetrahedral Finite Element matrix basic solvers.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Solvers * * * * * * * * * * * * * * * * * //

template<class Type>
lduSolverPerformance tetFemMatrix<Type>::solve
(
     const dictionary& solverControls
)
{
    if (debug)
    {
        Info<< "tetFemMatrix<Type>::solve(const dictionary&) : "
               "solving tetFemMatrix<Type>"
            << endl;
    }

    // Check the matrix
    if (debug > 1)
    {
        this->check();
    }

    BlockSolverPerformance<Type> solverPerfVec
    (
        "tetFemMatrix<Type>::solve",
        psi_.name()
    );

    // Add boundary source for gradient-type conditions
    addBoundarySourceDiag();

    // Store the boundary coefficients for insertion of boundary conditions
    storeBoundaryCoeffs();

    typename Type::labelType validComponents
    (
        pow
        (
            psi_.mesh()().solutionD(),
            pTraits<typename powProduct<Vector<label>, Type::rank>::type>::zero
        )
    );

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    // Cast into a non-const to solve.  HJ, 6/May/2016
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& psi =
        const_cast<GeometricField<Type, tetPolyPatchField, tetPointMesh>&>
        (
            psi_
        );

    for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1) continue;

        scalarField psiCmpt = psi_.internalField().component(cmpt);
        scalarField sourceCmpt = source_.component(cmpt);

        // Set component boundary conditions
        setComponentBoundaryConditions(cmpt, psiCmpt, sourceCmpt);

        // Add the coupling coefficients
        addCouplingCoeffs();

        addCouplingSource(sourceCmpt);

        // Prepare for coupled interface update
        FieldField<Field, scalar> coupledBouCoeffs
        (
            psi_.boundaryField().size()
        );

        FieldField<Field, scalar> coupledIntCoeffs
        (
            psi_.boundaryField().size()
        );

        forAll(psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

            coupledBouCoeffs.set
            (
                patchI,
                ptf.cutBouCoeffs(*this)
            );

            coupledIntCoeffs.set
            (
                patchI,
                ptf.cutIntCoeffs(*this)
            );
        }

        eliminateCouplingCoeffs();

        scalarField res(psi_.size(), 0);
        lduMatrix::residual
        (
            res,
            psiCmpt,
            sourceCmpt,
            coupledBouCoeffs,
            interfaces,
            cmpt
        );

        lduSolverPerformance solverPerf = lduSolver::New
        (
            psi_.name() + pTraits<Type>::componentNames[cmpt],
            *this,
            coupledBouCoeffs,
            coupledIntCoeffs,
            interfaces,
            solverControls
        )->solve(psiCmpt, sourceCmpt, cmpt);

        solverPerf.print();

        solverPerfVec.replace(cmpt, solverPerf);

        psi.internalField().replace(cmpt, psiCmpt);

        reconstructMatrix();
    }

    if (debug)
    {
        Info<< "tetFemMatrix<Type>::solve : correcting boundary conditions"
            << endl;
    }

    psi.correctBoundaryConditions();

    psi_.mesh().solutionDict().setSolverPerformance(psi_.name(), solverPerfVec);

    return solverPerfVec.max();
}


template<class Type>
lduSolverPerformance tetFemMatrix<Type>::solve()
{
    return solve(psi_.mesh().solutionDict().solver(psi_.name()));
}


// Return the matrix residual
template<class Type>
tmp<Field<Type> > tetFemMatrix<Type>::residual()
{
    tmp<Field<Type> > tres(psi_.size());

    // Store the boundary coefficients for insertion of boundary conditions
    storeBoundaryCoeffs();

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    // Loop over field components
    for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
    {
        scalarField PsiInternalCmpt = psi_.internalField().component(cmpt);
        scalarField sourceCmpt = source_.component(cmpt);

        setComponentBoundaryConditions(cmpt, PsiInternalCmpt, sourceCmpt);

        // Add the coupling coefficients
        addCouplingCoeffs();
        addCouplingSource(sourceCmpt);

        // Prepare for coupled interface update
        FieldField<Field, scalar> coupledBouCoeffs(psi_.boundaryField().size());

        forAll(psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];
            coupledBouCoeffs.set
            (
                patchI,
                new scalarField(ptf.cutBouCoeffs(*this))
            );
        }

        eliminateCouplingCoeffs();

        tres().replace
        (
            cmpt,
            lduMatrix::residual
            (
                psi_.internalField().component(cmpt),
                sourceCmpt,
                coupledBouCoeffs,
                interfaces,
                cmpt
            )
        );

        reconstructMatrix();
    }

    return tres;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
