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
    Finite-Area matrix basic solvers.

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set reference level for a component of the solution
// on a given patch face
template<class Type>
void faMatrix<Type>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction cmpt,
    const scalar value
)
{
    internalCoeffs_[patchi][facei].component(cmpt) +=
        diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]];

    boundaryCoeffs_[patchi][facei].component(cmpt) +=
        diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]]*value;
}


template<class Type>
lduSolverPerformance faMatrix<Type>::solve(const dictionary& solverControls)
{
    if (debug)
    {
        Info<< "faMatrix<Type>::solve(const dictionary&) : "
               "solving faMatrix<Type>"
            << endl;
    }

    BlockSolverPerformance<Type> solverPerfVec
    (
        "faMatrix<Type>::solve",
        psi_.name()
    );

    scalarField saveDiag = diag();

    Field<Type> source = source_;
    addBoundarySource(source);

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    // Cast into a non-const to solve.  HJ, 6/May/2016
    GeometricField<Type, faPatchField, areaMesh>& psi =
       const_cast<GeometricField<Type, faPatchField, areaMesh>&>(psi_);

    for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
    {
        // copy field and source

        scalarField psiCmpt = psi_.internalField().component(cmpt);
        addBoundaryDiag(diag(), cmpt);

        scalarField sourceCmpt = source.component(cmpt);

        FieldField<Field, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(cmpt)
        );

        FieldField<Field, scalar> intCoeffsCmpt
        (
            internalCoeffs_.component(cmpt)
        );

        // Use the initMatrixInterfaces and updateMatrixInterfaces to correct
        // bouCoeffsCmpt for the explicit part of the coupled boundary
        // conditions
        initMatrixInterfaces
        (
            bouCoeffsCmpt,
            interfaces,
            psiCmpt,
            sourceCmpt,
            cmpt
        );

        updateMatrixInterfaces
        (
            bouCoeffsCmpt,
            interfaces,
            psiCmpt,
            sourceCmpt,
            cmpt
        );

        lduSolverPerformance solverPerf;

        // Solver call
        solverPerf = lduSolver::New
        (
            psi_.name() + pTraits<Type>::componentNames[cmpt],
            *this,
            bouCoeffsCmpt,
            intCoeffsCmpt,
            interfaces,
            solverControls
        )->solve(psiCmpt, sourceCmpt, cmpt);

        solverPerf.print();

        solverPerfVec.replace(cmpt, solverPerf);

        psi.internalField().replace(cmpt, psiCmpt);
        diag() = saveDiag;
    }

    psi.correctBoundaryConditions();

    psi_.mesh().solutionDict().setSolverPerformance(psi_.name(), solverPerfVec);

    return solverPerfVec.max();
}


template<class Type>
lduSolverPerformance faMatrix<Type>::faSolver::solve()
{
    return solvei
    (
        faMat_.psi().mesh().solutionDict().solverDict
        (
            faMat_.psi().name()
        )
    );
}


template<class Type>
lduSolverPerformance faMatrix<Type>::solve()
{
    return solve
    (
        this->psi().mesh().solutionDict().solverDict
        (
            this->psi().name()
        )
    );
}


// Return the matrix residual
template<class Type>
tmp<Field<Type> > faMatrix<Type>::residual() const
{
    // Bug fix: Creating a tmp out of a const reference will change the field
    // HJ, 15/Apr/2011 and Petr Vita, Mark Olesen, 10/Sep/2020
    tmp<Field<Type> > tres(new Field<Type>(source_));
    Field<Type>& res = tres();

    addBoundarySource(res);

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    // Loop over field components
    for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
    {
        scalarField psiCmpt = psi_.internalField().component(cmpt);

        scalarField boundaryDiagCmpt(psi_.size(), 0.0);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);

        FieldField<Field, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(cmpt)
        );

        res.replace
        (
            cmpt,
            lduMatrix::residual
            (
                psiCmpt,
                res.component(cmpt) - boundaryDiagCmpt*psiCmpt,
                bouCoeffsCmpt,
                interfaces,
                cmpt
            )
        );
    }

    return tres;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
