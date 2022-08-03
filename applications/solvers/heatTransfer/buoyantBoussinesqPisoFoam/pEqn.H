{
    volScalarField rUA("rUA", 1.0/UEqn.A());
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

    U = rUA*UEqn.H();

    phi = fvc::interpolate(U) & mesh.Sf();
    adjustPhi(phi, U, p_rgh);

    surfaceScalarField buoyancyPhi(rUAf*ghf*fvc::snGrad(rhok)*mesh.magSf());
    phi -= buoyancyPhi;

    while (piso.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rUAf, p_rgh) == fvc::div(phi)
        );

        p_rghEqn.solve
        (
            mesh.solutionDict().solver(p_rgh.select(piso.finalInnerIter()))
        );

        if (piso.finalNonOrthogonalIter())
        {
            phi -= p_rghEqn.flux();
        }

        // Correct the momentum source with the pressure gradient flux
        U -= rUA*fvc::reconstruct((buoyancyPhi + p_rghEqn.flux())/rUAf);
        U.correctBoundaryConditions();
    }

#   include "continuityErrs.H"
}
