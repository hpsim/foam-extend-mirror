{
    // Calculate U from convection-diffusion matrix
    U = rAU*UEqn.H();

    phi = fvc::interpolate(U) & mesh.Sf();

    adjustPhi(phi, U, p);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p)
         ==
            fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve
        (
            mesh.solutionDict().solver(p.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

#   include "movingMeshContinuityErrs.H"

    U -= rAU*fvc::grad(p);
    U.correctBoundaryConditions();
}
