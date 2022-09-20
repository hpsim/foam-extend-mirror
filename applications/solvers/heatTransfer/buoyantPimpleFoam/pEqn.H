{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rhorUAf("(rho*(1|A(U)))", fvc::interpolate(rho*rUA));

    U = rUA*UEqn.H();
    U.boundaryField().updateCoupledPatchFields();

    surfaceScalarField phiU
    (
        "phiU",
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
        )
    );

    phi = phiU - ghf*fvc::snGrad(rho)*rhorUAf*mesh.magSf();

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::ddt(psi, p_rgh)
          + fvc::div(phi)
          - fvm::laplacian(rhorUAf, p_rgh)
        );

        p_rghEqn.solve
        (
            mesh.solutionDict().solver(p_rgh.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi += p_rghEqn.flux();
        }
    }

    U += rUA*fvc::reconstruct((phi - phiU)/rhorUAf);
    U.correctBoundaryConditions();

    // Update pressure
    p = p_rgh + rho*gh;
    p.correctBoundaryConditions();

    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

#   include "rhoEqn.H"
#   include "compressibleContinuityErrs.H"

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity

    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
        p.correctBoundaryConditions();

        rho = thermo.rho();

        p_rgh = p - rho*gh;
        p_rgh.correctBoundaryConditions();
    }
}
