    // Solve the Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff()
    );

    UEqn.relax();

    solve
    (
        UEqn
     ==
      - fvc::reconstruct(ghf*fvc::snGrad(rho)*mesh.magSf())
      - fvc::grad(p_rgh)
    );
