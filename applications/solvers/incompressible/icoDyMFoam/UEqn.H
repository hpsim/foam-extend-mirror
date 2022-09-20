    // Convection-diffusion matrix
    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      - fvm::laplacian(nu, U)
    );

    UEqn.relax();

    solve(UEqn == -fvc::grad(p));

    rAU = 1.0/UEqn.A();
