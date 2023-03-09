    // Solve the momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff()
    );

    UEqn.relax();

    solve(UEqn == -fvc::grad(p));
