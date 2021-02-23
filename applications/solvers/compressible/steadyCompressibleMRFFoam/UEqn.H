    // Added due to error during simulation start?
    U.storePrevIter();
    URel.storePrevIter();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff()
    );

    // MRF: add Coriolis force
    mrfZones.addCoriolis(rho, UEqn);

    UEqn.relax();

    solve(UEqn == -fvc::grad(p));
