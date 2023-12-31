{

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + fvm::div(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        DpDt
    );

    hEqn.solve();

    thermo.correct();

    // Recalculate density
    rho = thermo.rho();
    rho.correctBoundaryConditions();
}
