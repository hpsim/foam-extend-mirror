{
    kappaEff = rho*Cp*(turbulence->nu()/Pr + turbulence->nut()/Prt);

    rhoCpPhi = rho*Cp*phi;

    fvScalarMatrix TEqn
    (
        fvm::div(rhoCpPhi, T)
      - fvm::laplacian(kappaEff, T)
    );

    if (thermalModelFluid.active())
    {
        Info<< "Adding thermal source: "
            << gMax(thermalModelFluid.S()().internalField())
            << endl;

        TEqn -= thermalModelFluid.S();
    }

    TEqn.relax();
    TEqn.solve();

    rhok = 1.0 - beta*(T - TRef);
}
