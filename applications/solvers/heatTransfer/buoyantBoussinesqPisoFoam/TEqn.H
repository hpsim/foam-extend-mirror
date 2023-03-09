{
    kappaEff = thermalModelFluid.rho()*thermalModelFluid.C()*
        (turbulence->nu()/Pr + turbulence->nut()/Prt);

    rhoCpPhi =
        fvc::interpolate(thermalModelFluid.rho()*thermalModelFluid.C())*phi;

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(rhoCpPhi, T)
      - fvm::laplacian(kappaEff, T)
    );

    if (thermalModelFluid.active())
    {
        TEqn -= thermalModelFluid.S();
    }

    TEqn.solve();

    rhok = 1.0 - beta*(T - TRef);
}
