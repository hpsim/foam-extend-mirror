{
    kappaEff = turbulence->nu()/Pr + turbulence->nut()/Prt;

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
    );

    if (thermalModelFluid.active())
    {
        TEqn -=
            thermalModelFluid.S()/
            (thermalModelFluid.rho()*thermalModelFluid.C());
    }

    TEqn.solve();

    rhok = 1.0 - beta*(T - TRef);
}
