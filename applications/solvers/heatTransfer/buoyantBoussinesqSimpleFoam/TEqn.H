{
    volScalarField kappaEff
    (
        "kappaEff",
        turbulence->nu()/Pr + turbulence->nut()/Prt
    );

    fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      + fvm::SuSp(-fvc::div(phi), T)
      - fvm::laplacian(kappaEff, T)
    );

    if (thermalModelFluid.active())
    {
        Info<< "Adding thermal source: "
            << gMax(thermalModelFluid.S()().internalField())
            << endl;

        TEqn -=
            thermalModelFluid.S()/
            (thermalModelFluid.rho()*thermalModelFluid.C());
    }

    TEqn.relax();
    TEqn.solve();

    rhok = 1.0 - beta*(T - TRef);
}
