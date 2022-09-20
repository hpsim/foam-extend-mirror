    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff()
     ==
        rho.dimensionedInternalField()*g
      + coalParcels.SU()
      + limestoneParcels.SU()
      + massSource.Su()*momentumSource.Su()
    );

    UEqn.relax();

    pZones.addResistance(UEqn);

    solve(UEqn == -fvc::grad(p));
