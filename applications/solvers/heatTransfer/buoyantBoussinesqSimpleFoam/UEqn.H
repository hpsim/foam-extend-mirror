    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
      + fvm::SuSp(-fvc::div(phi), U)
    );

    pZones.addResistance(UEqn());

    UEqn().relax();

    solve
    (
        UEqn()
      ==
        fvc::reconstruct
        (
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        )
    );
