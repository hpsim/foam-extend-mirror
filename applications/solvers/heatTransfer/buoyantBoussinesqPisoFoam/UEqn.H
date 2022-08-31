    // Solve the momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    pZones.addResistance(UEqn);

    solve
    (
        UEqn
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
