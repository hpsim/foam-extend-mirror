{
    fvScalarMatrix eEqn
    (
        fvm::ddt(rho, e)
      + fvm::div(phi, e)
      - fvm::laplacian(turbulence->alphaEff(), e)
     ==
      - p*fvc::div(phi/fvc::interpolate(rho))
        // viscous heating?
    );

    eEqn.relax
    (
        mesh.solutionDict().equationRelaxationFactor(e.select(pimple.finalIter()))
    );

    eEqn.solve();

    // Bound the energy using TMin and TMax
    {
        dimensionedScalar Tstd("Tstd", dimTemperature, specie::Tstd());

        volScalarField Cv = thermo.Cv();
        volScalarField R = thermo.Cp() - Cv;

        e = Foam::min(e, TMax*Cv + R*Tstd);
        e = Foam::max(e, TMin*Cv + R*Tstd);
        e.correctBoundaryConditions();
    }

    thermo.correct();

    // Recalculate compressibility
    psis = thermo.psi()/thermo.Cp()*thermo.Cv();

    // Recalculate density
    rho = thermo.rho();
    rho.correctBoundaryConditions();
}
