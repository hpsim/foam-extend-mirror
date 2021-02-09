{
    // Solve the rothalpy equation

    // Create relative velocity
    URel == U;
    mrfZones.relativeVelocity(URel);

    // Bound the relative velocity to preserve the i to h conversion bound
    // HJ, 22/Mar/2017
    volScalarField magURel = mag(URel);

    if (max(magURel) > UMax)
    {
        volScalarField URelLimiter = pos(magURel - UMax)*UMax/(magURel + smallU)
            + neg(magURel - UMax);
        URelLimiter.max(scalar(0));
        URelLimiter.min(scalar(1));

        URel *= URelLimiter;
        URel.correctBoundaryConditions();
    }

    // Create rotational velocity (= omega x r)
    URot = U - URel;

    T.storePrevIter();

    i == h + 0.5*magSqr(U) - mag(UTheta)*mag(URot);
    i.correctBoundaryConditions();

    fvScalarMatrix iEqn
    (
        fvm::div(phi, i)
      + fvm::SuSp(-fvc::div(phi), i)
      - fvm::laplacian(turbulence->alphaEff(), i)
     ==
        // Viscous heating: note sign (devRhoReff has a minus in it)
      - (turbulence->devRhoReff() && fvc::grad(URel))
    );

    iEqn.relax();
    iEqn.solve();

    // Calculate enthalpy out of rothalpy
    h == i - 0.5*magSqr(U) + mag(UTheta)*mag(URot);
    h.correctBoundaryConditions();

    thermo.correct();
    psis = thermo.psi()/thermo.Cp()*thermo.Cv();
}
