{
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
    Urot == U - URel;

    fvScalarMatrix iEqn
    (
        fvm::div(phi, i)
      + fvm::SuSp(-fvc::div(phi), i)
      - fvm::laplacian(turbulence->alphaEff(), i)
     ==
      // Viscous heating: note sign (devRhoReff has a minus in it)
       - (turbulence->devRhoReff() && fvc::grad(U))
    );

    iEqn.relax();

    iEqn.solve();

    // From rothalpy, calculate enthalpy after solution of rothalpy equation
    h = i + 0.5*(magSqr(Urot) - magSqr(URel));
    h.correctBoundaryConditions();

    // Update thermo for new h
    thermo.correct();
    psis = thermo.psi()/thermo.Cp()*thermo.Cv();
}
