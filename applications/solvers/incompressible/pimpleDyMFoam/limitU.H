{
    // Bound the velocity
    volScalarField magU = mag(U);

    Info<< "mag(U) = " << max(magU).value()
        << " (" << gMax(magU.internalField()) << ")";
    
    if (max(magU) > UMax)
    {
        volScalarField Ulimiter =
            pos(magU - UMax)*UMax/(magU + smallU)
          + neg(magU - UMax);

        Ulimiter.max(scalar(0));
        Ulimiter.min(scalar(1));

        U *= Ulimiter;
        U.correctBoundaryConditions();

        Info<< " limiting: " << gMax(mag(U)().internalField());
    }

    Info<< endl;
}
