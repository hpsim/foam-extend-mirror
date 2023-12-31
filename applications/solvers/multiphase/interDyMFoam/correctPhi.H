{
#   include "continuityErrs.H"

    // Moved pcorrtypes because dynamic load balancing changes the number
    // of (processor) patches.  HJ, 23/Feb/2022
    wordList pcorrTypes
    (
        pd.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    for (label i = 0; i < pd.boundaryField().size(); i++)
    {
        if (pd.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", pd.dimensions(), 0.0),
        pcorrTypes
    );

    dimensionedScalar rAUf
    (
        "(1|A(U))",
        dimTime/rho.dimensions(),
        runTime.deltaT().value()
    );

    phi = (fvc::interpolate(U) & mesh.Sf());

    adjustPhi(phi, U, pcorr);

    mesh.schemesDict().setFluxRequired(pcorr.name());

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pdRefCell, pdRefValue);
        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }

#   include "continuityErrs.H"
#   include "CourantNo.H"

    // Recalculate rhoPhi from rho
    rhoPhi = fvc::interpolate(rho)*phi;
}
