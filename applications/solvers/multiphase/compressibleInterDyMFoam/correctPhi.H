{
    if (mesh.changing())
    {
        forAll(U.boundaryField(), patchi)
        {
            if (U.boundaryField()[patchi].fixesValue())
            {
                U.boundaryField()[patchi].initEvaluate();
            }
        }

        forAll(U.boundaryField(), patchi)
        {
            if (U.boundaryField()[patchi].fixesValue())
            {
                U.boundaryField()[patchi].evaluate();

                phi.boundaryField()[patchi] =
                U.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
            }
        }
    }

    #include "continuityErrs.H"

    // Moved pcorrtypes because dynamic load balancing changes the number
    // of (processor) patches.  HJ, 23/Feb/2022
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    for (label i=0; i<p.boundaryField().size(); i++)
    {
        if (p.boundaryField()[i].fixesValue())
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
        dimensionedScalar("pcorr", p.dimensions(), 0.0),
        pcorrTypes
    );

    dimensionedScalar rAUf("(1|A(U))", dimTime/rho.dimensions(), 1.0);

    adjustPhi(phi, U, pcorr);

    mesh.schemesDict().setFluxRequired(pcorr.name());

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi) - divU
        );

        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }

#   include "continuityErrs.H"
}
