tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.schemesDict().divScheme("div(phi,Yi_h)")
    )
);

{
    label inertIndex = -1;
    volScalarField Yt("Yt", 0.0*Y[0]);

    forAll(Y, i)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];
            solve
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
              ==
                coalParcels.Srho(i)
              + kappa*chemistry.RR(i)().dimensionedInternalField()
              + massSource.Su()*speciesSource.Su(i),
                mesh.solutionDict().solver("Yi")
            );

            Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
