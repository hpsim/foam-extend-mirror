    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading  thermal properties" << endl;
    thermalModel solidThermo(T);

    Info<< "Reading  diffusivity k\n" << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    k.correctBoundaryConditions();

    Info<< "Constructing rhoCp\n" << endl;
    volScalarField rhoCp
    (
        IOobject
        (
            "rhoCpsolid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidThermo.rho()*solidThermo.C()
    );
    rhoCp.oldTime();
    rhoCp.correctBoundaryConditions();
