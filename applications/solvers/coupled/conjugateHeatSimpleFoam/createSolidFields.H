    Info<< "Reading field Tsolid\n" << endl;
    volScalarField Tsolid
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            solidMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        solidMesh
    );

    Info<< "Reading solid thermal properties" << endl;
    thermalModel solidThermalModel(Tsolid);

    Info<< "Reading solid diffusivity k\n" << endl;
    volScalarField kSolid
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            solidMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        solidMesh
    );
