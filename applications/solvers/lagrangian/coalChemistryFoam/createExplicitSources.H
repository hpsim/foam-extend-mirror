    Info<< "\nConstructing explicit mass source" << endl;
    scalarTimeActivatedExplicitSourceList massSource
    (
        "mass",
        mesh,
        dimMass/dimTime/dimVolume,
        "mass"
    );

    Info<< "\nConstructing explicit momentum source" << endl;
    vectorTimeActivatedExplicitSourceList momentumSource
    (
        "momentum",
        mesh,
        dimVelocity,
        "U"
    );

    Info<< "\nConstructing explicit enthalpy source" << endl;
    scalarTimeActivatedExplicitSourceList enthalpySource
    (
        "energy",
        mesh,
        dimEnergy/dimMass,
        "hs"
    );

    Info<< "\nConstructing explicit species source" << endl;
    scalarTimeActivatedExplicitSourceList speciesSource
    (
        "species",
        mesh,
        dimless,
        composition.species()
    );
