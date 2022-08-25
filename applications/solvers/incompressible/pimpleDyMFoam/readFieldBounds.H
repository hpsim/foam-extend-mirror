    // Read field bounds
    dictionary fieldBounds = mesh.solutionDict().subDict("fieldBounds");

    // Pressure bounds
    dimensionedScalar pMin("pMin", dimPressure, -GREAT);
    dimensionedScalar pMax("pMax", dimPressure, GREAT);

    fieldBounds.lookup("p") >> pMin.value() >> pMax.value();

    // Velocity bound
    dimensionedScalar UMax("UMax", dimVelocity, GREAT);

    fieldBounds.lookup(U.name()) >> UMax.value();
    dimensionedScalar smallU("smallU", dimVelocity, 1e-10);
