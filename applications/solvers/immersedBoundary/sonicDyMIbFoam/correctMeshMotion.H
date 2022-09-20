{
    // Motion inside PIMPLE loop: simple IB update
    bool meshChanged = mesh.update();

#   include "updateIbPatchFields.H"
#   include "updateIbMasks.H"
#   include "volContinuity.H"

    // Mesh motion update
    if (meshChanged)
    {
        T.max(TMin);
        p.max(pMin);
        e == max(e, thermo.Cv()*TMin);

        thermo.correct();
        rho = thermo.rho();

        if (correctPhi)
        {
// #           include "correctPhi.H"
        }
    }

    if (checkMeshCourantNo)
    {
#       include "immersedBoundaryCompressibleCourantNo.H"
    }

    // Correct velocity boundary conditions after mesh motion, with relative
    // fluxes
    U.correctBoundaryConditions();
}
