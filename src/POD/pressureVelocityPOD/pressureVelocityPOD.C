/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    pressureVelocityPOD

\*---------------------------------------------------------------------------*/

#include "pressureVelocityPOD.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pressureVelocityPOD, 0);

    addToRunTimeSelectionTable
    (
        PODODE,
        pressureVelocityPOD,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pressureVelocityPOD::calcValidTimes() const
{
    if (validTimesPtr_)
    {
        FatalErrorInFunction
            << "Valid times already calculated"
            << abort(FatalError);
    }

    // Get times list
    Time& runTime = const_cast<Time&>(mesh().time());

    // Remember time index to restore it after the scan
    label origTimeIndex = runTime.timeIndex();

    instantList Times = runTime.times();

    validTimesPtr_ = new instantList(Times.size());
    instantList& vt = *validTimesPtr_;

    label nValidTimes = 0;

    forAll (Times, i)
    {
        // Check if zero field should be used
        if (Times[i].equal(0)  && !useZeroField_)
        {
            Info << "Skipping time " << Times[i] << endl;

            continue;
        }

        runTime.setTime(Times[i], i);

        IOobject UHeader
        (
            UName_,
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        IOobject pHeader
        (
            pName_,
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ
        );

        IOobject phiHeader
        (
            phiName_,
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ
        );

        if (UHeader.headerOk() && pHeader.headerOk() && phiHeader.headerOk())
        {
            // Record time as valid
            vt[nValidTimes] = Times[i];
            nValidTimes++;
        }
    }

    // Reset time index to initial state
    runTime.setTime(Times[origTimeIndex], origTimeIndex);

    // Check snapshots
    if (nValidTimes < 2)
    {
        FatalErrorInFunction
            << "Insufficient number of snapshots (validTimes): "
            << nValidTimes
            << abort(FatalError);
    }

    // Reset the list of times
    Info<< "Number of valid snapshots: " << nValidTimes << endl;
    vt.setSize(nValidTimes);
}


const Foam::instantList& Foam::pressureVelocityPOD::validTimes() const
{
    if (!validTimesPtr_)
    {
        calcValidTimes();
    }

    return *validTimesPtr_;
}


void Foam::pressureVelocityPOD::calcOrthoBase() const
{
    if
    (
        UBasePtr_
     || pBasePtr_
     || phiBasePtr_
     || reconUPtr_
     || reconPPtr_
    )
    {
        FatalErrorInFunction
            << "Orthogonal base already calculated"
            << abort(FatalError);
    }

    // Create ortho-normal base

    // Read accuracy of the velocity base
    scalar UAccuracy = readScalar(dict().lookup("UAccuracy"));

    // Read accuracy of the pressure base
    scalar pAccuracy = readScalar(dict().lookup("pAccuracy"));

    // Get times list
    Time& runTime = const_cast<Time&>(mesh().time());

    // Remember time index to restore it after the scan
    label origTimeIndex = runTime.timeIndex();

    const instantList& valTimes = validTimes();

    // Create a list of snapshots
    PtrList<volVectorField> UFields(valTimes.size());

    PtrList<volScalarField> pFields(valTimes.size());

    PtrList<surfaceScalarField> phiFields(valTimes.size());

    forAll (valTimes, timeI)
    {
        runTime.setTime(valTimes[timeI], timeI);

        Info<< "Reading snapshots from time = " << runTime.timeName() << endl;

        // Get velocity snapshot
        UFields.set
        (
            timeI,
            new volVectorField
            (
                IOobject
                (
                    UName_,
                    runTime.timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh()
            )
        );

        UFields[timeI].rename(UName_ + name(timeI));

        // Get pressure snapshot
        pFields.set
        (
            timeI,
            new volScalarField
            (
                IOobject
                (
                    pName_,
                    runTime.timeName(),
                    mesh(),
                    IOobject::MUST_READ
                ),
                this->mesh()
            )
        );

        pFields[timeI].rename(pName_ + name(timeI));

        // Get flux snapshot
        phiFields.set
        (
            timeI,
            new surfaceScalarField
            (
                IOobject
                (
                    phiName_,
                    runTime.timeName(),
                    mesh(),
                    IOobject::MUST_READ
                ),
                this->mesh()
            )
        );

        phiFields[timeI].rename(phiName_ + name(timeI));

        // Read the first time into reconU and reconP for
        // boundary conditions and general reconstruction settings
        if (!reconUPtr_)
        {
            Info<< "Reading " << "recon" << UName_ << endl;
            reconUPtr_ =
                new volVectorField
                (
                    "recon" + UName_,
                    UFields[timeI]
                );
        }

        if (!reconPPtr_)
        {
            Info<< "Reading " << "recon" << pName_ << endl;
            reconPPtr_ =
                new volScalarField
                (
                    "recon" + pName_,
                    pFields[timeI]
                );
        }
    }

    // Reset time index to initial state
    runTime.setTime(valTimes[0], 0);

    // Create velocity ortho-normal base
    UBasePtr_ = new vectorPODOrthoNormalBase(UFields, UAccuracy);

    // Create pressure ortho-normal base
    pBasePtr_ = new scalarPODOrthoNormalBase(pFields, pAccuracy);

    // Collect flux field base, consistent with velocity decomposition
    phiBasePtr_ = new PtrList<surfaceScalarField>(UBasePtr_->baseSize());
    UBasePtr_->getOrthoBase(phiFields, *phiBasePtr_);

    // Check orthogonality and magnitude of snapshots
    if (debug)
    {
        Info<< "Checking orthogonality of velocity base" << endl;
        UBasePtr_->checkBase();

        Info<< "Checking orthogonality of pressure base" << endl;
        pBasePtr_->checkBase();
    }


    // Write snapshots
    if (debug)
    {
        // Note: use UBasePtr_ and pBasePtr_ raw, as they have just been
        // calculated

        Info<< "Write reconstructed snapshots: check" << endl;

        forAll (valTimes, timeI)
        {
            runTime.setTime(valTimes[timeI], timeI);

            Info<< "Time = " << runTime.timeName() << endl;

            // Reconstruct velocity
            {
                const volVectorField& reconU = *reconUPtr_;

                volVectorField directReconU
                (
                    "directReconU",
                    reconU
                );

                directReconU =
                    dimensionedVector
                    (
                        "zero",
                        reconU.dimensions(),
                        vector::zero
                    );

                vectorField& UIn = directReconU.internalField();

                for (label uI = 0; uI < UBasePtr_->baseSize(); uI++)
                {
                    // Add to velocity
                    UIn +=
                        UBasePtr_->interpolationCoeffs()[timeI][uI]*
                        UBasePtr_->orthoField(uI);
                }

                // Internal field is set.  Correct boundary conditions
                directReconU.correctBoundaryConditions();
                directReconU.write();
            }

            // Reconstruct pressure
            {
                const volScalarField& reconP = *reconPPtr_;

                volScalarField directReconP
                (
                    "directReconP",
                    reconP
                );

                directReconP =
                    dimensionedScalar("zero", reconP.dimensions(), scalar(0));

                scalarField& pIn = directReconP.internalField();

                for (label pI = 0; pI < pBasePtr_->baseSize(); pI++)
                {
                    // Add to pressure
                    pIn +=
                        pBasePtr_->interpolationCoeffs()[timeI][pI]*
                        pBasePtr_->orthoField(pI);
                }

                directReconP.correctBoundaryConditions();
                directReconP.write();
            }
        }
    }

    // Reset time index to initial state
    runTime.setTime(runTime.times()[origTimeIndex], origTimeIndex);
}


void Foam::pressureVelocityPOD::calcDerivativeCoeffs() const
{
    if
    (
        convectionDerivativePtr_
     || diffusionDerivativePtr_
     || pressureGradDerivativePtr_
     || pressureLaplaceDerivativePtr_
     || pressureSourceDerivativePtr_
     || lagrangeDerPtr_
     || lagrangeSrcPtr_
    )
    {
        FatalErrorInFunction
            << "Derivative matrices already calculated"
            << abort(FatalError);
    }

    // Calculate coefficients for differential equation

    // Get velocity base
    const vectorPODOrthoNormalBase& Ub = UBase();

    // Get pressure base
    const scalarPODOrthoNormalBase& pb = pBase();

    // Flux fields (using velocity base)
    const PtrList<surfaceScalarField>& phib = phiBase();

    // Create derivative matrices

    // Convection derivative
    convectionDerivativePtr_ = new PtrList<scalarSquareMatrix>(Ub.baseSize());

    PtrList<scalarSquareMatrix>& convectionDerivative =
        *convectionDerivativePtr_;

    forAll (convectionDerivative, i)
    {
        convectionDerivative.set
        (
            i,
            new scalarSquareMatrix(Ub.baseSize(), scalar(0))
        );
    }

    // Diffusion derivative
    diffusionDerivativePtr_ =
        new scalarSquareMatrix(Ub.baseSize(), scalar(0));

    scalarSquareMatrix& diffusionDerivative = *diffusionDerivativePtr_;

    // Pressure gradient derivative
    pressureGradDerivativePtr_ =
        new scalarRectangularMatrix(Ub.baseSize(), pb.baseSize(), scalar(0));

    scalarSquareMatrix& pressureGradDerivative = *diffusionDerivativePtr_;

    // Pressure laplacian derivative
    pressureLaplaceDerivativePtr_ =
        new scalarSquareMatrix(pb.baseSize(), scalar(0));

    scalarSquareMatrix& pressureLaplaceDerivative =
        *pressureLaplaceDerivativePtr_;

    // Pressure source derivative
    pressureSourceDerivativePtr_ =
        new scalarRectangularMatrix(pb.baseSize(), Ub.baseSize(), scalar(0));

    scalarRectangularMatrix& pressureSourceDerivative =
        *pressureSourceDerivativePtr_;


    // Lagrange multiplier derivative
    lagrangeDerPtr_ = new scalarField(Ub.baseSize(), scalar(0));
    scalarField& lagrangeDer = *lagrangeDerPtr_;

    // Lagrange multiplier source
    lagrangeSrcPtr_ = new scalarField(Ub.baseSize(), scalar(0));
    scalarField& lagrangeSrc = *lagrangeSrcPtr_;

    // Get solution field for boundary conditions.  Force raw access without
    // checking
    const volVectorField& U = *reconUPtr_;

    // i = equation number (normal basis in dot-product)
    // j = summation over all bases
    // k = second summation for convection term

    // Read laplacian factor 1/deltaT
    dimensionedScalar deltaT
    (
        "deltaT",
        dimTime,
        readScalar(dict().lookup("laplaceFactor"))
    );

    // Calculate derivative by moving equation terms to rhs

    // Calculate volume for scaling - THIS IS WRONG!!! HJ, HERE!!!
    const scalar V = 0.2*gSum(mesh().V().field());
    // const scalar V = 1;

    // Derivatives assembly loop: Momentum equation
    for (label uI = 0; uI < Ub.baseSize(); uI++)
    {
        const volVectorField& snapUI = Ub.orthoField(uI);

        for (label uJ = 0; uJ < Ub.baseSize(); uJ++)
        {
            const volVectorField& snapUJ = Ub.orthoField(uJ);

            // Convection derivative
            for (label uK = 0; uK < Ub.baseSize(); uK++)
            {
                const surfaceScalarField& snapPhiK = phib[uK];

                convectionDerivative[uI][uJ][uK] =
                   -POD::projection
                    (
                        fvc::div
                        (
                            snapPhiK, snapUJ,
                            "div(" + phiName_ + "," + UName_ + ")"
                        ),
                        snapUI
                    );
            }

            // Diffusion derivative
            diffusionDerivative[uI][uJ] =
                POD::projection
                (
                    fvc::laplacian
                    (
                        nu_, snapUJ,
                        "laplacian(" + nu_.name() + "," + UName_ + ")"
                    ),
                    snapUI
                );

            // Lagrange multiplier is calculated on boundaries where
            // reconU fixes value
            forAll (U.boundaryField(), patchI)
            {
                if (U.boundaryField()[patchI].fixesValue())
                {
                    // Note pre-multiplication by beta and signs
                    // of derivative and source.  Su-Sp treatment
                    lagrangeDer[uI] += -beta_*
                        POD::projection
                        (
                            snapUJ.boundaryField()[patchI],
                            snapUI.boundaryField()[patchI]
                        );

                    lagrangeSrc[uI] += beta_*
                        POD::projection
                        (
                            U.boundaryField()[patchI],
                            snapUI.boundaryField()[patchI]
                        );
                }
            }
        }

        // Pressure gradient derivative
        for (label pJ = 0; pJ < pb.baseSize(); pJ++)
        {
            const volScalarField& snapPJ = pb.orthoField(pJ);

            pressureGradDerivative[uI][pJ] =
                - POD::projection(fvc::grad(snapPJ, "grad(p)"), snapUI)*V;
        }
    }

    // Derivatives assembly loop: Pressure equation
    for (label pI = 0; pI < pb.baseSize(); pI++)
    {
        const volScalarField& snapPI = pb.orthoField(pI);

        // Pressure laplacian
        for (label pJ = 0; pJ < pb.baseSize(); pJ++)
        {
            const volScalarField& snapPJ = pb.orthoField(pJ);

            pressureLaplaceDerivative[pI][pJ] =
                POD::projection
                (
                    fvc::laplacian
                    (
                        snapPJ,
                        "laplacian(" + pName_ + ")"
                    ),
                    snapPI
                );
        }

        volVectorField gradSnapPI = fvc::grad(snapPI);

        // Divergence of velocity in the pressure equation
        for (label uJ = 0; uJ < Ub.baseSize(); uJ++)
        {
            const volVectorField& snapUJ = Ub.orthoField(uJ);

            // Note: pressure laplacian appears on both sides
            // of the equation and all terms apart from gradSnapPI?
            pressureSourceDerivative[pI][uJ] =
                POD::projection
                (
                    fvc::div(snapUJ + deltaT*gradSnapPI, "div(U)"),
                    snapPI
                );
        }
    }

    // Renormalize pressure laplacian
    // Note: ODE equation format is d alpha / dt = dydx;
    for (label pI = 0; pI < pb.baseSize(); pI++)
    {
        const scalar diagCoeff = pressureLaplaceDerivative[pI][pI];

        // Pressure laplacian
        for (label pJ = 0; pJ < pb.baseSize(); pJ++)
        {
            pressureLaplaceDerivative[pI][pJ] /= diagCoeff;
        }

        // Divergence of velocity in the pressure equation
        for (label uJ = 0; uJ < Ub.baseSize(); uJ++)
        {
            pressureSourceDerivative[pI][uJ] /= diagCoeff;
        }

        // Remove zero coefficient
        pressureLaplaceDerivative[pI][pI] = 0;
    }

    // Divide pressure laplace coefficient by diagonal and remove zero


    // Info<< "convectionDerivative: " << convectionDerivative << nl
    //     << "diffusionDerivative: " << diffusionDerivative << nl
    //     << "lagrangeDer: " << lagrangeDer << nl
    //     << "lagrangeSrc: " << lagrangeSrc << endl;
}


void Foam::pressureVelocityPOD::updateFields() const
{
    if (fieldUpdateTimeIndex_ < mesh().time().timeIndex())
    {
        // Field update required.  Record update time index
        fieldUpdateTimeIndex_ = mesh().time().timeIndex();

        if (!reconUPtr_ || !reconPPtr_)
        {
            FatalErrorInFunction
                << "Reconstructed fields not allocated"
                << abort(FatalError);
        }

        Info<< "coeffs_: " << coeffs_ << endl;

        // Reconstruct velocity
        {
            volVectorField& reconU = *reconUPtr_;

            // Reset entire field, including boundary conditions
            reconU ==
                dimensionedVector("zero", reconU.dimensions(), vector::zero);

            const vectorPODOrthoNormalBase& Ub = UBase();

            for (label uI = 0; uI < Ub.baseSize(); uI++)
            {
                reconU == reconU + coeffs_[uI]*Ub.orthoField(uI);
            }

            // Internal field is set.  Correct boundary conditions
            reconU.correctBoundaryConditions();

            // Hack test, HJ HERE!!!
            const scalar fluxZero =
                gSum(reconU.boundaryField()[0] & mesh().Sf().boundaryField()[0]);

            Info<< "Testing flux in: "
                << fluxZero << " velocity "
                << fluxZero/gSum(mesh().magSf().boundaryField()[0])
                << endl;
        }

        // Reconstruct pressure
        {
            volScalarField& reconP = *reconPPtr_;

            // Reset entire field, including boundary conditions
            reconP = dimensionedScalar("zero", reconP.dimensions(), scalar(0));

            const scalarPODOrthoNormalBase& pb = pBase();

            const label pOffset = UBase().baseSize();

            for (label pI = 0; pI < pb.baseSize(); pI++)
            {
                // Update pressure
                reconP == reconP + coeffs_[pI + pOffset]*pb.orthoField(pI);
            }

            // Internal field is set.  Correct boundary conditions
            reconP.correctBoundaryConditions();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pressureVelocityPOD::pressureVelocityPOD
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PODODE(mesh, dict),
    UName_(dict.lookup("U")),
    pName_(dict.lookup("p")),
    phiName_(dict.lookup("phi")),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            this->mesh().time().constant(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nu_(transportProperties_.lookup("nu")),
    beta_(readScalar(dict.lookup("beta"))),
    useZeroField_(dict.lookup("useZeroField")),
    coeffs_(),
    validTimesPtr_(nullptr),
    convectionDerivativePtr_(nullptr),
    diffusionDerivativePtr_(nullptr),
    pressureGradDerivativePtr_(nullptr),
    pressureLaplaceDerivativePtr_(nullptr),
    pressureSourceDerivativePtr_(nullptr),
    lagrangeDerPtr_(nullptr),
    lagrangeSrcPtr_(nullptr),
    UBasePtr_(nullptr),
    pBasePtr_(nullptr),
    phiBasePtr_(nullptr),
    reconUPtr_(nullptr),
    reconPPtr_(nullptr),
    fieldUpdateTimeIndex_(-1)
{
    // Check beta
    if (beta_ < 0)
    {
        FatalErrorInFunction
            << "Negative beta: " << beta_
            << abort(FatalError);
    }

    // Grab coefficients from the first snapshot of the ortho-normal base
    coeffs_.setSize(UBase().baseSize() + pBase().baseSize());

    const scalarRectangularMatrix& UBaseCoeffs =
        UBase().interpolationCoeffs();

    // Insert velocity coefficients first
    for (label uI = 0; uI < UBase().baseSize(); uI++)
    {
        coeffs_[uI] = UBaseCoeffs[0][uI];
    }

    // Insert pressure coefficients with offset
    const scalarRectangularMatrix& pBaseCoeffs =
        pBase().interpolationCoeffs();

    const label pOffset = UBase().baseSize();

    for (label pI = 0; pI < pBase().baseSize(); pI++)
    {
        coeffs_[pI + pOffset] = pBaseCoeffs[0][pI];
    }

    Info<< "Zero coeffs: " << coeffs_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pressureVelocityPOD::~pressureVelocityPOD()
{
    deleteDemandDrivenData(validTimesPtr_);
    deleteDemandDrivenData(convectionDerivativePtr_);
    deleteDemandDrivenData(diffusionDerivativePtr_);
    deleteDemandDrivenData(pressureGradDerivativePtr_);
    deleteDemandDrivenData(pressureLaplaceDerivativePtr_);
    deleteDemandDrivenData(pressureSourceDerivativePtr_);
    deleteDemandDrivenData(lagrangeDerPtr_);
    deleteDemandDrivenData(lagrangeSrcPtr_);

    deleteDemandDrivenData(UBasePtr_);
    deleteDemandDrivenData(pBasePtr_);
    deleteDemandDrivenData(phiBasePtr_);

    deleteDemandDrivenData(reconUPtr_);
    deleteDemandDrivenData(reconPPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::pressureVelocityPOD::nEqns() const
{
    return coeffs().size();
}


Foam::scalarField& Foam::pressureVelocityPOD::coeffs()
{
    return coeffs_;
}


const Foam::scalarField& Foam::pressureVelocityPOD::coeffs() const
{
    return coeffs_;
}


void Foam::pressureVelocityPOD::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    if
    (
        !convectionDerivativePtr_
     || !diffusionDerivativePtr_
     || !pressureGradDerivativePtr_
     || !pressureLaplaceDerivativePtr_
     || !pressureSourceDerivativePtr_
     || !lagrangeDerPtr_
     || !lagrangeSrcPtr_
    )
    {
        calcDerivativeCoeffs();
    }

    // Notes:
    // x = time
    // y = current values of coefficients

    const label pOffset = UBase().baseSize();

    // Clear derivatives matrix
    dydx = 0;

    // Insert momentum equation
    {
        // Convection derivative
        const PtrList<scalarSquareMatrix>& convectionDerivative =
            *convectionDerivativePtr_;

        // Diffusion derivative
        const scalarSquareMatrix& diffusionDerivative =
            *diffusionDerivativePtr_;

        // Pressure gradient derivative
        const scalarRectangularMatrix& pressureGradDerivative =
            *pressureGradDerivativePtr_;

        // Lagrange multiplier derivative, boundary conditions
        const scalarField& lagrangeDer = *lagrangeDerPtr_;

        // Lagrange multiplier source, boundary conditions
        const scalarField& lagrangeSrc = *lagrangeSrcPtr_;

        for (label uI = 0; uI < UBase().baseSize(); uI++)
        {
            dydx[uI] = 0;
            dydx[uI] = lagrangeSrc[uI];

            for (label uJ = 0; uJ < UBase().baseSize(); uJ++)
            {
                // Convection
                for (label uK = 0; uK < UBase().baseSize(); uK++)
                {
                    dydx[uI] += convectionDerivative[uI][uJ][uK]*y[uJ]*y[uK];
                }

                // Diffusion
                dydx[uI] +=
                    diffusionDerivative[uI][uJ]*y[uJ]
                  + lagrangeDer[uI]*y[uJ];
            }

            // Pressure gradient
            // for (label pJ = 0; pJ < pBase().baseSize(); pJ++)
            // {
            //     dydx[uI] += pressureGradDerivative[uI][pJ]*y[pJ + pOffset];
            // }
        }
    }

    // Insert pressure equation.  Note the manipulation because of the
    // lack of ddt term
    {
        // Pressure laplacian derivative
        const scalarSquareMatrix& pressureLaplaceDerivative =
            *pressureLaplaceDerivativePtr_;

        // Pressure source derivative
        const scalarRectangularMatrix& pressureSourceDerivative =
            *pressureSourceDerivativePtr_;

        for (label pI = 0; pI < pBase().baseSize(); pI++)
        {
            dydx[pI + pOffset] = 0;

            // Pressure laplacian
            for (label pJ = 0; pJ < pBase().baseSize(); pJ++)
            {
                dydx[pI + pOffset] +=
                    pressureLaplaceDerivative[pI][pJ]*y[pJ + pOffset];
            }

            // Pressure source
            for (label uJ = 0; uJ < UBase().baseSize(); uJ++)
            {
                dydx[pI + pOffset] +=
                    pressureSourceDerivative[pI][uJ]*y[uJ];
            }
        }
    }
}


void Foam::pressureVelocityPOD::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    if
    (
        !convectionDerivativePtr_
     || !diffusionDerivativePtr_
     || !pressureGradDerivativePtr_
     || !pressureLaplaceDerivativePtr_
     || !pressureSourceDerivativePtr_
     || !lagrangeDerPtr_
     || !lagrangeSrcPtr_
    )
    {
        calcDerivativeCoeffs();
    }

    // Must calculate derivatives
    derivatives(x, y, dfdx);

    // Clear jacobian matrix
    forAll (y, i)
    {
        forAll (y, j)
        {
            dfdy[i][j] = 0;
        }
    }

    // Notes:
    // x = time
    // y = current values of coefficients

    const label pOffset = UBase().baseSize();

    // Insert momentum equation
    {
        // Convection derivative
        const PtrList<scalarSquareMatrix>& convectionDerivative =
            *convectionDerivativePtr_;

        // Diffusion derivative
        const scalarSquareMatrix& diffusionDerivative =
            *diffusionDerivativePtr_;

        // // Pressure gradient derivative
        // const scalarRectangularMatrix& pressureGradDerivative =
        //     *pressureGradDerivativePtr_;

        // Lagrange multiplier derivative, boundary conditions
        const scalarField& lagrangeDer = *lagrangeDerPtr_;

        for (label uI = 0; uI < UBase().baseSize(); uI++)
        {
            for (label uJ = 0; uJ < UBase().baseSize(); uJ++)
            {
                // Convection
                for (label uK = 0; uK < UBase().baseSize(); uK++)
                {
                    dfdy[uI][uJ] += convectionDerivative[uI][uJ][uK]*y[uK];
                }

                // Diffusion
                dfdy[uI][uJ] +=
                    diffusionDerivative[uI][uJ]
                  + lagrangeDer[uI];
            }

            // // Pressure gradient
            // for (label pJ = 0; pJ < pBase().baseSize(); pJ++)
            // {
            //     dfdy[uI][pJ + pOffset] = pressureGradDerivative[uI][pJ];
            // }
        }
    }
    //HJ: INCOMPLETE!!!
    // Insert pressure equation.  Note the manipulation because of the
    // lack of ddt term
    // {
    //     // Pressure laplacian derivative
    //     const scalarSquareMatrix& pressureLaplaceDerivative =
    //         *pressureLaplaceDerivativePtr_;

    //     // Pressure source derivative
    //     const scalarRectangularMatrix& pressureSourceDerivative =
    //         *pressureSourceDerivativePtr_;

    //     // Pressure laplacian
    //     for (label pI = 0; pI < pBase().baseSize(); pI++)
    //     {
    //         for (label pJ = 0; pJ < pBase().baseSize(); pJ++)
    //         {
    //             dfdy[pI + pOffset][pJ + pOffset] =
    //                 pressureLaplaceDerivative[pI][pJ];
    //         }
    //     }

    //     // Pressure source
    // }
}


void Foam::pressureVelocityPOD::update(const scalar delta)
{
    // Warning: lots of cost here: updating fields for top-level
    // function objects.  Can be removed from production code
    // HJ, 5/Aug/2020
    updateFields();
}


const Foam::vectorPODOrthoNormalBase& Foam::pressureVelocityPOD::UBase() const
{
    if (!UBasePtr_)
    {
        calcOrthoBase();
    }

    return *UBasePtr_;
}


const Foam::scalarPODOrthoNormalBase& Foam::pressureVelocityPOD::pBase() const
{
    if (!pBasePtr_)
    {
        calcOrthoBase();
    }

    return *pBasePtr_;
}


const Foam::PtrList<Foam::surfaceScalarField>&
Foam::pressureVelocityPOD::phiBase() const
{
    if (!phiBasePtr_)
    {
        calcOrthoBase();
    }

    return *phiBasePtr_;
}


const Foam::volVectorField& Foam::pressureVelocityPOD::reconU() const
{
    updateFields();

    return *reconUPtr_;
}


const Foam::volScalarField& Foam::pressureVelocityPOD::reconP() const
{
    updateFields();

    return *reconPPtr_;
}


void Foam::pressureVelocityPOD::writeSnapshots() const
{
    const vectorPODOrthoNormalBase& Ub = UBase();
    const scalarPODOrthoNormalBase& pb = pBase();

    Info<< "Writing POD base for Time = " << mesh().time().timeName() << endl;

    for (label uI = 0; uI < Ub.baseSize(); uI++)
    {
        // Write snapshots
        volVectorField snapU
        (
            UName_ + "POD" + name(uI),
            Ub.orthoField(uI)
        );
        snapU.write();
    }

    for (label pI = 0; pI < pb.baseSize(); pI++)
    {
        volScalarField snapP
        (
            pName_ + "POD" + name(pI),
            pb.orthoField(pI)
        );
        snapP.write();
    }
}


void Foam::pressureVelocityPOD::write() const
{
    // Recalculate field and force a write
    updateFields();
    reconU().write();
    reconP().write();
}


// ************************************************************************* //
