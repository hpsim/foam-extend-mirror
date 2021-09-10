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
    incompressibleFlowPOD

\*---------------------------------------------------------------------------*/

#include "incompressibleFlowPOD.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleFlowPOD, 0);

    addToRunTimeSelectionTable
    (
        PODODE,
        incompressibleFlowPOD,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::incompressibleFlowPOD::calcValidTimes() const
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


const Foam::instantList& Foam::incompressibleFlowPOD::validTimes() const
{
    if (!validTimesPtr_)
    {
        calcValidTimes();
    }

    return *validTimesPtr_;
}


void Foam::incompressibleFlowPOD::calcOrthoBase() const
{
    if
    (
        orthoBasePtr_
     || UBasePtr_
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
    scalar accuracy = readScalar(dict().lookup("accuracy"));

    // Get times list
    Time& runTime = const_cast<Time&>(mesh().time());

    // Remember time index to restore it after the scan
    label origTimeIndex = runTime.timeIndex();

    const instantList& valTimes = validTimes();

    // Create a list of snapshots
    PtrList<volVector4Field> UpFields(valTimes.size());

    PtrList<volVectorField> UFields(valTimes.size());

    PtrList<volScalarField> pFields(valTimes.size());

    PtrList<surfaceScalarField> phiFields(valTimes.size());

    forAll (valTimes, i)
    {
        runTime.setTime(valTimes[i], i);

        Info<< "Reading snapshots from time = " << runTime.timeName() << endl;

        // Get velocity snapshot
        UFields.set
        (
            i,
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

        UFields[i].rename(UName_ + name(i));

        // Get pressure snapshot
        pFields.set
        (
            i,
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

        pFields[i].rename(pName_ + name(i));

        // Get flux snapshot
        phiFields.set
        (
            i,
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

        phiFields[i].rename(phiName_ + name(i));

        UpFields.set
        (
            i,
            new volVector4Field
            (
                IOobject
                (
                    "Up" + name(i),
                    runTime.timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedVector4("zero", dimless, vector4::zero)
            )
        );

        // Insert U and p into shapshot
        UpFields[i].replace(0, UFields[i].component(vector::X));
        UpFields[i].replace(1, UFields[i].component(vector::Y));
        UpFields[i].replace(2, UFields[i].component(vector::Z));

        UpFields[i].replace(3, pFields[i]);

        // Read the first time into reconU and reconP for
        // boundary conditions and general reconstruction settings
        if (!reconUPtr_)
        {
            Info<< "Reading " << "recon" << UName_ << endl;
            reconUPtr_ =
                new volVectorField
                (
                    "recon" + UName_,
                    UFields[i]
                );
        }

        if (!reconPPtr_)
        {
            Info<< "Reading " << "recon" << pName_ << endl;
            reconPPtr_ =
                new volScalarField
                (
                    "recon" + pName_,
                    pFields[i]
                );
        }
    }

    // Reset time index to initial state
    runTime.setTime(valTimes[0], 0);

    // Create ortho-normal base for Up
    orthoBasePtr_ = new PODOrthoNormalBase<vector4>(UpFields, accuracy);

    // Check orthogonality and magnitude of snapshots
    if (debug)
    {
        orthoBasePtr_->checkBase();
    }

    // Collect velocity field base
    UBasePtr_ = new PtrList<volVectorField>(orthoBasePtr_->baseSize());
    PtrList<volVectorField>& UBase = *UBasePtr_;

    // Collect pressure field base
    pBasePtr_ = new PtrList<volScalarField>(orthoBasePtr_->baseSize());
    PtrList<volScalarField>& pBase = *pBasePtr_;

    for (label obpI = 0; obpI < orthoBasePtr_->baseSize(); obpI++)
    {
        // Set U and p pointers
        // Reconsider what to do with fluxes...
        UBase.set
        (
            obpI,
            new volVectorField
            (
                IOobject
                (
                    UFields[0].name() + "POD" + name(obpI),
                    UFields[0].time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedVector("zero", dimless, vector::zero)
            )
        );

        UBase[obpI].replace
        (
            vector::X,
            orthoBasePtr_->orthoField(obpI).component(0)
        );

        UBase[obpI].replace
        (
            vector::Y,
            orthoBasePtr_->orthoField(obpI).component(1)
        );

        UBase[obpI].replace
        (
            vector::Z,
            orthoBasePtr_->orthoField(obpI).component(2)
        );

        pBase.set
        (
            obpI,
            new volScalarField
            (
                IOobject
                (
                    pFields[0].name() + "POD" + name(obpI),
                    pFields[0].time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("zero", dimless, scalar(0))
            )
        );

        pBase[obpI] == orthoBasePtr_->orthoField(obpI).component(3);
    }

    // Collect flux field base
    phiBasePtr_ = new PtrList<surfaceScalarField>(orthoBasePtr_->baseSize());
    orthoBasePtr_->calcOrthoBase(phiFields, *phiBasePtr_);

    // Write snapshots
    if (debug)
    {
        Info<< "Write reconstructed snapshots: check" << endl;

        forAll (valTimes, timeI)
        {
            runTime.setTime(valTimes[timeI], timeI);

            Info<< "Time = " << runTime.timeName() << endl;

            const volVectorField& reconU = *reconUPtr_;
            const volScalarField& reconP = *reconPPtr_;

            volVectorField directReconU
            (
                "directReconU",
                reconU
            );

            volScalarField directReconP
            (
                "directReconP",
                reconP
            );

            // Reset entire field, including boundary conditions
            directReconU ==
                dimensionedVector
                (
                    "zero",
                    reconU.dimensions(),
                    vector::zero
                );

            // Reset entire field, including boundary conditions
            directReconP ==
                dimensionedScalar("zero", reconP.dimensions(), scalar(0));

            for (label obpI = 0; obpI < orthoBasePtr_->baseSize(); obpI++)
            {
                // Update velocity
                directReconU ==
                    directReconU
                  + orthoBasePtr_->interpolationCoeffs()[timeI][obpI]*
                    UBasePtr_->operator[](obpI);

                // Update pressure
                directReconP ==
                    directReconP
                  + orthoBasePtr_->interpolationCoeffs()[timeI][obpI]*
                    pBasePtr_->operator[](obpI);
            }

            // Internal field is set.  Correct boundary conditions
            directReconU.correctBoundaryConditions();
            directReconU.write();

            directReconP.correctBoundaryConditions();
            directReconP.write();
        }
    }

    // Reset time index to initial state
    runTime.setTime(runTime.times()[origTimeIndex], origTimeIndex);
}


void Foam::incompressibleFlowPOD::calcDerivativeCoeffs() const
{
    if
    (
        convectionDerivativePtr_
     || derivativePtr_
     || lagrangeDerPtr_
     || lagrangeSrcPtr_
    )
    {
        FatalErrorInFunction
            << "Derivative matrix already calculated"
            << abort(FatalError);
    }

    // Calculate coefficients for differential equation

    // Get bases
    const PODOrthoNormalBase<vector4>& b = orthoBase();

    // Velocity field base
    const PtrList<volVectorField>& Ub = UBase();

    // Pressure field base
    const PtrList<volScalarField>& pb = pBase();

    // Flux field base
    const PtrList<surfaceScalarField>& phib = phiBase();

    // Create derivative matrices

    // Convection derivative
    convectionDerivativePtr_ = new PtrList<scalarSquareMatrix>(b.baseSize());
    PtrList<scalarSquareMatrix>& convectionDerivative =
        *convectionDerivativePtr_;

    forAll (convectionDerivative, i)
    {
        convectionDerivative.set
        (
            i,
            new scalarSquareMatrix(b.baseSize(), scalar(0))
        );
    }

    // Diffusion and pressure derivative
    derivativePtr_ = new scalarSquareMatrix(b.baseSize(), scalar(0));
    scalarSquareMatrix& derivative = *derivativePtr_;

    // Lagrange multiplier derivative
    lagrangeDerPtr_ = new scalarField(b.baseSize(), scalar(0));
    scalarField& lagrangeDer = *lagrangeDerPtr_;

    // Lagrange multiplier source
    lagrangeSrcPtr_ = new scalarField(b.baseSize(), scalar(0));
    // scalarField& lagrangeSrc = *lagrangeSrcPtr_;

    // Get solution field for boundary conditions.  Force raw access without
    // checking
    const volVectorField& U = *reconUPtr_;

    // i = equation number (normal basis in dot-product)
    // j = summation over all bases
    // k = second summation for convection term

    // Calculate derivative by moving equation terms to rhs

    // Calculate volume for scaling - THIS IS WRONG!!! HJ, HERE!!!
    const scalar V = 0.5*gSum(mesh().V().field());

    // Derivatives assembly loop
    for (label i = 0; i < b.baseSize(); i++)
    {
        const volVectorField& snapUI = Ub[i];

        for (label j = 0; j < b.baseSize(); j++)
        {
            const volVectorField& snapUJ = Ub[j];

            const volScalarField& snapPJ = pb[j];

            // Note: volume scaling on all terms!

            // Convection derivative
            for (label k = 0; k < b.baseSize(); k++)
            {
                // const surfaceScalarField& snapPhiK = phib[k];

                const volVectorField& snapUK = Ub[k];

                // convectionDerivative[i][j][k] =
                //    -POD::projection
                //     (
                //         fvc::div
                //         (
                //             snapPhiK, snapUJ,
                //             "div(" + phiName_ + "," + UName_ + ")"
                //         ),
                //         snapUI
                //     )*V;

                convectionDerivative[i][j][k] =
                   -POD::projection
                    (
                        (snapUK & fvc::grad(snapUJ)),
                        snapUI
                    )*V;
            }

            derivative[i][j] =
                POD::projection
                (
                    fvc::laplacian
                    (
                        nu_, snapUJ,
                        "laplacian(" + nu_.name() + "," + UName_ + ")"
                    ),
                    snapUI
                )*V
                // Pressure derivative
              - POD::projection(fvc::grad(snapPJ, "grad(p)"), snapUI)*V;

            // Lagrange multiplier is calculated on boundaries where
            // reconU fixes value
            // Changed form of enforcement of boundary conditions
            // using ddt(bc) = 0
            // HJ, 14/Jul/2021
            forAll (U.boundaryField(), patchI)
            {
                if (U.boundaryField()[patchI].fixesValue())
                {
                    // Note pre-multiplication by beta and signs
                    // of derivative and source.  Su-Sp treatment
                    lagrangeDer[i] += beta_*
                        POD::projection
                        (
                            snapUJ.boundaryField()[patchI],
                            snapUI.boundaryField()[patchI]
                        );

                    // lagrangeSrc[i] += beta_*
                    //     POD::projection
                    //     (
                    //         U.boundaryField()[patchI],
                    //         snapUI.boundaryField()[patchI]
                    //     );
                }
            }
        }
    }

    // Assemble temporal coefficient
    // Note: under normal circumstances, the temporal matrix will be the
    // Cronecker delta.  However, if U is taken out of the larger POD base,
    // this is no longer the case.
    //
    // Note 2: changing the way the drift in the boundary condition is specified
    // via the ddt(b.c.) = 0 condition.  This changes the diagonal matrix
    // of the ddt term and all relevant matrices are re-scaled
    // HJ, 13/Jul/2021
    {
        scalarField diagScale(b.baseSize(), scalar(0));

        for (label i = 0; i < b.baseSize(); i++)
        {
            const volVectorField& snapUI = Ub[i];

            // Add the sqr(snapUI) in case it is not zero
            // Since the POD decomposition is the velocity field itself,
            // this is not strictly needed.  However, for other cases, this
            // is needed.  HJ, 13/Jul/2021
            diagScale[i] = POD::projection(snapUI, snapUI) + lagrangeDer[i];

            Info<< "diagScale " << i << POD::projection(snapUI, snapUI)
                << endl;
        }

        Info<< "diagScale: " << diagScale << endl;

        // Re-scale the derivatives
        for (label i = 0; i < b.baseSize(); i++)
        {
            const scalar curScale = diagScale[i];

            for (label j = 0; j < b.baseSize(); j++)
            {
                derivative[i][j] /= curScale;

                for (label k = 0; k < b.baseSize(); k++)
                {
                    convectionDerivative[i][j][k] /= curScale;
                }
            }
        }
    }

    Info
     // << "convectionDerivative: " << convectionDerivative << nl
     // << "derivative: " << derivative << nl
        << "lagrangeDer: " << lagrangeDer << endl;
    //     << "lagrangeSrc: " << lagrangeSrc << endl;
}


void Foam::incompressibleFlowPOD::updateFields() const
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

        volVectorField& reconU = *reconUPtr_;
        volScalarField& reconP = *reconPPtr_;

        // Reset entire field, including boundary conditions
        reconU == dimensionedVector("zero", reconU.dimensions(), vector::zero);
        reconP == dimensionedScalar("zero", reconP.dimensions(), scalar(0));

        // Velocity field base
        const PtrList<volVectorField>& Ub = UBase();

        // Pressure field base
        const PtrList<volScalarField>& pb = pBase();

        forAll (coeffs_, i)
        {
            // Update velocity
            reconU == reconU + coeffs_[i]*Ub[i];

            // Update pressure
            reconP == reconP + coeffs_[i]*pb[i];
        }

        // Internal field is set.  Correct boundary conditions
        reconU.correctBoundaryConditions();
        reconP.correctBoundaryConditions();

            // Hack test, HJ HERE!!!
            const scalar fluxZero =
                gSum(reconU.boundaryField()[0] & mesh().Sf().boundaryField()[0]);

            Info<< "Testing flux in: "
                << fluxZero << " velocity "
                << fluxZero/gSum(mesh().magSf().boundaryField()[0])
                << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::incompressibleFlowPOD::incompressibleFlowPOD
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
    derivativePtr_(nullptr),
    lagrangeDerPtr_(nullptr),
    lagrangeSrcPtr_(nullptr),
    orthoBasePtr_(nullptr),
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
    coeffs_.setSize(orthoBase().baseSize());

    const scalarRectangularMatrix& orthoBaseCoeffs =
        orthoBase().interpolationCoeffs();

    forAll (coeffs_, i)
    {
        coeffs_[i] = orthoBaseCoeffs[0][i];
    }
    Info<< "Zero coeffs: " << coeffs_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::incompressibleFlowPOD::~incompressibleFlowPOD()
{
    deleteDemandDrivenData(validTimesPtr_);
    deleteDemandDrivenData(convectionDerivativePtr_);
    deleteDemandDrivenData(derivativePtr_);
    deleteDemandDrivenData(lagrangeDerPtr_);
    deleteDemandDrivenData(lagrangeSrcPtr_);

    deleteDemandDrivenData(orthoBasePtr_);
    deleteDemandDrivenData(UBasePtr_);
    deleteDemandDrivenData(pBasePtr_);
    deleteDemandDrivenData(phiBasePtr_);

    deleteDemandDrivenData(reconUPtr_);
    deleteDemandDrivenData(reconPPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::incompressibleFlowPOD::nEqns() const
{
    return coeffs().size();
}


Foam::scalarField& Foam::incompressibleFlowPOD::coeffs()
{
    return coeffs_;
}


const Foam::scalarField& Foam::incompressibleFlowPOD::coeffs() const
{
    return coeffs_;
}


void Foam::incompressibleFlowPOD::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    if
    (
        !convectionDerivativePtr_
     || !derivativePtr_
     || !lagrangeDerPtr_
     || !lagrangeSrcPtr_
    )
    {
        calcDerivativeCoeffs();
    }

    // Convection derivative
    const PtrList<scalarSquareMatrix>& convectionDerivative =
        *convectionDerivativePtr_;

    // Diffusion and pressure derivative
    const scalarSquareMatrix& derivative = *derivativePtr_;

    // Clear derivatives matrix
    dydx = 0;

    forAll (dydx, i)
    {
        // dydx[i] = lagrangeSrc[i];

        forAll (y, j)
        {
            dydx[i] += derivative[i][j]*y[j];

            forAll (y, k)
            {
                dydx[i] += convectionDerivative[i][j][k]*y[j]*y[k];
            }
        }
    }
}


void Foam::incompressibleFlowPOD::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    if (!convectionDerivativePtr_ || !derivativePtr_)
    {
        calcDerivativeCoeffs();
    }

    // Convection derivative
    const PtrList<scalarSquareMatrix>& convectionDerivative =
        *convectionDerivativePtr_;

    // Diffusion and pressure derivative
    const scalarSquareMatrix& derivative = *derivativePtr_;

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

    forAll (y, i)
    {
        forAll (y, j)
        {
            // Convection
            forAll (y, k)
            {
                dfdy[i][j] += convectionDerivative[i][j][k]*y[k];
            }

            // Diffusion
            dfdy[i][j] += derivative[i][j];
        }
    }
}


void Foam::incompressibleFlowPOD::update(const scalar delta)
{
    // Warning: lots of cost here: updating fields for top-level
    // function objects.  Can be removed from production code
    // HJ, 5/Aug/2020
    updateFields();
}


const Foam::PODOrthoNormalBase<Foam::vector4>&
Foam::incompressibleFlowPOD::orthoBase() const
{
    if (!orthoBasePtr_)
    {
        calcOrthoBase();
    }

    return *orthoBasePtr_;
}


const Foam::PtrList<Foam::volVectorField>&
Foam::incompressibleFlowPOD::UBase() const
{
    if (!UBasePtr_)
    {
        calcOrthoBase();
    }

    return *UBasePtr_;
}


const Foam::PtrList<Foam::volScalarField>&
Foam::incompressibleFlowPOD::pBase() const
{
    if (!pBasePtr_)
    {
        calcOrthoBase();
    }

    return *pBasePtr_;
}


const Foam::PtrList<Foam::surfaceScalarField>&
Foam::incompressibleFlowPOD::phiBase() const
{
    if (!phiBasePtr_)
    {
        calcOrthoBase();
    }

    return *phiBasePtr_;
}


const Foam::volVectorField& Foam::incompressibleFlowPOD::reconU() const
{
    updateFields();

    return *reconUPtr_;
}


const Foam::volScalarField& Foam::incompressibleFlowPOD::reconP() const
{
    updateFields();

    return *reconPPtr_;
}


void Foam::incompressibleFlowPOD::writeSnapshots() const
{
    const PODOrthoNormalBase<vector4>& b = orthoBase();
    const PtrList<volVectorField>& Ub = UBase();
    const PtrList<volScalarField>& pb = pBase();

    Info<< "Writing POD base for Time = " << mesh().time().timeName() << endl;

    for (label baseI = 0; baseI < b.baseSize(); baseI++)
    {
        // Write snapshots
        volVectorField snapU
        (
            UName_ + "POD" + name(baseI),
            Ub[baseI]
        );
        snapU.write();

        volScalarField snapP
        (
            pName_ + "POD" + name(baseI),
            pb[baseI]
        );
        snapP.write();
    }
}


void Foam::incompressibleFlowPOD::write() const
{
    // Recalculate field and force a write
    updateFields();
    reconU().write();
    reconP().write();
}


// ************************************************************************* //
