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
    scalarTransportPOD

\*---------------------------------------------------------------------------*/

#include "scalarTransportPOD.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(scalarTransportPOD, 1);

    addToRunTimeSelectionTable
    (
        PODODE,
        scalarTransportPOD,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::scalarTransportPOD::calcValidTimes() const
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

        // Both field and flux need to be read for the snapshot to be valid
        // Check both

        // Field header
        IOobject fieldHeader
        (
            fieldName_,
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ
        );

        // Flux header
        IOobject phiHeader
        (
            phiName_,
            runTime.timeName(),
            this->mesh(),
            IOobject::MUST_READ
        );

        // Check field and phi header and report separately
        if (fieldHeader.headerOk() && !phiHeader.headerOk())
        {
            Info<< "    Field present but not phi: skipping"
                << endl;
        }

        if (fieldHeader.headerOk() && phiHeader.headerOk())
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


const Foam::instantList& Foam::scalarTransportPOD::validTimes() const
{
    if (!validTimesPtr_)
    {
        calcValidTimes();
    }

    return *validTimesPtr_;
}


void Foam::scalarTransportPOD::calcOrthoBase() const
{
    if (orthoBasePtr_)
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
    PtrList<volScalarField> fields(valTimes.size());

    forAll (valTimes, i)
    {
        runTime.setTime(valTimes[i], i);

        Info<< "    Reading " << fieldName_ << " from time "
            << runTime.timeName() << endl;

        fields.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    fieldName_,
                    runTime.timeName(),
                    mesh(),
                    IOobject::MUST_READ
                ),
                this->mesh()
            )
        );

        // Read the reconstructed field pointer
        if (!reconFieldPtr_)
        {
            Info<< "Setting up " << "recon" << fieldName_
                << " from " << fields[i].name() << endl;
            reconFieldPtr_ =
                new volScalarField
                (
                    IOobject
                    (
                        "recon" + fieldName_,
                        runTime.timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fields[i]
                );
        }
    }

    if (!reconFieldPtr_)
    {
        FatalErrorInFunction
            << "recon" << fieldName_ << " not read"
            << abort(FatalError);
    }

    // Reset time index to initial state
    runTime.setTime(valTimes[0], 0);

    // Create ortho-normal base for transported variable
    orthoBasePtr_ = new scalarPODOrthoNormalBase(fields, accuracy);

    // Check orthogonality and magnitude of snapshots
    if (debug)
    {
        orthoBasePtr_->checkBase();
    }

    if (debug)
    {
        Info<< "Write reconstructed snapshots: check" << endl;
        forAll (valTimes, i)
        {
            runTime.setTime(valTimes[i], i);

            Info<< "Time = " << runTime.timeName() << endl;

            // Field preserves boundary conditions but is reset to zero
            // for accumulation of ortho base values
            volScalarField directRecon
            (
                IOobject
                (
                    "directRecon" + fieldName_,
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                *reconFieldPtr_
            );

            // Reset entire field, including boundary conditions
            directRecon ==
                dimensionedScalar("zero", directRecon.dimensions(), 0);

            // Note: use raw pointer access to orthoBase,
            // as it has just been calculated
            for (label obpI = 0; obpI < orthoBasePtr_->baseSize(); obpI++)
            {
                directRecon == directRecon
                  + orthoBasePtr_->interpolationCoeffs()[i][obpI]*
                    orthoBasePtr_->orthoField(obpI);
            }

            // Internal field is set.  Correct boundary conditions
            directRecon.correctBoundaryConditions();
            directRecon.write();
        }
    }

    // Reset time index to initial state
    runTime.setTime(runTime.times()[origTimeIndex], origTimeIndex);
}


void Foam::scalarTransportPOD::calcDerivativeCoeffs() const
{
    if (derivativePtr_ || lagrangeDerPtr_ || lagrangeSrcPtr_)
    {
        FatalErrorInFunction
            << "Derivative matrix already calculated"
            << abort(FatalError);
    }

    Time& runTime = const_cast<Time&>(this->mesh().time());

    // Remember time index to restore it
    label origTimeIndex = runTime.timeIndex();

    runTime.setTime(validTimes()[0], 0);

    // Read diffusivity

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            this->mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading diffusivity D\n" << endl;

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );

    // Read first available flux field.  Note: flux is fixed for
    // scalar transport

    // Flux field
    surfaceScalarField phi
    (
        IOobject
        (
            phiName_,
            runTime.timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        this->mesh()
    );

    // Reset time index to initial state
    runTime.setTime(runTime.times()[origTimeIndex], origTimeIndex);

    // Create derivative matrix

    const scalarPODOrthoNormalBase& b = orthoBase();

    // Derivative
    derivativePtr_ = new scalarSquareMatrix(b.baseSize(), scalar(0));
    scalarSquareMatrix& derivative = *derivativePtr_;

    // Lagrange multiplier derivative
    lagrangeDerPtr_ = new scalarSquareMatrix(b.baseSize(), scalar(0));
    scalarSquareMatrix& lagrangeDer = *lagrangeDerPtr_;

    // Lagrange multiplier source.  Currently disabled
    lagrangeSrcPtr_ = new scalarField(b.baseSize(), scalar(0));
    scalarField& lagrangeSrc = *lagrangeSrcPtr_;

    const volScalarField& field = *reconFieldPtr_;

    for (label i = 0; i < b.baseSize(); i++)
    {
        const volScalarField& snapI = b.orthoField(i);

        for (label j = 0; j < b.baseSize(); j++)
        {
            const volScalarField& snapJ = b.orthoField(j);

            // Calculate derivative by moving equation terms to rhs
            // Note: both forms of Laplacian work.
            // Oscillations occur when too many snapshots are used.
            // HJ, 20/Jan/2021
            derivative[i][j] =
                POD::projection
                (
                    fvc::laplacian
                    (
                        DT, snapJ,
                        "laplacian(" + DT.name() + "," + fieldName_ + ")"
                    ),
                    snapI
                )
                // -DT.value()*
                // POD::projection
                // (
                //     fvc::grad(snapJ),
                //     fvc::grad(snapI)
                // )
              - POD::projection
                (
                    fvc::div
                    (
                        phi, snapJ,
                        "div(" + phiName_ + "," + fieldName_ + ")"
                    ),
                    snapI
                );

            // Lagrange multiplier is calculated on boundaries where
            // reconField fixes value
            // Changed form of enforcement of boundary conditions
            // using ddt(bc) = 0
            // HJ, 14/Jul/2021
            forAll (field.boundaryField(), patchI)
            {
                if (field.boundaryField()[patchI].fixesValue())
                {
                    // Note pre-multiplication by beta and signs
                    // of derivative and source.  New formulation
                    // HJ, 14/Jul/2021

                    // Accumulated across multiple patches
                    lagrangeDer[i][j] += beta_*
                        POD::projection
                        (
                            snapJ.boundaryField()[patchI],
                            snapI.boundaryField()[patchI]
                        );
                }
            }
        }

        // Calculate Lagrange source
        forAll (field.boundaryField(), patchI)
        {
            if (field.boundaryField()[patchI].fixesValue())
            {
                // Accumulated across multiple patches
                lagrangeSrc[i] += beta_*
                    POD::projection
                    (
                        field.boundaryField()[patchI],
                        snapI.boundaryField()[patchI]
                    );
            }
        }
    }

    // Assemble temporal coefficient
    // Note: under normal circumstances, the temporal matrix will be the
    // Cronecker delta.  However, if the field is taken out of the larger
    // POD base, this is no longer the case.
    //
    // Note 2: changing the way the drift in the boundary condition is specified
    // via the ddt(b.c.) = 0 condition.  This changes the diagonal matrix
    // of the ddt term and all relevant matrices are re-scaled
    // HJ, 13/Jul/2021
    {
        scalarField diagScale(b.baseSize(), scalar(0));

        for (label i = 0; i < b.baseSize(); i++)
        {
            const volScalarField& snapI = b.orthoField(i);

            // Add the sqr(snapUI) in case it is not zero
            // Since the POD decomposition is the velocity field itself,
            // this is not strictly needed.  However, for other cases, this
            // is needed.  HJ, 13/Jul/2021
            diagScale[i] = POD::projection(snapI, snapI) + lagrangeDer[i][i];
        }

        // Re-scale the derivatives
        for (label i = 0; i < b.baseSize(); i++)
        {
            const scalar curScale = diagScale[i];

            for (label j = 0; j < b.baseSize(); j++)
            {
                derivative[i][j] /= curScale;

                // Rescale Lagrange derivatives
                lagrangeDer[i][j] /= curScale;
            }

            // Kill diagonal Lagrange derivative
            lagrangeDer[i][i] = 0;

            // Rescale Lagrange source
            lagrangeSrc[i] /= curScale;
        }
    }
}


void Foam::scalarTransportPOD::updateFields() const
{
    if (fieldUpdateTimeIndex_ < mesh().time().timeIndex())
    {
        // Field update required.  Record update time index
        fieldUpdateTimeIndex_ = mesh().time().timeIndex();

        if (!reconFieldPtr_)
        {
            FatalErrorInFunction
                << "Reconstructed field not allocated"
                << abort(FatalError);
        }

        volScalarField& recon = *reconFieldPtr_;

        const scalarPODOrthoNormalBase& b = orthoBase();

        // Reset entire field, including boundary conditions
        recon == dimensionedScalar("zero", b.orthoField(0).dimensions(), 0);

        forAll (coeffs_, i)
        {
            // Update field
            recon == recon + coeffs_[i]*b.orthoField(i);
        }

        recon.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::scalarTransportPOD::scalarTransportPOD
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PODODE(mesh, dict),
    fieldName_(dict.lookup("field")),
    phiName_(dict.lookup("flux")),
    beta_(readScalar(dict.lookup("beta"))),
    useZeroField_(dict.lookup("useZeroField")),
    driftCorrection_(dict.lookup("driftCorrection")),
    coeffs_(),
    validTimesPtr_(nullptr),
    derivativePtr_(nullptr),
    lagrangeDerPtr_(nullptr),
    lagrangeSrcPtr_(nullptr),
    orthoBasePtr_(nullptr),
    reconFieldPtr_(nullptr),
    fieldUpdateTimeIndex_(-1)
{
    // Check beta
    if (beta_ < 0)
    {
        FatalErrorInFunction
            << "Negative beta: " << beta_
            << ".  Only positive values are allowed"
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

Foam::scalarTransportPOD::~scalarTransportPOD()
{
    deleteDemandDrivenData(validTimesPtr_);
    deleteDemandDrivenData(derivativePtr_);
    deleteDemandDrivenData(lagrangeDerPtr_);
    deleteDemandDrivenData(lagrangeSrcPtr_);

    deleteDemandDrivenData(orthoBasePtr_);

    deleteDemandDrivenData(reconFieldPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::scalarTransportPOD::nEqns() const
{
    return coeffs().size();
}


Foam::scalarField& Foam::scalarTransportPOD::coeffs()
{
    return coeffs_;
}


const Foam::scalarField& Foam::scalarTransportPOD::coeffs() const
{
    return coeffs_;
}


void Foam::scalarTransportPOD::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    if (!derivativePtr_ || !lagrangeDerPtr_ || !lagrangeSrcPtr_)
    {
        calcDerivativeCoeffs();
    }

    // Equation derivative
    const scalarSquareMatrix& derivative = *derivativePtr_;

    // Lagrange multiplier derivative
    const scalarSquareMatrix& lagrangeDer = *lagrangeDerPtr_;

    // Lagrange multiplier source
    const scalarField& lagrangeSrc = *lagrangeSrcPtr_;

    forAll (dydx, i)
    {
        dydx[i] = lagrangeSrc[i];

        forAll (y, j)
        {
            dydx[i] += (derivative[i][j] - lagrangeDer[i][j])*y[j];
        }
    }
}


void Foam::scalarTransportPOD::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    derivatives(x, y, dfdx);

    // Clear jacobian matrix
    // Note: this needs to be done component-by-components
    forAll (y, i)
    {
        forAll (y, j)
        {
            dfdy[i][j] = 0;
        }
    }
}


const Foam::scalarPODOrthoNormalBase&
Foam::scalarTransportPOD::orthoBase() const
{
    if (!orthoBasePtr_)
    {
        calcOrthoBase();
    }

    return *orthoBasePtr_;
}


const Foam::volScalarField& Foam::scalarTransportPOD::reconField() const
{
    updateFields();

    return *reconFieldPtr_;
}


void Foam::scalarTransportPOD::update(const scalar delta)
{
    if (driftCorrection_)
    {
        // Check (and rescale) Dirichlet boundary condition

        // Get ortho-normal base
        const scalarPODOrthoNormalBase& b = orthoBase();

        const volScalarField& field = *reconFieldPtr_;

        forAll (field.boundaryField(), patchI)
        {
            if (field.boundaryField()[patchI].fixesValue())
            {
                scalarField reconBC
                (
                    field.boundaryField()[patchI].size(),
                    scalar(0)
                );

                forAll (coeffs_, i)
                {
                    reconBC +=
                        coeffs_[i]*b.orthoField(i).boundaryField()[patchI];
                }

                Info<< "Patch " << patchI
                    << " fieldBC = "
                    << refCast<const scalarField>(field.boundaryField()[patchI])
                    << " Coeffs: " << coeffs_
                    << " reconBC = " << reconBC
                    << endl;

                // Adjust leading coefficient for drift

                // Field average
                const scalar avgFieldBC =
                    average(field.boundaryField()[patchI]);

                // Zeroth ortho base average
                const scalar avgBase0BC =
                    average(b.orthoField(0).boundaryField()[patchI]);

                if
                (
                    !field.boundaryField()[patchI].empty()
                 && mag(avgBase0BC) > SMALL
                )
                {
                    scalar delta = (avgFieldBC - average(reconBC))/avgBase0BC;

                    Info<< "Correction on patch " << patchI
                        << " coeff[0] = " << coeffs_[0]
                        << " delta: " << delta
                        << endl;

                    coeffs_[0] += delta;
                }
            }
        }
    }
}

void Foam::scalarTransportPOD::writeSnapshots() const
{
    const scalarPODOrthoNormalBase& b = orthoBase();

    Info<< "Writing POD base for Time = " << mesh().time().timeName() << endl;

    for (label i = 0; i < b.baseSize(); i++)
    {
        b.orthoField(i).write();
    }
}


void Foam::scalarTransportPOD::write() const
{
    // Recalculate field and force a write
    updateFields();

    reconField().write();
}


// ************************************************************************* //
