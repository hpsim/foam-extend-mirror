/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
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
    thermalModel

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "thermalModel.H"
#include "volFields.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalModel, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalModel::thermalModel(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "thermalProperties",
            T.time().constant(),
            T.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    T_(T),
    lawPtr_(thermalLaw::New("law", T_, subDict("thermal"))),
    source_(subDict("thermal"), T)
{
    {
        PtrList<entry> entries(subDict("thermal").lookup("gaps"));
        gapPtr_.setSize(entries.size());

        forAll (gapPtr_, gapI)
        {
            gapPtr_.set
            (
                gapI,
                thermalGap::New
                (
                    entries[gapI].keyword(),
                    T,
                    entries[gapI].dict()
                )
            );
        }
    }

    Info<< "SOURCE ACTIVE?: " << source_.active() << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermalModel::modifyResistance
(
    surfaceScalarField& kf
) const
{
    forAll(gapPtr_, gapI)
    {
        gapPtr_[gapI].modifyResistance(kf);
    }
}


tmp<fvScalarMatrix> thermalModel::laplacian(volScalarField& T)
{
    const word kScheme("laplacian(k,T)");

    surfaceScalarField kf = fvc::interpolate(lawPtr_->k());

    modifyResistance(kf);

    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix(fvm::laplacian(kf, T, kScheme))
    );
}


bool thermalModel::read()
{
    if (regIOobject::read())
    {
        lawPtr_ = thermalLaw::New("law", T_, subDict("thermal"));

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
