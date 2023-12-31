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

\*---------------------------------------------------------------------------*/

#include "IDDESDelta.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IDDESDelta, 0);
    addToRunTimeSelectionTable(LESdelta, IDDESDelta, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::IDDESDelta::calcDelta()
{
    label nD = mesh().nGeometricD();

    // initialise hwn as wall distance
    volScalarField hwn = wallDist(mesh()).y();

    scalar deltamaxTmp = 0.;

    const cellList& cells = mesh().cells();

    forAll(cells,cellI)
    {
        scalar deltaminTmp = 1.e10;
        const labelList& cFaces = mesh().cells()[cellI];
        const point& centrevector = mesh().cellCentres()[cellI];

        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];
            const point& facevector = mesh().faceCentres()[faceI];
            scalar tmp = mag(facevector - centrevector);

            if (tmp > deltamaxTmp)
            {
                deltamaxTmp = tmp;
            }
            if (tmp < deltaminTmp)
            {
                deltaminTmp = tmp;
            }
        }
        hwn[cellI] = 2.0*deltaminTmp;
    }

    dimensionedScalar deltamax("deltamax",dimLength,2.0*deltamaxTmp);

    if (nD == 3)
    {
        delta_.internalField() =
            deltaCoeff_
           *min
            (
                max(max(cw_*wallDist(mesh()).y(), cw_*deltamax), hwn),
                deltamax
            );
    }
    else if (nD == 2)
    {
        WarningIn("IDDESDelta::calcDelta()")
            << "Case is 2D, LES is not strictly applicable\n"
            << endl;

        delta_.internalField() =
            deltaCoeff_
           *min
            (
                max(max(cw_*wallDist(mesh()).y(), cw_*deltamax), hwn),
                deltamax
            );
    }
    else
    {
        FatalErrorIn("IDDESDelta::calcDelta()")
            << "Case is not 3D or 2D, LES is not strictly applicable"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IDDESDelta::IDDESDelta
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dd
)
:
    LESdelta(name, mesh),
    deltaCoeff_(readScalar(dd.subDict(type() + "Coeffs").lookup("deltaCoeff"))),
    cw_(0)
{
    dd.subDict(type() + "Coeffs").readIfPresent("cw", cw_);
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IDDESDelta::read(const dictionary& dd)
{
    dd.subDict(type() + "Coeffs").lookup("deltaCoeff") >> deltaCoeff_;
    calcDelta();
}


void Foam::IDDESDelta::correct()
{
    if (mesh_.changing())
    {
        calcDelta();
    }
}


// ************************************************************************* //
