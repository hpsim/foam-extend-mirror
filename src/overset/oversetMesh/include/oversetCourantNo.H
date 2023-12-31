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

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers with overset
    mesh correction.

\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;
scalar velMag = 0.0;

if (mesh.nInternalFaces())
{
    const oversetMesh& om = oversetMesh::New(mesh);

    const surfaceScalarField magPhi = mag(om.sGamma()*phi);

    const surfaceScalarField SfUfbyDelta =
        mesh.surfaceInterpolation::deltaCoeffs()*magPhi;

    CoNum = max(SfUfbyDelta/mesh.magSf())
        .value()*runTime.deltaT().value();

    meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf()))
        .value()*runTime.deltaT().value();

    velMag = max(magPhi/mesh.magSf()).value();
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum
    << " velocity magnitude: " << velMag
    << endl;

// ************************************************************************* //
