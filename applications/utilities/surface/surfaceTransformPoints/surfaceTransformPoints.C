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

Description
    Transform (scale/rotate) a surface. Like transformPoints but then for
    surfaces.

    The rollPitchYaw option takes three angles (degrees):
    - roll (rotation about x) followed by
    - pitch (rotation about y) followed by
    - yaw (rotation about z)

    or -rotateAlongVector (vector and angle in degrees)

    The yawPitchRoll does yaw followed by pitch followed by roll.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "argList.H"
#include "OFstream.H"
#include "IFstream.H"
#include "boundBox.H"
#include "transformField.H"
#include "Pair.H"
#include "quaternion.H"
#include "RodriguesRotation.H"

using namespace Foam;
using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");
    argList::validOptions.insert("translate", "vector");
    argList::validOptions.insert("rotate", "(vector vector)");
    argList::validOptions.insert("rotateAlongVector", "(vector angleInDegree)");
    argList::validOptions.insert("scale", "vector");
    argList::validOptions.insert("rollPitchYaw", "(roll pitch yaw)");
    argList::validOptions.insert("yawPitchRoll", "(yaw pitch roll)");
    argList args(argc, argv);

    fileName surfFileName(args.additionalArgs()[0]);

    Info<< "Reading surf from " << surfFileName << " ..." << endl;

    fileName outFileName(args.additionalArgs()[1]);

    Info<< "Writing surf to " << outFileName << " ..." << endl;

    if (args.options().empty())
    {
        FatalErrorIn(args.executable())
            << "No options supplied, please use one or more of "
               "-translate, -rotate or -scale options."
            << exit(FatalError);
    }

    triSurface surf1(surfFileName);

    pointField points(surf1.points());

    // Translation options

    if (args.optionFound("translate"))
    {
        vector transVector(args.optionLookup("translate")());

        Info<< "Translating points by " << transVector << endl;

        points += transVector;
    }

    // Rotation options

    if (args.optionFound("rotate"))
    {
        Pair<vector> n1n2(args.optionLookup("rotate")());
        n1n2[0] /= mag(n1n2[0]);
        n1n2[1] /= mag(n1n2[1]);

        tensor T = rotationTensor(n1n2[0], n1n2[1]);

        Info<< "Rotating points by " << T << endl;

        points = transform(T, points);
    }
    else if (args.optionFound("rollPitchYaw"))
    {
        vector v(args.optionLookup("rollPitchYaw")());

        Info<< "Rotating points by" << nl
            << "    roll  " << v.x() << nl
            << "    pitch " << v.y() << nl
            << "    yaw   " << v.z() << endl;


        // Convert to radians
        v *= pi/180.0;

        quaternion R(v.x(), v.y(), v.z());

        Info<< "Rotating points by quaternion " << R << endl;
        points = transform(R, points);
    }
    else if (args.optionFound("yawPitchRoll"))
    {
        vector v(args.optionLookup("yawPitchRoll")());

        Info<< "Rotating points by" << nl
            << "    yaw   " << v.x() << nl
            << "    pitch " << v.y() << nl
            << "    roll  " << v.z() << endl;

        // Convert to radians
        v *= pi/180.0;

        scalar yaw = v.x();
        scalar pitch = v.y();
        scalar roll = v.z();

        quaternion R = quaternion(vector(0, 0, 1), yaw);
        R *= quaternion(vector(0, 1, 0), pitch);
        R *= quaternion(vector(1, 0, 0), roll);

        Info<< "Rotating points by quaternion " << R << endl;
        points = transform(R, points);
    }
    else if (args.optionFound("rotateAlongVector"))
    {
        vector rotationAxis;
        scalar rotationAngle;

        args.optionLookup("rotateAlongVector")()
            >> rotationAxis
            >> rotationAngle;

        tensor T = RodriguesRotation(rotationAxis, rotationAngle);

        Info << "Rotating points by " << T << endl;
        points = transform(T, points);
    }

    // Scale options

    if (args.optionFound("scale"))
    {
        vector scaleVector(args.optionLookup("scale")());

        Info<< "Scaling points by " << scaleVector << endl;

        points.replace(vector::X, scaleVector.x()*points.component(vector::X));
        points.replace(vector::Y, scaleVector.y()*points.component(vector::Y));
        points.replace(vector::Z, scaleVector.z()*points.component(vector::Z));
    }

    triSurface surf2(surf1, surf1.patches(), points);

    surf2.write(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
