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
    Sample field data with a choice of interpolation schemes, sampling options
    and write formats.

    Keywords:

    @param setFormat : set output format, choice of \n
      - xmgr
      - jplot
      - gnuplot
      - raw

    @param surfaceFormat : surface output format, choice of \n
      - null        : suppress output
      - foamFile    : separate points, faces and values file
      - dx          : DX scalar or vector format
      - vtk         : VTK ascii format
      - raw         : x y z value format for use with e.g. gnuplot 'splot'.
      - obj         : Wavefron stl. Does not contain values!
      - stl         : ascii stl. Does not contain values!

    @param interpolationScheme : interpolation scheme, choice of \n
      - cell          : use cell-centre value; constant over cells (default)
      - cellPoint     : use cell-centre and vertex values
      - cellPointFace : use cell-centre, vertex and face values. \n
        -# vertex values determined from neighbouring cell-centre values
        -# face values determined using the current face interpolation scheme
           for the field (linear, limitedLinear, etc.)

    @param fields : list of fields to sample

    @param sets : list of sets to sample, choice of \n
      - uniform             evenly distributed points on line
      - face                one point per face intersection
      - midPoint            one point per cell, inbetween two face intersections
      - midPointAndFace     combination of face and midPoint

      - curve               specified points, not nessecary on line, uses
                            tracking
      - cloud               specified points, uses findCell

        Option axis: how to write point coordinate. Choice of
          - x/y/z: x/y/z coordinate only
          - xyz: three columns
            (probably does not make sense for anything but raw)
          - distance: distance from start of sampling line (if uses line)
            or distance from first specified sampling point

        Type specific options:
            uniform, face, midPoint, midPointAndFace : start and end coordinate
            uniform: extra number of sampling points
            curve, cloud: list of coordinates

    @param surfaces : list of surfaces to sample, choice of \n
      - plane : values on plane defined by point, normal.
      - patch : values on patch.

Notes
    Runs in parallel

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOsampledSets.H"
#include "IOsampledSurfaces.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"

    IOsampledSets sSets
    (
        sampledSets::typeName,
        mesh,
        word("sampleDict"),
        IOobject::MUST_READ,
        true
    );

    IOsampledSurfaces sSurfs
    (
        sampledSurfaces::typeName,
        mesh,
        word("sampleDict"),
        IOobject::MUST_READ,
        true
    );

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        // Handle geometry/topology changes
        polyMesh::readUpdateState state = mesh.readUpdate();

        sSets.readUpdate(state);
        sSurfs.readUpdate(state);

        sSets.write();
        sSurfs.write();

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
