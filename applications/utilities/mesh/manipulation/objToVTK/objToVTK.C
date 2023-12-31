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
    Read obj line (not surface!) file and convert into vtk.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"
#include <fstream>
#include <sstream>
#include "IStringStream.H"
#include "point.H"
#include "DynamicList.H"


using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

string getLine(std::ifstream& is)
{
    string line;
    do
    {
        std::getline(is, line);
    }
    while (line.size() && line[0] == '#');

    return line;
}


// Read space-separated vertices (with optional '/' arguments)
labelList parseVertices(const string& line)
{
    dynamicLabelList verts;

    // Assume 'l' is followed by space.
    string::size_type endNum = 1;

    do
    {
        string::size_type startNum = line.find_first_not_of(' ', endNum);

        if (startNum == string::npos)
        {
            break;
        }

        endNum = line.find(' ', startNum);

        string vertexSpec;
        if (endNum != string::npos)
        {
            vertexSpec = line.substr(startNum, endNum-startNum);
        }
        else
        {
            vertexSpec = line.substr(startNum, line.size() - startNum);
        }

        string::size_type slashPos = vertexSpec.find('/');

        label vertI = 0;
        if (slashPos != string::npos)
        {
            IStringStream intStream(vertexSpec.substr(0, slashPos));

            intStream >> vertI;
        }
        else
        {
            IStringStream intStream(vertexSpec);

            intStream >> vertI;
        }
        verts.append(vertI - 1);
    }
    while (true);

    return verts.shrink();
}


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("OBJ file");
    argList::validArgs.append("output VTK file");
    argList args(argc, argv);

    fileName objName(args.additionalArgs()[0]);
    fileName outName(args.additionalArgs()[1]);

    std::ifstream OBJfile(objName.c_str());

    if (!OBJfile.good())
    {
        FatalErrorIn(args.executable())
            << "Cannot read file " << objName << exit(FatalError);
    }

    // Points and lines
    DynamicList<point> points;
    DynamicList<labelList> polyLines;
    DynamicList<labelList> polygons;

    bool hasWarned = false;

    label lineNo = 0;
    while (OBJfile.good())
    {
        string line = getLine(OBJfile);
        lineNo++;

        // Read first word
        IStringStream lineStream(line);
        word cmd;
        lineStream >> cmd;

        if (cmd == "v")
        {
            scalar x, y, z;

            lineStream >> x >> y >> z;

            points.append(point(x, y, z));
        }
        else if (cmd == "l")
        {
            polyLines.append(parseVertices(line));
        }
        else if (cmd == "f")
        {
            polygons.append(parseVertices(line));
        }
        else if (cmd != "")
        {
            if (!hasWarned)
            {
                hasWarned = true;

                WarningIn(args.executable())
                    << "Unrecognized OBJ command " << cmd << nl
                    << "In line " << lineStream.str()
                    << " at linenumber " << lineNo << nl
                    << "Only recognized commands are 'v' and 'l'.\n"
                    << "If this is a surface command use surfaceConvert instead"
                    << " to convert to a file format that can be read by VTK"
                    << endl;
            }
        }
    }


    //
    // Write as vtk 'polydata' file
    //


    OFstream outFile(outName);

    outFile
        << "# vtk DataFile Version 2.0\n"
        << objName << nl
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS " << points.size() << " float\n";

    forAll(points, i)
    {
        const point& pt = points[i];

        outFile << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }

    label nItems = 0;
    forAll(polyLines, polyI)
    {
        nItems += polyLines[polyI].size() + 1;
    }

    outFile
        << "LINES " << polyLines.size() << ' ' << nItems << nl;

    forAll(polyLines, polyI)
    {
        const labelList& line = polyLines[polyI];

        outFile << line.size();

        forAll(line, i)
        {
            outFile << ' ' << line[i];
        }
        outFile << nl;
    }


    nItems = 0;
    forAll(polygons, polyI)
    {
        nItems += polygons[polyI].size() + 1;
    }

    outFile
        << "POLYGONS " << polygons.size() << ' ' << nItems << nl;

    forAll(polygons, polyI)
    {
        const labelList& line = polygons[polyI];

        outFile << line.size();

        forAll(line, i)
        {
            outFile << ' ' << line[i];
        }
        outFile << nl;
    }


    outFile
        << "POINT_DATA " << points.size() << nl
        << "SCALARS pointID float 1\n"
        << "LOOKUP_TABLE default\n";

    forAll(points, i)
    {
        outFile << i;

        if ((i % 10) == 1)
        {
            outFile << nl;
        }
        else
        {
            outFile << ' ';
        }
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
