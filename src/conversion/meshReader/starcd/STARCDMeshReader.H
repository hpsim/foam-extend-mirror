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
    Foam::meshReaders::STARCD

Description
    Read pro-STAR vrt/cel/bnd files.
    The protected data in meshReader are filled.

    Starting with pro-STAR version 4, the files have become easier to read.
    - vertices are space-delimited.
    - the cell format is logical.
    - trimmed and degenerate cells are saved as polyhedral.
    - the boundaries corresponds to cells and their faces.

SourceFiles
    STARCDMeshReader.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDMeshReader_H
#define STARCDMeshReader_H

# include "meshReader.H"
# include "boundaryRegion.H"
# include "cellShape.H"
# include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace meshReaders
{

/*---------------------------------------------------------------------------*\
                    Class meshReaders::STARCD Declaration
\*---------------------------------------------------------------------------*/

class STARCD
:
    public meshReader
{

protected:

    // Protected Data

        static const char* defaultBoundaryName;
        static const char* defaultSolidBoundaryName;

        //- Face addressing from pro-STAR faces -> foam faces
        static const int starToFoamFaceAddr[4][6];

        //- Cell shapes
        cellShapeList cellShapes_;

        //- Point labels (imported Point numbering not necessarily contiguous)
        labelList mapToFoamPointId_;

        //- Cell labels (imported Cell numbering not necessarily contiguous)
        labelList mapToFoamCellId_;

        //- Boundary region data
        boundaryRegion boundaryRegion_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        STARCD(const STARCD&);

        //- Disallow default bitwise assignment
        void operator=(const STARCD&);

        //- Read the mesh from the file(s)
        virtual bool readGeometry(const scalar scaleFactor = 1.0);

        //- Read points from file
        void readPoints(const fileName&, const scalar scaleFactor);

        //- Read cell connectivities from file
        virtual void readCells(const fileName&);

        //- Remove unused points
        void cullPoints();

        //- Read boundary (cell/face) definitions
        void readBoundary(const fileName&);

        //- Read auxiliary data from constant/{boundaryRegion,cellTable}
        void readAux(const objectRegistry&);

        //- Read and discard to newline
        static void readToNewline(IFstream&);

        //- Read header
        static bool readHeader(IFstream&, word fileSignature);


protected:

        enum cellType
        {
            starcdFluidType   = 1,
            starcdSolidType   = 2,
            starcdBaffleType  = 3,
            starcdShellType   = 4,
            starcdLineType    = 5,
            starcdPointType   = 6
        };

        enum shapeType
        {
            starcdPoint = 1,
            starcdLine  = 2,
            starcdShell = 3,
            starcdHex   = 11,
            starcdPrism = 12,
            starcdTet   = 13,
            starcdPyr   = 14,
            starcdPoly  = 255
        };


public:

    // Static data

        //- Keep solids (default false)
        static bool keepSolids;


    // Constructors

        //- Construct from case name
        STARCD
        (
            const fileName& prefix,
            const objectRegistry&,
            const scalar scaleFactor = 1.0
        );


    //- Destructor
    virtual ~STARCD();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meshReaders
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
