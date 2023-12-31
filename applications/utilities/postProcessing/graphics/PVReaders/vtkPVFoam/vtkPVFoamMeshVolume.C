#include "vtkPVFoam.H"
#include "vtkPVFoamReader.h"

// FOAM includes
#include "fvMesh.H"
#include "cellModeller.H"
#include "vtkPVFoamPoints.H"
#include "Swap.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkIdTypeArray.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vtkUnstructuredGrid* Foam::vtkPVFoam::volumeVTKMesh
(
    const fvMesh& mesh,
    polyDecomp& decompInfo
)
{
    const cellModel& tet = *(cellModeller::lookup("tet"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& wedge = *(cellModeller::lookup("wedge"));
    const cellModel& tetWedge = *(cellModeller::lookup("tetWedge"));
    const cellModel& hex = *(cellModeller::lookup("hex"));

    vtkUnstructuredGrid* vtkmesh = vtkUnstructuredGrid::New();

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::volumeVTKMesh" << endl;
        printMemory();
    }

    const cellShapeList& cellShapes = mesh.cellShapes();

    // Number of additional points needed by the decomposition of polyhedra
    label nAddPoints = 0;

    // Number of additional cells generated by the decomposition of polyhedra
    label nAddCells = 0;

    // face owner is needed to determine the face orientation
    const labelList& owner = mesh.faceOwner();

    labelList& superCells = decompInfo.superCells();
    labelList& addPointCellLabels = decompInfo.addPointCellLabels();

    // Scan for cells which need to be decomposed and count additional points
    // and cells
    if (!reader_->GetUseVTKPolyhedron())
    {
        if (debug)
        {
            Info<< "... scanning for polyhedra" << endl;
        }

        forAll(cellShapes, celli)
        {
            const cellModel& model = cellShapes[celli].model();

            if
            (
                model != hex
             && model != wedge
             && model != prism
             && model != pyr
             && model != tet
             && model != tetWedge
            )
            {
                const cell& cFaces = mesh.cells()[celli];

                forAll(cFaces, cFacei)
                {
                    const face& f = mesh.faces()[cFaces[cFacei]];

                    label nQuads = 0;
                    label nTris = 0;
                    f.nTrianglesQuads(mesh.points(), nTris, nQuads);

                    nAddCells += nQuads + nTris;
                }

                nAddCells--;
                nAddPoints++;
            }
        }
    }

    // Set size of additional point addressing array
    // (from added point to original cell)
    addPointCellLabels.setSize(nAddPoints);

    // Set size of additional cells mapping array
    // (from added cell to original cell)

    if (debug)
    {
        Info<<" mesh nCells     = " << mesh.nCells() << nl
            <<"      nPoints    = " << mesh.nPoints() << nl
            <<"      nAddCells  = " << nAddCells << nl
            <<"      nAddPoints = " << nAddPoints << endl;
    }

    superCells.setSize(mesh.nCells() + nAddCells);

    if (debug)
    {
        Info<< "... converting points" << endl;
    }

    // Convert OpenFOAM mesh vertices to VTK
    vtkPoints* vtkpoints = vtkPoints::New();
    vtkpoints->Allocate(mesh.nPoints() + nAddPoints);

    const Foam::pointField& points = mesh.points();

    forAll(points, i)
    {
        vtkPVFoamInsertNextPoint(vtkpoints, points[i]);
    }


    if (debug)
    {
        Info<< "... converting cells" << endl;
    }

    vtkmesh->Allocate(mesh.nCells() + nAddCells);

    // Set counters for additional points and additional cells
    label addPointi = 0, addCelli = 0;

    // Create storage for points - needed for mapping from OpenFOAM to VTK
    // data types - max 'order' = hex = 8 points
    vtkIdType nodeIds[8];

    // face-stream for a polyhedral cell
    // [numFace0Pts, id1, id2, id3, numFace1Pts, id1, id2, id3, ...]
    DynamicList<vtkIdType> faceStream(256);

    forAll(cellShapes, celli)
    {
        const cellShape& cellShape = cellShapes[celli];
        const cellModel& cellModel = cellShape.model();

        superCells[addCelli++] = celli;

        if (cellModel == tet)
        {
            for (int j = 0; j < 4; j++)
            {
                nodeIds[j] = cellShape[j];
            }
            vtkmesh->InsertNextCell
            (
                VTK_TETRA,
                4,
                nodeIds
            );
        }
        else if (cellModel == pyr)
        {
            for (int j = 0; j < 5; j++)
            {
                nodeIds[j] = cellShape[j];
            }
            vtkmesh->InsertNextCell
            (
                VTK_PYRAMID,
                5,
                nodeIds
            );
        }
        else if (cellModel == prism)
        {
            // VTK has a different node order for VTK_WEDGE
            // their triangles point outwards!
            nodeIds[0] = cellShape[0];
            nodeIds[1] = cellShape[2];
            nodeIds[2] = cellShape[1];
            nodeIds[3] = cellShape[3];
            nodeIds[4] = cellShape[5];
            nodeIds[5] = cellShape[4];

            vtkmesh->InsertNextCell
            (
                VTK_WEDGE,
                6,
                nodeIds
            );
        }
        else if (cellModel == tetWedge && !reader_->GetUseVTKPolyhedron())
        {
            // Treat as squeezed prism (VTK_WEDGE)

            nodeIds[0] = cellShape[0];
            nodeIds[1] = cellShape[2];
            nodeIds[2] = cellShape[1];
            nodeIds[3] = cellShape[3];
            nodeIds[4] = cellShape[4];
            nodeIds[5] = cellShape[3];

            vtkmesh->InsertNextCell
            (
                VTK_WEDGE,
                6,
                nodeIds
            );
        }
        else if (cellModel == wedge)
        {
            // Treat as squeezed hex

            nodeIds[0] = cellShape[0];
            nodeIds[1] = cellShape[1];
            nodeIds[2] = cellShape[2];
            nodeIds[3] = cellShape[2];
            nodeIds[4] = cellShape[3];
            nodeIds[5] = cellShape[4];
            nodeIds[6] = cellShape[5];
            nodeIds[7] = cellShape[6];

            vtkmesh->InsertNextCell
            (
                VTK_HEXAHEDRON,
                8,
                nodeIds
            );
        }
        else if (cellModel == hex)
        {
            for (int j = 0; j < 8; j++)
            {
                nodeIds[j] = cellShape[j];
            }
            vtkmesh->InsertNextCell
            (
                VTK_HEXAHEDRON,
                8,
                nodeIds
            );
        }
        else if (reader_->GetUseVTKPolyhedron())
        {
            // Polyhedral cell - use VTK_POLYHEDRON
            const labelList& cFaces = mesh.cells()[celli];

#ifdef HAS_VTK_POLYHEDRON
            vtkIdType nFaces = cFaces.size();
            vtkIdType nLabels = nFaces;

            // count size for face stream
            forAll(cFaces, cFacei)
            {
                const face& f = mesh.faces()[cFaces[cFacei]];
                nLabels += f.size();
            }

            // build face-stream
            // [numFace0Pts, id1, id2, id3, numFace1Pts, id1, id2, id3, ...]
            // point Ids are global
            faceStream.clear();
            faceStream.reserve(nLabels + nFaces);

            forAll(cFaces, cFacei)
            {
                const face& f = mesh.faces()[cFaces[cFacei]];
                const bool isOwner = (owner[cFaces[cFacei]] == celli);
                const label nFacePoints = f.size();

                // number of labels for this face
                faceStream.append(nFacePoints);

                if (isOwner)
                {
                    forAll(f, fp)
                    {
                        faceStream.append(f[fp]);
                    }
                }
                else
                {
                    // fairly immaterial if we reverse the list
                    // or use face::reverseFace()
                    forAllReverse(f, fp)
                    {
                        faceStream.append(f[fp]);
                    }
                }
            }

            vtkmesh->InsertNextCell(VTK_POLYHEDRON, nFaces, faceStream.data());
#else
            // this is a horrible substitute
            // but avoids crashes when there is no vtkPolyhedron support

            // establish unique node ids used
            HashSet<vtkIdType, Hash<label> > hashUniqId(2*256);

            forAll(cFaces, cFacei)
            {
                const face& f = mesh.faces()[cFaces[cFacei]];

                forAll(f, fp)
                {
                    hashUniqId.insert(f[fp]);
                }
            }

            // use face stream to store unique node ids:
            faceStream = hashUniqId.sortedToc();

            vtkmesh->InsertNextCell
            (
                VTK_CONVEX_POINT_SET,
                vtkIdType(faceStream.size()),
                faceStream.data()
            );
#endif
        }
        else
        {
            // Polyhedral cell. Decompose into tets + prisms.

            // Mapping from additional point to cell
            addPointCellLabels[addPointi] = celli;

            // The new vertex from the cell-centre
            const label newVertexLabel = mesh.nPoints() + addPointi;
            vtkPVFoamInsertNextPoint(vtkpoints, mesh.C()[celli]);

            // Whether to insert cell in place of original or not.
            bool substituteCell = true;

            const labelList& cFaces = mesh.cells()[celli];
            forAll(cFaces, cFacei)
            {
                const face& f = mesh.faces()[cFaces[cFacei]];
                const bool isOwner = (owner[cFaces[cFacei]] == celli);

                // Number of triangles and quads in decomposition
                label nTris = 0;
                label nQuads = 0;
                f.nTrianglesQuads(mesh.points(), nTris, nQuads);

                // Do actual decomposition into triFcs and quadFcs.
                faceList triFcs(nTris);
                faceList quadFcs(nQuads);
                label trii = 0;
                label quadi = 0;
                f.trianglesQuads(mesh.points(), trii, quadi, triFcs, quadFcs);

                forAll(quadFcs, quadI)
                {
                    if (substituteCell)
                    {
                        substituteCell = false;
                    }
                    else
                    {
                        superCells[addCelli++] = celli;
                    }

                    const face& quad = quadFcs[quadI];

                    // Ensure we have the correct orientation for the
                    // base of the primitive cell shape.
                    // If the cell is face owner, the orientation needs to be
                    // flipped.
                    // At the moment, VTK doesn't actually seem to care if
                    // negative cells are defined, but we'll do it anyhow
                    // (for safety).
                    if (isOwner)
                    {
                        nodeIds[0] = quad[3];
                        nodeIds[1] = quad[2];
                        nodeIds[2] = quad[1];
                        nodeIds[3] = quad[0];
                    }
                    else
                    {
                        nodeIds[0] = quad[0];
                        nodeIds[1] = quad[1];
                        nodeIds[2] = quad[2];
                        nodeIds[3] = quad[3];
                    }
                    nodeIds[4] = newVertexLabel;
                    vtkmesh->InsertNextCell
                    (
                        VTK_PYRAMID,
                        5,
                        nodeIds
                    );
                }

                forAll(triFcs, triI)
                {
                    if (substituteCell)
                    {
                        substituteCell = false;
                    }
                    else
                    {
                        superCells[addCelli++] = celli;
                    }

                    const face& tri = triFcs[triI];

                    // See note above about the orientation.
                    if (isOwner)
                    {
                        nodeIds[0] = tri[2];
                        nodeIds[1] = tri[1];
                        nodeIds[2] = tri[0];
                    }
                    else
                    {
                        nodeIds[0] = tri[0];
                        nodeIds[1] = tri[1];
                        nodeIds[2] = tri[2];
                    }
                    nodeIds[3] = newVertexLabel;

                    vtkmesh->InsertNextCell
                    (
                        VTK_TETRA,
                        4,
                        nodeIds
                    );
                }
            }

            addPointi++;
        }
    }

    vtkmesh->SetPoints(vtkpoints);
    vtkpoints->Delete();

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::volumeVTKMesh" << endl;
        printMemory();
    }

    return vtkmesh;
}


// ************************************************************************* //
