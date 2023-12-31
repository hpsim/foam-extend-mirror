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
    Foam::vtkPVFoam

Description
    Provides a reader interface for FOAM to VTK interaction.

SourceFiles
    vtkPVFoam.C
    vtkPVFoam.H
    vtkPVFoamI.H
    vtkPVFoamFields.C
    vtkPVFoamMesh.C
    vtkPVFoamMeshLagrangian.C
    vtkPVFoamMeshPatch.C
    vtkPVFoamMeshSet.C
    vtkPVFoamMeshVolume.C
    vtkPVFoamMeshZone.C
    vtkPVFoamFaceField.H
    vtkPVFoamLagrangianFields.H
    vtkPVFoamPatchField.H
    vtkPVFoamPointFields.H
    vtkPVFoamPoints.H
    vtkPVFoamUpdateInfo.C
    vtkPVFoamUpdateInfoFields.H
    vtkPVFoamUtilities.C
    vtkPVFoamVolFields.H
    vtkPVFoamAddToSelection.H

    // Needed by VTK:
    vtkDataArrayTemplateImplicit.txx

\*---------------------------------------------------------------------------*/

#ifndef vtkPVFoam_H
#define vtkPVFoam_H

// do not include legacy strstream headers
#ifndef  VTK_EXCLUDE_STRSTREAM_HEADERS
# define VTK_EXCLUDE_STRSTREAM_HEADERS
#endif

#include "className.H"
#include "fileName.H"
#include "stringList.H"
#include "wordList.H"
#include "primitivePatch.H"
#include "PrimitivePatchInterpolationTemplate.H"
#include "volPointInterpolation.H"

#undef VTKPVFOAM_DUALPORT

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkDataArraySelection;
class vtkDataSet;
class vtkPoints;
class vtkPVFoamReader;
class vtkRenderer;
class vtkTextActor;
class vtkMultiBlockDataSet;
class vtkPolyData;
class vtkUnstructuredGrid;
class vtkIndent;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Foam class forward declarations
class argList;
class Time;
class fvMesh;
class IOobjectList;
class polyPatch;
class faceSet;
class pointSet;

template<class Type> class IOField;
template<class Type> class List;

/*---------------------------------------------------------------------------*\
                        Class vtkPVFoam Declaration
\*---------------------------------------------------------------------------*/

class vtkPVFoam
{
    // Private classes

        //- Bookkeeping for GUI checklists and the multi-block organization
        class partInfo
        {
            const char *name_;
            int block_;
            int start_;
            int size_;

        public:

            partInfo(const char *name, const int blockNo=0)
            :
                name_(name),
                block_(blockNo),
                start_(-1),
                size_(0)
            {}

            //- Return the block holding these datasets
            int block() const
            {
                return block_;
            }

            //- Assign block number, return previous value
            int block(int blockNo)
            {
                int prev = block_;
                block_ = blockNo;
                return prev;
            }

            const char* name() const
            {
                return name_;
            }

            int start() const
            {
                return start_;
            }

            int end() const
            {
                return start_ + size_;
            }

            int size() const
            {
                return size_;
            }

            bool empty() const
            {
                return !size_;
            }

            void reset()
            {
                start_ = -1;
                size_ = 0;
            }

            //- Assign new start and reset the size
            void operator=(const int i)
            {
                start_ = i;
                size_ = 0;
            }

            //- Increment the size
            void operator+=(const int n)
            {
                size_ += n;
            }
        };

        //- bookkeeping for polyhedral cell decomposition
        //  hide in extra pointMap (cellSet/cellZone) for now
        class polyDecomp
        {
            labelList superCells_;
            labelList addPointCellLabels_;
            labelList pointMap_;

        public:

            polyDecomp()
            {}

            //- Label of original cell for decomposed cells
            labelList& superCells()
            {
                return superCells_;
            }

            //- Label of original cell for decomposed cells
            const labelList& superCells() const
            {
                return superCells_;
            }

            //- Cell-centre labels for additional points of decomposed cells
            labelList& addPointCellLabels()
            {
                return addPointCellLabels_;
            }

            //- Cell-centre labels for additional points of decomposed cells
            const labelList& addPointCellLabels() const
            {
                return addPointCellLabels_;
            }

            //- Point labels for subsetted meshes
            labelList& pointMap()
            {
                return pointMap_;
            }

            //- Point labels for subsetted meshes
            const labelList& pointMap() const
            {
                return pointMap_;
            }


            //- Clear
            void clear()
            {
                superCells_.clear();
                addPointCellLabels_.clear();
                pointMap_.clear();
            }
        };


    // Private Data

        //- Access to the controlling vtkPVFoamReader
        vtkPVFoamReader* reader_;

        //- Foam time control
        autoPtr<Time> dbPtr_;

        //- Foam mesh
        fvMesh* meshPtr_;

        //- The mesh region
        word meshRegion_;

        //- The mesh directory for the region
        fileName meshDir_;

        //- The time index
        int timeIndex_;

        //- Track changes in mesh geometry
        bool meshChanged_;

        //- Track changes in fields
        bool fieldsChanged_;

        //- Selected geometrical parts (internalMesh, patches, ...)
        boolList partStatus_;

        //- Datasets corresponding to selected geometrical pieces
        //  a negative number indicates that no vtkmesh exists for this piece
        labelList partDataset_;

        //- First instance and size of various mesh parts
        //  used to index into partStatus_ and partDataset_
        partInfo partInfoVolume_;
        partInfo partInfoPatches_;
        partInfo partInfoLagrangian_;
        partInfo partInfoCellZones_;
        partInfo partInfoFaceZones_;
        partInfo partInfoPointZones_;
        partInfo partInfoCellSets_;
        partInfo partInfoFaceSets_;
        partInfo partInfoPointSets_;

        //- Decomposed cells information (mesh regions)
        //  TODO: regions
        List<polyDecomp> regionPolyDecomp_;

        //- Decomposed cells information (cellZone meshes)
        List<polyDecomp> zonePolyDecomp_;

        //- Decomposed cells information (cellSet meshes)
        List<polyDecomp> csetPolyDecomp_;

        //- List of patch names for rendering to window
        List<vtkTextActor*> patchTextActorsPtrs_;

    // Private Member Functions

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static void AddToBlock
        (
            vtkMultiBlockDataSet* output,
            vtkDataSet* dataset,
            const partInfo&,
            const label datasetNo,
            const string& datasetName
        );

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static vtkDataSet* GetDataSetFromBlock
        (
            vtkMultiBlockDataSet* output,
            const partInfo&,
            const label datasetNo
        );

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static label GetNumberOfDataSets
        (
            vtkMultiBlockDataSet* output,
            const partInfo&
        );

        //- Reset data counters
        void resetCounters();

        // Update information helper functions

            //- Update the mesh parts selected in the GUI
            void updateMeshPartsStatus();

            //- Internal mesh info
            void updateInfoInternalMesh();

            //- Lagrangian info
            void updateInfoLagrangian();

            //- Patch info
            void updateInfoPatches();

            //- Set info
            void updateInfoSets();

            //- Zone info
            void updateInfoZones();

            //- Read zone names for zoneType from file
            wordList readZoneNames(const word& zoneType);

            //- Add objects of Type to paraview array selection
            template<class Type>
            label addToSelection
            (
                vtkDataArraySelection*,
                const IOobjectList&,
                const string& suffix=string::null
            );

            //- Field info
            template<template<class> class patchType, class meshType>
            void updateInfoFields(vtkDataArraySelection*);

            //- Lagrangian field info
            void updateInfoLagrangianFields();


        // Update helper functions

            //- Foam mesh
            void updateFoamMesh();

            //- Reduce memory footprint after conversion
            void reduceMemory();

            //- Volume fields
            void updateVolFields(vtkMultiBlockDataSet*);

            //- Point fields
            void updatePointFields(vtkMultiBlockDataSet*);

            //- Lagrangian fields
            void updateLagrangianFields(vtkMultiBlockDataSet*);


        // Mesh conversion functions

            //- Volume mesh
            void convertMeshVolume(vtkMultiBlockDataSet*, int& blockNo);

            //- Lagrangian mesh
            void convertMeshLagrangian(vtkMultiBlockDataSet*, int& blockNo);

            //- Patch meshes
            void convertMeshPatches(vtkMultiBlockDataSet*, int& blockNo);

            //- Cell zone meshes
            void convertMeshCellZones(vtkMultiBlockDataSet*, int& blockNo);

            //- Face zone meshes
            void convertMeshFaceZones(vtkMultiBlockDataSet*, int& blockNo);

            //- Point zone meshes
            void convertMeshPointZones(vtkMultiBlockDataSet*, int& blockNo);

            //- Cell set meshes
            void convertMeshCellSets(vtkMultiBlockDataSet*, int& blockNo);

            //- Face set meshes
            void convertMeshFaceSets(vtkMultiBlockDataSet*, int& blockNo);

            //- Point set meshes
            void convertMeshPointSets(vtkMultiBlockDataSet*, int& blockNo);


        // Add mesh functions

            //- Add internal mesh/cell set meshes
            vtkUnstructuredGrid* volumeVTKMesh(const fvMesh&, polyDecomp&);

            //- Add Lagrangian mesh
            vtkPolyData* lagrangianVTKMesh
            (
                const fvMesh&,
                const word& cloudName
            );

            //- Add patch mesh
            vtkPolyData* patchVTKMesh(const polyPatch&);

            //- Add face zone mesh
            vtkPolyData* faceZoneVTKMesh
            (
                const fvMesh&,
                const labelList& faceLabels
            );

            //- Add point zone
            vtkPolyData* pointZoneVTKMesh
            (
                const fvMesh&,
                const labelList& pointLabels
            );

            //- Add face set mesh
            vtkPolyData* faceSetVTKMesh
            (
                const fvMesh&,
                const faceSet&
            );

            //- Add point mesh
            vtkPolyData* pointSetVTKMesh
            (
                const fvMesh&,
                const pointSet&
            );

        // Field conversion functions

            //- Convert volume fields
            void convertVolFields(vtkMultiBlockDataSet*);

            //- Convert point fields
            void convertPointFields(vtkMultiBlockDataSet*);

            //- Convert Lagrangian fields
            void convertLagrangianFields(vtkMultiBlockDataSet*);


        //- Add the fields in the selected time directory to the selection
        //  lists
        template<class GeoField>
        label addObjectsToSelection
        (
            vtkDataArraySelection*,
            const IOobjectList&,
            const string& suffix=string::null
        );


        // Convert Foam fields

            //- Volume fields - all types
            template<class Type>
            void convertVolFields
            (
                const fvMesh&,
                const PtrList<PrimitivePatchInterpolation<primitivePatch> >&,
                const IOobjectList&,
                vtkMultiBlockDataSet* output
            );

            //- Volume field - all selected parts
            template<class Type>
            void convertVolFieldBlock
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                autoPtr<GeometricField<Type, pointPatchField, pointMesh> >&,
                vtkMultiBlockDataSet* output,
                const partInfo& selector,
                const List<polyDecomp>& decompLst
            );

            //- Volume field
            template<class Type>
            void convertVolField
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                vtkMultiBlockDataSet* output,
                const partInfo&,
                const label datasetNo,
                const polyDecomp&
            );

            //- Patch field
            template<class Type>
            void convertPatchField
            (
                const word& name,
                const Field<Type>&,
                vtkMultiBlockDataSet* output,
                const partInfo&,
                const label datasetNo
            );

            //- face set/zone field
            template<class Type>
            void convertFaceField
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                vtkMultiBlockDataSet* output,
                const partInfo&,
                const label datasetNo,
                const fvMesh&,
                const labelList& faceLabels
            );

            //- face set/zone field
            template<class Type>
            void convertFaceField
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                vtkMultiBlockDataSet* output,
                const partInfo&,
                const label datasetNo,
                const fvMesh&,
                const faceSet&
            );

            //- Lagrangian fields - all types
            template<class Type>
            void convertLagrangianFields
            (
                const IOobjectList&,
                vtkMultiBlockDataSet* output,
                const label datasetNo
            );

            //- Lagrangian field
            template<class Type>
            void convertLagrangianField
            (
                const IOField<Type>&,
                vtkMultiBlockDataSet* output,
                const partInfo&,
                const label datasetNo
            );

            //- Point fields - all types
            template<class Type>
            void convertPointFields
            (
                const fvMesh&,
                const pointMesh&,
                const IOobjectList&,
                vtkMultiBlockDataSet* output
            );

            //- Point field - all selected parts
            template<class Type>
            void convertPointFieldBlock
            (
                const GeometricField<Type, pointPatchField, pointMesh>&,
                vtkMultiBlockDataSet* output,
                const partInfo& selector,
                const List<polyDecomp>&
            );

            //- Point fields
            template<class Type>
            void convertPointField
            (
                const GeometricField<Type, pointPatchField, pointMesh>&,
                const GeometricField<Type, fvPatchField, volMesh>&,
                vtkMultiBlockDataSet* output,
                const partInfo&,
                const label datasetNo,
                const polyDecomp&
            );

            //- Patch point field
            template<class Type>
            void convertPatchPointField
            (
                const word& name,
                const Field<Type>&,
                vtkMultiBlockDataSet* output,
                const partInfo&,
                const label datasetNo
            );


       // GUI selection helper functions

            //- Extract up to the first non-word characters
            inline static word getFirstWord(const char*);

            //- Only keep what is listed in hashSet
            static void pruneObjectList
            (
                IOobjectList&,
                const wordHashSet&
            );

            //- Retrieve the current selections
            static wordHashSet getSelected(vtkDataArraySelection*);

            //- Retrieve a sub-list of the current selections
            static wordHashSet getSelected
            (
                vtkDataArraySelection*,
                const partInfo&
            );

            //- Retrieve the current selections
            static stringList getSelectedArrayEntries(vtkDataArraySelection*);

            //- Retrieve a sub-list of the current selections
            static stringList getSelectedArrayEntries
            (
                vtkDataArraySelection*,
                const partInfo&
            );

            //- Set selection(s)
            static void setSelectedArrayEntries
            (
                vtkDataArraySelection*,
                const stringList&
            );

            //- Get the first word from the mesh parts selection
            word getPartName(int);


        //- Disallow default bitwise copy construct
        vtkPVFoam(const vtkPVFoam&);

        //- Disallow default bitwise assignment
        void operator=(const vtkPVFoam&);


public:

    //- Static data members

        ClassName("vtkPVFoam");


    // Constructors

        //- Construct from components
        vtkPVFoam
        (
            const char* const FileName,
            vtkPVFoamReader* reader
        );


    //- Destructor

        ~vtkPVFoam();


    // Member Functions

        //- Update
        void updateInfo();

        void Update
        (
            vtkMultiBlockDataSet* output,
            vtkMultiBlockDataSet* lagrangianOutput
        );

        //- Clean any storage
        void CleanUp();

        //- Allocate and return a list of selected times
        //  returns the count via the parameter
        double* findTimes(int& nTimeSteps);

        //- Add patch names to the display
        void addPatchNames(vtkRenderer* renderer);

        //- Remove patch names from the display
        void removePatchNames(vtkRenderer* renderer);

        //- set the runTime to the first plausible request time,
        //  returns the timeIndex
        //  sets to "constant" on error
        int setTime(int count, const double requestTimes[]);


        //- The current time index
        int timeIndex() const
        {
           return timeIndex_;
        }


     // Access

        //- Debug information
        void PrintSelf(ostream&, vtkIndent) const;

        //- Simple memory used debugging information
        static void printMemory();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#    include "vtkPVFoamI.H"

#endif

// ************************************************************************* //
