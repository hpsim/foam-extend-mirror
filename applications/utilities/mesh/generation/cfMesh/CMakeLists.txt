# --------------------------------------------------------------------------
#   =========                 |
#   \\      /  F ield         | foam-extend: Open Source CFD
#    \\    /   O peration     | Version:     5.0
#     \\  /    A nd           | Web:         http://www.foam-extend.org
#      \\/     M anipulation  | For copyright notice see file Copyright
# --------------------------------------------------------------------------
# License
#     This file is part of foam-extend.
#
#     foam-extend is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     foam-extend is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
#
# Description
#     CMakeLists.txt file for libraries and applications
#
# Author
#     Henrik Rusche, Wikki GmbH, 2017. All rights reserved
#
#
# --------------------------------------------------------------------------

add_subdirectory(FMSToVTK)
add_subdirectory(scaleMesh)
add_subdirectory(cartesian2DMesh)
add_subdirectory(extrudeEdgesInto2DSurface)
add_subdirectory(scaleSurfaceMesh)
add_subdirectory(surfaceFeatureEdges)
add_subdirectory(patchesToSubsets)
add_subdirectory(FLMAToSurface)
add_subdirectory(mergeSurfacePatches)
add_subdirectory(generateBoundaryLayers)
add_subdirectory(surfaceToFMS)
add_subdirectory(improveSymmetryPlanes)
add_subdirectory(pMesh)
add_subdirectory(tetMesh)
add_subdirectory(importSurfaceAsSubset)
add_subdirectory(checkSurfaceMesh)
add_subdirectory(subsetToPatch)
add_subdirectory(improveMeshQuality)
add_subdirectory(preparePar)
add_subdirectory(surfaceGenerateBoundingBox)
add_subdirectory(removeSurfaceFacets)
add_subdirectory(FMSToSurface)
add_subdirectory(cartesianMesh)
add_subdirectory(copySurfaceParts)
add_subdirectory(meshToFPMA)