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

add_subdirectory(engine)
add_subdirectory(solidMechanics)
add_subdirectory(multiphase)
add_subdirectory(lagrangian)
add_subdirectory(financial)
add_subdirectory(discreteMethods)
add_subdirectory(coupled)
add_subdirectory(compressible)
add_subdirectory(multiSolver)
add_subdirectory(incompressible)
add_subdirectory(equationReaderDemo)
add_subdirectory(heatTransfer)
add_subdirectory(electromagnetics)
add_subdirectory(combustion)
add_subdirectory(finiteArea)
add_subdirectory(DNS)
add_subdirectory(basic)
add_subdirectory(immersedBoundary)
add_subdirectory(surfaceTracking)
add_subdirectory(viscoelastic)
add_subdirectory(overset)

