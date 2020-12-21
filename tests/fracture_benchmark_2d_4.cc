// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief Test for the immisicible VCVF discretization with only a single phase
 */
#include "config.h"
#include "problems/benchmark2dproblem4.hh"

#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/multidomain/utils/start.hh>


const auto FILE_NAME2D = "./data/benchmark4_2.txt";
const auto FILE_NAME1D = "./data/benchmark4_1.txt";
const auto FILE_NAME0D = "./data/benchmark4_0.txt";

const auto FILE_NAME_MORTAR1D = "./data/benchmark4_mortar_1.txt";
const auto FILE_NAME_MORTAR0D = "./data/benchmark4_mortar_0.txt";
const auto FILE_NAME_MAPPING2D1D = "./data/benchmark4_mapping_1.txt";
const auto FILE_NAME_MAPPING1D0D = "./data/benchmark4_mapping_0.txt";

BEGIN_PROPERTIES
NEW_TYPE_TAG(BenchmarkProblem, INHERITS_FROM(ImmiscibleSinglePhaseModel, Benchmark4Problem));
END_PROPERTIES

#include "benchmark2d.hh"

int main(int argc, char **argv)
{
    typedef TTAG(MultiDimModel) MixedDimModelTypeTag;
    start<MixedDimModelTypeTag>(argc, argv);
}