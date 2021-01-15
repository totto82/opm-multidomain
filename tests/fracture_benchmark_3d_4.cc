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
#include "problems/benchmark3dproblem4.hh"

#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/multidomain/utils/multidomainstart.hh>


constexpr auto FILE_NAME3D = "./data/benchmark_3d_case_4_3.txt";
constexpr auto FILE_NAME2D = "./data/benchmark_3d_case_4_2.txt";
constexpr auto FILE_NAME1D = "./data/benchmark_3d_case_4_1.txt";
constexpr auto FILE_NAME0D = "";

constexpr auto FILE_NAME_MORTAR2D = "./data/benchmark_3d_case_4_mortar_2.txt";
constexpr auto FILE_NAME_MORTAR1D = "./data/benchmark_3d_case_4_mortar_1.txt";
constexpr auto FILE_NAME_MORTAR0D = "";

constexpr auto FILE_NAME_MAPPING3D2D = "./data/benchmark_3d_case_4_mapping_2.txt";
constexpr auto FILE_NAME_MAPPING2D1D = "./data/benchmark_3d_case_4_mapping_1.txt";
constexpr auto FILE_NAME_MAPPING1D0D = "";

// Create new type tags
namespace Opm::Properties {

namespace TTag {
struct BenchmarkProblem {
    using InheritsFrom = std::tuple<Benchmark3d4Problem>;
};
}  // end namespace TTag
}  // end namespace Opm::Properties

#include "benchmark3d.hh"

// Define problem specific properties
namespace Opm::Properties {
// Set end time of problem
template <class TypeTag>
struct EndTime<TypeTag, TTag::MultiDimModel> {
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 2000;
};

// Set time step size of problem
template <class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::MultiDimModel> {
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 50;
};

// Set max time step size of problem
template <class TypeTag>
struct MaxTimeStepSize<TypeTag, TTag::MultiDimModel> {
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 50;
};

// Set active domains of problem
template <class TypeTag>
struct SubTypeTag<TypeTag, TTag::MultiDimModel> {
    using type = Opm::MultiDomainProperties<TTag::Domain3D, TTag::Domain2D, TTag::Domain1D>;
};

// Set active couplers of problem
template <class TypeTag>
struct CouplerTypeTag<TypeTag, TTag::MultiDimModel> {
    using type = Opm::MultiCouplerProperties<TTag::Coupler32, TTag::Coupler21>;
};
}  // end namespace Opm::Properties

int main(int argc, char **argv)
{
    using MixedDimModelTypeTag = Opm::Properties::TTag::MultiDimModel;
    Opm::multidomainStart<MixedDimModelTypeTag>(argc, argv);
}
