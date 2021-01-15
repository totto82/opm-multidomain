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

#ifndef OPM_FRACTURE_BENCHMARK_2D_HH
#define OPM_FRACTURE_BENCHMARK_2D_HH

#include "benchmarkproperties.hh"

#include <dune/common/indices.hh>


namespace Opm::Properties {
// // Create new type tags
// namespace TTag {
// struct Domain2D { using InheritsFrom = std::tuple<BenchmarkProblem>; };

// struct Domain1D { using InheritsFrom = std::tuple<Domain2D>; };

// struct Domain0D { using InheritsFrom = std::tuple<Domain1D>; };

// struct Coupler21 {using InheritsFrom = std::tuple<DarcyCoupler>; };

// struct Coupler10 {using InheritsFrom = std::tuple<Coupler21>; };

// struct MultiDimModel {using InheritsFrom = std::tuple<MultiDomainBaseModel>; };
// }  // end namespace TTag

// // Set the properties of the 2d domain
// template <class TypeTag>
// struct Vanguard<TypeTag, TTag::Domain2D> { using type = Opm::UnstructuredGridVanguard<TypeTag>; };

// // use the element centered finite volume spatial discretization
// template<class TypeTag>
// struct SpatialDiscretizationSplice<TypeTag, TTag::Domain2D> { using type = TTag::EcfvDiscretization; };

// // use automatic differentiation for this simulator
// template<class TypeTag>
// struct LocalLinearizerSplice<TypeTag, TTag::Domain2D> { using type = TTag::AutoDiffLocalLinearizer; };

// // Set the grid dimension
// template<class TypeTag>
// struct GridDim<TypeTag, TTag::Domain2D> { static constexpr int value = 2; };

// // Set the domain dimension
// template<class TypeTag>
// struct DomainDim<TypeTag, TTag::Domain2D> { static constexpr int value = 2; };

// // Set the grid file name
// template<class TypeTag>
// struct GridFile<TypeTag, TTag::Domain2D> { static constexpr auto value = FILE_NAME2D; };

// // Set the output file surfix
// template<class TypeTag>
// struct NameSurfix<TypeTag, TTag::Domain2D> { static constexpr auto value = "2D"; };

// // Do not use volumetric residual
// template<class TypeTag>
// struct UseVolumetricResidual<TypeTag, TTag::Domain2D> { static constexpr bool value = false; };

// //////////////////////////////////////////////////////////////////////////////////////
// // Set the properties of the 1d domain
// // Set the grid dimension
// template<class TypeTag>
// struct GridDim<TypeTag, TTag::Domain1D> { static constexpr int value = 1; };

// // Set the domain dimension
// template<class TypeTag>
// struct DomainDim<TypeTag, TTag::Domain1D> { static constexpr int value = 1; };

// // Set the grid file name
// template<class TypeTag>
// struct GridFile<TypeTag, TTag::Domain1D> { static constexpr auto value = FILE_NAME1D; };

// // Set the output file surfix
// template<class TypeTag>
// struct NameSurfix<TypeTag, TTag::Domain1D> { static constexpr auto value = "1D"; };


// //////////////////////////////////////////////////////////////////////////////////////
// // Set the properties of the 0d domain
// // Set the grid dimension
// template<class TypeTag>
// struct GridDim<TypeTag, TTag::Domain0D> { static constexpr int value = 1; };

// // Set the domain dimension
// template<class TypeTag>
// struct DomainDim<TypeTag, TTag::Domain0D> { static constexpr int value = 0; };

// // Set the grid file name
// template<class TypeTag>
// struct GridFile<TypeTag, TTag::Domain0D> { static constexpr auto value = FILE_NAME0D; };

// // Do not write output for 0d domain
// template<class TypeTag>
// struct EnableVtkOutput<TypeTag, TTag::Domain0D> { static constexpr bool value = false; };

// //////////////////////////////////////////////////////////////////////////////////////
// // Set the properties of the 2d-1d Coupler
// template <class TypeTag>
// struct Vanguard<TypeTag, TTag::Coupler21> { using type = Opm::UnstructuredGridVanguard<TypeTag>; };

// // Set the grid dimension
// template<class TypeTag>
// struct GridDim<TypeTag, TTag::Coupler21> { static constexpr int value = 1; };

// // Set leaf view
// template<class TypeTag>
// struct MortarView<TypeTag, TTag::Coupler21>
// { using type = typename GetPropType<TTag::Domain1D, Properties::Grid>::LeafGridView; };

// // Set the grid file name
// template<class TypeTag>
// struct GridFile<TypeTag, TTag::Coupler21> { static constexpr auto value = FILE_NAME_MORTAR1D; };

// // Set the file containing the mapping between domains
// template<class TypeTag>
// struct MappingFile<TypeTag, TTag::Coupler21> { static constexpr auto value = FILE_NAME_MAPPING2D1D; };

// template<class TypeTag>
// struct Grid<TypeTag, TTag::Coupler21>
// { using type = GetPropType<TTag::Domain1D, Properties::Grid>; };

// template<class TypeTag>
// struct Scalar<TypeTag, TTag::Coupler21>
// { using type = GetPropType<TTag::Domain1D, Properties::Scalar>; };

// template<class TypeTag>
// struct CouplingMapper<TypeTag, TTag::Coupler21>
// { using type = Opm::FaceElementMapper<TypeTag>; };

// template<class TypeTag>
// struct DomainI<TypeTag, TTag::Coupler21>
// { using type = Dune::index_constant<0>; };

// template<class TypeTag>
// struct DomainJ<TypeTag, TTag::Coupler21>
// { using type = Dune::index_constant<1>; };

// template<class TypeTag>
// struct SubTypeTag<TypeTag, TTag::Coupler21>
// {
//     using type = Opm::MultiDomainProperties<TTag::Domain2D, TTag::Domain1D>;
// };

// //////////////////////////////////////////////////////////////////////////////////////
// // Set the properties of the 1d-0d Coupler
// // Set the grid dimension
// template<class TypeTag>
// struct GridDim<TypeTag, TTag::Coupler10> { static constexpr int value = 0; };

// // Set the grid file name
// template<class TypeTag>
// struct GridFile<TypeTag, TTag::Coupler10> { static constexpr auto value = FILE_NAME_MORTAR0D; };

// // Set the file containing the mapping between domains
// template<class TypeTag>
// struct MappingFile<TypeTag, TTag::Coupler10> { static constexpr auto value = FILE_NAME_MAPPING1D0D; };

// template<class TypeTag>
// struct DomainI<TypeTag, TTag::Coupler10>
// { using type = Dune::index_constant<1>; };

// template<class TypeTag>
// struct DomainJ<TypeTag, TTag::Coupler10>
// { using type = Dune::index_constant<2>; };

// template<class TypeTag>
// struct SubTypeTag<TypeTag, TTag::Coupler10>
// {
//     using type = Opm::MultiDomainProperties<TTag::Domain1D, TTag::Domain0D>;
// };


// //////////////////////////////////////////////////////////////////////////////////////
// // Set the properties of the multi-dimensional problem
// // Set the linearizer
// template<class TypeTag>
// struct Linearizer<TypeTag, TTag::MultiDimModel>
// {
//     using type = Opm::MultiDomainLinearizer< TypeTag >;
// };

// Set end time of simulation
template<class TypeTag>
struct EndTime<TypeTag, TTag::MultiDimModel> {
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 1;
};

// Set initial time step of simulation
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::MultiDimModel> {
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 1;
};

// Set active domain type tags
template<class TypeTag>
struct SubTypeTag<TypeTag, TTag::MultiDimModel>
{
    using type = Opm::MultiDomainProperties<TTag::Domain2D, TTag::Domain1D, TTag::Domain0D>;
};
// Set active copuler type tags
template<class TypeTag>
struct CouplerTypeTag<TypeTag, TTag::MultiDimModel>
{
    using type = Opm::MultiCouplerProperties<TTag::Coupler21, TTag::Coupler10>;
};

// Define domain copuling ids:
// 2d -> 1d
template <class TypeTag>
struct DomainI<TypeTag, TTag::Coupler21> { using type = Dune::index_constant<0>; };
template <class TypeTag>
struct DomainJ<TypeTag, TTag::Coupler21> { using type = Dune::index_constant<1>; };

// 1d -> 0d
template <class TypeTag>
struct DomainI<TypeTag, TTag::Coupler10> { using type = Dune::index_constant<1>; };
template <class TypeTag>
struct DomainJ<TypeTag, TTag::Coupler10> { using type = Dune::index_constant<2>; };

} // end namespace Opm::Properties

#endif // OPM_FRACTURE_BENCHMARK_2D_HH