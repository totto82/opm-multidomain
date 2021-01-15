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
 * \brief Two-phase test for the immiscible model which uses the
 *        vertex-centered finite volume discretization
 */
#include "config.h"

#include "problems/co2StorageProblem.hh"

#include <opm/grid/polyhedralgrid.hh>

#include <opm/models/immiscible/immisciblemodel.hh>

#include <opm/multidomain/dualporositycoupler.hh>
#include <opm/multidomain/multidomainlinearizer.hh>
#include <opm/multidomain/matrixconverter.hh>
#include <opm/multidomain/multidomainmodel.hh>
#include "opm/multidomain/utils/multidomainstart.hh"

constexpr auto fileMatrixDomain = "./data/dualPorosity.txt";


namespace Opm::Properties {

//! The generic type tag for problems using the dual porosity coupler
namespace TTag {
struct MatrixDomain { using InheritsFrom = std::tuple<Co2StorageProblem, ImmiscibleTwoPhaseModel>; };
struct FractureDomain { using InheritsFrom = std::tuple<MatrixDomain>; };
struct CouplerMF { using InheritsFrom = std::tuple<DualPorosityCoupler>; };
struct DualPorModel { using InheritsFrom = std::tuple<MultiDomainBaseModel>; };
} // end namespace TTag

// Set vanguard
template <class TypeTag>
struct Vanguard<TypeTag, TTag::MatrixDomain> { using type = Opm::UnstructuredGridVanguard<TypeTag>; };

// use the element centered finite volume spatial discretization
template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::MatrixDomain> { using type = TTag::EcfvDiscretization; };

// use automatic differentiation for this simulator
template<class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::MatrixDomain> { using type = TTag::AutoDiffLocalLinearizer; };

template<class TypeTag>
struct NameSurfix<TypeTag, TTag::MatrixDomain> { static constexpr auto value = "matrix"; };

// Define grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::MatrixDomain> { using type = Dune::PolyhedralGrid<2, 2>; };

template<class TypeTag>
struct DomainDim<TypeTag, TTag::MatrixDomain> { static constexpr int value = 2; };

template<class TypeTag>
struct GridFile<TypeTag, TTag::MatrixDomain> { static constexpr auto value = fileMatrixDomain; };

/////////////////////////////////////////////////////////////////////////////////////////////////
// Fracture domain Properties
template<class TypeTag>
struct NameSurfix<TypeTag, TTag::FractureDomain> { static constexpr auto value = "fracture"; };
template<class TypeTag>
struct FractureDomain<TypeTag, TTag::FractureDomain> { static constexpr bool value = true; };

/////////////////////////////////////////////////////////////////////////////////////////////////
// Set Coupler properties
// Vanguard
template <class TypeTag>
struct Vanguard<TypeTag, TTag::CouplerMF> { using type = Opm::UnstructuredGridVanguard<TypeTag>; };

// Set leaf view
template<class TypeTag>
struct MortarView<TypeTag, TTag::CouplerMF>
{ using type = typename GetPropType<TTag::MatrixDomain, Properties::Grid>::LeafGridView; };

// Set the grid file name
template<class TypeTag>
struct GridFile<TypeTag, TTag::CouplerMF> { static constexpr auto value = fileMatrixDomain; };

template<class TypeTag>
struct MappingFile<TypeTag, TTag::CouplerMF> { static constexpr auto value = ""; };

template<class TypeTag>
struct Grid<TypeTag, TTag::CouplerMF>
{ using type = GetPropType<TTag::MatrixDomain, Properties::Grid>; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::CouplerMF>
{ using type = GetPropType<TTag::MatrixDomain, Properties::Scalar>; };

template<class TypeTag>
struct CouplingMapper<TypeTag, TTag::CouplerMF>
{ using type = Opm::ElementElementMapper<TypeTag>; };

template<class TypeTag>
struct SubTypeTag<TypeTag, TTag::CouplerMF>
{ using type = Opm::MultiDomainProperties<TTag::MatrixDomain, TTag::FractureDomain>; };

// Define domain copuling ids:
template <class TypeTag>
struct DomainI<TypeTag, TTag::CouplerMF> { using type = Dune::index_constant<0>; };
template <class TypeTag>
struct DomainJ<TypeTag, TTag::CouplerMF> { using type = Dune::index_constant<1>; };
//////////////////////////////////////////////////////////////////////////////////////
// Set the properties of the multi-dimensional problem
// Set the linearizer
template<class TypeTag>
struct Linearizer<TypeTag, TTag::DualPorModel>
{
    using type = Opm::MultiDomainLinearizer< TypeTag >;
};
// Set end time of simulation
template<class TypeTag>
struct EndTime<TypeTag, TTag::DualPorModel> {
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 2000 * 60 * 60 * 24;
};

// Set initial time step of simulation
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::DualPorModel> {
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 60 * 60;
};

// Set max time step of simulation
template<class TypeTag>
struct MaxTimeStepSize<TypeTag, TTag::DualPorModel> {
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 20 * 60 * 60 * 24;
};

// Set active domain type tags
template<class TypeTag>
struct SubTypeTag<TypeTag, TTag::DualPorModel>
{
    using type = Opm::MultiDomainProperties<TTag::MatrixDomain, TTag::FractureDomain>;
};
// Set active copuler type tags
template<class TypeTag>
struct CouplerTypeTag<TypeTag, TTag::DualPorModel>
{
    using type = Opm::MultiCouplerProperties<TTag::CouplerMF>;
};
} //end namespace Opm::Properties


int main(int argc, char **argv)
{
    using DualPorosityProblem = Opm::Properties::TTag::DualPorModel;
    Opm::multidomainStart<DualPorosityProblem>(argc, argv);
}