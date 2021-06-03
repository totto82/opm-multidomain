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
 * \brief Contains type tag specifications for a multi-dimensional problem
 */

#ifndef OPM_MULTI_DIM_PROPERTIES_HH
#define OPM_MULTI_DIM_PROPERTIES_HH

#include "config.h"

#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/io/unstructuredgridvanguard.hh>
#include <opm/multidomain/multidomaincoupler.hh>
#include <opm/multidomain/multidomainmodel.hh>
#include <opm/multidomain/multidomainlinearizer.hh>
#include <opm/grid/polyhedralgrid.hh>


namespace Opm::Properties {
// Create new type tags
namespace TTag {

struct MultiDimModel {using InheritsFrom = std::tuple<MultiDomainBaseModel>; };

struct Domain3D { using InheritsFrom = std::tuple<ImplicitModel>; };

struct Domain2D { using InheritsFrom = std::tuple<Domain3D>; };

struct Domain1D { using InheritsFrom = std::tuple<Domain2D>; };

struct Domain0D { using InheritsFrom = std::tuple<Domain1D>; };

struct Coupler32 {using InheritsFrom = std::tuple<DarcyCoupler>; };

struct Coupler21 {using InheritsFrom = std::tuple<Coupler32>; };

struct Coupler10 {using InheritsFrom = std::tuple<Coupler21>; };

}  // end namespace TTag

// Set the grid dimension
template<class TypeTag, class MyTypeTag>
struct NameSurfix { using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DomainDim { using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct WorldDim { using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct GridDim { using type = UndefinedProperty; };


template<class TypeTag>
struct WorldDim<TypeTag, TTag::Domain3D> { static constexpr int value = 3; };

// Set the properties of the 3d domain
template <class TypeTag>
struct Vanguard<TypeTag, TTag::Domain3D> { using type = Opm::UnstructuredGridVanguard<TypeTag>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Domain3D> {
private:
    enum{ gridDim = getPropValue<TypeTag, Properties::GridDim>() };
    enum{ worldDim = getPropValue<TypeTag, Properties::WorldDim>() };
public:
    using type = Dune::PolyhedralGrid<gridDim, worldDim>;
};

// use the element centered finite volume spatial discretization
template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::Domain3D> { using type = TTag::EcfvDiscretization; };

// use automatic differentiation for this simulator
template<class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::Domain3D> { using type = TTag::AutoDiffLocalLinearizer; };

// Set the grid dimension
template<class TypeTag>
struct GridDim<TypeTag, TTag::Domain3D> { static constexpr int value = 3; };

// Set the domain dimension
template<class TypeTag>
struct DomainDim<TypeTag, TTag::Domain3D> { static constexpr int value = 3; };

// Set the output file surfix (used for vtk)
template<class TypeTag>
struct NameSurfix<TypeTag, TTag::Domain3D> { static constexpr auto value = "3D"; };

// Do not use volumetric residual
template<class TypeTag>
struct UseVolumetricResidual<TypeTag, TTag::Domain3D> { static constexpr bool value = false; };

//////////////////////////////////////////////////////////////////////////////////////
// Set the properties of the 2d domain
// Set the grid dimension
template<class TypeTag>
struct GridDim<TypeTag, TTag::Domain2D> { static constexpr int value = 2; };

// Set the domain dimension
template<class TypeTag>
struct DomainDim<TypeTag, TTag::Domain2D> { static constexpr int value = 2; };

// Set the output file surfix
template<class TypeTag>
struct NameSurfix<TypeTag, TTag::Domain2D> { static constexpr auto value = "2D"; };

//////////////////////////////////////////////////////////////////////////////////////
// Set the properties of the 1d domain
// Set the grid dimension
template<class TypeTag>
struct GridDim<TypeTag, TTag::Domain1D> { static constexpr int value = 1; };

// Set the domain dimension
template<class TypeTag>
struct DomainDim<TypeTag, TTag::Domain1D> { static constexpr int value = 1; };

// Set the output file surfix
template<class TypeTag>
struct NameSurfix<TypeTag, TTag::Domain1D> { static constexpr auto value = "1D"; };


//////////////////////////////////////////////////////////////////////////////////////
// Set the properties of the 0d domain
// Set the grid dimension
template<class TypeTag>
struct GridDim<TypeTag, TTag::Domain0D> { static constexpr int value = 1; };

// Set the domain dimension
template<class TypeTag>
struct DomainDim<TypeTag, TTag::Domain0D> { static constexpr int value = 0; };

// Do not write output for 0d domain
template<class TypeTag>
struct EnableVtkOutput<TypeTag, TTag::Domain0D> { static constexpr bool value = false; };

//////////////////////////////////////////////////////////////////////////////////////
// Set the properties of the 3d-2d Coupler
template <class TypeTag>
struct Vanguard<TypeTag, TTag::Coupler32> { using type = Opm::UnstructuredGridVanguard<TypeTag>; };

// Set the grid dimension
template<class TypeTag>
struct GridDim<TypeTag, TTag::Coupler32> { static constexpr int value = 2; };

// Set leaf view
template<class TypeTag>
struct MortarView<TypeTag, TTag::Coupler32>
{ using type = typename GetPropType<TTag::Domain2D, Properties::Grid>::LeafGridView; };

template<class TypeTag>
struct Grid<TypeTag, TTag::Coupler32>
{ using type = GetPropType<TTag::Domain2D, Properties::Grid>; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::Coupler32>
{ using type = GetPropType<TTag::Domain2D, Properties::Scalar>; };

template<class TypeTag>
struct CouplingMapper<TypeTag, TTag::Coupler32>
{ using type = Opm::FaceElementMapper<TypeTag>; };


//////////////////////////////////////////////////////////////////////////////////////
// Set the properties of the 2d-1d Coupler
// Set the grid dimension
template<class TypeTag>
struct GridDim<TypeTag, TTag::Coupler21> { static constexpr int value = 1; };

template<class TypeTag>
struct MortarView<TypeTag, TTag::Coupler21>
{ using type = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView; };

template<class TypeTag>
struct Grid<TypeTag, TTag::Coupler21>
{ 
    using SubTypes = GetPropType<TypeTag, Properties::SubTypeTag>;
    using LowerDimType = typename SubTypes::template TypeTag<1>;
    using type = GetPropType<LowerDimType, Properties::Grid>; 
};


//////////////////////////////////////////////////////////////////////////////////////
// Set the properties of the 1d-0d Coupler
// Set the grid dimension
template<class TypeTag>
struct GridDim<TypeTag, TTag::Coupler10> { static constexpr int value = 0; };

template<class TypeTag>
struct MortarView<TypeTag, TTag::Coupler10>
{ using type = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView; };

template<class TypeTag>
struct Grid<TypeTag, TTag::Coupler10>
{
private:
    enum{ gridDim = getPropValue<TypeTag, Properties::GridDim>() };
    enum{ worldDim = getPropValue<TypeTag, Properties::WorldDim>() };
public:
    using type = Dune::PolyhedralGrid<gridDim, worldDim>;
};

//////////////////////////////////////////////////////////////////////////////////////
// Set the properties of the multi-dimensional problem
// Set the linearizer
template<class TypeTag>
struct Linearizer<TypeTag, TTag::MultiDimModel>
{
    using type = Opm::MultiDomainLinearizer< TypeTag >;
};


} // end namespace Opm::Properties

#endif // OPM_MULTI_DIM_PROPERTIES_HH