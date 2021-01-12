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

#include "config.h"

#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/io/unstructuredgridvanguard.hh>
#include <opm/multidomain/multidomaincoupler.hh>
#include <opm/multidomain/multidomainmodel.hh>
#include <opm/multidomain/multidomainproperties.hh>
#include <opm/multidomain/multidomainlinearizer.hh>

BEGIN_PROPERTIES
NEW_TYPE_TAG(Domain2D, INHERITS_FROM(BenchmarkProblem));
SET_TYPE_PROP(Domain2D, Vanguard, Opm::UnstructuredGridVanguard<TypeTag>);
SET_TAG_PROP(Domain2D, LocalLinearizerSplice, AutoDiffLocalLinearizer);
// use the element centered finite volume spatial discretization
SET_INT_PROP(Domain2D, GridDim, 2);
SET_INT_PROP(Domain2D, DomainDim, 2);
SET_TAG_PROP(Domain2D, SpatialDiscretizationSplice, EcfvDiscretization);
SET_STRING_PROP(Domain2D, GridFile, FILE_NAME2D);
SET_STRING_PROP(Domain2D, NameSurfix, "2D");
SET_BOOL_PROP(Domain2D, UseVolumetricResidual, 0);

NEW_TYPE_TAG(Domain1D, INHERITS_FROM(Domain2D));
SET_INT_PROP(Domain1D, GridDim, 1);
SET_INT_PROP(Domain1D, DomainDim, 1);
SET_STRING_PROP(Domain1D, GridFile, FILE_NAME1D);
SET_STRING_PROP(Domain1D, NameSurfix, "1D");

NEW_TYPE_TAG(Domain0D, INHERITS_FROM(Domain1D));
SET_INT_PROP(Domain0D, GridDim, 1);
SET_INT_PROP(Domain0D, DomainDim, 0);
SET_STRING_PROP(Domain0D, GridFile, FILE_NAME0D);
SET_BOOL_PROP(Domain0D, EnableVtkOutput, false);


// Add the Domain2D Domain1D coupler
NEW_TYPE_TAG(Coupler21, INHERITS_FROM(DarcyCoupler));
SET_TYPE_PROP(Coupler21, Vanguard, Opm::UnstructuredGridVanguard<TypeTag>);
SET_INT_PROP(Coupler21, GridDim, 1);
SET_INT_PROP(Coupler21, WorldDim, 2);
SET_STRING_PROP(Coupler21, GridFile, FILE_NAME_MORTAR1D);
SET_STRING_PROP(Coupler21, MappingFile, FILE_NAME_MAPPING2D1D);
SET_TYPE_PROP(Coupler21, Grid, GET_PROP_TYPE(TTAG(Domain1D), Grid));
SET_TYPE_PROP(Coupler21, Scalar, GET_PROP_TYPE(TTAG(Domain2D), Scalar));
SET_TYPE_PROP(Coupler21, CouplingMapper, Opm::FaceElementMapper<TypeTag>);
SET_TYPE_PROP(Coupler21, DomainI, Dune::index_constant<0>);
SET_TYPE_PROP(Coupler21, DomainJ, Dune::index_constant<1>);
SET_PROP(Coupler21, SubTypeTag)
{
    typedef TTAG(Domain2D) Domain2DTypeTag;
    typedef TTAG(Domain1D) Domain1DTypeTag;

public:
    typedef Opm::MultiDomainProperties<Domain2DTypeTag, Domain1DTypeTag> type;
};


NEW_TYPE_TAG(Coupler10, INHERITS_FROM(Coupler21));
SET_STRING_PROP(Coupler10, GridFile, FILE_NAME_MORTAR0D);
SET_STRING_PROP(Coupler10, MappingFile, FILE_NAME_MAPPING1D0D);
SET_INT_PROP(Coupler10, GridDim, 0);
SET_TYPE_PROP(Coupler10, DomainI, Dune::index_constant<1>);
SET_TYPE_PROP(Coupler10, DomainJ, Dune::index_constant<2>);
SET_PROP(Coupler10, SubTypeTag)
{
    typedef TTAG(Domain1D) Domain1DTypeTag;
    typedef TTAG(Domain0D) Domain0DTypeTag;

public:
    typedef Opm::MultiDomainProperties<Domain1DTypeTag, Domain0DTypeTag> type;
};

NEW_TYPE_TAG(MultiDimModel, INHERITS_FROM(MultiDomainBaseModel));
SET_TYPE_PROP(MultiDomain, Linearizer, Opm::MultiDomainLinearizer< TypeTag >);
SET_TYPE_PROP(MultiDimModel, MortarView, typename GET_PROP_TYPE(TTAG(Domain1D), Grid)::LeafGridView);
SET_SCALAR_PROP(MultiDimModel, EndTime, 1);
SET_SCALAR_PROP(MultiDimModel, InitialTimeStepSize, 1);
SET_PROP(MultiDimModel, SubTypeTag)
{
    typedef TTAG(Domain2D) Domain2DType;
    typedef TTAG(Domain1D) Domain1DType;
    typedef TTAG(Domain0D) Domain0DType;

public:
    typedef Opm::MultiDomainProperties<Domain2DType, Domain1DType, Domain0DType> type;
};
SET_PROP(MultiDimModel, CouplerTypeTag)
{
    typedef TTAG(Coupler21) Coupler0;
    typedef TTAG(Coupler10) Coupler1;

public:
    typedef Opm::MultiCouplerProperties<Coupler0, Coupler1> type;
};

END_PROPERTIES

#endif