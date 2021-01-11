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

#ifndef OPM_FRACTURE_BENCHMARK_3D_HH
#define OPM_FRACTURE_BENCHMARK_3D_HH

#include "config.h"

#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/multidomain/multidomaincoupler.hh>
#include <opm/multidomain/multidomainmodel.hh>


BEGIN_PROPERTIES
NEW_TYPE_TAG(Domain3d, INHERITS_FROM(BenchmarkProblem));
SET_TYPE_PROP(Domain3d, Vanguard, Opm::MultiDomainVanguard<TypeTag>);
SET_TAG_PROP(Domain3d, LocalLinearizerSplice, AutoDiffLocalLinearizer);
// use the element centered finite volume spatial discretization
SET_INT_PROP(Domain3d, GridDim, 3);
SET_INT_PROP(Domain3d, DomainDim, 3);
SET_TAG_PROP(Domain3d, SpatialDiscretizationSplice, EcfvDiscretization);
SET_STRING_PROP(Domain3d, GridFile, FILE_NAME3D);
SET_STRING_PROP(Domain3d, NameSurfix, "3D");
SET_BOOL_PROP(Domain3d, UseVolumetricResidual, 0);

NEW_TYPE_TAG(Domain2d, INHERITS_FROM(Domain3d));
SET_INT_PROP(Domain2d, GridDim, 2);
SET_INT_PROP(Domain2d, DomainDim, 2);
SET_STRING_PROP(Domain2d, GridFile, FILE_NAME2D);
SET_STRING_PROP(Domain2d, NameSurfix, "2D");

NEW_TYPE_TAG(Domain1d, INHERITS_FROM(Domain3d));
SET_INT_PROP(Domain1d, GridDim, 1);
SET_INT_PROP(Domain1d, DomainDim, 1);
SET_STRING_PROP(Domain1d, GridFile, FILE_NAME1D);
SET_STRING_PROP(Domain1d, NameSurfix, "1D");

NEW_TYPE_TAG(Domain0d, INHERITS_FROM(Domain3d));
SET_INT_PROP(Domain0d, GridDim, 1);
SET_INT_PROP(Domain0d, DomainDim, 0);
SET_STRING_PROP(Domain0d, GridFile, FILE_NAME0D);
SET_BOOL_PROP(Domain0d, EnableVtkOutput, false);

// Add the Domain3d Domain2d coupler
NEW_TYPE_TAG(Coupler32, INHERITS_FROM(DarcyCoupler));
SET_TYPE_PROP(Coupler32, Vanguard, Opm::MultiDomainVanguard<TypeTag>);
SET_INT_PROP(Coupler32, GridDim, 2);
SET_INT_PROP(Coupler32, WorldDim, 3);
SET_STRING_PROP(Coupler32, GridFile, FILE_NAME_MORTAR2D);
SET_STRING_PROP(Coupler32, MappingFile, FILE_NAME_MAPPING3D2D);
SET_TYPE_PROP(Coupler32, Scalar, GET_PROP_TYPE(TTAG(Domain3d), Scalar));
SET_TYPE_PROP(Coupler32, CouplingMapper, Opm::FaceElementMapper<TypeTag>);
SET_TYPE_PROP(Coupler32, DomainI, Dune::index_constant<0>);
SET_TYPE_PROP(Coupler32, DomainJ, Dune::index_constant<1>);
SET_PROP(Coupler32, SubTypeTag)
{
    typedef TTAG(Domain3d) Domain0TypeTag;
    typedef TTAG(Domain2d) Domain1TypeTag;

public:
    typedef Opm::MultiDomainProperties<Domain0TypeTag, Domain1TypeTag> type;
};

NEW_TYPE_TAG(Coupler21, INHERITS_FROM(Coupler32));
SET_INT_PROP(Coupler21, GridDim, 1);
SET_STRING_PROP(Coupler21, GridFile, FILE_NAME_MORTAR1D);
SET_STRING_PROP(Coupler21, MappingFile, FILE_NAME_MAPPING2D1D);
SET_TYPE_PROP(Coupler21, DomainI, Dune::index_constant<1>);
SET_TYPE_PROP(Coupler21, DomainJ, Dune::index_constant<2>);
SET_PROP(Coupler21, SubTypeTag)
{
    typedef TTAG(Domain2d) Domain0TypeTag;
    typedef TTAG(Domain1d) Domain1TypeTag;

public:
    typedef Opm::MultiDomainProperties<Domain0TypeTag, Domain1TypeTag> type;
};
NEW_TYPE_TAG(Coupler10, INHERITS_FROM(Coupler32));
SET_INT_PROP(Coupler10, GridDim, 1);
SET_STRING_PROP(Coupler10, GridFile, FILE_NAME_MORTAR0D);
SET_STRING_PROP(Coupler10, MappingFile, FILE_NAME_MAPPING1D0D);
SET_TYPE_PROP(Coupler10, DomainI, Dune::index_constant<2>);
SET_TYPE_PROP(Coupler10, DomainJ, Dune::index_constant<3>);
SET_PROP(Coupler10, SubTypeTag)
{
    typedef TTAG(Domain1d) Domain0TypeTag;
    typedef TTAG(Domain0d) Domain1TypeTag;

public:
    typedef Opm::MultiDomainProperties<Domain0TypeTag, Domain1TypeTag> type;
};

NEW_TYPE_TAG(MultiDimModel, INHERITS_FROM(MultiDomainBaseModel));
NEW_PROP_TAG(SubTypeTag);
NEW_PROP_TAG(MortarView);
NEW_PROP_TAG(MaxTimeStep);
SET_TYPE_PROP(MultiDimModel, MortarView, typename GET_PROP_TYPE(TTAG(Domain2d), Grid)::LeafGridView);
SET_TYPE_PROP(MultiDimModel, Scalar, double);
NEW_PROP_TAG(CouplerTypeTag);

END_PROPERTIES

#endif