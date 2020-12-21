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
#include "problems/benchmark2dproblem3.hh"

#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/multidomain/multidomaincoupler.hh>
#include <opm/multidomain/utils/start.hh>


const auto fileDomain2d = "./data/benchmark3_2.txt";
const auto fileDomain1d = "./data/benchmark3_1.txt";
const auto fileDomain0d = "./data/benchmark3_0.txt";

const auto fileMortar1d = "./data/benchmark3_mortar_1.txt";
const auto fileMortar0d = "./data/benchmark3_mortar_0.txt";
const auto fileMapping2d1d = "./data/benchmark3_mapping_1.txt";
const auto fileMapping1d0d = "./data/benchmark3_mapping_0.txt";

BEGIN_PROPERTIES
NEW_TYPE_TAG(Domain2d, INHERITS_FROM(ImmiscibleSinglePhaseModel, Benchmark3Problem));

SET_BOOL_PROP(Domain2d, LeftRight, 1);
SET_TYPE_PROP(Domain2d, Vanguard, Opm::MultiDomainVanguard<TypeTag>);
SET_TAG_PROP(Domain2d, LocalLinearizerSplice, AutoDiffLocalLinearizer);
// use the element centered finite volume spatial discretization
SET_INT_PROP(Domain2d, GridDim, 2);
SET_TAG_PROP(Domain2d, SpatialDiscretizationSplice, EcfvDiscretization);
SET_STRING_PROP(Domain2d, GridFile, fileDomain2d);
SET_STRING_PROP(Domain2d, NameSurfix, "2d");
SET_BOOL_PROP(Domain2d, UseVolumetricResidual, 0);

NEW_TYPE_TAG(Domain1d, INHERITS_FROM(Domain2d));
SET_INT_PROP(Domain1d, GridDim, 1);
SET_STRING_PROP(Domain1d, GridFile, fileDomain1d);
SET_STRING_PROP(Domain1d, NameSurfix, "1d");

NEW_TYPE_TAG(Domain0d, INHERITS_FROM(Domain1d));
SET_INT_PROP(Domain0d, GridDim, 1);
SET_STRING_PROP(Domain0d, GridFile, fileDomain0d);
SET_BOOL_PROP(Domain0d, EnableVtkOutput, false);


// Add the Domain2d Domain1d coupler
NEW_TYPE_TAG(Coupler21, INHERITS_FROM(DarcyCoupler));
SET_TYPE_PROP(Coupler21, Vanguard, Opm::MultiDomainVanguard<TypeTag>);
SET_INT_PROP(Coupler21, GridDim, 1);
SET_INT_PROP(Coupler21, WorldDim, 2);
SET_STRING_PROP(Coupler21, GridFile, fileMortar1d);
SET_STRING_PROP(Coupler21, MappingFile, fileMapping2d1d);
SET_TYPE_PROP(Coupler21, Scalar, GET_PROP_TYPE(TTAG(Domain2d), Scalar));
SET_TYPE_PROP(Coupler21, CouplingMapper, Opm::FaceElementMapper<TypeTag>);
SET_TYPE_PROP(Coupler21, DomainI, Dune::index_constant<0>);
SET_TYPE_PROP(Coupler21, DomainJ, Dune::index_constant<1>);
SET_PROP(Coupler21, SubTypeTag)
{
    typedef TTAG(Domain2d) Domain2dTypeTag;
    typedef TTAG(Domain1d) Domain1dTypeTag;

public:
    typedef Opm::MultiDomainProperties<Domain2dTypeTag, Domain1dTypeTag> type;
};


NEW_TYPE_TAG(Coupler10, INHERITS_FROM(Coupler21));
SET_STRING_PROP(Coupler10, GridFile, fileMortar0d);
SET_STRING_PROP(Coupler10, MappingFile, fileMapping1d0d);
SET_TYPE_PROP(Coupler10, DomainI, Dune::index_constant<1>);
SET_TYPE_PROP(Coupler10, DomainJ, Dune::index_constant<2>);
SET_PROP(Coupler10, SubTypeTag)
{
    typedef TTAG(Domain1d) Domain1dTypeTag;
    typedef TTAG(Domain0d) Domain0dTypeTag;

public:
    typedef Opm::MultiDomainProperties<Domain1dTypeTag, Domain0dTypeTag> type;
};

NEW_TYPE_TAG(MultiDimModel, INHERITS_FROM(MultiDomainBaseModel));
NEW_PROP_TAG(SubTypeTag);
NEW_PROP_TAG(MortarView);
SET_TYPE_PROP(MultiDimModel, MortarView, typename GET_PROP_TYPE(TTAG(Domain1d), Grid)::LeafGridView);
SET_SCALAR_PROP(MultiDimModel, EndTime, 1);
SET_SCALAR_PROP(MultiDimModel, InitialTimeStepSize, 1);
NEW_PROP_TAG(CouplerTypeTag);
SET_PROP(MultiDimModel, SubTypeTag)
{
    typedef TTAG(Domain2d) Domain2dType;
    typedef TTAG(Domain1d) Domain1dType;
    typedef TTAG(Domain0d) Domain0dType;

public:
    typedef Opm::MultiDomainProperties<Domain2dType, Domain1dType, Domain0dType> type;
};
SET_PROP(MultiDimModel, CouplerTypeTag)
{
    typedef TTAG(Coupler21) Coupler0;
    typedef TTAG(Coupler10) Coupler1;

public:
    typedef Opm::MultiCouplerProperties<Coupler0, Coupler1> type;
};

END_PROPERTIES

int main(int argc, char **argv)
{
    typedef TTAG(MultiDimModel) MixedDimModelTypeTag;
    start<MixedDimModelTypeTag>(argc, argv);
}