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
 * \copydoc Opm::DarcyCoupler
 */

#ifndef OPM_MULTI_DOMAIN_COUPLER_HH
#define OPM_MULTI_DOMAIN_COUPLER_HH

#include <dune/common/indices.hh>
#include <opm/material/densead/Math.hpp>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/multidomain/couplingelementcontext.hh>
#include <opm/multidomain/ecfvcouplingstencil.hh>
#include <opm/multidomain/multidomainmapper.hh>
#include <opm/multidomain/multidomainproperties.hh>
#include <stdexcept>

namespace Opm {
template <class TypeTag>
class DarcyCoupler;
}  // namespace Opm

namespace Opm::Properties {
//! The generic type tag for problems using the dual porosity coupler
namespace TTag {
struct DarcyCoupler {
    using InheritsFrom = std::tuple<ImmiscibleSinglePhaseModel, MultiDomain>;
};
}  // end namespace TTag

template <class TypeTag, class MyTypeTag>
struct MappingFile {
    using type = UndefinedProperty;
};

//! set the Coupler type
template <class TypeTag>
struct Coupler<TypeTag, TTag::DarcyCoupler> {
    using type = Opm::DarcyCoupler<TypeTag>;
};

//! set the Coupler mapper
template <class TypeTag>
struct CouplingMapper<TypeTag, TTag::DarcyCoupler> {
    using type = Opm::FaceElementMapper<TypeTag>;
};

//! set the Coupler Element context
template <class TypeTag>
struct CouplingElementContext<TypeTag, TTag::DarcyCoupler> {
    using type = Opm::CouplingElementContext<TypeTag>;
};

//! The simulator is the same as the coupler
template <class TypeTag>
struct Simulator<TypeTag, TTag::DarcyCoupler> {
    using type = Opm::DarcyCoupler<TypeTag>;
};

//! The discretization history
template <class TypeTag>
struct TimeDiscHistorySize<TypeTag, TTag::DarcyCoupler> {
    static constexpr int value = 0;
};

//! The mapping file must be set by the problem
template <class TypeTag>
struct MappingFile<TypeTag, TTag::DarcyCoupler> {
    using type = UndefinedProperty;
};

//! Set the copuling stencil
template <class TypeTag>
struct Stencil<TypeTag, TTag::DarcyCoupler> {
   private:
    using Scalar = GetPropType<TypeTag, Scalar>;
    using GridView = GetPropType<TypeTag, GridView>;
    using Mapper = GetPropType<TypeTag, CouplingMapper>;
    using SubTypeTag = GetPropType<TypeTag, SubTypeTag>;

   public:
    using type = Opm::EcfvMixedDimStencil<Scalar, GridView, Mapper, SubTypeTag>;
};
}  // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup MultiDomainModel
 *
 * \brief Couple two domains in a multidomain model
 *
 * A multidomain model consist of several subdomains which are initiated as
 * standard OPM models. This class defines the coupling between two subdomains
 * by a Darcy law. The class allows for two different types of coupling: the
 * first is between two domains of equal dimension through a shared interface.
 * The second is between a higher dimensional domain and a lower-dimensional
 * domain (fracture). The mortar grid is associated with the coupler.
 *
 * The Darcy coupler does to some extent act as a Simulator class for the
 * standard models. However, it is also responsible for discretization, linearization
 * and problem definition. This will in the feature be split into different classes.
 *
 * You should not expect the DarcyCoupler to work for anything else than immicible
 * models.
 *
 */
template <class TypeTag>
class DarcyCoupler {
    using Implementation = GetPropType<TypeTag, Properties::Coupler>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using DomainI = GetPropType<TypeTag, Properties::DomainI>;
    using DomainJ = GetPropType<TypeTag, Properties::DomainJ>;
    using MortarGrid = GetPropType<TypeTag, Properties::Grid>;
    using MortarView = GetPropType<TypeTag, Properties::GridView>;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
    using SubTypeTag = GetPropType<TypeTag, Properties::SubTypeTag>;
    using CouplingMapper = GetPropType<TypeTag, Properties::CouplingMapper>;
    using CouplingElementContext =
        GetPropType<TypeTag, Properties::CouplingElementContext>;

    template <std::size_t i>
    using Simulator = typename SubTypeTag::template Simulator<i>;
    template <std::size_t i>
    using Model = typename SubTypeTag::template Model<i>;
    template <std::size_t i>
    using Stencil = typename SubTypeTag::template Stencil<i>;
    template <std::size_t i>
    using Grid = typename SubTypeTag::template GridView<i>::Grid;
    template <std::size_t i>
    using ElementContext = typename SubTypeTag::template ElementContext<i>;
    template <std::size_t i>
    using MaterialLaw = typename SubTypeTag::template MaterialLaw<i>;

    using Evaluation = GetPropType<typename SubTypeTag::template SubDomain<0>::TypeTag,
                                   Properties::Evaluation>;
    using GlobalEqVector = typename SubTypeTag::GlobalEqVector;
    using JacobianMatrix = typename SubTypeTag::JacobianMatrix;

    enum {
        numPhases = getPropValue<typename SubTypeTag::template SubDomain<0>::TypeTag,
                                 Properties::NumPhases>(),
        EnableGravity_ = getPropValue<typename SubTypeTag::template SubDomain<0>::TypeTag,
                                      Properties::EnableGravity>()
    };
    enum { dimWorld = MortarGrid::dimensionworld, dim = MortarGrid::dimension };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using EvalDimVector = Dune::FieldVector<Evaluation, dimWorld>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using EvalDimMatrix = Dune::FieldMatrix<Evaluation, dimWorld, dimWorld>;
    using Toolbox = Opm::MathToolbox<Evaluation>;

   public:
    DarcyCoupler(Simulator<0>& simulator0, Simulator<1>& simulator1)
        : simulator0_{simulator0}, simulator1_{simulator1} {
        // We check that the two subproblems have the same number of phases.
        if (getPropValue<typename SubTypeTag::template SubDomain<0>::TypeTag,
                         Properties::NumPhases>() !=
            getPropValue<typename SubTypeTag::template SubDomain<1>::TypeTag,
                         Properties::NumPhases>())
            throw std::runtime_error(
                "Not implemented: Can only couple two models with the same number of "
                "faces ");
        // And the same type of evaluation.
        if (!(std::is_same<
                GetPropType<typename SubTypeTag::template SubDomain<0>::TypeTag,
                            Properties::Evaluation>,
                GetPropType<typename SubTypeTag::template SubDomain<1>::TypeTag,
                            Properties::Evaluation>>::value))
            throw std::runtime_error(
                "Not implemented: Can only couple two models with the same Evaluation");
        // And both have the same gravity
        if (getPropValue<typename SubTypeTag::template SubDomain<0>::TypeTag,
                           Properties::EnableGravity>() !=
            getPropValue<typename SubTypeTag::template SubDomain<1>::TypeTag,
                           Properties::EnableGravity>())
            throw std::runtime_error(
                "Not implemented: Can only couple two models with the same number of "
                "faces ");

        vanguard_.reset(new Vanguard(asImp_()));
        finalizeInit_();
    }

    static void registerParameters() {
        Vanguard::registerParameters();
        EWOMS_REGISTER_PARAM(TypeTag, std::string, MappingFile,
                             "The file name of the mapping file to load");
    }

    /*!
     * \brief Return the volume flux of a fluid phase at the mortar cell's integration
     * point \f$[m^3/s / m^2]\f$
     *
     * This is the fluid volume of a phase per second and per square meter of face
     * area.
     *
     * \param elemCtx The element context of the mortar cell
     */
    void volumeFlux(const CouplingElementContext& elemCtx) {
        const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);
        const auto& face = stencil.template interiorFace<0>(0);
        auto focusDofIdx = elemCtx.focusDofIndex();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            flux_[phaseIdx] = 0.0;
            Evaluation p0;  // Pressure in model 0
            Evaluation p1;  // Pressure in model 1
            // Only carry along the derivative from the model we focus on
            if (focusDofIdx == 0) {
                p0 = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0)
                         .fluidState()
                         .pressure(phaseIdx);
                p1 = Opm::getValue(
                    elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0)
                        .fluidState()
                        .pressure(phaseIdx));
            }
            else if (focusDofIdx == 1) {
                p0 = Opm::getValue(
                    elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0)
                        .fluidState()
                        .pressure(phaseIdx));
                p1 = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0)
                         .fluidState()
                         .pressure(phaseIdx);
            }
            else
                DUNE_THROW(Dune::NotImplemented,
                           "Can only couple two degrees of freedom");

            Evaluation pface;
            auto deltay = p1 - p0;

            // Find geometric centers
            const auto& interiorPos = stencil.template subControlVolume<0>(0).globalPos();
            const auto& exteriorPos = stencil.template subControlVolume<1>(0).globalPos();
            auto& faceCenter = face.integrationPos();
            // Calculate aperture
            Scalar extrusionFactor1 =
                elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0)
                    .extrusionFactor();
            Scalar aperture = std::pow(extrusionFactor1, 1.0 / (dimWorld - dim));
            // Make sure faceNormal points outwards
            auto faceNormal = face.normal();
            Scalar sign = 0.0;
            for (unsigned dimIdx = 0; dimIdx < faceNormal.size(); ++dimIdx)
                sign += (exteriorPos[dimIdx] - interiorPos[dimIdx]) * faceNormal[dimIdx];
            if (sign < 0.0) faceNormal *= -1.0;

            // Calculate vector form cells to face
            auto distVec0 = faceCenter;
            DimVector distVec1(0.0);
            for (unsigned dimIdx = 0; dimIdx < faceNormal.size(); ++dimIdx) {
                distVec0[dimIdx] -= interiorPos[dimIdx];
                distVec1[dimIdx] -= faceNormal[dimIdx] * aperture / 2.0;
            }
            if (EnableGravity_) {
                const auto& gIn = elemCtx.template problem<0>().gravity(elemCtx, 0, 0);
                const auto& gEx = elemCtx.template problem<1>().gravity(elemCtx, 0, 0);

                const auto& intQuantsIn = elemCtx.template intensiveQuantities<0>(0, 0);
                const auto& intQuantsEx = elemCtx.template intensiveQuantities<1>(0, 0);

                // calculate the hydrostatic pressure at the integration point of the face
                Evaluation pStatIn;
                if (std::is_same<Scalar, Evaluation>::value ||
                    0 == static_cast<int>(focusDofIdx)) {
                    const Evaluation& rhoIn = intQuantsIn.fluidState().density(phaseIdx);
                    pStatIn = -rhoIn * (gIn * interiorPos);
                }
                else {
                    Scalar rhoIn =
                        Toolbox::value(intQuantsIn.fluidState().density(phaseIdx));
                    pStatIn = -rhoIn * (gIn * interiorPos);
                }
                Evaluation pStatEx;
                if (std::is_same<Scalar, Evaluation>::value ||
                    1 == static_cast<int>(focusDofIdx)) {
                    const Evaluation& rhoEx = intQuantsEx.fluidState().density(phaseIdx);
                    pStatEx = -rhoEx * (gEx * (exteriorPos - distVec1));
                }
                else {
                    Scalar rhoEx =
                        Toolbox::value(intQuantsEx.fluidState().density(phaseIdx));
                    pStatEx = -rhoEx * (gEx * (exteriorPos - distVec1));
                }
                p0 += pStatIn;
                p1 += pStatEx;
            }

            Scalar distSquared = 0.0;
            for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
                Scalar tmp = exteriorPos[dimIdx] - interiorPos[dimIdx];
                distSquared += tmp * tmp;
            }
            // divide the gradient by the squared distance between the centers of the
            // sub-control volumes: the gradient is the normalized directional vector
            // between the two centers times the ratio of the difference of the values and
            // their distance, i.e., d/abs(d) * delta y / abs(d) = d*delta y / abs(d)^2.
            EvalDimVector quantityGrad;
            for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
                Scalar tmp = exteriorPos[dimIdx] - interiorPos[dimIdx];
                quantityGrad[dimIdx] =
                    (tmp / distSquared) * (p1 - p0);  // deltay * (tmp / distSquared);
            }

            Evaluation K0 = 0.0;
            Evaluation K1 = 0.0;
            const auto K0mat = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0)
                                   .intrinsicPermeability();
            const auto K1mat = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0)
                                   .intrinsicPermeability();
            for (unsigned dimIdx = 0; dimIdx < faceNormal.size(); dimIdx++) {
                K0 += Opm::abs(K0mat[dimIdx] * faceNormal);
                K1 += Opm::abs(K1mat[dimIdx] * faceNormal);
            }

            // determine the upstream and downstream DOFs
            short upstreamDofIdx;
            short downstreamDofIdx;
            Evaluation tmp = 0.0;
            for (unsigned dimIdx = 0; dimIdx < faceNormal.size(); ++dimIdx)
                tmp += quantityGrad[dimIdx] * faceNormal[dimIdx];

            if (tmp > 0) {
                upstreamDofIdx = 1;
                downstreamDofIdx = 0;
            }
            else {
                upstreamDofIdx = 0;
                downstreamDofIdx = 1;
            }
            // we only carry the derivatives along if the upstream DOF is the one which
            // we currently focus on
            Evaluation Kr0;
            Evaluation Kr1;
            const auto& materialParams0 = simulator0_.problem().materialLawParams(
                elemCtx.template subElementContext<0>(), 0, 0);
            const auto& materialParams1 = simulator1_.problem().materialLawParams(
                elemCtx.template subElementContext<1>(), 0, 0);
            if (upstreamDofIdx == 0) {
                Kr0 = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0)
                          .relativePermeability(phaseIdx);
                // const auto& fluidState = elemCtx.template intensiveQuantities<0>(0,
                // /*timeIdx=*/0).fluidState(); EvalDimVector saturations;
                // MaterialLaw<0>::saturations(saturations, materialParams0, fluidState);
                // if (phaseIdx == MaterialLaw<0>::Traits::wettingPhaseIdx) {
                //     Kr1 = MaterialLaw<1>::twoPhaseSatKrw(materialParams1,
                //     saturations[MaterialLaw<1>::Traits::wettingPhaseIdx]);
                //     //Kr0 = MaterialLaw<0>::twoPhaseSatKrw(materialParams0,
                //     saturations[MaterialLaw<0>::Traits::wettingPhaseIdx]);
                // } else if (phaseIdx == MaterialLaw<0>::Traits::nonWettingPhaseIdx) {
                //     Kr1 = MaterialLaw<1>::twoPhaseSatKrn(materialParams1,
                //     saturations[MaterialLaw<1>::Traits::wettingPhaseIdx]);
                //     //Kr0 = MaterialLaw<0>::twoPhaseSatKrn(materialParams0,
                //     saturations[MaterialLaw<0>::Traits::wettingPhaseIdx]);
                // } else
                //     throw std::runtime_error("asdg");
                Kr1 = Kr0;
                if (static_cast<int>(focusDofIdx) == 1) {
                    Kr1 = Opm::getValue(Kr1);
                    Kr0 = Opm::getValue(Kr0);
                }
            }
            else {
                Kr1 = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0)
                          .relativePermeability(phaseIdx);
                // const auto& fluidState = elemCtx.template intensiveQuantities<1>(0,
                // /*timeIdx=*/0).fluidState(); EvalDimVector saturations;
                // MaterialLaw<0>::saturations(saturations, materialParams0, fluidState);
                // if (phaseIdx == MaterialLaw<0>::Traits::wettingPhaseIdx) {
                //     //Kr1 = MaterialLaw<1>::twoPhaseSatKrw(materialParams1,
                //     saturations[MaterialLaw<1>::Traits::wettingPhaseIdx]); Kr0 =
                //     MaterialLaw<0>::twoPhaseSatKrw(materialParams0,
                //     saturations[MaterialLaw<0>::Traits::wettingPhaseIdx]);
                // } else if (phaseIdx == MaterialLaw<0>::Traits::nonWettingPhaseIdx) {
                //     //Kr1 = MaterialLaw<1>::twoPhaseSatKrn(materialParams1,
                //     saturations[MaterialLaw<1>::Traits::wettingPhaseIdx]); Kr0 =
                //     MaterialLaw<0>::twoPhaseSatKrn(materialParams0,
                //     saturations[MaterialLaw<0>::Traits::wettingPhaseIdx]);
                // } else
                //     throw std::runtime_error("asdg");

                // Kr0 = Opm::harmonicMean(Kr1, Kr0);
                Kr0 = Kr1;
                if (static_cast<int>(focusDofIdx) == 0) {
                    Kr1 = Opm::getValue(Kr1);
                    Kr0 = Opm::getValue(Kr0);
                }
            }
            if (dim == 0) {
                if (upstreamDofIdx == 0) {
                    Kr0 = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0)
                              .relativePermeability(phaseIdx);
                    Kr1 = Kr0;
                    if (static_cast<int>(focusDofIdx) == 1) {
                        Kr1 = Opm::getValue(Kr1);
                        Kr0 = Opm::getValue(Kr0);
                    }
                }
                else {
                    Kr1 = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0)
                              .relativePermeability(phaseIdx);
                    Kr0 = Kr1;
                    if (static_cast<int>(focusDofIdx) == 0) {
                        Kr1 = Opm::getValue(Kr1);
                        Kr0 = Opm::getValue(Kr0);
                    }
                }
            }
            Evaluation T0;
            Evaluation T1;
            Evaluation T;
            // If the two subdomains have the same dimension, we calculate the
            // permeability as the harmonic mean.
            if (false && Grid<0>::dimension == Grid<1>::dimension &&
                Grid<0>::dimension != 1) {
                // Entry-wise harmonic mean. this is almost certainly wrong if
                // you have off-main diagonal entries in your permeabilities!
                Scalar ndotDistIn = faceNormal * distVec0;
                Scalar ndotDistExt = faceNormal * distVec1;
                Scalar distSquaredIn = distVec0 * distVec0;
                Scalar distSquaredEx = distVec1 * distVec1;

                T0 = K0 * Kr0 * ndotDistIn / distSquaredIn;
                T1 = K1 * Kr1 * ndotDistIn / distSquaredIn;

                pface = (T0 * p0 + T1 * p1) / (T0 + T1);
                // K = Opm::harmonicMean(T0, T1);
#if false
                for (unsigned i = 0; i < dimWorld; ++i)
                    for (unsigned j = 0; j < dimWorld; ++j)
                        K[i][j] = Opm::harmonicMean(K0[i][j], K1[i][j]);
#endif
            }
            else if (true) {
                // The harmonic mean is not sufficient when the two subdomains are of
                // different dimension. We then have to scale the permeability by the
                // distance to the interface (which for the lower-dimensional domain is
                // equal the aperture / 2).
                {
                    Scalar ndotDistIn = faceNormal * distVec0;
                    Scalar ndotDistExt = -1.0 * (faceNormal * distVec1);
                    Scalar distSquaredIn = distVec0 * distVec0;
                    Scalar distSquaredEx = distVec1 * distVec1;

                    T0 = K0 * ndotDistIn / distSquaredIn;
                    T1 = K1 * ndotDistExt / distSquaredEx;
                    T = T0 * T1 / (T0 + T1);
                    if (Kr1 <= 0) {
                        pface = p0;
                    }
                    else if (Kr0 <= 0) {
                        pface = p1;
                    }
                    else {
                        pface = (T0 * p0 + T1 * p1) / (T0 + T1);
                    }
                }
            }

            // we only carry the derivatives along if the upstream DOF is the one which
            // we currently focus on
            Evaluation mobility;
            Evaluation rho;
            Evaluation viscosity;

            if (upstreamDofIdx == 0) {
                if (static_cast<int>(focusDofIdx) == 0) {
                    mobility = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0)
                                   .mobility(phaseIdx);
                    viscosity = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0)
                                    .fluidState()
                                    .viscosity(phaseIdx);
                    rho = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0)
                              .fluidState()
                              .density(phaseIdx);
                }
                else {
                    mobility = Opm::getValue(
                        elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0)
                            .mobility(phaseIdx));
                    rho = Opm::getValue(
                        elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0)
                            .fluidState()
                            .density(phaseIdx));
                    viscosity = Opm::getValue(
                        elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0)
                            .fluidState()
                            .viscosity(phaseIdx));
                }
            }
            else {
                if (static_cast<int>(focusDofIdx) == 1) {
                    mobility = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0)
                                   .mobility(phaseIdx);
                    rho = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0)
                              .fluidState()
                              .density(phaseIdx);
                    viscosity = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0)
                                    .fluidState()
                                    .viscosity(phaseIdx);
                }
                else {
                    mobility = Opm::getValue(
                        elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0)
                            .mobility(phaseIdx));
                    rho = Opm::getValue(
                        elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0)
                            .fluidState()
                            .density(phaseIdx));
                    viscosity = Opm::getValue(
                        elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0)
                            .fluidState()
                            .viscosity(phaseIdx));
                }
            }

            // Evaluation filterVelocity;
            // filterVelocity = -T * (p0 - p1);//T0 * (p0 - pface);
            // filterVelocity *= K;
            // K.mv(quantityGrad, filterVelocity);
            // filterVelocity *= -mobility;
            flux_[phaseIdx] += T * mobility * rho * (p0 - p1);
            // if (upstreamDofIdx == 0)
            //     flux_[phaseIdx] += -T0 / viscosity * rho * (pface - p0);
            // else
            //     flux_[phaseIdx] += T1 / viscosity * rho * (pface - p1);

            // for (unsigned i = 0; i < faceNormal.size(); ++i)
            //     flux_[phaseIdx] += rho * (filterVelocity[i] * faceNormal[i]);

            // Scale the face area by the higher dimensional extrusion factor
            Scalar alpha =
                face.area() * elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0)
                                  .extrusionFactor();
            Opm::Valgrind::CheckDefined(alpha);
            assert(alpha > 0.0);
            assert(Opm::isfinite(alpha));
            flux_[phaseIdx] *= alpha;
        }
    }

    /*!
     * \brief Linearize and calculate the residual and its jacobian
     *
     * The advective Jacobian and the residual is added to the coupling
     * Jacobian and residual.
     */
    void advectiveFluxCoupling() {
        const auto& idxSet = mortarView().indexSet();
        const auto& gridViewL = simulator0_.gridView();
        const auto& gridViewR = simulator1_.gridView();
        const auto simulators = std::forward_as_tuple(simulator0_, simulator1_);
        CouplingElementContext elemCtx(mortarView(), simulators, *map_);
        for (const auto& e : elements(mortarView())) {
            elemCtx.updateAll(e);

            // compute the local residual and its Jacobian
            unsigned numPrimaryDof = elemCtx.numPrimaryDof(/*timeIdx=*/0);
            for (unsigned focusDofIdx = 0; focusDofIdx < numPrimaryDof; ++focusDofIdx) {
                elemCtx.setFocusDofIndex(focusDofIdx);
                elemCtx.updateAllExtensiveQuantities();
                // calculate the local residual
                asImp_().volumeFlux(elemCtx);
                // Global here refers to global within each domain
                unsigned globI =
                    elemCtx.template globalSpaceIndex<0>(/*spaceIdx=*/0, /*timeIdx=*/0);
                unsigned globJ =
                    elemCtx.template globalSpaceIndex<1>(/*spaceIdx=*/0, /*timeIdx=*/0);
                for (unsigned eqIdx = 0; eqIdx < numPhases; ++eqIdx) {
                    assert(Opm::isfinite(flux_[eqIdx]));
                    if (focusDofIdx == 0) {
                        residual_[_0][globI][eqIdx] += Opm::getValue(flux_[eqIdx]);
                        for (unsigned pvIdx = 0; pvIdx < numPhases; ++pvIdx) {
                            jacobian_[_0][_0][globI][globI][eqIdx][pvIdx] +=
                                flux_[eqIdx].derivative(pvIdx);
                            jacobian_[_1][_0][globJ][globI][eqIdx][pvIdx] -=
                                flux_[eqIdx].derivative(pvIdx);
                        }
                    }
                    else if (focusDofIdx == 1) {
                        residual_[_1][globJ][eqIdx] -= Opm::getValue(flux_[eqIdx]);

                        for (unsigned pvIdx = 0; pvIdx < numPhases; ++pvIdx) {
                            jacobian_[_0][_1][globI][globJ][eqIdx][pvIdx] +=
                                flux_[eqIdx].derivative(pvIdx);
                            jacobian_[_1][_1][globJ][globJ][eqIdx][pvIdx] -=
                                flux_[eqIdx].derivative(pvIdx);
                        }
                    }
                    else
                        std::logic_error(
                            "NotImplemented: advectiveFluxCoupling can only couple two "
                            "dofs");
                }
            }
        }
    }

    /*!
     * \brief Linearize the full system of non-linear equations.
     *
     * This means the spatial subdomains plus all auxiliary equations.
     */
    void linearize() {
        if (firstIteration_) {
            initFirstIteration_();
        }
        resetSystem_();  // reset the global linear system of equations.

        // Linearize coupling
        advectiveFluxCoupling();
    }

    /*!
     * \brief Return reference to global Jacobian matrix backend.
     */
    JacobianMatrix& jacobian() { return jacobian_; }

    /*!
     * \brief Return reference to global residual vector.
     */
    GlobalEqVector& residual() { return residual_; }

    /*!
     * \brief Return the Dune::IndexConstant of coupling domain 0
     */
    template <std::size_t k, typename std::enable_if_t<(k == 0), int> = 0>
    const auto subDomainIndex() const {
        DomainI K;
        return K;
    }

    /*!
     * \brief Return the Dune::IndexConstant of coupling domain 1
     */
    template <std::size_t k, typename std::enable_if_t<(k == 1), int> = 0>
    const auto subDomainIndex() const {
        DomainJ K;
        return K;
    }

    /*!
     * \brief Return the mortar grid view for which the simulation is done
     */
    const MortarView& mortarView() const { return vanguard_->gridView(); }

    /*!
     * \brief Return the projection mapper of the coupling
     */
    const CouplingMapper& projectionMapper() const { return *map_; }

    /*!
     * \brief Return the first model of the coupling
     */
    template <std::size_t k, typename std::enable_if_t<(k == 0), int> = 0>
    const auto& model() const {
        return simulator0_.model();
    }

    /*!
     * \brief Return the second model of the coupling
     */
    template <std::size_t k, typename std::enable_if_t<(k == 1), int> = 0>
    const auto& model() const {
        return simulator1_.model();
    }

   protected:
    // reset the linear system of equations.
    void resetSystem_() {
        residual_ = 0.0;
        // zero all matrix entries
        jacobian_ = 0.0;
    }

    /*!
     * Set the jacobian build mode
     */
    void setJacobianBuildMode() {
        using namespace Dune::Hybrid;
        forEach(jacobian_, [](auto& jacRow) {
            forEach(jacRow, [](auto& jacBlock) {
                using BlockType = std::decay_t<decltype(jacBlock)>;
                if (jacBlock.buildMode() == BlockType::BuildMode::unknown)
                    jacBlock.setBuildMode(BlockType::BuildMode::random);
                else if (jacBlock.buildMode() != BlockType::BuildMode::random)
                    DUNE_THROW(Dune::NotImplemented,
                               "Only BCRS matrices with random build mode are supported "
                               "at the moment");
            });
        });
    }

    // Set the sparcity pattern for a Jacobian matrix
    template <std::size_t i, std::size_t j, class Set>
    void reserve_(const std::vector<Set>& sparsityPattern) {
        Dune::index_constant<i> I;
        Dune::index_constant<j> J;
        auto rows_ = model<i>().numTotalDof();
        auto cols_ = model<j>().numTotalDof();

        // make sure sparsityPattern is consistent with number of rows
        assert(rows_ == sparsityPattern.size());

        // allocate space for the rows of the matrix
        for (size_t dofIdx = 0; dofIdx < rows_; ++dofIdx)
            jacobian_[I][J].setrowsize(dofIdx, sparsityPattern[dofIdx].size());
        jacobian_[I][J].endrowsizes();

        // fill the rows with indices. each degree of freedom talks to
        // all of its neighbors. (it also talks to itself since
        // degrees of freedom are sometimes quite egocentric.)
        for (size_t dofIdx = 0; dofIdx < rows_; ++dofIdx) {
            auto nIt = sparsityPattern[dofIdx].begin();
            auto nEndIt = sparsityPattern[dofIdx].end();
            for (; nIt != nEndIt; ++nIt) jacobian_[I][J].addindex(dofIdx, *nIt);
        }
        jacobian_[I][J].endindices();
    }

    // Initialize the Jacobian matrix
    void createMatrix() {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(jacobian_)), [&](const auto domainI) {
            forEach(integralRange(Dune::Hybrid::size(jacobian_[domainI])),
                    [&](const auto domainJ) { setJacobianPattern_(domainI, domainJ); });
        });
    }

    void initFirstIteration_() {
        firstIteration_ = false;
        setJacobianBuildMode();
        createMatrix();
    }

    void finalizeInit_() {
        const std::string mappingFileName =
            getPropValue<TypeTag, Properties::MappingFile>();

        if (mappingFileName.size() > 0)
            map_.reset(new CouplingMapper(mappingFileName, mortarView(),
                                          simulator0_.gridView(),
                                          simulator1_.gridView()));
        else
            map_.reset(new CouplingMapper(mortarView(), simulator0_.gridView(),
                                          simulator1_.gridView()));
    }

    // Set the jacobian sparcity pattern
    template <std::size_t i, std::size_t j>
    void setJacobianPattern_(Dune::index_constant<i> I, Dune::index_constant<j> J) {
        typedef std::set<unsigned> NeighborSet;
        std::vector<NeighborSet> sparsityPattern(model<i>().numTotalDof());

        auto rows_ = model<i>().numTotalDof();
        auto cols_ = model<j>().numTotalDof();

        jacobian_[I][J].setSize(rows_, cols_);
        residual_[I].resize(rows_);
        if (i == j) {
            addDiagonalPattern_<i>(sparsityPattern);
        }
        else {
            addOffDiagonalPattern_<i, j>(sparsityPattern);
        }

        // allocate raw matrix

        // create matrix structure based on sparsity pattern
        reserve_<i, j>(sparsityPattern);
    }

    // Add the off diagonal part of the sparcity pattern
    template <std::size_t i, std::size_t j, class Set>
    void addOffDiagonalPattern_(std::vector<Set>& sparsityPattern) {
        if (0 == i && 1 == j) addCouplerPattern_(false, sparsityPattern);
        if (0 == j && 1 == i) addCouplerPattern_(true, sparsityPattern);

        return;
    }

    // Add the copuling pattern
    template <class Set>
    void addCouplerPattern_(bool swap, Set& sparsityPattern) {
        // for the main model, find out the global indices of the neighboring degrees of
        // freedom of each primary degree of freedom
        int numDofs;
        if (!swap)
            numDofs = model<0>().numTotalDof();
        else
            numDofs = model<1>().numTotalDof();

        const auto& idxSet = mortarView().indexSet();
        for (const auto& elem : elements(mortarView())) {
            const auto idx = idxSet.index(elem);
            // const auto face = model1().gridView().grid().entity(map_->[i][idx]);
            const auto element1 = map_->template toElement<0>(elem);
            unsigned elIdx1 = model<0>().gridView().indexSet().index(element1);

            const auto element2 = map_->template toElement<1>(elem);
            unsigned elIdx2 = model<1>().gridView().indexSet().index(element2);
            if (!swap)
                sparsityPattern[elIdx1].insert(elIdx2);
            else
                sparsityPattern[elIdx2].insert(elIdx1);
        }
    }

    // Add the diagonal sparsity patter. It is here assumed we have a cell-centered
    // final volume discretization (tpfa).
    template <std::size_t i, class Set>
    void addDiagonalPattern_(std::vector<Set>& sparsityPattern) const {
        // for the main model, find out the global indices of the neighboring degrees of
        // freedom of each primary degree of freedom
        Stencil<i> stencil(model<i>().gridView(), model<i>().dofMapper());

        auto elemIt = model<i>().gridView().template begin<0>();
        const auto elemEndIt = model<i>().gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            stencil.update(elem);

            for (unsigned primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof();
                 ++primaryDofIdx) {
                unsigned myIdx = stencil.globalSpaceIndex(primaryDofIdx);

                for (unsigned dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx) {
                    unsigned neighborIdx = stencil.globalSpaceIndex(dofIdx);
                    sparsityPattern[myIdx].insert(neighborIdx);
                }
            }
        }

        // add the additional neighbors and degrees of freedom caused by the auxiliary
        // equations
        size_t numAuxMod = model<i>().numAuxiliaryModules();
        for (unsigned auxModIdx = 0; auxModIdx < numAuxMod; ++auxModIdx)
            model<i>().auxiliaryModule(auxModIdx)->addNeighbors(sparsityPattern);
    }

    // private:
    Dune::index_constant<0> _0;
    Dune::index_constant<1> _1;

    Simulator<0>& simulator0_;
    Simulator<1>& simulator1_;
    std::unique_ptr<CouplingMapper> map_;
    JacobianMatrix jacobian_;
    GlobalEqVector residual_;
    bool verbose_{true};
    std::unique_ptr<Vanguard> vanguard_;

    Evaluation flux_[numPhases];

    bool firstIteration_{true};

   private:
    Implementation& asImp_() { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const {
        return *static_cast<const Implementation*>(this);
    }
};

}  // namespace Opm

#endif
