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
 * \copydoc Opm::DualPorosityCoupler
 */

#ifndef OPM_DUAL_POROSITY_COUPLER_HH
#define OPM_DUAL_POROSITY_COUPLER_HH

#include <opm/models/utils/propertysystem.hh>

#include <opm/material/densead/Math.hpp>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/multiDomain/couplingelementcontext.hh>
#include <opm/models/multiDomain/ecfvcouplingstencil.hh>
#include <opm/models/multiDomain/multidomainmapper.hh>
#include <opm/models/multiDomain/multidomaincoupler.hh>
#include <opm/models/multiDomain/multidomainproperties.hh>

#include <dune/common/indices.hh>

#include <stdexcept>

namespace Opm {
template <class TypeTag>
class DualPorosityCoupler;
} // namespace Opm

BEGIN_PROPERTIES
NEW_TYPE_TAG(DualPorosityCoupler, INHERITS_FROM(DarcyCoupler));

SET_TYPE_PROP(DualPorosityCoupler, Coupler, Opm::DualPorosityCoupler<TypeTag>);
SET_TYPE_PROP(DualPorosityCoupler, CouplingMapper, Opm::ElementElementMapper<TypeTag>);
SET_TYPE_PROP(DualPorosityCoupler, CouplingElementContext, Opm::CouplingElementContext<TypeTag>);
SET_TYPE_PROP(DualPorosityCoupler, Simulator, Opm::DualPorosityCoupler<TypeTag>);

SET_INT_PROP(DualPorosityCoupler, TimeDiscHistorySize, 0);
SET_STRING_PROP(DualPorosityCoupler, MappingFile, "");

SET_PROP(DualPorosityCoupler, Stencil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, CouplingMapper) Mapper;
    typedef typename GET_PROP_TYPE(TypeTag, SubTypeTag) SubTypeTag;

public:
    typedef Opm::EcfvDualPorosityStencil<Scalar, GridView, Mapper, SubTypeTag> type;
};

END_PROPERTIES

namespace Opm {

/*!
 * \ingroup MultiDomainModel
 *
 * \brief Couple two domains in a multidomain model
 * 
 * A multidomain model consist of several subdomains which are initiated as
 * standard OPM models. This class defines the coupling between two subdomains
 * by a dual porosity law. The class allows between the cells of two domains of
 * equal dimension
 * 
 * You should not expect the DualPorosityCoupler to work for anything else than immicible
 * models.
 * 
*/
template <class TypeTag>
class DualPorosityCoupler: public Opm::DarcyCoupler<TypeTag> {
    using ParentType = Opm::DarcyCoupler<TypeTag>;
    typedef typename GET_PROP_TYPE(TypeTag, SubTypeTag) SubTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, CouplingElementContext) CouplingElementContext;


    typedef typename GET_PROP_TYPE(typename SubTypeTag::template SubDomain<0>::TypeTag, Evaluation) Evaluation;

    enum {
        numPhases = GET_PROP_VALUE(typename SubTypeTag::template SubDomain<0>::TypeTag, NumPhases),
        EnableGravity_ = GET_PROP_VALUE(typename SubTypeTag::template SubDomain<0>::TypeTag, EnableGravity)
    };


public:
    using ParentType::ParentType;

    static void registerParameters()
    {
        ParentType::registerParameters();
    }
    /*!
     * \brief Return the volume flux of a fluid phase at the mortar cell's integration point
     *        \f$[m^3/s / m^2]\f$
     *
     * This is the fluid volume of a phase per second and per square meter of face
     * area.
     *
     * \param elemCtx The element context of the mortar cell
     */
    void volumeFlux(const CouplingElementContext& elemCtx)
    {
        const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);

        const auto& face = stencil.template interiorFace<0>(0);
        auto focusDofIdx = elemCtx.focusDofIndex();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            this->flux_[phaseIdx] = 0.0;
            Evaluation p0; // Pressure in model 0
            Evaluation p1; // Pressure in model 1

            // Only carry along the derivative from the model we focus on
            if (focusDofIdx == 0) {
                p0 = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).fluidState().pressure(phaseIdx);
                p1 = Opm::getValue(elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).fluidState().pressure(phaseIdx));
            } else if (focusDofIdx == 1) {
                p0 = Opm::getValue(elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).fluidState().pressure(phaseIdx));
                p1 = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).fluidState().pressure(phaseIdx);
            } else
                DUNE_THROW(Dune::NotImplemented, "Can only couple two degrees of freedom");

            Evaluation pface;
            auto deltay = p1 - p0;

            short upstreamDofIdx;
            short downstreamDofIdx;
            if (deltay > 0) {
                upstreamDofIdx = 1;
                downstreamDofIdx = 0;
            } else {
                upstreamDofIdx = 0;
                downstreamDofIdx = 1;
            }

            Evaluation mobility;
            Evaluation rho;
            if (upstreamDofIdx == 0) {
                if (static_cast<int>(focusDofIdx) == 0) {
                    mobility = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).mobility(phaseIdx);
                    rho = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).fluidState().density(phaseIdx);
                } else {
                    mobility = Opm::getValue(elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).mobility(phaseIdx));
                    rho = Opm::getValue(elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).fluidState().density(phaseIdx));
                }
            } else {
                if (static_cast<int>(focusDofIdx) == 1) {
                    mobility = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).mobility(phaseIdx);
                    rho = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).fluidState().density(phaseIdx);
                } else {
                    mobility = Opm::getValue(elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).mobility(phaseIdx));
                    rho = Opm::getValue(elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).fluidState().density(phaseIdx));
                }
            }

            //const auto Km = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).intrinsicPermeability().infinity_norm();
            
            double Km = 21.0*1e-3 * 9.869233 * 1e-13;
            const auto Kf = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).intrinsicPermeability().infinity_norm();
            const auto phi = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).porosity();
            Evaluation Sf;
            Evaluation Sm;
            Evaluation rhom;
            if (static_cast<int>(focusDofIdx) == 0) {
                Sf = Opm::getValue(elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).fluidState().saturation(1));
                Sm = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).fluidState().saturation(1);
                rhom = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).fluidState().density(1);
            }
            else{
                Sf = elemCtx.template intensiveQuantities<1>(0, /*timeIdx=*/0).fluidState().saturation(1);
                Sm = Opm::getValue(elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).fluidState().saturation(1));
                rhom = Opm::getValue(elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).fluidState().density(1));
            }
            if (Sf < 0.0)
                 Sf = 0.0;
            if (Sf > 1.0)
                 Sf = 1.0;
            double beta = 8.48 * std::pow(10, -8);
            auto F = (1 - Opm::exp(-Opm::sqrt(Kf / Km) * Sf)) / (1 - Opm::exp(-Opm::sqrt(Kf / Km)));
            auto Tn = phi * rhom * beta * F * (0.39 - Sm);
            if (phaseIdx==0)
                this->flux_[phaseIdx] += Tn;
            else
                this->flux_[phaseIdx] -= Tn;


            // Scale the face area by the higher dimensional extrusion factor
            auto alpha = elemCtx.template intensiveQuantities<0>(0, /*timeIdx=*/0).extrusionFactor();
            auto mortarArea = stencil.element(0).geometry().volume();
            Opm::Valgrind::CheckDefined(alpha);
            assert(alpha > 0.0);
            assert(Opm::isfinite(alpha));
            this->flux_[phaseIdx] *= mortarArea * alpha;
        }
    }

};
} // namespace Opm

#endif
