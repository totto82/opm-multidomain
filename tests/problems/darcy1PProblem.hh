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
 * \copydoc Opm::Darcy1PProblem
 */

#ifndef OPM_DARCY_1_P_PROBLEM_HH
#define OPM_DARCY_1_P_PROBLEM_HH


#include "opm/grid/polyhedralgrid.hh"
#include <opm/material/components/Unit.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/fluidsystems/LiquidPhase.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/models/utils/basicproperties.hh>
#include <opm/simulators/linalg/parallelistlbackend.hh>

#include <dune/common/version.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <sstream>
#include <string>

namespace Opm
{
template <class TypeTag>
class Darcy1PProblem;
}

namespace Opm::Properties {
// Create new type tags
namespace TTag {
struct Darcy1PBaseProblem { using InheritsFrom = std::tuple<ImmiscibleSinglePhaseModel, ImplicitModel>; };
} // end namespace TTag

template<class TypeTag, class MyTypeTag>
struct NameSurfix { using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DomainDim { using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct WorldDim { using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct GridDim { using type = UndefinedProperty; };

template<class TypeTag>
struct NameSurfix<TypeTag, TTag::Darcy1PBaseProblem> { static constexpr auto value = ""; };

template<class TypeTag>
struct WorldDim<TypeTag, TTag::Darcy1PBaseProblem> { static constexpr int value = 2; };

template<class TypeTag>
struct DomainDim<TypeTag, TTag::Darcy1PBaseProblem> { static constexpr int value = 2; };

template<class TypeTag>
struct Fluid<TypeTag, TTag::Darcy1PBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Opm::LiquidPhase<Scalar, Opm::Unit<Scalar>>;
};

template<class TypeTag>
struct Problem<TypeTag, TTag::Darcy1PBaseProblem> { using type = Opm::Darcy1PProblem<TypeTag>; };

template<class TypeTag>
struct Grid<TypeTag, TTag::Darcy1PBaseProblem> {
private:
    enum{ gridDim = getPropValue<TypeTag, Properties::GridDim>() };
    enum{ worldDim = getPropValue<TypeTag, Properties::WorldDim>() };
public:
    using type = Dune::PolyhedralGrid<gridDim, worldDim>;
};

template<class TypeTag>
struct EnableGravity<TypeTag, TTag::Darcy1PBaseProblem> { static constexpr bool value = false; };

template<class TypeTag>
struct EndTime<TypeTag, TTag::Darcy1PBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1;
};

template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::Darcy1PBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1;
};

template<class TypeTag>
struct VtkWriteSaturations<TypeTag, TTag::Darcy1PBaseProblem> { static constexpr bool value = false; };
template<class TypeTag>
struct VtkWriteMobilities<TypeTag, TTag::Darcy1PBaseProblem> { static constexpr bool value = false; };
template<class TypeTag>
struct VtkWriteRelativePermeabilities<TypeTag, TTag::Darcy1PBaseProblem> { static constexpr bool value = false; };

} // end namespace Opm::Properties


namespace Opm
{

/*!
 * \ingroup TestProblems
 *
 * \brief Test for the immisicible VCVF discretization with only a single phase
 *
 * This problem is inspired by groundwater flow. Don't expect it to be
 * realistic, though: For two dimensions, the domain size is 1m times
 * 1m. On the left and right of the domain, no-flow boundaries are
 * used, while at the top and bottom free flow boundaries with a
 * pressure of 2 bar and 1 bar are used. The center of the domain is
 * occupied by a rectangular lens of lower permeability.
 */
template <class TypeTag>
class Darcy1PProblem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    enum
    {
        numPhases = FluidSystem::numPhases,

        // Grid and world dimension
        dim = getPropValue<TypeTag, Properties::DomainDim>(),
        dimWorld = GridView::dimensionworld,

        // indices of the primary variables
        pressure0Idx = Indices::pressure0Idx
    };


    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;


public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    Darcy1PProblem(Simulator &simulator) : ParentType(simulator)
    {
        eps_ = 1e-10;
        initiated_ = true;
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();
        double K;
        if (dim == 2)
            K = 1;
        else if (dim == 1)
            K = 1e4;
        else if (dim == 0)
            K = 1e4;
        else
            throw std::logic_error("Unknown dimension");
        intrinsicPerm_ = this->toDimMatrix_(K);
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, std::string, NameSurfix, "Id to append to grid name");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "darcy1PProblem_" << getPropValue<TypeTag, Properties::NameSurfix>();
        return oss.str();
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        this->model().checkConservativeness();

        // Calculate storage terms
        EqVector storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0)
        {
            std::cout << "Storage: " << storage << std::endl
                      << std::flush;
        }
#endif // NDEBUG
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context &context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    {
        return 273.15 + 10;
    } // 10C

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context OPM_UNUSED,
                    unsigned spaceIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED) const
    {
        return 0.4;
    }
    template <class Context>
    Scalar extrusionFactor(const Context &context OPM_UNUSED,
                           unsigned spaceIdx OPM_UNUSED,
                           unsigned timeIdx OPM_UNUSED) const
    {
        return extrusionFactor();
    }

    Scalar extrusionFactor() const
    {
        return std::pow(1e-4, dimWorld - dim);
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context,
                                           unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        return intrinsicPerm_;
    }
    //! \}
    /*!
     * \name Boundary conditions
     */
    //! \{

    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);

        if (onRightBoundary_(globalPos))
        {
            Scalar pressure;
            Scalar T = temperature(context, spaceIdx, timeIdx);

            pressure = 1.0;
            Opm::ImmiscibleFluidState<Scalar, FluidSystem, /*storeEnthalpy=*/false> fs;
            fs.setSaturation(/*phaseIdx=*/0, 1.0);
            fs.setPressure(/*phaseIdx=*/0, pressure);
            fs.setTemperature(T);

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            paramCache.updateAll(fs);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                fs.setDensity(phaseIdx, FluidSystem::density(fs, paramCache, phaseIdx));
                fs.setViscosity(phaseIdx, FluidSystem::viscosity(fs, paramCache, phaseIdx));
            }

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (onLeftBoundary_(globalPos))
        {
            values.setMassRate(-1.0);
        }
        else
        {
            // no flow boundary
            values.setNoFlow();
        }
    }

    //! \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class PrimaryVariables, class Context>
    void initial(PrimaryVariables &values,
                 const Context &context OPM_UNUSED,
                 unsigned spaceIdx OPM_UNUSED,
                 unsigned timeIdx OPM_UNUSED) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        values[pressure0Idx] = (1 - pos[0]);
    }

    /*!
     * \copydoc FvBaseProblem::source
     */
    template <class Context>
    void source(RateVector &rate,
                const Context &context OPM_UNUSED,
                unsigned spaceIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED) const
    {
        rate = Scalar(0.0);
    }

    //! \}

protected:
    bool onLeftBoundary_(const GlobalPosition &pos) const
    {
        return pos[0] < eps_;
    }

    bool onRightBoundary_(const GlobalPosition &pos) const
    {
        return pos[0] > 1 - eps_;
    }
    
    bool onBottomBoundary_(const GlobalPosition &pos) const
    {
        return pos[1] < eps_;
    }

    bool onTopBoundary_(const GlobalPosition &pos) const
    {
        return pos[1] > 1 - eps_;
    }

    DimMatrix intrinsicPerm_;
    Scalar eps_;

    bool initiated_;
};

} // namespace Opm

#endif
