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
#ifndef OPM_BENCHMARK_4_PROBLEM_HH
#define OPM_BENCHMARK_4_PROBLEM_HH

#include "darcy1PProblem.hh"

namespace Opm
{
template <class TypeTag>
class Benchmark4Problem;
}

namespace Opm::Properties {
// Create new type tags
namespace TTag {
struct Benchmark4Problem { using InheritsFrom = std::tuple<Darcy1PBaseProblem>; };
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::Benchmark4Problem> { using type = Opm::Benchmark4Problem<TypeTag>; };

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
class Benchmark4Problem : public Opm::Darcy1PProblem<TypeTag>
{
    using ParentType = Opm::Darcy1PProblem<TypeTag>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;

    enum
    {
        numPhases = FluidSystem::numPhases,
        // Grid and world dimension
        dim = getPropValue<TypeTag, Properties::DomainDim>(),
        dimWorld = GridView::dimensionworld,
    };

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    using ParentType::ParentType;

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix intrinsicPermeability(const Context &context,
                                          unsigned spaceIdx,
                                          unsigned timeIdx) const
    {
        double K;
        if (dim == 2)
            K = 10e-14;
        else if (dim == 1)
            K = 10e-8;
        else if (dim == 0)
            K = 10e-8;
        else
            throw std::logic_error("Unknown dimension");
        return this->toDimMatrix_(K);
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
        return std::pow(1e-2, dimWorld - dim);
    }
    //! \}
    /*!
     * \name Boundary conditions
     */
    //! \{
    std::string name() const
    {
        return "benchmark4_" + ParentType::name();
    }

    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);
        if (onDirichletBoundary_(globalPos))
        {
            Scalar pressure;
            Scalar T = this->temperature(context, spaceIdx, timeIdx);

            if (onInjectionBoundary_(globalPos))
                pressure = 1013250;
            else
                pressure = 0.0;

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
        else
        {
            // no flow boundary
            values.setNoFlow();
        }
    }
    static void registerParameters()
    {
        ParentType::registerParameters();
    }

protected:
    bool onDirichletBoundary_(GlobalPosition pos) const {
        return this->onLeftBoundary_(pos) || this->onRightBoundary_(pos);
    }
    bool onInjectionBoundary_(GlobalPosition pos) const {
        return this->onLeftBoundary_(pos);
    }

    bool onLeftBoundary_(const GlobalPosition &pos) const
    {
        return pos[0] < eps_;
    }

    bool onRightBoundary_(const GlobalPosition &pos) const
    {
        return pos[0] > 700 - eps_;
    }
    
    bool onBottomBoundary_(const GlobalPosition &pos) const
    {
        return pos[1] < eps_;
    }

    bool onTopBoundary_(const GlobalPosition &pos) const
    {
        return pos[1] > 600 - eps_;
    }
private:
    double eps_{1e-8};
};

} // namespace Opm

#endif
