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
#ifndef OPM_BENCHMARK3d_1_PROBLEM_HH
#define OPM_BENCHMARK3d_1_PROBLEM_HH

#include "darcy2PProblem.hh"

namespace Opm
{
template <class TypeTag>
class Benchmark3d1Problem;
}
namespace Opm::Properties {
// Create new type tags
namespace TTag {
struct Benchmark3d1Problem { using InheritsFrom = std::tuple<Darcy2pBaseProblem>; };
} // end namespace TTag

template<class TypeTag, class MyTypeTag>
struct LeftRight { using type = UndefinedProperty; };

template<class TypeTag>
struct Problem<TypeTag, TTag::Benchmark3d1Problem> { using type = Opm::Benchmark3d1Problem<TypeTag>; };

template<class TypeTag>
struct WorldDim<TypeTag, TTag::Benchmark3d1Problem> { static constexpr int value = 3; };
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
class Benchmark3d1Problem : public Opm::Darcy2pProblem<TypeTag>
{
    using ParentType = Opm::Darcy2pProblem<TypeTag>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using WettingPhase = GetPropType<TypeTag, Properties::WettingPhase>;
    using NonwettingPhase = GetPropType<TypeTag, Properties::NonwettingPhase>;

    enum
    {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        wettingPhaseIdx = FluidSystem::wettingPhaseIdx,
        nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx,

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
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);
        double tol = 1e-3;
        float f;
        double K;
        double Kh = 1e-5;
        double Kl = 1e-6;
        double Kf = 1e-3 / extrusionFactor(context, spaceIdx, timeIdx);
        if (dim == 3)
        {
            if (inHighPermRegion_(globalPos))
                K = Kh;
            else
                K = Kl;
        }
        else if (dim == 2)
        {
            K = Kf;
        }
        else{
            throw std::runtime_error("NotImplemented: Intrinsic permeability is only defined for dim=2 or 3");
        }
        return this->toDimMatrix_(K);
    }

    template <class Context>
    Scalar extrusionFactor(const Context &context ,
                           unsigned spaceIdx ,
                           unsigned timeIdx ) const
    {
        return extrusionFactor();
    }

    Scalar extrusionFactor() const
    {
        return std::pow(1e-2, dimWorld - dim);
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context ,
                    unsigned spaceIdx ,
                    unsigned timeIdx ) const
    {
                const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);
        double tol = 1e-3;
        double phi;
        if (dim == 3)
        {
            if (inHighPermRegion_(globalPos))
                phi = 0.25;
            else
                phi = 0.2;
        }
        else if (dim == 2)
        {
            phi=0.4;
        }
        else{
            throw std::runtime_error("NotImplemented: Intrinsic permeability is only defined for dim=2 or 3");
        }
        return phi;
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
        if (onDirichletBoundary_(globalPos))
        {
            // Free flow boundary. we assume incompressible fluids
            Scalar pw, Sw;
            Scalar T = this->temperature(context, spaceIdx, timeIdx);
            Scalar densityW = WettingPhase::density(T, /*pressure=*/Scalar(1));
            Scalar densityN = NonwettingPhase::density(T, /*pressure=*/Scalar(1));

            if (onInjectionBoundary_(globalPos))
            {
                pw = 4.0;
                Sw = 1e-2;
            }
            else
            {
                pw = 1.0;
                Sw = 0.0;
            }

            // specify a full fluid state using pw and Sw
            const MaterialLawParams &matParams = this->materialLawParams(context, spaceIdx, timeIdx);

            Opm::ImmiscibleFluidState<Scalar, FluidSystem,
                                      /*storeEnthalpy=*/false>
                fs;
            fs.setSaturation(wettingPhaseIdx, Sw);
            fs.setSaturation(nonWettingPhaseIdx, 1 - Sw);
            fs.setTemperature(T);

            Scalar pC[numPhases];
            MaterialLaw::capillaryPressures(pC, matParams, fs);
            fs.setPressure(wettingPhaseIdx, pw);
            fs.setPressure(nonWettingPhaseIdx, pw + pC[nonWettingPhaseIdx] - pC[wettingPhaseIdx]);

            fs.setDensity(wettingPhaseIdx, densityW);
            fs.setDensity(nonWettingPhaseIdx, densityN);

            fs.setViscosity(wettingPhaseIdx, WettingPhase::viscosity(T, /*pressure=*/Scalar(1e5)));
            fs.setViscosity(nonWettingPhaseIdx, NonwettingPhase::viscosity(T, /*pressure=*/Scalar(1e5)));

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else
        {
            // no flow boundary
            values.setNoFlow();
        }
    }
    template <class PrimaryVariables, class Context>
    void initial(PrimaryVariables &values, const Context &context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

        Scalar pw = 0.0;

        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
        fs.setPressure(wettingPhaseIdx, pw);

        Scalar Sw = 0.0;
        fs.setSaturation(wettingPhaseIdx, Sw);
        fs.setSaturation(nonWettingPhaseIdx, 1 - Sw);

        fs.setTemperature(273);

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updatePhase(fs, wettingPhaseIdx);

        // calculate the capillary pressure
        const MaterialLawParams &matParams = this->materialLawParams(context, spaceIdx, timeIdx);
        Scalar pC[numPhases];
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        // make a full fluid state
        fs.setPressure(wettingPhaseIdx, pw);
        fs.setPressure(nonWettingPhaseIdx, pw + (pC[wettingPhaseIdx] - pC[nonWettingPhaseIdx]));

        // assign the primary variables
        values.assignNaive(fs);
    }

    std::string name() const
    {
        return "benchmark3d_1" + ParentType::name();
    }

    static void registerParameters()
    {
        ParentType::registerParameters();
    }

protected:
    bool onDirichletBoundary_(GlobalPosition pos) const
    {
        return onInjectionBoundary_(pos) || onOutflowBoundary_(pos);
    }

    bool onInjectionBoundary_(GlobalPosition pos) const
    {
        return this->onWestBoundary_(pos) && pos[2] > 90;
    }
    bool onOutflowBoundary_(GlobalPosition pos) const
    {
        return this->onSouthBoundary_(pos) && pos[2] < 10;
    }
    bool onWestBoundary_(const GlobalPosition &pos) const
    {
        return pos[0] < eps_;
    }

    bool onEastBoundary_(const GlobalPosition &pos) const
    {
        return pos[0] > 100 - eps_;
    }

    bool onSouthBoundary_(const GlobalPosition &pos) const
    {
        return pos[1] < eps_;
    }

    bool onNorthBoundary_(const GlobalPosition &pos) const
    {
        return pos[1] > 100 - eps_;
    }

    bool inHighPermRegion_(const GlobalPosition &pos) const{
        return pos[2] < 10;
    }

    Scalar eps_{1e-9};
};

} // namespace Opm

#endif
