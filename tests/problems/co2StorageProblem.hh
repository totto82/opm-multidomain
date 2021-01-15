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
 * \copydoc Opm::Co2StorageProblem
 */
#ifndef CO2_STORAGE_PROBLEM_HH
#define CO2_STORAGE_PROBLEM_HH

#include <opm/models/discretization/common/fvbaseadlocallinearizer.hh>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/immiscible/immiscibleproperties.hh>
#include <opm/models/io/vtkmultiphasemodule.hh>
#include <opm/models/io/unstructuredgridvanguard.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/multidomain/multidomainproperties.hh>

#include <opm/material/common/Unused.hpp>
#include <opm/material/components/Lnapl.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/Unit.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/fluidsystems/TwoPhaseImmiscibleFluidSystem.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include "cmath"
#include <iostream>
#include <sstream>
#include <string>

namespace Opm {
template <class TypeTag>
class Co2StorageProblem;
}

namespace Opm {
template <class Scalar>
class Water : public Opm::Unit<Scalar> {
public:
    template <class Evaluation>
    static Evaluation liquidViscosity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    {
        return 8 * std::pow(10.0, -4.0);
    }
    static Scalar molarMass()
    {
        return 1.0;
    }

    template <class Evaluation>
    static Evaluation liquidDensity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    {
        return 1012;
    }
};
template <class Scalar>
class CO2 {
public:
    template <class Evaluation>
    static Evaluation liquidViscosity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    {
        return 5 * std::pow(10.0, -5.0);
    }
    static Scalar molarMass()
    {
        return 1.0;
    }
    template <class Evaluation>
    static Evaluation liquidDensity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    {
        return 714;
    }
};
} //namespace opm


namespace Opm::Properties {

//! The generic type tag for problems using the dual porosity coupler
namespace TTag {
struct Co2StorageProblem {};
} // end namespace TTag


template<class TypeTag, class MyTypeTag>
struct DomainDim { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct NameSurfix { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct FractureDomain { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct SymmetricDomain { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct DomainXmax { using type = UndefinedProperty; };

template<class TypeTag>
struct FractureDomain<TypeTag, TTag::Co2StorageProblem> { static constexpr bool value = false; };
template<class TypeTag>
struct SymmetricDomain<TypeTag, TTag::Co2StorageProblem> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Co2StorageProblem> { using type = Opm::Co2StorageProblem<TypeTag>; };

// Set the wetting phase
template<class TypeTag>
struct WettingPhase<TypeTag, TTag::Co2StorageProblem> 
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Opm::LiquidPhase<Scalar, Water<Scalar>>;
};

// Set the non-wetting phase
template<class TypeTag>
struct NonwettingPhase<TypeTag, TTag::Co2StorageProblem> 
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Opm::LiquidPhase<Scalar, CO2<Scalar>>;
};

// Set the material Law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::Co2StorageProblem>
{
private:
    using FluidSystem = GetPropType<TypeTag, FluidSystem>;
    using Scalar = GetPropType<TypeTag, Scalar>;

    enum { wettingPhaseIdx = FluidSystem::wettingPhaseIdx };
    enum { nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx };

    using Traits =  Opm::TwoPhaseMaterialTraits<Scalar,
        /*wettingPhaseIdx=*/FluidSystem::wettingPhaseIdx,
        /*nonWettingPhaseIdx=*/FluidSystem::nonWettingPhaseIdx>;

    // define the material law which is parameterized by effective
    // saturations
    using EffectiveLaw = Opm::RegularizedBrooksCorey<Traits>;

public:
    // define the material law parameterized by absolute saturations
    using type = Opm::EffToAbsLaw<EffectiveLaw>;
};

template<class TypeTag>
struct NewtonWriteConvergence<TypeTag, TTag::Co2StorageProblem> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::Co2StorageProblem> { static constexpr bool value = false; };

template<class TypeTag>
struct EndTime<TypeTag, TTag::Co2StorageProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 2000 * 60 * 60 * 24;
};

template<class TypeTag>
struct MaxTimeStepSize<TypeTag, TTag::Co2StorageProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value =  20 * 60 * 60 * 24;
};

template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::Co2StorageProblem>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value =  60 * 60;
};


// By default, include the intrinsic permeability tensor to the VTK output files
template<class TypeTag>
struct VtkWriteIntrinsicPermeabilities<TypeTag, TTag::Co2StorageProblem> { static constexpr bool value = true; };

// enable the storage cache by default for this problem
template<class TypeTag>
struct EnableStorageCache<TypeTag, TTag::Co2StorageProblem> { static constexpr bool value = true; };

// enable the cache for intensive quantities by default for this problem
template<class TypeTag>
struct EnableIntensiveQuantityCache<TypeTag, TTag::Co2StorageProblem> { static constexpr bool value = true; };

}

namespace Opm {

/*!
 * \ingroup TestProblems
 *
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability. Note that
 * this problem is discretized using only two dimensions, so from the
 * point of view of the model, the depth of the domain is implicitly
 * assumed to be 1 m everywhere.
 *
 * On the top and the bottom of the domain no-flow boundary conditions
 * are used, while free-flow conditions apply on the left and right
 * boundaries; DNAPL is injected at the top boundary from 3m to 4m at
 * a rate of 0.04 kg/(s m^2).
 *
 * At the boundary on the left, a free-flow condition using the
 * hydrostatic pressure scaled by a factor of 1.125 is imposed, while
 * on the right, it is just the hydrostatic pressure. The DNAPL
 * saturation on both sides is zero.
 */
template <class TypeTag>
class Co2StorageProblem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using WettingPhase = GetPropType<TypeTag, Properties::WettingPhase>;
    using NonwettingPhase = GetPropType<TypeTag, Properties::NonwettingPhase>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    enum {
        isFracture = getPropValue<TypeTag, Properties::FractureDomain>(),
        isSymmetric = getPropValue<TypeTag, Properties::SymmetricDomain>()
    };
    
    enum
    {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        wettingPhaseIdx = FluidSystem::wettingPhaseIdx,
        nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx,

        // equation indices
        contiNEqIdx = Indices::conti0EqIdx + wettingPhaseIdx,

        // Grid and world dimension
        dim = getPropValue<TypeTag, Properties::DomainDim>(),
        dimWorld = GridView::dimensionworld,
    };


    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    Co2StorageProblem(Simulator& simulator)
        : ParentType(simulator)
    {
        init();
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void init()
    {
        //ParentType::finishInit();

        eps_ = 3e-6;
        FluidSystem::init();

        temperature_ = 273.15 + 20; // -> 20°C

        if (isFracture) {
            materialParams_.setEntryPressure(0.0 * 1000);
            // materialParams_.setResidualSaturation(0, 0.0);
            // materialParams_.setResidualSaturation(1, 0.0);
        } else {
            materialParams_.setEntryPressure(10 * 1000);
            materialParams_.setResidualSaturation(0, 0.42);
            materialParams_.setResidualSaturation(1, 0.0);
        }
        materialParams_.finalize();

        double milliDarcy = 1e-3 * 9.869233 * 1e-13;
        double Kx = 4053.0 * milliDarcy;//405.0 * milliDarcy;
        double Kz = 4053.0 * milliDarcy;

        if (isFracture) {
            K_ = this->toDimMatrix_(Kx);
            K_[dimWorld - 1][dimWorld - 1] = Kz;
        } else
            K_ = this->toDimMatrix_(21.0 * milliDarcy);

        if (dimWorld != 2) {
            throw std::logic_error("Not implemented for world dimension different than 2");
        } else {
            this->gravity_ = 0;
            this->gravity_[1] = -9.81;
        }
        numWellElements_ = 0;
        for (const auto& e : elements(this->simulator().gridView())) {
            if (isWell_(e.geometry().center()))
                numWellElements_++;
        }
    }
    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();
        EWOMS_REGISTER_PARAM(TypeTag, std::string, NameSurfix, "Id to append to grid name");
        EWOMS_REGISTER_PARAM(TypeTag, int, DomainDim, "Dimension of domain");
        EWOMS_REGISTER_PARAM(TypeTag, bool, FractureDomain, "Is the domain a fracture domain? If false, then matrix");
        EWOMS_REGISTER_PARAM(TypeTag, bool, SymmetricDomain, "If true, assume domain (and flow) is symmetric along line x=500 m");
    }

    /*!
     * \copydoc FvBaseProblem::briefDescription
     */
    static std::string briefDescription()
    {
        std::string thermal = "isothermal";
        bool enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>();
        if (enableEnergy)
            thermal = "non-isothermal";

        std::string deriv = "finite difference";
        using LLS = GetPropType<TypeTag, Properties::LocalLinearizerSplice>;
        bool useAutoDiff = std::is_same<LLS, Properties::TTag::AutoDiffLocalLinearizer>::value;
        if (useAutoDiff)
            deriv = "automatic differentiation";

        std::string disc = "vertex centered finite volume";
        using D = GetPropType<TypeTag, Properties::Discretization>;
        bool useEcfv = std::is_same<D, Opm::EcfvDiscretization<TypeTag>>::value;
        if (useEcfv)
            disc = "element centered finite volume";

        return std::string("") + "Ground remediation problem where CO2 infiltrates " + "an aquifer with an embedded low-permability lens. " + "This is the binary for the " + thermal + " variant using " + deriv + "and the " + disc + " discretization";
    }

    /*!
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context, unsigned spaceIdx,
        unsigned timeIdx) const
    {
        const GlobalPosition& globalPos = context.pos(spaceIdx, timeIdx);
        return K_;
    }
    template <class Context>
    Scalar extrusionFactor(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {

        Scalar aperture = 0.00962;
        if (isSymmetric && onSymmetryLine_(context.pos(spaceIdx, timeIdx)))
            aperture /= 2.0;
        return std::pow(aperture, dimWorld - dim);
    }

    Scalar extrusionFactor() const
    {
        throw std::runtime_error("need position to calculate extrusion factor");
    }
    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context OPM_UNUSED,
        unsigned spaceIdx OPM_UNUSED,
        unsigned timeIdx OPM_UNUSED) const
    {
        if (isFracture){
            if (dim<dimWorld)
                return 1.0;
            else
                return 0.01;
        }
        else
            return 0.19;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
        unsigned spaceIdx, unsigned timeIdx) const
    {
        return materialParams_;
    }
    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    const MaterialLawParams& materialLawParams() const
    {
        return materialParams_;
    }
    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
        unsigned spaceIdx OPM_UNUSED,
        unsigned timeIdx OPM_UNUSED) const
    {
        return temperature_;
    }

    //! \}

    /*!
     * \name Auxiliary methods
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "gravityDominatedFlow" << Model::name()
            << "_" << getPropValue<TypeTag, Properties::NameSurfix>();
        return oss.str();
    }

    /*!
     * \copydoc FvBaseProblem::beginTimeStep
     */
    void beginTimeStep()
    {
    }

    /*!
     * \copydoc FvBaseProblem::beginIteration
     */
    void beginIteration()
    {
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
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: " << storage << std::endl
                      << std::flush;
        }
#endif // NDEBUG
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
        const Context& context,
        unsigned spaceIdx,
        unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        const Scalar time = this->simulator().time();
        if (onLeftBoundary_(pos) || onRightBoundary_(pos)) {
            // free flow boundary. we assume incompressible fluids
            Scalar T = temperature(context, spaceIdx, timeIdx);
            Scalar pw, Sw;

            // set wetting phase pressure and saturation
            Scalar y = pos[1];
            Scalar depth = this->boundingBoxMax()[1] - y;
            pw = depth * WettingPhase::density(T, 1e5) * this->gravity().two_norm();
            //pw = 1e5;
            Sw = 1.0;

            // specify a full fluid state using pw and Sw
            const MaterialLawParams& matParams = this->materialLawParams(context, spaceIdx, timeIdx);

            Opm::ImmiscibleFluidState<Scalar, FluidSystem, /*storeEnthalpy=*/false> fs;

            // Upstream weight saturation based on wetting phase pressure
            const auto& insideIntQuants = context.intensiveQuantities(spaceIdx, timeIdx);
            auto pInside = Opm::getValue(insideIntQuants.fluidState().pressure(wettingPhaseIdx));

            if (pw > pInside) {
                fs.setSaturation(wettingPhaseIdx, Sw);
                fs.setSaturation(nonWettingPhaseIdx, 1 - Sw);
            } else {
                const auto SwInside = Opm::getValue(insideIntQuants.fluidState().saturation(wettingPhaseIdx));
                fs.setSaturation(wettingPhaseIdx, SwInside);
                fs.setSaturation(nonWettingPhaseIdx, 1 - SwInside);
            }
            // fs.setSaturation(wettingPhaseIdx, Sw);
            // fs.setSaturation(nonWettingPhaseIdx, 1 - Sw);

            fs.setTemperature(T);

            Scalar pC[numPhases];
            MaterialLaw::capillaryPressures(pC, matParams, fs);
            fs.setPressure(wettingPhaseIdx, pw);
            fs.setPressure(nonWettingPhaseIdx, pw + (pC[nonWettingPhaseIdx] - pC[wettingPhaseIdx]));

            Scalar densityW = WettingPhase::density(temperature_, fs.pressure(wettingPhaseIdx));
            Scalar densityN = NonwettingPhase::density(temperature_, fs.pressure(nonWettingPhaseIdx));

            fs.setDensity(wettingPhaseIdx, densityW);
            fs.setDensity(nonWettingPhaseIdx, densityN);

            fs.setViscosity(wettingPhaseIdx, WettingPhase::viscosity(temperature_, fs.pressure(wettingPhaseIdx)));
            fs.setViscosity(nonWettingPhaseIdx, NonwettingPhase::viscosity(temperature_, fs.pressure(nonWettingPhaseIdx)));

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        } else {
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
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        const Scalar y = pos[1];
        Scalar depth = this->boundingBoxMax()[1] - y;

        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
        fs.setPressure(wettingPhaseIdx, /*pressure=*/1e5);

        Scalar Sw = 1.0;
        // if (pos[1] > 50 +1e-6){
        //     Sw = 1.0;
        // }
        // if (isFracture){
        //     Sw = 0.0;
        // }
        fs.setSaturation(wettingPhaseIdx, Sw);
        fs.setSaturation(nonWettingPhaseIdx, 1 - Sw);

        fs.setTemperature(temperature_);

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updatePhase(fs, wettingPhaseIdx);
        Scalar densityW = FluidSystem::density(fs, paramCache, wettingPhaseIdx);

        // hydrostatic pressure (assuming incompressibility)
        Scalar pw = densityW * this->gravity().two_norm() * depth;

        // calculate the capillary pressure
        const MaterialLawParams& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
        Scalar pC[numPhases];
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        // make a full fluid state
        fs.setPressure(wettingPhaseIdx, pw);
        fs.setPressure(nonWettingPhaseIdx, pw + (pC[nonWettingPhaseIdx] - pC[wettingPhaseIdx]));

        // assign the primary variables
        values.assignNaive(fs);
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
        const Context& context,
        unsigned spaceIdx,
        unsigned timeIdx) const
    {

        rate = Scalar(0.0);
        if (numWellElements_ > 0) {
            const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
            auto rho = context.intensiveQuantities(spaceIdx, timeIdx).fluidState().density(1);
            const Scalar time = this->simulator().time();
            auto volume = context.dofVolume(spaceIdx, timeIdx) * extrusionFactor(context, spaceIdx, timeIdx);
            if (time < 1000 * day_ && isWell_(pos))
                rate[1] = 1502.41 * 1000/ (1000 * day_) / numWellElements_ / volume;
        }
        if (isSymmetric)
            rate /= 2.0;
    }

    //! \}
protected:
    Dune::FieldVector<double, 2> boundingBoxMax() const
    {
        Dune::FieldVector<double, 2> max { 1000, 300 };
        return max;
    }
    Dune::FieldVector<double, 2> boundingBoxMin() const
    {
        Dune::FieldVector<double, 2> min { 0, 0};
        return min;
    }

    bool onLeftBoundary_(const GlobalPosition& pos) const
    {
        return pos[0] < this->boundingBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition& pos) const
    {
        return pos[0] > this->boundingBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    {
        return pos[1] < this->boundingBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    {
        return pos[1] > this->boundingBoxMax()[1] - eps_;
    }
    bool onSymmetryLine_(const GlobalPosition& pos) const
    {
        return (std::abs(pos[0] - 500.0) < eps_) && isFracture && isSymmetric;
    }

    bool isWell_(const GlobalPosition& pos) const
    {
        return (std::abs(pos[0] - 500.0) < eps_) && isFracture && (dim > 0);
    }

private:
    DimMatrix K_;
    MaterialLawParams materialParams_;

    Scalar temperature_;
    Scalar eps_;
    Scalar year_ = 60 * 60 * 24 * 365;
    Scalar day_ = 60 * 60 * 24;
    unsigned numWellElements_;
};

} // namespace Opm

#endif
