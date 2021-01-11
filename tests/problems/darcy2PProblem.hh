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
 * \copydoc Opm::Darcy2pProblem
 */
#ifndef EWOMS_DARCY_2p_PROBLEM_HH
#define EWOMS_DARCY_2p_PROBLEM_HH

#include <opm/models/immiscible/immiscibleproperties.hh>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/io/multidomainvanguard.hh>

#include <opm/material/components/Unit.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidsystems/TwoPhaseImmiscibleFluidSystem.hpp>

#include <dune/grid/yaspgrid.hh>

namespace Opm
{
template <class TypeTag>
class Darcy2pProblem;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(Darcy2pBaseProblem, INHERITS_FROM(MultiDomainVanguard));

NEW_PROP_TAG(NameSurfix);
NEW_PROP_TAG(DomainDim);

// Set the problem property
SET_TYPE_PROP(Darcy2pBaseProblem, Problem, Opm::Darcy2pProblem<TypeTag>);

// Set the wetting phase
SET_PROP(Darcy2pBaseProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::Unit<Scalar>> type;
};

// Set the non-wetting phase
SET_PROP(Darcy2pBaseProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::Unit<Scalar>> type;
};

// Set the material Law
SET_PROP(Darcy2pBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum
    {
        wettingPhaseIdx = FluidSystem::wettingPhaseIdx
    };
    enum
    {
        nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx
    };
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::wettingPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::nonWettingPhaseIdx>
        Traits;

    // define the material law which is parameterized by effective
    // saturations
    typedef Opm::LinearMaterial<Traits> EffectiveLaw;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::EffToAbsLaw<EffectiveLaw> type;
};

// Dissable gravity
SET_BOOL_PROP(Darcy2pBaseProblem, EnableGravity, false);

// The default for the initial time step size of the simulation. This might be overridden
// in a multidomain simulation.
SET_SCALAR_PROP(Darcy2pBaseProblem, InitialTimeStepSize, 0.01);

// By default, include the intrinsic permeability tensor to the VTK output files
SET_BOOL_PROP(Darcy2pBaseProblem, VtkWriteIntrinsicPermeabilities, true);

// enable the storage cache by default for this problem
SET_BOOL_PROP(Darcy2pBaseProblem, EnableStorageCache, true);

// enable the cache for intensive quantities by default for this problem
SET_BOOL_PROP(Darcy2pBaseProblem, EnableIntensiveQuantityCache, true);


SET_STRING_PROP(Darcy2pBaseProblem, NameSurfix, "");
SET_INT_PROP(Darcy2pBaseProblem, WorldDim, 3);
SET_INT_PROP(Darcy2pBaseProblem, DomainDim, 3);

END_PROPERTIES

namespace Opm
{

/*!
 * \ingroup TestProblems
 *
 * \brief An incompressible two phase problem.
 */
template <class TypeTag>
class Darcy2pProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    enum
    {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        wettingPhaseIdx = FluidSystem::wettingPhaseIdx,
        nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx,

        // equation indices
        contiNEqIdx = Indices::conti0EqIdx + nonWettingPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    Darcy2pProblem(Simulator &simulator)
        : ParentType(simulator)
    {
        initiated_ = true;
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 3e-6;
        FluidSystem::init();

        temperature_ = 273.15 + 20; // -> 20°C

        // residual saturations
        materialParams_.setResidualSaturation(wettingPhaseIdx, 0.0);
        materialParams_.setResidualSaturation(nonWettingPhaseIdx, 0.0);

        materialParams_.finalize();

        K_ = this->toDimMatrix_(1.0);

        if (dimWorld == 3)
        {
            this->gravity_ = 0;
            this->gravity_[1] = 9.8;
        }
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();
    }
    bool recycleFirstIterationStorage() const
    {
        return false;
    }
    /*!
     * \copydoc FvBaseProblem::briefDescription
     */
    static std::string briefDescription()
    {
        std::string thermal = "isothermal";
        bool enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy);
        if (enableEnergy)
            thermal = "non-isothermal";

        std::string disc = "vertex centered finite volume";
        typedef typename GET_PROP_TYPE(TypeTag, Discretization) D;
        bool useEcfv = std::is_same<D, Opm::EcfvDiscretization<TypeTag>>::value;
        if (useEcfv)
            disc = "element centered finite volume";

        return std::string("") +
               "Ground remediation problem where a dense oil infiltrates " +
               "an aquifer. " +
               "This is the binary for the " + thermal + " variant using the" + disc + " discretization";
    }

    /*!
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        return K_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context OPM_UNUSED,
                    unsigned spaceIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED) const
    {
        return 1.0;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams &materialLawParams(const Context &context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        return materialParams_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context &context OPM_UNUSED,
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
        oss << "darcy2PProblem_" << GET_PROP_VALUE(TypeTag, NameSurfix);
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
        if (this->gridView().comm().rank() == 0)
        {
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
    void boundary(BoundaryRateVector &values,
                  const Context &context,
                  unsigned spaceIdx,
                  unsigned timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos) || onRightBoundary_(pos))
        {
            // free flow boundary. we assume incompressible fluids
            Scalar densityW = WettingPhase::density(temperature_, /*pressure=*/Scalar(1));
            Scalar densityN = NonwettingPhase::density(temperature_, /*pressure=*/Scalar(1));

            Scalar T = temperature(context, spaceIdx, timeIdx);
            Scalar pw, Sw;

            // set wetting phase pressure and saturation
            if (onLeftBoundary_(pos))
            {
                // hydrostatic pressure + 1
                pw = 700.0;
                Sw = 0.0;
            }
            else
            {
                // hydrostatic pressure
                pw = 0.0;
                Sw = 1.0;
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

            fs.setViscosity(wettingPhaseIdx, WettingPhase::viscosity(temperature_, /*pressure=*/Scalar(1e5)));
            fs.setViscosity(nonWettingPhaseIdx, NonwettingPhase::viscosity(temperature_, /*pressure=*/Scalar(1e5)));

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
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
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

        Scalar pw = 700 * (1 - 1.0 / 700 * pos[0]);

        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
        fs.setPressure(wettingPhaseIdx, pw);

        Scalar Sw = 1.0;
        fs.setSaturation(wettingPhaseIdx, Sw);
        fs.setSaturation(nonWettingPhaseIdx, 1 - Sw);

        fs.setTemperature(temperature_);

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

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
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
    template <class Context>
    Scalar extrusionFactor(const Context &context OPM_UNUSED,
                           unsigned spaceIdx OPM_UNUSED,
                           unsigned timeIdx OPM_UNUSED) const
    {
        if (dim == dimWorld - 1)
            return 1e-4; // fracture apperture
        else
            return 1.0;
    }

private:
    bool onLeftBoundary_(const GlobalPosition &pos) const
    {
        return pos[0] < 0 + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &pos) const
    {
        return pos[0] > 1 - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &pos) const
    {
        return pos[1] < this->boundingBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &pos) const
    {
        return pos[1] > this->boundingBoxMax()[1] - eps_;
    }

    DimMatrix K_;
    MaterialLawParams materialParams_;

    Scalar temperature_;
    Scalar eps_;
    bool initiated_{false};
};

} // namespace Opm

#endif
