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

BEGIN_PROPERTIES

NEW_TYPE_TAG(Darcy1PBaseProblem, INHERITS_FROM(ImmiscibleSinglePhaseModel));

NEW_PROP_TAG(NameSurfix);
NEW_PROP_TAG(DomainDim);
NEW_PROP_TAG(WorldDim);
NEW_PROP_TAG(GridDim);

SET_STRING_PROP(Darcy1PBaseProblem, NameSurfix, "");
SET_INT_PROP(Darcy1PBaseProblem, WorldDim, 2);
SET_INT_PROP(Darcy1PBaseProblem, DomainDim, 2);

SET_PROP(Darcy1PBaseProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::Unit<Scalar>> type;
};

SET_TYPE_PROP(Darcy1PBaseProblem, Problem,
              Opm::Darcy1PProblem<TypeTag>);

// Define grid type
SET_TYPE_PROP(Darcy1PBaseProblem, Grid, Dune::PolyhedralGrid<GET_PROP_VALUE(TypeTag, GridDim), GET_PROP_VALUE(TypeTag, WorldDim)>);

// Dissable gravity
SET_BOOL_PROP(Darcy1PBaseProblem, EnableGravity, false);

// The default for the end time of the simulation
SET_SCALAR_PROP(Darcy1PBaseProblem, EndTime, 1);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(Darcy1PBaseProblem, InitialTimeStepSize, 1);

// Use the conjugated gradient linear solver with the default preconditioner (i.e.,
// ILU-0) from dune-istl
SET_TAG_PROP(Darcy1PBaseProblem, LinearSolverSplice, ParallelIstlLinearSolver);
SET_TYPE_PROP(Darcy1PBaseProblem, LinearSolverWrapper,
              Opm::Linear::SolverWrapperConjugatedGradients<TypeTag>);
// For debugging
SET_BOOL_PROP(Darcy1PBaseProblem, VtkWriteSaturations, true);
SET_BOOL_PROP(Darcy1PBaseProblem, VtkWriteMobilities, true);
SET_BOOL_PROP(Darcy1PBaseProblem, VtkWriteRelativePermeabilities, true);

END_PROPERTIES

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
class Darcy1PProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum
    {
        numPhases = FluidSystem::numPhases,

        // Grid and world dimension
        dim = GET_PROP_VALUE(TypeTag, DomainDim), 
        dimWorld = GridView::dimensionworld,

        // indices of the primary variables
        pressure0Idx = Indices::pressure0Idx
    };

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

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
        oss << "darcy1PProblem_" << GET_PROP_VALUE(TypeTag, NameSurfix);
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
    template <class Context>
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
