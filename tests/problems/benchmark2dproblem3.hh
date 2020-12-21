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
#ifndef OPM_BENCHMARK_3_PROBLEM_HH
#define OPM_BENCHMARK_3_PROBLEM_HH

#include "darcy1PProblem.hh"

namespace Opm
{
template <class TypeTag>
class Benchmark3Problem;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(Benchmark3Problem, INHERITS_FROM(Darcy1PBaseProblem));
NEW_PROP_TAG(LeftRight);
SET_BOOL_PROP(Benchmark3Problem, LeftRight, 1);

SET_TYPE_PROP(Benchmark3Problem, Problem,
              Opm::Benchmark3Problem<TypeTag>);

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
class Benchmark3Problem : public Opm::Darcy1PProblem<TypeTag>
{
    typedef Opm::Darcy1PProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;

    enum
    {
        numPhases = FluidSystem::numPhases,
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

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
        double Kh = 1e4;
        double Kl = 1e-4;
        if (dim == 2)
            K = 1;
        else if (dim == 1)
        {
            if (std::abs((globalPos[0] - 0.15) / (0.4 - 0.15) - (globalPos[1] - 0.9167) / (0.5 - 0.9167)) < tol)
                K = Kl;
            else if (std::abs((globalPos[0] - 0.65) / (0.849723 - 0.65) - (globalPos[1] - 0.8333) / (0.167625 - 0.8333)) < tol)
                K = Kl;
            else K = Kh;
        }
        else if (dim == 0)
        {
            if (std::abs((globalPos[0] - 0.152174)) < tol && (globalPos[1] - 0.203478) < tol)
                K = Kh;
            else if (std::abs((globalPos[0] - 0.37326)) < tol && (globalPos[1] - 0.958111) < tol)
                K = Kh;
            else
                K = Kl;
        }
        else
            throw std::runtime_error("Benchmark3problem can only handle 0, 1 and 2 dimensions");

        return this->toDimMatrix_(K);
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
            Scalar pressure;
            Scalar T = this->temperature(context, spaceIdx, timeIdx);

            if (onInjectionBoundary_(globalPos))
                pressure = 4.0;
            else
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
        else
        {
            // no flow boundary
            values.setNoFlow();
        }
    }

    std::string name() const
    {
        std::string leftRight;

        if (GET_PROP_VALUE(TypeTag, LeftRight))
            leftRight = "left_to_right_";
        else
            leftRight = "top_to_bottom_";
        return "benchmark3_" + leftRight + ParentType::name();
    }

    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, LeftRight,
                             "Flow from left -> right if true. Else, top -> bottom");
    }

protected:
    bool onDirichletBoundary_(GlobalPosition pos) const
    {
        if (GET_PROP_VALUE(TypeTag, LeftRight))
            return this->onLeftBoundary_(pos) || this->onRightBoundary_(pos);
        else
            return this->onBottomBoundary_(pos) || this->onTopBoundary_(pos);
    }
    bool onInjectionBoundary_(GlobalPosition pos) const
    {
        if (GET_PROP_VALUE(TypeTag, LeftRight))
            return this->onLeftBoundary_(pos);
        else
            return this->onTopBoundary_(pos);
    }
};

} // namespace Opm

#endif
