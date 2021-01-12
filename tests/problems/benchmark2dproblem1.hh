
#ifndef OPM_BENCHMARK_1_PROBLEM_HH
#define OPM_BENCHMARK_1_PROBLEM_HH

#include "darcy1PProblem.hh"

#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/immiscible/immisciblemodel.hh>

namespace Opm
{
template <class TypeTag>
class Benchmark1Problem;
}

BEGIN_PROPERTIES
NEW_TYPE_TAG(Benchmark1Problem, INHERITS_FROM(Darcy1PBaseProblem));

SET_TYPE_PROP(Benchmark1Problem, Problem,
              Opm::Benchmark1Problem<TypeTag>);

END_PROPERTIES

namespace Opm
{

template <class TypeTag>
class Benchmark1Problem : public Opm::Darcy1PProblem<TypeTag>
{
    typedef Opm::Darcy1PProblem<TypeTag> ParentType;

public:
    using ParentType::ParentType;

    std::string name() const
    {
        return "benchmark1" + ParentType::name();
    }
};
} // namespace Opm

#endif