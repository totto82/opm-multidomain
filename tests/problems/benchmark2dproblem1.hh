
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

namespace Opm::Properties {

// Create new type tags
namespace TTag {
struct Benchmark1Problem { using InheritsFrom = std::tuple<Darcy1PBaseProblem>; };
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::Benchmark1Problem> { using type = Opm::Benchmark1Problem<TypeTag>; };
} // end namespace Opm::Properties


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