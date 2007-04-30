#include "Random.h"

#include <cstdlib>

namespace mg
{

Random::Random()
{
    std::srand(std::time(0));
}

Random::Random(Index seed)
{
    std::srand(seed);
}

Precision Random::apply(Precision, Precision) const
{
	return static_cast<Precision>(std::rand())/RAND_MAX;
}

}
