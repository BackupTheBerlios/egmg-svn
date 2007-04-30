#ifndef RANDOM_H_
#define RANDOM_H_

#include "Function.h"

namespace mg
{

class Random: public Function
{

public:
	Random();
	Random(Index seed);

protected:
	virtual Precision apply(Precision x, Precision y) const;
};

}

#endif /*RANDOM_H_*/
