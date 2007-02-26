#ifndef CONSTANT_H_
#define CONSTANT_H_

#include "Function.h"

namespace mg
{

class Constant: public Function
{
const Precision constant_;

public:
	Constant();
	Constant(Precision constant);

protected:
	virtual Precision apply(Precision x, Precision y) const;
};

}

#endif
