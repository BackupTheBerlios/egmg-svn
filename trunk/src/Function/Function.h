#ifndef FUNCTION_H_
#define FUNCTION_H_

#include "../general/parameters.h"

namespace mg
{

class Function
{
	Precision operator() (Precision x, Precision y) const;
protected:
	virtual Precision apply(Precision x, Precision y) const =0;
};

}

#endif
