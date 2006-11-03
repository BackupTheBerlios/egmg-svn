#include "Function.h"

namespace mg
{

Precision Function::operator() (Precision x, Precision y) const
{
	return apply(x, y);
}

Function::~Function()
{}

}
