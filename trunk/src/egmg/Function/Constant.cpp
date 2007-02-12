#include "Constant.h"

namespace mg
{

Constant::Constant(Precision constant): constant_(constant)
{}

Precision mg::Constant::apply(Precision, Precision) const
{
	return constant_;
}

}
