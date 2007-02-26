#include "Constant.h"

namespace mg
{

Constant::Constant() : constant_(1.0)
{}

Constant::Constant(Precision constant): constant_(constant)
{}

Precision Constant::apply(Precision, Precision) const
{
	return constant_;
}

}
