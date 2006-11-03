#include "TestRightSide1.h"

namespace mg
{

Precision mg::TestRightSide1::apply(Precision x, Precision y) const
{
    return x*x*x*x*y*y*y;
}

}
