#include "TestFunction2.h"

namespace mg
{

Precision mg::TestFunction2::apply(Precision x, Precision y) const
{
    return c_*cos(2*Pi*x)*y*y;
}

}
