#include "TestFunction2.h"

namespace mg
{

Precision mg::TestFunction2::apply(Precision x, Precision y) const
{
    return c_*cos(4*Pi*x)*cos(4*Pi*y);
}

}
