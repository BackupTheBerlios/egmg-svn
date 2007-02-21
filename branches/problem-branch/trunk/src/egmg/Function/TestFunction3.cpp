#include "TestFunction3.h"

namespace mg
{

Precision mg::TestFunction3::apply(Precision x, Precision y) const
{
    return c_*sin(2*Pi*x)*sin(2*Pi*y);
}

}
