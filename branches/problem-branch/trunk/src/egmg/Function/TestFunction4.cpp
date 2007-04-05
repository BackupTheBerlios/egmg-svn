#include "TestFunction4.h"

namespace mg
{

Precision mg::TestFunction4::apply(Precision x, Precision y) const
{
    return -cos(2*Pi*x)*(2-4*Pi*Pi*y*y-(x*x+y*y)*y*y);
}

}
