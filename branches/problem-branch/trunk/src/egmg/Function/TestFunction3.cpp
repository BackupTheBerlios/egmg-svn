#include "TestFunction3.h"

namespace mg
{

Precision mg::TestFunction3::apply(Precision x, Precision y) const
{
    return x*x+y*y;
}

}
