#include "TestFunction1.h"

namespace mg
{

Precision mg::TestFunction1::apply(Precision x, Precision y) const
{
    return -6*(x*x*y*(x*x+2*y*y));
}

}
