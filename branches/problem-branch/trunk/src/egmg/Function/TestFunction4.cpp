#include "TestFunction4.h"

namespace mg
{

Precision mg::TestFunction4::apply(Precision x, Precision y) const
{
    return 16.0*Pi*Pi*sin(2*Pi*x+1)+60*Pi*Pi*cos(2*Pi*(2*x+y));
}

}
