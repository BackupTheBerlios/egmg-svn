/** \file putrhs.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the function putrhs
 */
#include "putrhs.h"

namespace mg
{
    void putrhs(
        NumericArray& fv,
        const Index nx,
        const Index ny,
        const function2D f)
    {
        const Precision hx=1.0/nx;
        const Precision hy=1.0/ny;
        for (Index sy=0; sy<=ny; sy++)
            for (Index sx=0; sx<=nx; sx++)
                fv[sy*(nx+1)+sx]=f(sx*hx,sy*hy);
    }
}
