/** \file putrhs.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the function putrhs
 */
#include "putrhs.h"

namespace mg
{
    void putrhs(
        NumericArray& fv,
        const size_t nx,
        const size_t ny,
        const function2D f)
    {
        const Precision hx=1.0/nx;
        const Precision hy=1.0/ny;
        for (size_t sy=0; sy<=ny; sy++)
            for (size_t sx=0; sx<=nx; sx++)
                fv[sy*(nx+1)+sx]=f(sx*hx,sy*hy);
    }
}
