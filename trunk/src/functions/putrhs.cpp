/** \file putrhs.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the function putrhs
 */
#include "putrhs.h"

namespace mg
{
    void putrhs(
        std::valarray<Precision>& fv,
        const size_t nx,
        const size_t ny,
        const function2D f)
    {
        const Precision hx=1.0/nx;
        const Precision hy=1.0/ny;
        for (size_t j=0; j<=ny; j++)
            for (size_t i=0; i<=nx; i++)
                fv[j*(nx+1)+i]=f(i*hx,j*hy);
    }
}
