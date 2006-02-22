/** \file putbval.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the function putbval
 */

#include "putbval.h"

namespace mg
{
    void putbval(
        std::valarray<Precision>& u,
        const size_t nx,
        const size_t ny,
        const function2D g)
    {
        const Precision hx = 1.0/nx;
        const Precision hy = 1.0/ny;
        for (size_t i=0; i<=nx; i++)
        {
            u[i]=g(i*hx,0);                     //top border
            u[u.size()-1-i]=g((nx-i)*hx,1);     //bottom border
        }
        for (size_t i=0; i<=ny; i++)
        {
            u[i*(nx+1)]=g(0,i*hy);              //left border
            u[(i+1)*(nx+1)-1]=g(1,i*hy);        //right border
        }
    }
}
