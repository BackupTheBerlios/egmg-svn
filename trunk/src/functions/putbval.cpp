/** \file putbval.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the function putbval
 */

#include "putbval.h"

namespace mg
{
    void putbval(
        NumericArray& u,
        const Index nx,
        const Index ny,
        const function2D g)
    {
        const Precision hx = 1.0/nx;
        const Precision hy = 1.0/ny;
        for (Index sx=0; sx<=nx; sx++)
        {
            u[sx]=g(sx*hx,0);                     //top border
            u[u.size()-1-sx]=g((nx-sx)*hx,1);     //bottom border
        }
        for (Index sy=0; sy<=ny; sy++)
        {
            u[sy*(nx+1)]=g(0,sy*hy);              //left border
            u[(sy+1)*(nx+1)-1]=g(1,sy*hy);        //right border
        }
    }
}
