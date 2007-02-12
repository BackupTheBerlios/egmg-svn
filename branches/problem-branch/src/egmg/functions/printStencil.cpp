/** \file printStencil.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief contains the implementation of the function printStencil
 */
#include "printStencil.h"
#include <iomanip>
#include "../Stencil/Stencil.h"

namespace mg
{

void printStencil(
    const NumericArray& L,
    const PositionArray& jX,
    const PositionArray& jY,
    std::ostream& out)
{
    for ( mg::Integer sy=jY.max(); sy>=jY.min(); --sy )
    {
        for ( mg::Integer sx=jX.min(); sx<=jX.max(); ++sx )
        {
            mg::Index k=0;
            for ( k=0; k<jX.size(); ++k )
            {
                if ( jY[k] == sy && jX[k] == sx )
                {
                    if ( sx == 0 && sy == 0 )
                    {
                        out<<"[";
                    }
                    else
                    {
                        out<<" ";
                    }
                    out<<std::setw(10)<<std::setprecision(6)<<L[k];
                    if ( sx == 0 && sy == 0 )
                    {
                        out<<"]";
                    }
                    else
                    {
                        out<<" ";
                    }
                    break;
                }   
            }
            if ( k == jX.size() && sx == 0 && sy == 0 )
                out<<" [        ] ";
            else if ( k == jX.size() )
                out<<"     --     ";
            out<<" ";
        }
        out<<"\n\n";
    }
}

void printAllStencils(
    const Stencil& stencil,
    const Index nx,
    const Index ny,
    std::ostream& out,
    Precision scale)
{
    if (scale==0)
        scale=nx*ny;
    for (Index sx = 1; sx<nx; ++sx)
        for (Index sy = 1; sy<ny; ++sy)
        {
            out<<sx<<" "<<sy<<std::endl;
            printStencil(   
                    stencil.getL(C,sx,sy,nx,ny)/scale,
                    stencil.getJx(C,nx,ny),
                    stencil.getJy(C,nx,ny),
                    out);
        }
}
}
