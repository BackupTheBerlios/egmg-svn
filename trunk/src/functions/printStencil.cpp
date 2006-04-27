/** \file printStencil.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief contains the implementation of the function printStencil
 */
#include "printStencil.h"
#include <iomanip>

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
                    out<<std::setw(8)<<std::setprecision(2)<<L[k];
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
            if ( k == jX.size() )
                out<<"          ";
            out<<" ";
        }
        out<<"\n\n";
    }
}

}
