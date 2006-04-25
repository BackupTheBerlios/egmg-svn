/** \file generatePositionArrays.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains implementation of the function generatePositionArrays.
 */
#include "generatePositionArrays.h"
#include <stack>

namespace mg
{

void generatePositionArrays(
    PositionArray& jX,
    PositionArray& jY,
    const Position pos,
    const Index size,
    const Index sizeToBorder)
{
    //ASSERT( sizeToBorder > 1 )
    if ( C == pos )
    {
        jX.resize((2*size+1)*(2*size+1));
        jY.resize((2*size+1)*(2*size+1));
        jX=0;
        jY=0;
    }
    else if ( W == pos || N == pos || E == pos || S == pos )
    {
        jX.resize((size+sizeToBorder+1)*(2*size+1));
        jY.resize((size+sizeToBorder+1)*(2*size+1));
        jX=0;
        jY=0;
    }
    else
    {
        jX.resize((size+sizeToBorder+1)*(size+sizeToBorder+1));
        jY.resize((size+sizeToBorder+1)*(size+sizeToBorder+1));
        jX=0;
        jY=0;
    }
    std::stack<Integer> xPositions;
    std::stack<Integer> yPositions;
    
    Integer intSize = size;
    Integer intSizeToBorder = sizeToBorder;
    for ( Integer k=intSize; k>=1; --k)
    {
        xPositions.push(k);
        if ( k > intSizeToBorder && ( W == pos || SW == pos || NW == pos ||
                                      E == pos || NE == pos || SE == pos ) )
        {
            continue;
        }
        xPositions.push(k);
    }
    for ( Integer k=intSize; k>=1; --k)
    {
        yPositions.push(k);
        if ( k > intSizeToBorder && ( S == pos || SW == pos || SE == pos ||
                                      N == pos || NE == pos || NW == pos ) )
        {
            continue;
        }
        yPositions.push(k);
    }
    jX[0]=0;
    jY[0]=0;
    for ( Index i=1; i<=4*sizeToBorder; ++i )
    {
        if ( i%4 == 0 )
        {
            jX[i] = 0;
            jY[i] = -1*yPositions.top();
            yPositions.pop();
        }
        else if ( i%4 == 1 )
        {
            jX[i] = -1*xPositions.top();
            xPositions.pop();
            jY[i] = 0;
        }
        else if ( i%4 == 2 )
        {
            jX[i] = 0;
            jY[i] = yPositions.top();
            yPositions.pop();
        }
        else if ( i%4 == 3 )
        {
            jX[i] = xPositions.top();
            xPositions.pop();
            jY[i] = 0;
        }
    }
    if ( pos == C )
    {
        for ( Index i=4*sizeToBorder+1; i<4*size+1; ++i )
        {
            if ( i%4 == 0 )
            {
                jX[i] = 0;
                jY[i] = -1*yPositions.top();
                yPositions.pop();
            }
            else if ( i%4 == 1 )
            {
                jX[i] = -1*xPositions.top();
                xPositions.pop();
                jY[i] = 0;
            }
            else if ( i%4 == 2 )
            {
                jX[i] = 0;
                jY[i] = yPositions.top();
                yPositions.pop();
            }
            else if ( i%4 == 3 )
            {
                jX[i] = xPositions.top();
                xPositions.pop();
                jY[i] = 0;
            }
        }
    }
    else if ( pos == W )
    {
        for ( Index i=4*sizeToBorder+1; i<=3*size+sizeToBorder; ++i )
        {
            if ( i%3 == 2 )
            {
                jY[i] = yPositions.top();
                yPositions.pop();
                jX[i] = 0;
            }
            else if ( i%3 == 0 )
            {
                jX[i] = xPositions.top();
                xPositions.pop();
                jY[i] = 0;
            }
            else if ( i%3 == 1 )
            {
                jX[i] = 0;
                jY[i] = -1*yPositions.top();
                yPositions.pop();
            }
        }
    }
    else if ( pos == N )
    {
        for ( Index i=4*sizeToBorder+1; i<=3*size+sizeToBorder; ++i )
        {
            if ( i%3 == 2 )
            {
                jX[i] = -1*xPositions.top();
                xPositions.pop();
                jY[i] = 0;
            }
            else if ( i%3 == 0 )
            {
                jX[i] = xPositions.top();
                xPositions.pop();
                jY[i] = 0;
            }
            else if ( i%3 == 1 )
            {
                jX[i] = 0;
                jY[i] = -1*yPositions.top();
                yPositions.pop();
            }
        }
    }
    else if ( pos == E )
    {
        for ( Index i=4*sizeToBorder+1; i<=3*size+sizeToBorder; ++i )
        {
            if ( i%3 == 2 )
            {
                jX[i] = -1*xPositions.top();
                xPositions.pop();
                jY[i] = 0;
            }
            else if ( i%3 == 0 )
            {
                jY[i] = yPositions.top();
                yPositions.pop();
                jX[i] = 0;
            }
            else if ( i%3 == 1 )
            {
                jX[i] = 0;
                jY[i] = -1*yPositions.top();
                yPositions.pop();
            }
        }
    }
    else if ( pos == S )
    {
        for ( Index i=4*sizeToBorder+1; i<=3*size+sizeToBorder; ++i )
        {
            if ( i%3 == 2 )
            {
                jX[i] = -1*xPositions.top();
                xPositions.pop();
                jY[i] = 0;
            }
            else if ( i%3 == 0 )
            {
                jX[i] = 0;
                jY[i] = yPositions.top();
                yPositions.pop();
            }
            else if ( i%3 == 1 )
            {
                jX[i] = xPositions.top();
                xPositions.pop();
                jY[i] = 0;
            }
        }
    }
    else if ( pos == NW )
    {
        for ( Index i=4*sizeToBorder+1; i<=2*size+2*sizeToBorder; ++i )
        {
            if ( i%2 == 1 )
            {
                jX[i] = xPositions.top();
                xPositions.pop();
                jY[i] = 0;
            }
            else if ( i%2 == 0 )
            {
                jX[i] = 0;
                jY[i] = -1*yPositions.top();
                yPositions.pop();
            }
        }
    }
    else if ( pos == NE )
    {
        for ( Index i=4*sizeToBorder+1; i<=2*size+2*sizeToBorder; ++i )
        {
            if ( i%2 == 1 )
            {
                jX[i] = -1*xPositions.top();
                xPositions.pop();
                jY[i] = 0;
            }
            else if ( i%2 == 0 )
            {
                jX[i] = 0;
                jY[i] = -1*yPositions.top();
                yPositions.pop();
            }
        }
    }
    else if ( pos == SE )
    {
        for ( Index i=4*sizeToBorder+1; i<=2*size+2*sizeToBorder; ++i )
        {
            if ( i%2 == 1 )
            {
                jX[i] = -1*xPositions.top();
                xPositions.pop();
                jY[i] = 0;
            }
            else if ( i%2 == 0 )
            {
                jX[i] = 0;
                jY[i] = yPositions.top();
                yPositions.pop();
            }
        }
    }
    else if ( pos == SW )
    {
        for ( Index i=4*sizeToBorder+1; i<=2*size+2*sizeToBorder; ++i )
        {
            if ( i%2 == 1 )
            {
                jX[i] = 0;
                jY[i] = yPositions.top();
                yPositions.pop();
            }
            else if ( i%2 == 0 )
            {
                jX[i] = xPositions.top();
                xPositions.pop();
                jY[i] = 0;
            }
        }
    }
    //ASSERT( xPositions.empty() && yPositions.empty() );
    for ( Integer k=-intSize; k<=-1; ++k )
    {
        if ( k < -1*intSizeToBorder && ( W == pos || SW == pos || NW == pos ) )
        {
            continue;
        }
        xPositions.push(k);
    }
    for ( Integer k=-intSize; k<=-1; ++k )
    {
        if ( k < -1*intSizeToBorder && ( S == pos || SW == pos || SE == pos ) )
        {
            continue;
        }
        yPositions.push(k);
    }
    for ( Integer k=1; k<=intSize; ++k )
    {
        if ( k > intSizeToBorder && ( N == pos || NW == pos || NE == pos ) )
        {
            continue;
        }
        yPositions.push(k);   
    }
    for ( Integer k=1; k<=intSize; ++k )
    {
        if ( k > intSizeToBorder && ( E == pos || NE == pos || SE == pos ) )
        {
            continue;
        }
        xPositions.push(k);
    }
    Index i=2*size+2*sizeToBorder+1;
    if ( C == pos )
    {
        i=4*size+1;
    }
    else if ( W == pos || N == pos || E == pos || S == pos )
    {
        i=3*size+sizeToBorder+1;
    }
    while( i<jX.size() )
    {
        std::stack<Integer> xP = xPositions;
        Integer yPosition = yPositions.top();
        yPositions.pop();
        while( !xP.empty() )
        {
            jX[i] = xP.top();
            xP.pop();
            jY[i] = yPosition;
            ++i;
            if ( i>=jX.size() )
                break;
        }
        //ASSERT( !yPositions.empty() );    
    }
}

}
