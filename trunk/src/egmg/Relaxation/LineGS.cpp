/** \file LineGS.cpp
 * \author Andre Oeckerath
 * \brief LineGS.cpp contains the implementation of the class LineGS.
 * \see LineGS.h
 * \todo clean up ninepointxline, ninepointyline, xline and yline by
 *       doing redudant things in seperate functions
 */

#include <iostream>
#include "LineGS.h"

namespace mg
{
void LineGS::ninePointX(
    NumericArray &u,
    const NumericArray &f,
    const Stencil &stencil,
    const Index nx,
    const Index ny) const
{
    if(nx > 2)
    {
        solveNinePointXLine(u,f,u,stencil,nx,ny,1,1,SW,S,SE);
        
        for(Index sy=2; sy<ny-1 ; sy++)
            solveNinePointXLine(u,f,u,stencil,nx,ny,1,sy,W,C,E);

        solveNinePointXLine(u,f,u,stencil,nx,ny,1,ny-1,NW,N,NE);
    }
    else
    {
        gsLexicographic_.relax(u,f,stencil,nx,ny);
    }
}

void LineGS::ninePointY(
    NumericArray &u,
    const NumericArray &f,
    const Stencil &stencil,
    const Index nx,
    const Index ny) const
{
    if(ny > 2)
    {
        solveNinePointYLine(u, f, u, stencil, nx, ny, 1, 1, SW, W, NW);

        for(Index sx=2; sx<nx-1 ; sx++)
            solveNinePointYLine(u, f, u, stencil, nx, ny, sx, 1, S, C, N);

        solveNinePointYLine(u, f, u, stencil, nx, ny, nx-1, 1, SE, E, NE);
    }

    else
    {
        gsLexicographic_.relax(u,f,stencil,nx,ny);
    }
}

void LineGS::fullX(
    NumericArray &u,
    const NumericArray &f,
    const Stencil &stencil,
    const Index nx,
    const Index ny) const

{
    if ((ny > 4) && (nx > 4))
    {
        solveFullXLine(u, f, u, stencil, nx, ny, 1, 1, SW, S, SE);

        for(Index sy=2; sy<ny-1; sy++)
            solveFullXLine(u, f, u, stencil, nx, ny, 1, sy, W, C, E);

        solveFullXLine(u, f, u, stencil, nx, ny, 1, nx-1, NW, N, NE);
    }
    else
    {
        for(int k=0; k<2; k++)
        {
            gsLexicographic_.relax(u,f,stencil,nx,ny);
        }
    }
}

void LineGS::fullY(
    NumericArray &u,
    const NumericArray &f,
    const Stencil &stencil,
    const Index nx,
    const Index ny) const
{
    if ((ny > 4) && (nx > 4))
    {
        solveFullYLine(u, f, u, stencil, nx, ny, 1, 1, SW, W, NW);

        for(Index sx=2; sx<nx-1; sx++)
            solveFullYLine(u, f, u, stencil, nx, ny, sx, 1, S, C, N);

        solveFullYLine(u, f, u, stencil, nx, ny, nx-1, 1, SE, E, NE);
    }
    else
    {
        for(int k=0; k<2; k++)
        {
            gsLexicographic_.relax(u,f,stencil,nx,ny);
        }
    }
}
}
