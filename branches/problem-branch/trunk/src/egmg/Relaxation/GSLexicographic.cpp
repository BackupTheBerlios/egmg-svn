/** \file GSLexicographic.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementaion of the class GSLexicographic.
 * \see GSLexicographic.h
 */
 

#include "GSLexicographic.h"

namespace mg
{
void GSLexicographic::relax(
    Problem& problem) const
{
    Precision factor = 1.0;
    Stencil& stencil = problem.getStencil();
    NumericArray& u = problem.getSolution();
    const NumericArray& f = problem.getRightHandSide();
    Point llc = problem.getLowerLeftCorner();
    Point fp = problem.getFirstPoint();
    Point lp = problem.getLastPoint();
    const Index nx = problem.getNx();
    const Index ny = problem.getNy();
    if (stencil.size() < 2)
    {
        for (Index sy=fp.sy;sy<=lp.sy;sy++)
            for (Index sx=fp.sx;sx<=lp.sx;sx++)
            {
                factor = 1.0/stencil.getCenter(C,llc.sx+sx,llc.sy+sy,nx,ny);
                u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
                        -stencil.apply(u,C,sx,sy,nx,ny));
            }
    }
    else
    {
        factor = 1.0/stencil.getCenter(SW,llc.sx+fp.sx,llc.sy+fp.sy,nx,ny);
        u[1*(nx+1)+1]+=factor*(f[1*(nx+1)+1]
                -stencil.apply(u,SW,llc.sx+fp.sx,llc.sy+fp.sy,nx,ny));
        for (Index sx=fp.sx+1;sx<=lp.sx-1;sx++)
        {
            factor = 1.0/stencil.getCenter(S,sx,llc.sy+fp.sy,nx,ny);
            u[1*(nx+1)+sx]+=factor*(f[1*(nx+1)+sx]
                -stencil.apply(u,S,sx,llc.sy+fp.sy,nx,ny));
        }
        factor = 1.0/stencil.getCenter(SE,llc.sx+lp.sx,llc.sy+fp.sy,nx,ny);
        u[1*(nx+1)+(nx-1)]+=factor*(f[1*(nx+1)+(nx-1)]
                -stencil.apply(u,SE,llc.sx+lp.sx,llc.sy+fp.sy,nx,ny));
        //everything up to north west corner
        for (Index sy=fp.sy+1;sy<=lp.sy-1;sy++)
        {
            factor = 1.0/stencil.getCenter(W,llc.sx+fp.sx,sy,nx,ny);
            u[sy*(nx+1)+1]+=factor*(f[sy*(nx+1)+1]
                -stencil.apply(u,W,llc.sx+fp.sx,sy,nx,ny));
            for (Index sx=fp.sx+1;sx<=lp.sx-1;sx++)
            {
                factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
                u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
                    -stencil.apply(u,C,sx,sy,nx,ny));
            }
            factor = 1.0/stencil.getCenter(E,llc.sx+lp.sx,sy,nx,ny);
            u[sy*(nx+1)+(nx-1)]+=factor*(f[sy*(nx+1)+(nx-1)]
                -stencil.apply(u,E,llc.sx+lp.sx,sy,nx,ny));
        }
        factor = 1.0/stencil.getCenter(NW,llc.sx+fp.sx,llc.sy+lp.sy,nx,ny);
        u[(ny-1)*(nx+1)+1]+=factor*(f[(ny-1)*(nx+1)+1]
                -stencil.apply(u,NW,llc.sx+fp.sx,llc.sy+lp.sy,nx,ny));
        for (Index sx=fp.sx+1;sx<=lp.sx-1;sx++)
        {
            factor = 1.0/stencil.getCenter(N,sx,llc.sy+lp.sy,nx,ny);
            u[(nx-1)*(nx+1)+sx]+=factor*(f[(nx-1)*(nx+1)+sx]
                -stencil.apply(u,N,sx,llc.sy+lp.sy,nx,ny));
        }
        factor = 1.0/stencil.getCenter(NE,llc.sx+lp.sx,llc.sy+lp.sy,nx,ny);
        u[(nx-1)*(nx+1)+(nx-1)]+=factor*(f[(nx-1)*(nx+1)+(nx-1)]
            -stencil.apply(u,NE,llc.sx+lp.sx,llc.sy+lp.sy,nx,ny));
    }
    problem.applyBoundaryConstraint();
}
}
