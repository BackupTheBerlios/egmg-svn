/** \file GSLexicographic.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementaion of the class GSLexicographic.
 * \see GSLexicographic.h
 */
 

#include "GSLexicographic.h"
#include "../general/DiscreteFunction.h"

namespace mg
{
void GSLexicographic::relax(
    Problem& problem) const
{
    problem.applyBoundaryConstraint();
    Precision factor = 1.0;
    Stencil& stencil = problem.getStencil();
    DiscreteFunction& u = problem.getSolution();
    const DiscreteFunction& f = problem.getRightHandSide();
    const IndexPair fp = problem.getFirstPoint();
    const IndexPair lp = problem.getLastPoint();
    const Index nx = problem.getNx();
    const Index ny = problem.getNy();
    const Precision hx = problem.getHx();
    const Precision hy = problem.getHy();
    const Point origin = problem.getOrigin();
    if (stencil.size() < 2)
    {
        for (Index sy=fp.sy;sy<=lp.sy;sy++)
            for (Index sx=fp.sx;sx<=lp.sx;sx++)
            {
                factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny,hx,hy,origin);
                u(sx,sy)+=factor*(f(sx,sy)
                        -stencil.apply(u,C,sx,sy));
            }
    }
    else
    {
        /*factor = 1.0/stencil.getCenter(SW,fp.sx-llc.sx,fp.sy-llc.sy,nx,ny);
        u[fp.sx*(nx+3)+fp.sy]+=factor*(f[1*(nx+3)+1]
                -stencil.apply(u,SW,fp.sx-llc.sx,fp.sy-llc.sy,nx,ny));
        for (Index sx=fp.sx+1;sx<=lp.sx-1;sx++)
        {
            factor = 1.0/stencil.getCenter(S,sx-llc.sx,fp.sy-llc.sy,nx,ny);
            u[1*(nx+3)+sx]+=factor*(f[1*(nx+3)+sx]
                -stencil.apply(u,S,sx-llc.sx,fp.sy-llc.sy,nx,ny));
        }
        factor = 1.0/stencil.getCenter(SE,lp.sx-llc.sx,fp.sy-llc.sy,nx,ny);
        u[1*(nx+3)+(nx-1)]+=factor*(f[1*(nx+3)+(nx-1)]
                -stencil.apply(u,SE,lp.sx-llc.sx,fp.sy-llc.sy,nx,ny));
        //everything up to north west corner
        for (Index sy=fp.sy+1;sy<=lp.sy-1;sy++)
        {
            factor = 1.0/stencil.getCenter(W,fp.sx-llc.sx,sy-llc.sy,nx,ny);
            u[sy*(nx+3)+1]+=factor*(f[sy*(nx+3)+1]
                -stencil.apply(u,W,fp.sx-llc.sx,sy-llc.sy,nx,ny));
            for (Index sx=fp.sx+1;sx<=lp.sx-1;sx++)
            {
                factor = 1.0/stencil.getCenter(C,sx-llc.sx,sy-llc.sy,nx,ny);
                u[sy*(nx+3)+sx]+=factor*(f[sy*(nx+3)+sx]
                    -stencil.apply(u,C,sx-llc.sx,sy-llc.sy,nx,ny));
            }
            factor = 1.0/stencil.getCenter(E,lp.sx-llc.sx,sy-llc.sy,nx,ny);
            u[sy*(nx+3)+(nx-1)]+=factor*(f[sy*(nx+3)+(nx-1)]
                -stencil.apply(u,E,lp.sx-llc.sx,sy-llc.sy,nx,ny));
        }
        factor = 1.0/stencil.getCenter(NW,fp.sx-llc.sx,lp.sy-llc.sy,nx,ny);
        u[(ny-1)*(nx+3)+1]+=factor*(f[(ny-1)*(nx+3)+1]
                -stencil.apply(u,NW,fp.sx-llc.sx,lp.sy-llc.sy,nx,ny));
        for (Index sx=fp.sx+1;sx<=lp.sx-1;sx++)
        {
            factor = 1.0/stencil.getCenter(N,sx-llc.sx,lp.sy-llc.sy,nx,ny);
            u[(nx-1)*(nx+3)+sx]+=factor*(f[(nx-1)*(nx+3)+sx]
                -stencil.apply(u,N,sx-llc.sx,lp.sy-llc.sy,nx,ny));
        }
        factor = 1.0/stencil.getCenter(NE,lp.sx-llc.sx,lp.sy-llc.sy,nx,ny);
        u[(nx-1)*(nx+3)+(nx-1)]+=factor*(f[(nx-1)*(nx+3)+(nx-1)]
            -stencil.apply(u,NE,lp.sx-llc.sx,lp.sy-llc.sy,nx,ny));*/
    }
    problem.applyBoundaryConstraint();
}
}
