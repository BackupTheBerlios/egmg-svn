/** \file LineRelaxation.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief LineRelaxation.cpp contains the impl. of the class LineRelaxation.
 * \see LineRelaxation.h
 */
#include<iostream>
#include "LineRelaxation.h"

namespace mg
{
void LineRelaxation::xLRSolver(
    NumericArray& u,
    const size_t sy,
    const size_t nx,
    NumericArray& rhs,
    NumericArray& ndiagL,
    NumericArray& diagR,
    const NumericArray& ndiagR) const
{
    for(size_t sx=1; sx<nx-1; sx++) 
    {
        ndiagL[sx-1]=ndiagL[sx-1]/diagR[sx-1];  
        diagR[sx]-=ndiagL[sx-1]*ndiagR[sx-1]; 
        rhs[sy*(nx+1)+1+sx]=rhs[sy*(nx+1)+1+sx]-ndiagL[sx-1]*rhs[sy*(nx+1)+sx];  
    }
    // solve the linear system of equations R u = rhs
    u[sy*(nx+1)+(nx-1)]=rhs[sy*(nx+1)+(nx-1)]/diagR[nx-2];
    for(size_t sx=nx-2; sx>0; sx--)
    {
        u[sy*(nx+1)+sx]=1/diagR[sx-1]*(rhs[sy*(nx+1)+sx]
                        -ndiagR[sx-1]*u[sy*(nx+1)+sx+1]);
    }
}
void LineRelaxation::yLRSolver(
    NumericArray& u,
    const size_t sx,
    const size_t nx,
    const size_t ny,
    NumericArray& rhs,
    NumericArray& ndiagL,
    NumericArray& diagR,
    const NumericArray& ndiagR) const
{
    // LR-decomposition + transformation of the rhs vector
    for(size_t sy=1; sy<ny-1; sy++) 
    {
        ndiagL[sy-1]=ndiagL[sy-1]/diagR[sy-1]; 
        diagR[sy]-=ndiagL[sy-1] * ndiagR[sy-1];
        rhs[(sy+1)*(nx+1)+sx]=rhs[(sy+1)*(nx+1)+sx]
                                        -ndiagL[sy-1]*rhs[sy*(nx+1)+sx];  
    }
    // solve the linear system of equations R u = rhs
    u[sx+(nx+1)*(ny-1)]=rhs[sx+(nx+1)*(ny-1)]/diagR[ny-2];
    for(size_t sy=ny-2; sy>0; sy--)
    {
        u[sx+sy*(nx+1)]=1/diagR[sy-1]*(rhs[sx+sy*(nx+1)]
                                        -ndiagR[sy-1]*u[sx+(sy+1)*(nx+1)]);
    }
}
void LineRelaxation::xLRSolver(
    NumericArray& u,
    const size_t sy,
    const size_t nx,
    NumericArray& rhs,
    NumericArray& ndiagL1,
    NumericArray& ndiagL2,
    NumericArray& diagR,
    NumericArray& ndiagR1,
    const NumericArray& ndiagR2) const
{
    // LR-decomposition + transformation of the rhs
    for(size_t sx=1; sx<nx-2; sx++)  
    {
        ndiagL1[sx-1]=ndiagL1[sx-1]/diagR[sx-1];
        diagR[sx]-=ndiagL1[sx-1]*ndiagR1[sx-1];
        ndiagR1[sx]-=ndiagL1[sx-1]*ndiagR2[sx-1];
        rhs[sy*(nx+1)+1+sx]=rhs[sy*(nx+1)+1+sx]-ndiagL1[sx-1]*rhs[sy*(nx+1)+sx]; 
        ndiagL2[sx-1]=ndiagL2[sx-1]/diagR[sx-1];
        ndiagL1[sx]-=ndiagL2[sx-1]*ndiagR1[sx-1];
        diagR[sx+1]-=ndiagL2[sx-1]*ndiagR2[sx-1];
        rhs[sy*(nx+1)+1+sx+1]=rhs[sy*(nx+1)+1+sx+1]
                                               -ndiagL2[sx-1]*rhs[sy*(nx+1)+sx];
    }
    ndiagL1[nx-2-1]=ndiagL1[nx-2-1]/diagR[nx-2-1];
    diagR[nx-2]-=ndiagL1[nx-2-1]*ndiagR1[nx-2-1];
    rhs[sy*(nx+1)+1+nx-2]=rhs[sy*(nx+1)+1+nx-2]
                                           -ndiagL1[nx-2-1]*rhs[sy*(nx+1)+nx-2];
    // solve the linear system of equations R u = rhs
    u[sy*(nx+1)+(nx-1)]=rhs[sy*(nx+1)+(nx-1)]/diagR[nx-2];
    for(size_t sx=nx-2; sx>1; sx--)
    {
        u[sy*(nx+1)+sx]=1/diagR[sx-1]*(rhs[sy*(nx+1)+sx]
                                              -ndiagR1[sx-1]*u[sy*(nx+1)+sx+1]);
        rhs[sy*(nx+1)+sx-1]-=ndiagR2[sx-2]*u[sy*(nx+1)+sx+1];
    }
    u[sy*(nx+1)+1]=1/diagR[0]*(rhs[sy*(nx+1)+1]-ndiagR1[0]*u[sy*(nx+1)+1+1]);
}
void LineRelaxation::yLRSolver(
    NumericArray& u,
    const size_t sx,
    const size_t nx,
    const size_t ny,
    NumericArray& rhs,
    NumericArray& ndiagL1,
    NumericArray& ndiagL2,
    NumericArray& diagR,
    NumericArray& ndiagR1,
    const NumericArray& ndiagR2) const
{
    // LR-decomposition + transformation of the rhs
    for(size_t sy=1; sy<ny-2; sy++)  
    {
        ndiagL1[sy-1]=ndiagL1[sy-1]/diagR[sy-1];
        diagR[sy]-=ndiagL1[sy-1]*ndiagR1[sy-1];
        ndiagR1[sy]-=ndiagL1[sy-1]*ndiagR2[sy-1];
        rhs[(sy+1)*(nx+1)+sx]=rhs[(sy+1)*(nx+1)+sx]
                                               -ndiagL1[sy-1]*rhs[sy*(nx+1)+sx]; 
        ndiagL2[sy-1]=ndiagL2[sy-1]/diagR[sy-1];
        ndiagL1[sy]-=ndiagL2[sy-1]*ndiagR1[sy-1];
        diagR[sy+1]-=ndiagL2[sy-1]*ndiagR2[sy-1];
        rhs[(sy+2)*(nx+1)+sx]=rhs[(sy+2)*(nx+1)+sx]
                                               -ndiagL2[sy-1]*rhs[sy*(nx+1)+sx];
    }
    ndiagL1[ny-2-1]=ndiagL1[ny-2-1]/diagR[ny-2-1];
    diagR[ny-2]-=ndiagL1[ny-2-1]*ndiagR1[ny-2-1];
    rhs[(ny-1)*(nx+1)+sx]=rhs[(ny-1)*(nx+1)+sx]
                                        -ndiagL1[ny-2-1]*rhs[(ny-2)*(nx+1)+sx];
    // solve the linear system of equations R u = rhs
    u[sx+(nx+1)*(ny-1)]=rhs[sx+(nx+1)*(ny-1)]/diagR[ny-2];
    for(size_t sy=ny-2; sy>1; sy--)
    {
        u[sx+sy*(nx+1)]=1/diagR[sy-1]*(rhs[sx+sy*(nx+1)]
                                            -ndiagR1[sy-1]*u[sx+(sy+1)*(nx+1)]);
        rhs[(sy-1)*(nx+1)+sx]-=ndiagR2[sy-2]*u[(sy+1)*(nx+1)+sx];
    }
    u[(nx+1)+sx]=1/diagR[0]*(rhs[(nx+1)+sx]-ndiagR1[0]*u[2*(nx+1)+sx]);
}
}
