/** \file ZebraLineGS.cpp
 * \author Andre Oeckerath
 * \brief ZebraLineGS.cpp contains the implementation of the class ZebraLineGS.
 * \see ZebraLineGS.h
 * \todo clean up ninepointxzebra, ninepointyzebra, xzebra and yzebra by
 *       doing redudant things in seperate functions
 */
#include<iostream>
#include "ZebraLineGS.h"

namespace mg
{       
void ZebraLineGS::relax(
    NumericArray &u,
    const NumericArray &f, 
    const Stencil &stencil,
    const Index nx,
    const Index ny) const
{
    // valarrays needed for LR-decomposition of a tridiagonal matrix
    NumericArray rhs(0.0,u.size());
    switch (stencil.size())
    {
    case 1:  // stencil of size 1
        switch (direction_)
        {
        case ALTDIR:
            ninepointxzebra(u, f, rhs, stencil, nx, ny);
            ninepointyzebra(u, f, rhs, stencil, nx, ny);
            break;
        case XDIR:
            ninepointxzebra(u, f, rhs, stencil, nx, ny);
            break;
        case YDIR:
            ninepointyzebra(u, f, rhs, stencil, nx, ny);
            break;
        default:
            std::cerr << "Error in direction of the line relaxation!\n";
            break;
        }
        break;          
    case 2:  //stencil of size 2
        switch (direction_)
        {
        case ALTDIR:
            xzebra(u, f, rhs, stencil, nx, ny);
            yzebra(u, f, rhs, stencil, nx, ny);
            break;
        case XDIR:
            xzebra(u, f, rhs, stencil, nx, ny);
            break;
        case YDIR:
            yzebra(u, f, rhs, stencil, nx, ny);
            break;
        default:
            std::cerr << "Error in direction of the line relaxation!\n";
            break;
        }
        break;
    default:
        std::cerr << "Stencil is too big (size>2)!" << std::endl;
        break;
    }
}
void ZebraLineGS::ninepointxzebra(
    NumericArray &u,
    const NumericArray &f, 
    NumericArray &rhs,
    const Stencil &stencil,
    const Index nx, 
    const Index ny) const                      
{ 
    //valarrays needed for saving the tridiagonal matrix A of linear system
    //A u = rhs
    NumericArray diagR(nx-1);
    NumericArray ndiagR(nx-2);
    NumericArray ndiagL(nx-2);
    if(stencil.isConstant())
    {
        // get const operator L
        const NumericArray operatorL = stencil.getL(C,2,2,nx,ny);
        const PositionArray jX = stencil.getJx(C,nx,ny);
        const PositionArray jY = stencil.getJy(C,nx,ny);
            
        // for each line: correction of the rhs given by 
        // rhs = fv - [L[n]  0  L[s]]^t * u and elimination of the 
        // boundary condition in first and last inner point
        // odd lines
        for(Index sy=1; sy<ny ; sy+=2) 
        {
            rhs[1+sy*(nx+1)] = f[1+sy*(nx+1)]
                    -operatorL[N]*u[1+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[1+(sy+jY[S])*(nx+1)]
                    -operatorL[W]*u[sy*(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[1+sy*(nx+1)]-=
                                operatorL[i]*u[1+jX[i]+(sy+jY[i])*(nx+1)];
            }
            for(Index sx=2; sx<nx-1; sx++)  
            {
                rhs[sx+sy*(nx+1)] = f[sx+sy*(nx+1)] 
                    -operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[sx+(sy+jY[S])*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+sy*(nx+1)]-=
                                operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                }
            }
            rhs[(nx-1)+sy*(nx+1)] = f[nx-1+sy*(nx+1)]
                    -operatorL[N]*u[nx-1+(sy+jY[N])*(nx+1)] 
                    -operatorL[S]*u[nx-1+(sy+jY[S])*(nx+1)]
                    -operatorL[E]*u[nx+sy*(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[nx-1+sy*(nx+1)]-=
                                operatorL[i]*u[nx-1+jX[i]+(sy+jY[i])*(nx+1)];
            }
            // set tridiagonalmatrix for solving A u = rhs
            // A[i][i] = L[c]; A[i-1][i] = L[w]; A[i+1][i] = L[e]
            diagR = operatorL[C];
            ndiagR = operatorL[E];
            ndiagL = operatorL[W];
            // solve the tridiagonal system with LR-decomposition
            xLRSolver(u,sy,nx,rhs,ndiagL,diagR,ndiagR);
        }
        // same for even lines 
        for(Index sy=2; sy<ny ; sy+=2)
        {
            rhs[1+sy*(nx+1)] = f[1+sy*(nx+1)]
                -operatorL[N]*u[1+(sy+jY[N])*(nx+1)]
                -operatorL[S]*u[1+(sy+jY[S])*(nx+1)] 
                -operatorL[W]*u[sy*(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[1+sy*(nx+1)]-=operatorL[i]*u[1+jX[i]+(sy+jY[i])*(nx+1)];
            }
            for(Index sx=2; sx<nx-1; sx++) 
            {
                rhs[sx+sy*(nx+1)] = f[sx+sy*(nx+1)]
                    -operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[sx+(sy+jY[S])*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                }
            }
            rhs[(nx-1)+sy*(nx+1)] = f[nx-1+sy*(nx+1)]
                -operatorL[N]*u[nx-1+(sy+jY[N])*(nx+1)]
                -operatorL[S]*u[nx-1+(sy+jY[S])*(nx+1)] 
                -operatorL[E]*u[nx+sy*(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[nx-1+sy*(nx+1)]-=
                                   operatorL[i]*u[nx-1+jX[i]+(sy+jY[i])*(nx+1)];
            }
            diagR = operatorL[C];
            ndiagR = operatorL[E];
            ndiagL = operatorL[W];
            xLRSolver(u,sy,nx,rhs,ndiagL,diagR,ndiagR);
        }
    }
    else
    {
        //Stencil ist not constant, so L needs to be evaluated in each grid
        //point no other change in the algorithm  
        NumericArray operatorL = stencil.getL(C,2,2,nx,ny);
        PositionArray jX = stencil.getJx(C,nx,ny);
        PositionArray jY = stencil.getJy(C,nx,ny);
        if(nx > 2)
        {
            operatorL=stencil.getL(SW,1,1,nx,ny);
            diagR[0]=operatorL[C];
            ndiagR[0]=operatorL[E];
            rhs[1+nx+1]=f[1+nx+1]
                -operatorL[N]*u[1+(1+jY[N])*(nx+1)]
                -operatorL[S]*u[1+(1+jY[S])*(nx+1)] 
                -operatorL[W]*u[nx+1];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[1+nx+1]-=operatorL[i]*u[1+jX[i]+(1+jY[i])*(nx+1)];
            }
            for(Index sx=2; sx<nx-1; sx++)
            {
                operatorL=stencil.getL(S,sx,1,nx,ny);
                diagR[sx-1]=operatorL[C];
                ndiagR[sx-1]=operatorL[E];
                ndiagL[sx-2]=operatorL[W];
                rhs[sx+nx+1]=f[sx+nx+1]
                    -operatorL[N]*u[sx+(1+jY[N])*(nx+1)]
                    -operatorL[S]*u[sx+(1+jY[S])*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+nx+1]-=operatorL[i]*u[sx+jX[i]+(1+jY[i])*(nx+1)];
                }
            }
            operatorL=stencil.getL(SE,nx-1,1,nx,ny);
            diagR[nx-2]=operatorL[C];
            ndiagL[nx-3]=operatorL[W];
            rhs[(nx-1)+(nx+1)]=f[nx-1+(nx+1)] 
                -operatorL[N]*u[nx-1+(1+jY[N])*(nx+1)]
                -operatorL[S]*u[nx-1+(1+jY[S])*(nx+1)] 
                -operatorL[E]*u[nx+(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[nx-1+(nx+1)]-=operatorL[i]*u[nx-1+jX[i]+(1+jY[i])*(nx+1)];
            }
            xLRSolver(u,1,nx,rhs,ndiagL,diagR,ndiagR);
            for(Index sy=3; sy<ny-1 ; sy+=2)
            {
                operatorL=stencil.getL(W,1,sy,nx,ny);
                diagR[0]=operatorL[C];
                ndiagR[0]=operatorL[E];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorL[N]*u[1+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[1+(sy+jY[S])*(nx+1)] 
                    -operatorL[W]*u[sy*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[1+sy*(nx+1)]-=operatorL[i]*u[1+jX[i]+(sy+jY[i])*(nx+1)];
                }
                for(Index sx=2; sx<nx-1; sx++)
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    diagR[sx-1]=operatorL[C];
                    ndiagR[sx-1]=operatorL[E];
                    ndiagL[sx-2]=operatorL[W];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
                        -operatorL[S]*u[sx+(sy+jY[S])*(nx+1)];
                    for(Index i=5; i<operatorL.size(); i++)
                    {
                        rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                    }
                }
                operatorL=stencil.getL(E,nx-1,sy,nx,ny);
                diagR[nx-2]=operatorL[C];
                ndiagL[nx-3]=operatorL[W];
                rhs[(nx-1)+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorL[N]*u[nx-1+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[nx-1+(sy+jY[S])*(nx+1)] 
                    -operatorL[E]*u[nx+sy*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                   operatorL[i]*u[nx-1+jX[i]+(sy+jY[i])*(nx+1)];
                }
                xLRSolver(u,sy,nx,rhs,ndiagL,diagR,ndiagR);
            }
            operatorL=stencil.getL(NW,1,ny-1,nx,ny);
            diagR[0]=operatorL[C];
            ndiagR[0]=operatorL[E];
            rhs[1+(ny-1)*(nx+1)] = f[1+(ny-1)*(nx+1)]
                -operatorL[N]*u[1+(ny-1+jY[N])*(nx+1)]
                -operatorL[S]*u[1+(ny-1+jY[S])*(nx+1)] 
                -operatorL[W]*u[(ny-1)*(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[1+(ny-1)*(nx+1)]-=
                                    operatorL[i]*u[1+jX[i]+(ny-1+jY[i])*(nx+1)];
            }
            for(Index sx=2; sx<nx-1; sx++)
            {
                operatorL=stencil.getL(N,sx,ny-1,nx,ny);
                diagR[sx-1]=operatorL[C];
                ndiagR[sx-1]=operatorL[E];
                ndiagL[sx-2]=operatorL[W];
                rhs[sx+(ny-1)*(nx+1)] = f[sx+(ny-1)*(nx+1)]
                    -operatorL[N]*u[sx+(ny-1+jY[N])*(nx+1)]
                    -operatorL[S]*u[sx+(ny-1+jY[S])*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+(ny-1)*(nx+1)]-=
                                   operatorL[i]*u[sx+jX[i]+(ny-1+jY[i])*(nx+1)];
                }
            }
            operatorL=stencil.getL(NE,nx-1,ny-1,nx,ny);
            diagR[nx-2]=operatorL[C];
            ndiagL[nx-3]=operatorL[W];
            rhs[(nx-1)+(ny-1)*(nx+1)]=f[nx-1+(ny-1)*(nx+1)]
                -operatorL[N]*u[nx-1+(ny-1+jY[N])*(nx+1)]
                -operatorL[S]*u[nx-1+(ny-1+jY[S])*(nx+1)] 
                -operatorL[E]*u[nx+(ny-1)*(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[nx-1+(ny-1)*(nx+1)]-=
                                operatorL[i]*u[nx-1+jX[i]+(ny-1+jY[i])*(nx+1)];
            }
            xLRSolver(u,(ny-1),nx,rhs,ndiagL,diagR,ndiagR);
            //even lines
            for(Index sy=2; sy<ny ; sy+=2) 
            {                   
                operatorL=stencil.getL(W,1,sy,nx,ny);
                diagR[0]=operatorL[C];
                ndiagR[0]=operatorL[E];
                rhs[1+sy*(nx+1)] = f[1+sy*(nx+1)]
                    -operatorL[N]*u[1+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[1+(sy+jY[S])*(nx+1)] 
                    -operatorL[W]*u[sy*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[1+sy*(nx+1)]-=operatorL[i]*u[1+jX[i]+(sy+jY[i])*(nx+1)];
                }
                // L im Zentrum im Punkt (j/i)
                for(Index sx=2; sx<nx-1; sx++)  
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    diagR[sx-1]=operatorL[C];
                    ndiagR[sx-1]=operatorL[E];
                    ndiagL[sx-2]=operatorL[W];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
                        -operatorL[S]*u[sx+(sy+jY[S])*(nx+1)];
                    for(Index i=5; i<operatorL.size(); i++)
                    {
                        rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                    }
                }
                operatorL=stencil.getL(E,nx-1,sy,nx,ny);
                diagR[nx-2]=operatorL[C];
                ndiagL[nx-3]=operatorL[W];
                rhs[(nx-1)+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorL[N]*u[nx-1+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[nx-1+(sy+jY[S])*(nx+1)] 
                    -operatorL[E]*u[nx+sy*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                   operatorL[i]*u[nx-1+jX[i]+(sy+jY[i])*(nx+1)];
                }
                xLRSolver(u,sy,nx,rhs,ndiagL,diagR,ndiagR);
            }
        }       
        else //nx,ny to small
        {
            gsRedBlack_.relax(u,f,stencil,nx,ny);
        }
    }                       
}
void ZebraLineGS::ninepointyzebra(
    NumericArray &u,
    const NumericArray &f, 
    NumericArray &rhs,
    const Stencil &stencil,
    const Index nx, 
    const Index ny) const
{ 
    //valarrays needed for saving the tridiagonal matrix A of linear system
    //A u = rhs
    NumericArray diagR(ny-1);
    NumericArray ndiagR(ny-2);
    NumericArray ndiagL(ny-2);
    if(stencil.isConstant())
    {
        // get const operator L
        const NumericArray operatorL=stencil.getL(C,2,2,nx,ny);
        const PositionArray jX = stencil.getJx(C,nx,ny);
        const PositionArray jY = stencil.getJy(C,nx,ny);
        // for each odd line correction of the rhs given by
        // rhs = fv + [1  0  1] * u and elimination of the 
        // boundary condition in first and last inner point  
        for(Index sx=1; sx<nx ; sx+=2)
        {
            rhs[sx+(nx+1)]=f[sx+(nx+1)]
                -operatorL[W]*u[sx+jX[W]+(nx+1)]
                -operatorL[E]*u[sx+jX[E]+(nx+1)]
                -operatorL[S]*u[sx+(1+jY[S])*(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[sx+nx+1]-=operatorL[i]*u[sx+jX[i]+(1+jY[i])*(nx+1)];
            }
            for(Index sy=2; sy<ny-1; sy++)
            {
                rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+sy*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                }
            }
            rhs[sx+(ny-1)*(nx+1)]=f[sx+(ny-1)*(nx+1)]
                -operatorL[W]*u[sx+jX[W]+(ny-1)*(nx+1)]
                -operatorL[E]*u[sx+jX[E]+(ny-1)*(nx+1)] 
                -operatorL[N]*u[sx+(ny-1+jY[N])*(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[sx+(ny-1)*(nx+1)]-=
                                   operatorL[i]*u[sx+jX[i]+(ny-1+jY[i])*(nx+1)];
            }
            // set tridiagonalmatrix for solving A u = rhs
            // A[i][i] = L[c]; A[i-1][i] = L[s]; A[i+1][i] = L[n]
            diagR = operatorL[C];
            ndiagR = operatorL[N];
            ndiagL = operatorL[S];
            yLRSolver(u,sx,nx,ny,rhs,ndiagL,diagR,ndiagR);
        }
        // same for each even line
        for(Index sx=2; sx<nx ; sx+=2) 
        {
            rhs[sx+(nx+1)]=f[sx+(nx+1)]
                -operatorL[W]*u[sx+jX[W]+(nx+1)]
                -operatorL[E]*u[sx+jX[E]+(nx+1)] 
                -operatorL[S]*u[sx+(1+jY[S])*(nx+1)]; 
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[sx+nx+1]-=operatorL[i]*u[sx+jX[i]+(1+jY[i])*(nx+1)];
            }
            for(Index sy=2; sy<ny-1; sy++) 
            {
                rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+sy*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                }
            }
            rhs[sx+(ny-1)*(nx+1)]=f[sx+(ny-1)*(nx+1)]
                -operatorL[W]*u[sx+jX[W]+(ny-1)*(nx+1)] 
                -operatorL[E]*u[sx+jX[E]+(ny-1)*(nx+1)]
                -operatorL[N]*u[sx+(ny-1+jY[N])*(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[sx+(ny-1)*(nx+1)]-=
                                   operatorL[i]*u[sx+jX[i]+(ny-1+jY[i])*(nx+1)];
            }
            diagR=operatorL[C];
            ndiagR=operatorL[N];
            ndiagL=operatorL[S];
            yLRSolver(u,sx,nx,ny,rhs,ndiagL,diagR,ndiagR);
        }
    }
    else
    {
        //Stencil ist not constant, so L needs to be evaluated in each grid
        //point no other change in the algorithm  
        NumericArray operatorL = stencil.getL(C,2,2,nx,ny);
        PositionArray jX = stencil.getJx(C,nx,ny);
        PositionArray jY = stencil.getJy(C,nx,ny);
        if(ny > 2)
        {
            operatorL=stencil.getL(SW,1,1,nx,ny);
            diagR[0]=operatorL[C];
            ndiagR[0]=operatorL[N];
            rhs[1+(nx+1)]=f[1+(nx+1)]
                -operatorL[W]*u[1+jX[W]+(nx+1)]
                -operatorL[E]*u[1+jX[E]+(nx+1)] 
                -operatorL[S]*u[1+(1+jY[S])*(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                    rhs[1+nx+1]-=operatorL[i]*u[1+jX[i]+(1+jY[i])*(nx+1)];
            }
            for(Index sy=2; sy<ny-1; sy++) 
            {
                operatorL = stencil.getL(W,1,sy,nx,ny);
                diagR[sy-1]=operatorL[C];
                ndiagR[sy-1]=operatorL[N];
                ndiagL[sy-2]=operatorL[S];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorL[W]*u[1+jX[W]+sy*(nx+1)]
                    -operatorL[E]*u[1+jX[E]+sy*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[1+sy*(nx+1)]-=operatorL[i]*u[1+jX[i]+(sy+jY[i])*(nx+1)];
                }        
            }
            operatorL=stencil.getL(NW,1,ny-1,nx,ny);
            diagR[ny-2]=operatorL[C];
            ndiagL[ny-3]=operatorL[S];
            rhs[1+(ny-1)*(nx+1)]=f[1+(ny-1)*(nx+1)]
                -operatorL[W]*u[1+jX[W]+(ny-1)*(nx+1)] 
                -operatorL[E]*u[1+jX[E]+(ny-1)*(nx+1)]
                -operatorL[N]*u[1+(ny-1+jY[N])*(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[1+(ny-1)*(nx+1)]-=
                                    operatorL[i]*u[1+jX[i]+(ny-1+jY[i])*(nx+1)];
            }
            yLRSolver(u,1,nx,ny,rhs,ndiagL,diagR,ndiagR);
            for(Index sx=3; sx<nx-1 ; sx+=2)
            {
                operatorL=stencil.getL(S,sx,1,nx,ny);
                diagR[0]=operatorL[C];
                ndiagR[0]=operatorL[N];
                rhs[sx+(nx+1)]=f[sx+(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+(nx+1)] 
                    -operatorL[S]*u[sx+(1+jY[S])*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+nx+1]-=operatorL[i]*u[sx+jX[i]+(1+jY[i])*(nx+1)];
                }
                for(Index sy=2; sy<ny-1; sy++) 
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    diagR[sy-1]=operatorL[C];
                    ndiagR[sy-1]=operatorL[N];
                    ndiagL[sy-2]=operatorL[S];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
                        -operatorL[E]*u[sx+jX[E]+sy*(nx+1)];
                    for(Index i=5; i<operatorL.size(); i++)
                    {
                        rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                    }
                }
                operatorL=stencil.getL(N,sx,ny-1,nx,ny);
                diagR[ny-2]=operatorL[C];
                ndiagL[ny-3]=operatorL[S];
                rhs[sx+(ny-1)*(nx+1)]=f[sx+(ny-1)*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+(ny-1)*(nx+1)] 
                    -operatorL[E]*u[sx+jX[E]+(ny-1)*(nx+1)]
                    -operatorL[N]*u[sx+(ny-1+jY[N])*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+(ny-1)*(nx+1)]-=
                                   operatorL[i]*u[sx+jX[i]+(ny-1+jY[i])*(nx+1)];
                }
                yLRSolver(u,sx,nx,ny,rhs,ndiagL,diagR,ndiagL);
            }
            operatorL=stencil.getL(SE,nx-1,1,nx,ny);
            diagR[0]=operatorL[C];
            ndiagR[0]=operatorL[N];
            rhs[nx-1+(nx+1)]=f[nx-1+(nx+1)]
                -operatorL[W]*u[nx-1+jX[W]+(nx+1)]
                -operatorL[E]*u[nx-1+jX[E]+(nx+1)] 
                -operatorL[S]*u[nx-1+(1+jY[S])*(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[nx-1+nx+1]-=operatorL[i]*u[nx-1+jX[i]+(1+jY[i])*(nx+1)];
            }
            for(Index sy=2; sy<ny-1; sy++) 
            {
                operatorL=stencil.getL(E,nx-1,sy,nx,ny);
                diagR[sy-1]=operatorL[C];
                ndiagR[sy-1]=operatorL[N];
                ndiagL[sy-2]=operatorL[S];
                rhs[nx-1+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorL[W]*u[nx-1+jX[W]+sy*(nx+1)]
                    -operatorL[E]*u[nx-1+jX[E]+sy*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                   operatorL[i]*u[nx-1+jX[i]+(sy+jY[i])*(nx+1)];
                }
            }
            operatorL=stencil.getL(NE,nx-1,ny-1,nx,ny);
            diagR[ny-2]=operatorL[C];
            ndiagL[ny-3]=operatorL[S];
            rhs[nx-1+(ny-1)*(nx+1)]=f[nx-1+(ny-1)*(nx+1)]
                -operatorL[W]*u[nx-1+jX[W]+(ny-1)*(nx+1)] 
                -operatorL[E]*u[nx-1+jX[E]+(ny-1)*(nx+1)]
                -operatorL[N]*u[nx-1+(ny-1+jY[N])*(nx+1)];
            for(Index i=5; i<operatorL.size(); i++)
            {
                rhs[nx-1+(ny-1)*(nx+1)]-=
                                operatorL[i]*u[nx-1+jX[i]+(ny-1+jY[i])*(nx+1)];
            }
            yLRSolver(u,nx-1,nx,ny,rhs,ndiagL,diagR,ndiagR);
            for(Index sx=2; sx<nx ; sx+=2)
            {
                operatorL=stencil.getL(S,sx,1,nx,ny);
                diagR[0]=operatorL[C];
                ndiagR[0]=operatorL[N];
                rhs[sx+(nx+1)]=f[sx+(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+(nx+1)] 
                    -operatorL[S]*u[sx+(1+jY[S])*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+nx+1]-=operatorL[i]*u[sx+jX[i]+(1+jY[i])*(nx+1)];
                }
                for(Index sy=2; sy<ny-1; sy++) 
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    diagR[sy-1]=operatorL[C];
                    ndiagR[sy-1]=operatorL[N];                                     
                    ndiagL[sy-2]=operatorL[S];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
                        -operatorL[E]*u[sx+jX[E]+sy*(nx+1)]; 
                    for(Index i=5; i<operatorL.size(); i++)
                    {
                        rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                    }
                }
                operatorL=stencil.getL(N,sx,ny-1,nx,ny);
                diagR[nx-2]=operatorL[C];
                ndiagL[nx-3]=operatorL[S];
                rhs[sx+(ny-1)*(nx+1)]=f[sx+(ny-1)*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+(ny-1)*(nx+1)] 
                    -operatorL[E]*u[sx+jX[E]+(ny-1)*(nx+1)]
                    -operatorL[N]*u[sx+(ny-1+jY[N])*(nx+1)];
                for(Index i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+(ny-1)*(nx+1)]-=
                                   operatorL[i]*u[sx+jX[i]+(ny-1+jY[i])*(nx+1)];
                }
                yLRSolver(u,sx,nx,ny,rhs,ndiagL,diagR,ndiagR);
            }
        }
        else //nx,ny to small
        {
            gsRedBlack_.relax(u,f,stencil,nx,ny);
        }
    }                       
}
void ZebraLineGS::xzebra(
    NumericArray &u,
    const NumericArray &f, 
    NumericArray &rhs,
    const Stencil &stencil,
    const Index nx, 
    const Index ny) const
{
    if((ny > 4) && (nx > 4))
    {   
        NumericArray diagR(0.0,nx-1);
        NumericArray ndiagR1(0.0,nx-2);
        NumericArray ndiagL1(0.0,nx-2);
        NumericArray ndiagR2(0.0,nx-3);
        NumericArray ndiagL2(0.0,nx-3);
        if(stencil.isConstant())
        {
            // get const operator L
            const NumericArray operatorL=stencil.getL(C,2,2,nx,ny);
            const PositionArray jX=stencil.getJx(C,nx,ny);
            const PositionArray jY=stencil.getJy(C,nx,ny);
            NumericArray operatorLB=stencil.getL(S,2,1,nx,ny);
            PositionArray jXB=stencil.getJx(S,nx,ny);
            PositionArray jYB=stencil.getJy(S,nx,ny);        
            NumericArray operatorLC=stencil.getL(SW,1,1,nx,ny);
            PositionArray jXC=stencil.getJx(SW,nx,ny);
            PositionArray jYC=stencil.getJy(SW,nx,ny);
            // set right hand side for line 1                   
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[E];
            ndiagR2[0]=operatorLC[NE];
            rhs[1+nx+1]=f[1+nx+1]
                -operatorLC[N]*u[1+(1+jYC[N])*(nx+1)]
                -operatorLC[S]*u[1+(1+jYC[S])*(nx+1)]
                -operatorLC[W]*u[nx+1]
                -operatorLC[NW]*u[1+(1+jYC[NW])*(nx+1)];  
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[1+nx+1]-=operatorLC[i]*u[1+jXC[i]+(1+jYC[i])*(nx+1)];
            }
            ndiagL1[0]=operatorLB[W];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[E];                    
            ndiagR2[1]=operatorLB[SE];
            rhs[2+nx+1]=f[2+nx+1]
                -operatorLB[N]*u[2+(1+jYB[N])*(nx+1)]
                -operatorLB[S]*u[2+(1+jYB[S])*(nx+1)]
                -operatorLB[NE]*u[2+(1+jYB[NE])*(nx+1)]
                -operatorLB[NW]*u[2+jXB[NW]+nx+1];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[2+nx+1]-=operatorLB[i]*u[2+jXB[i]+(1+jYB[i])*(nx+1)];
            }
            for(Index sx=3; sx<nx-2; sx++)  
            {
                ndiagL2[sx-3]=operatorLB[NW];
                ndiagL1[sx-2]=operatorLB[W];
                diagR[sx-1]=operatorLB[C];
                ndiagR1[sx-1]=operatorLB[E];                  
                ndiagR2[sx-1]=operatorLB[SE];
                rhs[sx+nx+1]=f[sx+nx+1]
                    -operatorLB[N]*u[sx+(1+jYB[N])*(nx+1)]
                    -operatorLB[S]*u[sx+(1+jYB[S])*(nx+1)]
                    -operatorLB[NE]*u[sx+(1+jYB[NE])*(nx+1)];
                for(Index i=8; i<operatorLB.size(); i++)
                {
                    rhs[sx+nx+1]-=operatorLB[i]*u[sx+jXB[i]+(1+jYB[i])*(nx+1)];
                }
            }
            ndiagL2[nx-5]=operatorLB[NW];
            ndiagL1[nx-4]=operatorLB[W];
            diagR[nx-3]=operatorLB[C];
            ndiagR1[nx-3]=operatorLB[E];
            rhs[nx-2+nx+1]=f[nx-2+nx+1]
                -operatorLB[N]*u[nx-2+(1+jYB[N])*(nx+1)]
                -operatorLB[S]*u[nx-2+(1+jYB[S])*(nx+1)]
                -operatorLB[NE]*u[nx-2+(1+jYB[NE])*(nx+1)]
                -operatorLB[SE]*u[nx-2+jXB[SE]+nx+1];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[nx-2+nx+1]-=operatorLB[i]*u[nx-2+jXB[i]+(1+jYB[i])*(nx+1)];
            }
            operatorLC=stencil.getL(SE,nx-1,1,nx,ny);
            jXC=stencil.getJx(SE,nx,ny);
            jYC=stencil.getJy(SE,nx,ny);
            ndiagL2[nx-4]=operatorLC[NW];
            ndiagL1[nx-3]=operatorLC[W];
            diagR[nx-2]=operatorLC[C];
            rhs[(nx-1)+(nx+1)]=f[nx-1+nx+1]
                -operatorLC[N]*u[nx-1+(1+jYC[N])*(nx+1)] 
                -operatorLC[S]*u[nx-1+(1+jYC[S])*(nx+1)]
                -operatorLC[E]*u[nx+nx+1]
                -operatorLC[NE]*u[nx-1+(1+jYC[NE])*(nx+1)];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[nx-1+(nx+1)]-=
                                operatorLC[i]*u[nx-1+jXC[i]+(1+jYC[i])*(nx+1)];
            }
            xLRSolver(u,1,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            for(Index sy=3; sy<ny-2; sy+=2)
            {                   
                operatorLB=stencil.getL(W,1,sy,nx,ny);
                jXB=stencil.getJx(W,nx,ny);
                jYB=stencil.getJy(W,nx,ny);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[E];
                ndiagR2[0]=operatorLB[NE];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorLB[N]*u[1+(sy+jYB[N])*(nx+1)]
                    -operatorLB[S]*u[1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[W]*u[sy*(nx+1)]
                    -operatorLB[NW]*u[1+(sy+jYB[NW])*(nx+1)]
                    -operatorLB[SE]*u[1+(sy+jYB[SE])*(nx+1)];
                for(Index i=8; i<operatorLB.size(); i++)
                {
                    rhs[1+sy*(nx+1)]-=
                                   operatorLB[i]*u[1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
                ndiagL1[0]=operatorL[W];
                diagR[1]=operatorL[C];
                ndiagR1[1]=operatorL[E];                  
                ndiagR2[1]=operatorL[SE];
                rhs[2+sy*(nx+1)]=f[2+sy*(nx+1)]
                    -operatorL[N]*u[2+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[2+(sy+jY[S])*(nx+1)]
                    -operatorL[NE]*u[2+(sy+jY[NE])*(nx+1)]
                    -operatorL[SW]*u[2+(sy+jY[SW])*(nx+1)]
                    -operatorL[NW]*u[2+jX[NW]+sy*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[2+sy*(nx+1)]-=operatorL[i]*u[2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                for(Index sx=3; sx<nx-2; sx++)  
                {
                    ndiagL2[sx-3]=operatorL[NW];
                    ndiagL1[sx-2]=operatorL[W];
                    diagR[sx-1]=operatorL[C];
                    ndiagR1[sx-1]=operatorL[E];                    
                    ndiagR2[sx-1]=operatorL[SE];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
                        -operatorL[S]*u[sx+(sy+jY[S])*(nx+1)]
                        -operatorL[NE]*u[sx+(sy+jY[NE])*(nx+1)]
                        -operatorL[SW]*u[sx+(sy+jY[SW])*(nx+1)];
                    for(Index i=9; i<operatorL.size(); i++)
                    {
                        rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                    }
                }
                ndiagL2[nx-5]=operatorL[NW];
                ndiagL1[nx-4]=operatorL[W];
                diagR[nx-3]=operatorL[C];
                ndiagR1[nx-3]=operatorL[E];
                rhs[nx-2+sy*(nx+1)]=f[nx-2+sy*(nx+1)]
                    -operatorL[N]*u[nx-2+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[nx-2+(sy+jY[S])*(nx+1)]
                    -operatorL[NE]*u[nx-2+(sy+jY[NE])*(nx+1)]
                    -operatorL[SW]*u[nx-2+(sy+jY[SW])*(nx+1)]
                    -operatorL[SE]*u[nx-2+jX[SE]+sy*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[nx-2+sy*(nx+1)]-=
                                   operatorL[i]*u[nx-2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                operatorLB=stencil.getL(E,nx-1,sy,nx,ny);
                jXB=stencil.getJx(E,nx,ny);
                jYB=stencil.getJy(E,nx,ny);
                ndiagL2[nx-4]=operatorLB[NW];
                ndiagL1[nx-3]=operatorLB[W];
                diagR[nx-2]=operatorLB[C];
                rhs[(nx-1)+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorLB[N]*u[nx-1+(sy+jYB[N])*(nx+1)] 
                    -operatorLB[S]*u[nx-1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[E]*u[nx+sy*(nx+1)]
                    -operatorLB[NE]*u[nx-1+(sy+jYB[NE])*(nx+1)]
                    -operatorLB[SE]*u[nx-1+(sy+jYB[SE])*(nx+1)];
                for(Index i=9; i<operatorLB.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                operatorLB[i]*u[nx-1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
                xLRSolver(u,sy,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);    
            }
            //relax top line Zeile
            //set rhs in top line
            operatorLC=stencil.getL(NW,1,ny-1,nx,ny);
            jXC=stencil.getJx(NW,nx,ny);
            jYC=stencil.getJy(NW,nx,ny);
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[E];
            ndiagR2[0]=operatorLC[NW];
            rhs[1+(ny-1)*(nx+1)]=f[1+(ny-1)*(nx+1)]
                -operatorLC[N]*u[1+(ny-1+jYC[N])*(nx+1)]
                -operatorLC[S]*u[1+(ny-1+jYC[S])*(nx+1)]
                -operatorLC[W]*u[(ny-1)*(nx+1)]
                -operatorLC[NE]*u[1+(ny-1+jYC[NE])*(nx+1)];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[1+(ny-1)*(nx+1)]-=
                                operatorLC[i]*u[1+jXC[i]+(ny-1+jYC[i])*(nx+1)];
            }
            operatorLB=stencil.getL(N,2,ny-1,nx,ny);
            jXB=stencil.getJx(N,nx,ny);
            jYB=stencil.getJy(N,nx,ny);
            ndiagL1[0]=operatorLB[W];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[E];                    
            ndiagR2[1]=operatorLB[NE];
            rhs[2+(ny-1)*(nx+1)]=f[2+(ny-1)*(nx+1)]
                -operatorLB[N]*u[2+(ny-1+jYB[N])*(nx+1)]
                -operatorLB[S]*u[2+(ny-1+jYB[S])*(nx+1)]
                -operatorLB[SE]*u[2+(ny-1+jYB[SE])*(nx+1)]
                -operatorLB[NW]*u[2+jXB[NW]+(ny-1)*(nx+1)];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[2+(ny-1)*(nx+1)]-=
                                operatorLB[i]*u[2+jXB[i]+(ny-1+jYB[i])*(nx+1)];
            }
            for(Index sx=3; sx<nx-2; sx++)  
            {
                ndiagL2[sx-3]=operatorLB[NW];
                ndiagL1[sx-2]=operatorLB[W];
                diagR[sx-1]=operatorLB[C];
                ndiagR1[sx-1]=operatorLB[E];                  
                ndiagR2[sx-1]=operatorLB[NE];
                rhs[sx+(ny-1)*(nx+1)]=f[sx+(ny-1)*(nx+1)]
                    -operatorLB[N]*u[sx+(ny-1+jYB[N])*(nx+1)]
                    -operatorLB[S]*u[sx+(ny-1+jYB[S])*(nx+1)]
                    -operatorLB[SE]*u[sx+(ny-1+jYB[SE])*(nx+1)];
                for(Index i=8; i<operatorLB.size(); i++)
                {
                    rhs[sx+(ny-1)*(nx+1)]-=
                                operatorLB[i]*u[sx+jXB[i]+(ny-1+jYB[i])*(nx+1)];
                }
            }
            ndiagL2[nx-5]=operatorLB[NW];
            ndiagL1[nx-4]=operatorLB[W];
            diagR[nx-3]=operatorLB[C];
            ndiagR1[nx-3]=operatorLB[E];
            rhs[nx-2+(ny-1)*(nx+1)]=f[nx-2+(ny-1)*(nx+1)]
                -operatorLB[N]*u[nx-2+(ny-1+jYB[N])*(nx+1)]
                -operatorLB[S]*u[nx-2+(ny-1+jYB[S])*(nx+1)]
                -operatorLB[SE]*u[nx-2+(ny-1+jYB[SE])*(nx+1)]
                -operatorLB[NE]*u[nx-2+jXB[NE]+(ny-1)*(nx+1)];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[nx-2+(ny-1)*(nx+1)]-=
                            operatorLB[i]*u[nx-2+jXB[i]+(ny-1+jYB[i])*(nx+1)];
            }
            operatorLC=stencil.getL(NE,nx-1,ny-1,nx,ny);
            jXC=stencil.getJx(NE,nx,ny);
            jYC=stencil.getJy(NE,nx,ny);
            ndiagL2[nx-4]=operatorLC[NW];
            ndiagL1[nx-3]=operatorLC[W];
            diagR[nx-2]=operatorLC[C];
            rhs[(nx-1)+(ny-1)*(nx+1)]=f[nx-1+(ny-1)*(nx+1)]
                -operatorLC[N]*u[nx-1+(ny-1+jYC[N])*(nx+1)]
                -operatorLC[S]*u[nx-1+(ny-1+jYC[S])*(nx+1)]
                -operatorLC[E]*u[nx+(ny-1)*(nx+1)]
                -operatorLC[NE]*u[nx-1+(ny-1+jYC[NE])*(nx+1)];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[nx-1+(ny-1)*(nx+1)]-=
                            operatorLC[i]*u[nx-1+jXC[i]+(ny-1+jYC[i])*(nx+1)];
            }
            xLRSolver(u,ny-1,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            // relax even inner lines
            for(Index sy=2; sy < ny-1; sy+=2)
            {
                // set rhs           
                operatorLB=stencil.getL(W,1,sy,nx,ny);
                jXB=stencil.getJx(W,nx,ny);
                jYB=stencil.getJy(W,nx,ny);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[E];
                ndiagR2[0]=operatorLB[NE];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorLB[N]*u[1+(sy+jYB[N])*(nx+1)]
                    -operatorLB[S]*u[1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[W]*u[sy*(nx+1)]
                    -operatorLB[NW]*u[1+(sy+jYB[NW])*(nx+1)]
                    -operatorLB[SE]*u[1+(sy+jYB[SE])*(nx+1)];
                for(Index i=8; i<operatorLB.size(); i++)
                {
                    rhs[1+sy*(nx+1)]-=
                                operatorLB[i]*u[1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
                ndiagL1[0]=operatorL[W];
                diagR[1]=operatorL[C];
                ndiagR1[1]=operatorL[E];                  
                ndiagR2[1]=operatorL[SE];
                rhs[2+sy*(nx+1)]=f[2+sy*(nx+1)]
                    -operatorL[N]*u[2+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[2+(sy+jY[S])*(nx+1)]
                    -operatorL[NE]*u[2+(sy+jY[NE])*(nx+1)]
                    -operatorL[SW]*u[2+(sy+jY[SW])*(nx+1)]
                    -operatorL[NW]*u[2+jX[NW]+sy*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[2+sy*(nx+1)]-=operatorL[i]*u[2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                for(Index sx=3; sx<nx-2; sx++)  
                {
                    ndiagL2[sx-3]=operatorL[NW];
                    ndiagL1[sx-2]=operatorL[W];
                    diagR[sx-1]=operatorL[C];
                    ndiagR1[sx-1]=operatorL[E];                    
                    ndiagR2[sx-1]=operatorL[SE];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
                        -operatorL[S]*u[sx+(sy+jY[S])*(nx+1)]
                        -operatorL[NE]*u[sx+(sy+jY[NE])*(nx+1)]
                        -operatorL[SW]*u[sx+(sy+jY[SW])*(nx+1)];
                    for(Index i=9; i<operatorL.size(); i++)
                    {
                        rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                    }
                }
                ndiagL2[nx-5]=operatorL[NW];
                ndiagL1[nx-4]=operatorL[W];
                diagR[nx-3]=operatorL[C];
                ndiagR1[nx-3]=operatorL[E];
                rhs[nx-2+sy*(nx+1)]=f[nx-2+sy*(nx+1)]
                    -operatorL[N]*u[nx-2+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[nx-2+(sy+jY[S])*(nx+1)]
                    -operatorL[NE]*u[nx-2+(sy+jY[NE])*(nx+1)]
                    -operatorL[SW]*u[nx-2+(sy+jY[SW])*(nx+1)]
                    -operatorL[SE]*u[nx-2+jX[SE]+sy*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[nx-2+sy*(nx+1)]-=
                                operatorL[i]*u[nx-2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                operatorLB=stencil.getL(E,nx-1,sy,nx,ny);
                jXB=stencil.getJx(E,nx,ny);
                jYB=stencil.getJy(E,nx,ny);
                ndiagL2[nx-4]=operatorLB[NW];
                ndiagL1[nx-3]=operatorLB[W];
                diagR[nx-2]=operatorLB[C];
                rhs[(nx-1)+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorLB[N]*u[nx-1+(sy+jYB[N])*(nx+1)] 
                    -operatorLB[S]*u[nx-1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[E]*u[nx+sy*(nx+1)]
                    -operatorLB[NE]*u[nx-1+(sy+jYB[NE])*(nx+1)]
                    -operatorLB[SE]*u[nx-1+(sy+jYB[SE])*(nx+1)];
                for(Index i=9; i<operatorLB.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                operatorLB[i]*u[nx-1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
                xLRSolver(u,sy,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);    
            }
        }
        else // stencil not constant
        {
            NumericArray operatorL=stencil.getL(C,2,2,nx,ny);
            PositionArray jX=stencil.getJx(C,nx,ny);
            PositionArray jY=stencil.getJy(C,nx,ny);
            NumericArray operatorLB=stencil.getL(S,2,1,nx,ny);
            PositionArray jXB=stencil.getJx(S,nx,ny);
            PositionArray jYB=stencil.getJy(S,nx,ny);
            NumericArray operatorLC=stencil.getL(SW,1,1,nx,ny);
            PositionArray jXC=stencil.getJx(SW,nx,ny);
            PositionArray jYC=stencil.getJy(SW,nx,ny);
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[E];
            ndiagR2[0]=operatorLC[NE];
            rhs[1+nx+1]=f[1+nx+1]
                -operatorLC[N]*u[1+(1+jYC[N])*(nx+1)]
                -operatorLC[S]*u[1+(1+jYC[S])*(nx+1)]
                -operatorLC[W]*u[nx+1]
                -operatorLC[NW]*u[1+(1+jYC[NW])*(nx+1)];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[1+nx+1]-=operatorLC[i]*u[1+jXC[i]+(1+jYC[i])*(nx+1)];
            }
            operatorLB=stencil.getL(S,2,1,nx,ny);
            jXB=stencil.getJx(S,nx,ny);
            jYB=stencil.getJy(S,nx,ny);
            ndiagL1[0]=operatorLB[W];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[E];                    
            ndiagR2[1]=operatorLB[SE];
            rhs[2+nx+1]=f[2+nx+1]
                -operatorLB[N]*u[2+(1+jYB[N])*(nx+1)]
                -operatorLB[S]*u[2+(1+jYB[S])*(nx+1)]
                -operatorLB[NE]*u[2+(1+jYB[NE])*(nx+1)]
                -operatorLB[NW]*u[2+jXB[NW]+nx+1];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[2+nx+1]-=operatorLB[i]*u[2+jXB[i]+(1+jYB[i])*(nx+1)];
            }
            for(Index sx=3; sx<nx-2; sx++)  
            {
                operatorLB=stencil.getL(S,sx,1,nx,ny);
                ndiagL2[sx-3]=operatorLB[NW];
                ndiagL1[sx-2]=operatorLB[W];
                diagR[sx-1]=operatorLB[C];
                ndiagR1[sx-1]=operatorLB[E];                  
                ndiagR2[sx-1]=operatorLB[SE];
                rhs[sx+nx+1]=f[sx+nx+1]
                    -operatorLB[N]*u[sx+(1+jYB[N])*(nx+1)]
                    -operatorLB[S]*u[sx+(1+jYB[S])*(nx+1)]
                    -operatorLB[NE]*u[sx+(1+jYB[NE])*(nx+1)];
                for(Index i=8; i<operatorLB.size(); i++)
                {
                    rhs[sx+nx+1]-=operatorLB[i]*u[sx+jXB[i]+(1+jYB[i])*(nx+1)];
                }
            }
            operatorLB=stencil.getL(S,nx-2,1,nx,ny);
            ndiagL2[nx-5]=operatorLB[NW];
            ndiagL1[nx-4]=operatorLB[W];
            diagR[nx-3]=operatorLB[C];
            ndiagR1[nx-3]=operatorLB[E];
            rhs[nx-2+nx+1]=f[nx-2+nx+1]
                -operatorLB[N]*u[nx-2+(1+jYB[N])*(nx+1)]
                -operatorLB[S]*u[nx-2+(1+jYB[S])*(nx+1)]
                -operatorLB[NE]*u[nx-2+(1+jYB[NE])*(nx+1)]
                -operatorLB[SE]*u[nx-2+jXB[SE]+nx+1];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[nx-2+nx+1]-=operatorLB[i]*u[nx-2+jXB[i]+(1+jYB[i])*(nx+1)];
            }
            operatorLC=stencil.getL(SE,nx-1,1,nx,ny);
            jXC=stencil.getJx(SE,nx,ny);
            jYC=stencil.getJy(SE,nx,ny);
            ndiagL2[nx-4]=operatorLC[NW];
            ndiagL1[nx-3]=operatorLC[W];
            diagR[nx-2]=operatorLC[C];
            rhs[(nx-1)+(nx+1)]=f[nx-1+nx+1]
                -operatorLC[N]*u[nx-1+(1+jYC[N])*(nx+1)] 
                -operatorLC[S]*u[nx-1+(1+jYC[S])*(nx+1)]
                -operatorLC[E]*u[nx+nx+1]
                -operatorLC[NE]*u[nx-1+(1+jYC[NE])*(nx+1)];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[nx-1+(nx+1)]-=
                                operatorLC[i]*u[nx-1+jXC[i]+(1+jYC[i])*(nx+1)];
            }
            xLRSolver(u,1,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            // do odd inner lines
            for(Index sy=3; sy < ny-2; sy+=2)
            {
                // set rhs                   
                operatorLB=stencil.getL(W,1,sy,nx,ny);
                jXB=stencil.getJx(W,nx,ny);
                jYB=stencil.getJy(W,nx,ny);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[E];
                ndiagR2[0]=operatorLB[NE];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorLB[N]*u[1+(sy+jYB[N])*(nx+1)]
                    -operatorLB[S]*u[1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[W]*u[sy*(nx+1)]
                    -operatorLB[NW]*u[1+(sy+jYB[NW])*(nx+1)]
                    -operatorLB[SE]*u[1+(sy+jYB[SE])*(nx+1)];
                for(Index i=8; i<operatorLB.size(); i++)
                {
                    rhs[1+sy*(nx+1)]-=
                                operatorLB[i]*u[1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
                operatorL=stencil.getL(C,2,sy,nx,ny);
                ndiagL1[0]=operatorL[W];
                diagR[1]=operatorL[C];
                ndiagR1[1]=operatorL[E];                  
                ndiagR2[1]=operatorL[SE];
                rhs[2+sy*(nx+1)]=f[2+sy*(nx+1)]
                    -operatorL[N]*u[2+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[2+(sy+jY[S])*(nx+1)]
                    -operatorL[NE]*u[2+(sy+jY[NE])*(nx+1)]
                    -operatorL[SW]*u[2+(sy+jY[SW])*(nx+1)]
                    -operatorL[NW]*u[2+jX[NW]+sy*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[2+sy*(nx+1)]-=operatorL[i]*u[2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                for(Index sx=3; sx<nx-2; sx++)  
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    ndiagL2[sx-3]=operatorL[NW];
                    ndiagL1[sx-2]=operatorL[W];
                    diagR[sx-1]=operatorL[C];
                    ndiagR1[sx-1]=operatorL[E];                    
                    ndiagR2[sx-1]=operatorL[SE];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
                        -operatorL[S]*u[sx+(sy+jY[S])*(nx+1)]
                        -operatorL[NE]*u[sx+(sy+jY[NE])*(nx+1)]
                        -operatorL[SW]*u[sx+(sy+jY[SW])*(nx+1)];
                    for(Index i=9; i<operatorL.size(); i++)
                    {
                        rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                    }
                }
                operatorL=stencil.getL(C,nx-2,sy,nx,ny);
                ndiagL2[nx-5]=operatorL[NW];
                ndiagL1[nx-4]=operatorL[W];
                diagR[nx-3]=operatorL[C];
                ndiagR1[nx-3]=operatorL[E];
                rhs[nx-2+sy*(nx+1)]=f[nx-2+sy*(nx+1)]
                    -operatorL[N]*u[nx-2+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[nx-2+(sy+jY[S])*(nx+1)]
                    -operatorL[NE]*u[nx-2+(sy+jY[NE])*(nx+1)]
                    -operatorL[SW]*u[nx-2+(sy+jY[SW])*(nx+1)]
                    -operatorL[SE]*u[nx-2+jX[SE]+sy*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[nx-2+sy*(nx+1)]-=
                                operatorL[i]*u[nx-2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                operatorLB=stencil.getL(E,nx-1,sy,nx,ny);
                jXB=stencil.getJx(E,nx,ny);
                jYB=stencil.getJy(E,nx,ny);
                ndiagL2[nx-4]=operatorLB[NW];
                ndiagL1[nx-3]=operatorLB[W];
                diagR[nx-2]=operatorLB[C];
                rhs[(nx-1)+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorLB[N]*u[nx-1+(sy+jYB[N])*(nx+1)] 
                    -operatorLB[S]*u[nx-1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[E]*u[nx+sy*(nx+1)]
                    -operatorLB[NE]*u[nx-1+(sy+jYB[NE])*(nx+1)]
                    -operatorLB[SE]*u[nx-1+(sy+jYB[SE])*(nx+1)];
                for(Index i=9; i<operatorLB.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                operatorLB[i]*u[nx-1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
                xLRSolver(u,sy,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            }
            //relax top line
            // set rhs in top line                
            operatorLC=stencil.getL(NW,1,ny-1,nx,ny);
            jXC=stencil.getJx(NW,nx,ny);
            jYC=stencil.getJy(NW,nx,ny);
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[E];
            ndiagR2[0]=operatorLC[NW];
            rhs[1+(ny-1)*(nx+1)]=f[1+(ny-1)*(nx+1)]
                -operatorLC[N]*u[1+(ny-1+jYC[N])*(nx+1)]
                -operatorLC[S]*u[1+(ny-1+jYC[S])*(nx+1)]
                -operatorLC[W]*u[(ny-1)*(nx+1)]
                -operatorLC[NE]*u[1+(ny-1+jYC[NE])*(nx+1)];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[1+(ny-1)*(nx+1)]-=
                            operatorLC[i]*u[1+jXC[i]+(ny-1+jYC[i])*(nx+1)];
            }
            operatorLB=stencil.getL(N,2,ny-1,nx,ny);
            jXB=stencil.getJx(N,nx,ny);
            jYB=stencil.getJy(N,nx,ny);
            ndiagL1[0]=operatorLB[W];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[E];                    
            ndiagR2[1]=operatorLB[NE];
            rhs[2+(ny-1)*(nx+1)]=f[2+(ny-1)*(nx+1)]
                -operatorLB[N]*u[2+(ny-1+jYB[N])*(nx+1)]
                -operatorLB[S]*u[2+(ny-1+jYB[S])*(nx+1)]
                -operatorLB[SE]*u[2+(ny-1+jYB[SE])*(nx+1)]
                -operatorLB[NW]*u[2+jXB[NW]+(ny-1)*(nx+1)];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[2+(ny-1)*(nx+1)]-=
                                operatorLB[i]*u[2+jXB[i]+(ny-1+jYB[i])*(nx+1)];
            }
            for(Index sx=3; sx<nx-2; sx++)  
            {
                operatorLB=stencil.getL(N,sx,ny-1,nx,ny);
                ndiagL2[sx-3]=operatorLB[NW];
                ndiagL1[sx-2]=operatorLB[W];
                diagR[sx-1]=operatorLB[C];
                ndiagR1[sx-1]=operatorLB[E];                  
                ndiagR2[sx-1]=operatorLB[NE];
                rhs[sx+(ny-1)*(nx+1)]=f[sx+(ny-1)*(nx+1)]
                    -operatorLB[N]*u[sx+(ny-1+jYB[N])*(nx+1)]
                    -operatorLB[S]*u[sx+(ny-1+jYB[S])*(nx+1)]
                    -operatorLB[SE]*u[sx+(ny-1+jYB[SE])*(nx+1)];
                for(Index i=8; i<operatorLB.size(); i++)
                {
                    rhs[sx+(ny-1)*(nx+1)]-=
                                operatorLB[i]*u[sx+jXB[i]+(ny-1+jYB[i])*(nx+1)];
                }
            }
            operatorLB=stencil.getL(N,nx-2,ny-1,nx,ny);
            ndiagL2[nx-5]=operatorLB[NW];
            ndiagL1[nx-4]=operatorLB[W];
            diagR[nx-3]=operatorLB[C];
            ndiagR1[nx-3]=operatorLB[E];
            rhs[nx-2+(ny-1)*(nx+1)]=f[nx-2+(ny-1)*(nx+1)]
                -operatorLB[N]*u[nx-2+(ny-1+jYB[N])*(nx+1)]
                -operatorLB[S]*u[nx-2+(ny-1+jYB[S])*(nx+1)]
                -operatorLB[SE]*u[nx-2+(ny-1+jYB[SE])*(nx+1)]
                -operatorLB[NE]*u[nx-2+jXB[NE]+(ny-1)*(nx+1)];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[nx-2+(ny-1)*(nx+1)]-=
                            operatorLB[i]*u[nx-2+jXB[i]+(ny-1+jYB[i])*(nx+1)];
            }
            operatorLC=stencil.getL(NE,nx-1,ny-1,nx,ny);
            jXC=stencil.getJx(NE,nx,ny);
            jYC=stencil.getJy(NE,nx,ny);
            ndiagL2[nx-4]=operatorLC[NW];
            ndiagL1[nx-3]=operatorLC[W];
            diagR[nx-2]=operatorLC[C];
            rhs[(nx-1)+(ny-1)*(nx+1)]=f[nx-1+(ny-1)*(nx+1)]
                -operatorLC[N]*u[nx-1+(ny-1+jYC[N])*(nx+1)]
                -operatorLC[S]*u[nx-1+(ny-1+jYC[S])*(nx+1)]
                -operatorLC[E]*u[nx+(ny-1)*(nx+1)]
                -operatorLC[NE]*u[nx-1+(ny-1+jYC[NE])*(nx+1)];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[nx-1+(ny-1)*(nx+1)]-=
                            operatorLC[i]*u[nx-1+jXC[i]+(ny-1+jYC[i])*(nx+1)];
            }
            xLRSolver(u,ny-1,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            // relax even inner lines
            for(Index sy=2; sy<ny-1; sy+=2)
            {
                // set rhs                   
                operatorLB=stencil.getL(W,1,sy,nx,ny);
                jXB=stencil.getJx(W,nx,ny);
                jYB=stencil.getJy(W,nx,ny);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[E];
                ndiagR2[0]=operatorLB[NE];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorLB[N]*u[1+(sy+jYB[N])*(nx+1)]
                    -operatorLB[S]*u[1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[W]*u[sy*(nx+1)]
                    -operatorLB[NW]*u[1+(sy+jYB[NW])*(nx+1)]
                    -operatorLB[SE]*u[1+(sy+jYB[SE])*(nx+1)];
                for(Index i=8; i<operatorLB.size(); i++)
                {
                    rhs[1+sy*(nx+1)]-=
                                operatorLB[i]*u[1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
                operatorL=stencil.getL(C,2,sy,nx,ny);
                ndiagL1[0]=operatorL[W];
                diagR[1]=operatorL[C];
                ndiagR1[1]=operatorL[E];                  
                ndiagR2[1]=operatorL[SE];
                rhs[2+sy*(nx+1)]=f[2+sy*(nx+1)]
                    -operatorL[N]*u[2+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[2+(sy+jY[S])*(nx+1)]
                    -operatorL[NE]*u[2+(sy+jY[NE])*(nx+1)]
                    -operatorL[SW]*u[2+(sy+jY[SW])*(nx+1)]
                    -operatorL[NW]*u[2+jX[NW]+sy*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[2+sy*(nx+1)]-=operatorL[i]*u[2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                for(Index sx=3; sx<nx-2; sx++)  
                { 
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    ndiagL2[sx-3]=operatorL[NW];
                    ndiagL1[sx-2]=operatorL[W];
                    diagR[sx-1]=operatorL[C];
                    ndiagR1[sx-1]=operatorL[E];                    
                    ndiagR2[sx-1]=operatorL[SE];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
                        -operatorL[S]*u[sx+(sy+jY[S])*(nx+1)]
                        -operatorL[NE]*u[sx+(sy+jY[NE])*(nx+1)]
                        -operatorL[SW]*u[sx+(sy+jY[SW])*(nx+1)];
                    for(Index i=9; i<operatorL.size(); i++)
                    {
                        rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                    }
                }
                operatorL=stencil.getL(C,nx-2,sy,nx,ny);
                ndiagL2[nx-5]=operatorL[NW];
                ndiagL1[nx-4]=operatorL[W];
                diagR[nx-3]=operatorL[C];
                ndiagR1[nx-3]=operatorL[E];
                rhs[nx-2+sy*(nx+1)]=f[nx-2+sy*(nx+1)]
                    -operatorL[N]*u[nx-2+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[nx-2+(sy+jY[S])*(nx+1)]
                    -operatorL[NE]*u[nx-2+(sy+jY[NE])*(nx+1)]
                    -operatorL[SW]*u[nx-2+(sy+jY[SW])*(nx+1)]
                    -operatorL[SE]*u[nx-2+jX[SE]+sy*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[nx-2+sy*(nx+1)]-=
                                operatorL[i]*u[nx-2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                operatorLB=stencil.getL(E,nx-1,sy,nx,ny);
                jXB=stencil.getJx(E,nx,ny);
                jYB=stencil.getJy(E,nx,ny);
                ndiagL2[nx-4]=operatorLB[NW];
                ndiagL1[nx-3]=operatorLB[W];
                diagR[nx-2]=operatorLB[C];
                rhs[(nx-1)+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorLB[N]*u[nx-1+(sy+jYB[N])*(nx+1)] 
                    -operatorLB[S]*u[nx-1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[E]*u[nx+sy*(nx+1)]
                    -operatorLB[NE]*u[nx-1+(sy+jYB[NE])*(nx+1)]
                    -operatorLB[SE]*u[nx-1+(sy+jYB[SE])*(nx+1)];
                for(Index i=9; i<operatorLB.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                operatorLB[i]*u[nx-1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
                xLRSolver(u,sy,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);    
            }
        }
    }
    else //nx,ny to small
    {
        for(int i=0; i<2; i++)
        {
            gsRedBlack_.relax(u,f,stencil,nx,ny);
        }
    }
}
void ZebraLineGS::yzebra(
    NumericArray &u,
    const NumericArray &f, 
    NumericArray &rhs,
    const Stencil &stencil,
    const Index nx, 
    const Index ny) const
{
    if((ny > 4) && (nx > 4))
    {
        NumericArray diagR(0.0,ny-1);
        NumericArray ndiagR1(0.0,ny-2);
        NumericArray ndiagL1(0.0,ny-2);
        NumericArray ndiagR2(0.0,ny-3);
        NumericArray ndiagL2(0.0,ny-3);
        if(stencil.isConstant())
        {
            // get const operator L
            const NumericArray operatorL=stencil.getL(C,2,2,nx,ny);
            const PositionArray jX=stencil.getJx(C,nx,ny);
            const PositionArray jY=stencil.getJy(C,nx,ny);
            NumericArray operatorLB=stencil.getL(W,1,2,nx,ny);
            PositionArray jXB=stencil.getJx(W,nx,ny);
            PositionArray jYB=stencil.getJy(W,nx,ny);
            NumericArray operatorLC=stencil.getL(SW,1,1,nx,ny);
            PositionArray jXC=stencil.getJx(SW,nx,ny);
            PositionArray jYC=stencil.getJy(SW,nx,ny);
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[N];
            ndiagR2[0]=operatorLC[NW];
            rhs[1+nx+1]=f[1+nx+1]
                -operatorLC[S]*u[1+(1+jYC[S])*(nx+1)]
                -operatorLC[W]*u[nx+1]
                -operatorLC[E]*u[2+nx+1]
                -operatorLC[NE]*u[1+jXC[NE]+nx+1];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[1+nx+1]-=operatorLC[i]*u[1+jXC[i]+(1+jYC[i])*(nx+1)];
            }
            operatorLB=stencil.getL(W,1,2,nx,ny);
            jXB=stencil.getJx(W,nx,ny);
            jYB=stencil.getJy(W,nx,ny);
            ndiagL1[0]=operatorLB[S];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[N];                    
            ndiagR2[1]=operatorLB[NW];
            rhs[1+2*(nx+1)]=f[1+2*(nx+1)]
                -operatorLB[W]*u[2*(nx+1)]
                -operatorLB[E]*u[2+2*(nx+1)]
                -operatorLB[NE]*u[1+jXB[NE]+2*(nx+1)]
                -operatorLB[SE]*u[1+(2+jYB[SE])*(nx+1)];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[1+2*(nx+1)]-=operatorLB[i]*u[1+jXB[i]+(2+jYB[i])*(nx+1)];
            }
            for(Index sy=3; sy<ny-2; sy++)  
            {
                ndiagL2[sy-3]=operatorLB[SE];
                ndiagL1[sy-2]=operatorLB[S];
                diagR[sy-1]=operatorLB[C];
                ndiagR1[sy-1]=operatorLB[N];                  
                ndiagR2[sy-1]=operatorLB[NW];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorLB[W]*u[1+jXB[W]+sy*(nx+1)]
                    -operatorLB[E]*u[1+jXB[E]+sy*(nx+1)]
                    -operatorLB[NE]*u[1+jXB[NE]+sy*(nx+1)];
                for(Index i=8; i<operatorLB.size(); i++)
                {
                    rhs[1+sy*(nx+1)]-=
                                operatorLB[i]*u[1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
            }
            ndiagL2[ny-5]=operatorLB[SE];
            ndiagL1[ny-4]=operatorLB[S];
            diagR[ny-3]=operatorLB[C];
            ndiagR1[ny-3]=operatorLB[N];
            rhs[1+(ny-2)*(nx+1)]=f[1+(ny-2)*(nx+1)]
                -operatorLB[W]*u[(ny-2)*(nx+1)]
                -operatorLB[E]*u[2+(ny-2)*(nx+1)]
                -operatorLB[NE]*u[1+jXB[NE]+(ny-2)*(nx+1)]
                -operatorLB[NW]*u[1+(ny-2+jYB[NW])*(nx+1)];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[1+(ny-2)*(nx+1)]-=
                                operatorLB[i]*u[1+jXB[i]+(ny-2+jYB[i])*(nx+1)];
            }
            operatorLC=stencil.getL(NW,1,ny-1,nx,ny);
            jXC=stencil.getJx(NW,nx,ny);
            jYC=stencil.getJy(NW,nx,ny);
            ndiagL2[ny-4]=operatorLC[NE];
            ndiagL1[ny-3]=operatorLC[S];
            diagR[ny-2]=operatorLC[C];
            rhs[1+(ny-1)*(nx+1)]=f[1+(ny-1)*(nx+1)]
                -operatorLC[W]*u[(ny-1)*(nx+1)]
                -operatorLC[E]*u[2+(ny-1)*(nx+1)]
                -operatorLC[NW]*u[1+jXC[NW]+(ny-1)*(nx+1)]
                -operatorLC[N]*u[1+(ny-1+jYC[N])*(nx+1)];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[1+(ny-1)*(nx+1)]-=
                            operatorLC[i]*u[1+jXC[i]+(ny-1+jYC[i])*(nx+1)];
            }
            yLRSolver(u,1,nx,ny,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
////////////////////////////////////////////////////////////////////////////////
            // process all odd inner columns
            for(Index sx=3; sx < nx-2; sx+=2)
            {
                // set rhs
                operatorLB=stencil.getL(S,sx,1,nx,ny);
                jXB=stencil.getJx(S,nx,ny);
                jYB=stencil.getJy(S,nx,ny);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[N];
                ndiagR2[0]=operatorLB[NE];
                rhs[sx+nx+1]=f[sx+nx+1]
                    -operatorLB[S]*u[sx+(1+jYB[S])*(nx+1)]
                    -operatorLB[W]*u[sx+jXB[W]+nx+1] 
                    -operatorLB[E]*u[sx+jXB[E]+nx+1]
                    -operatorLB[NW]*u[sx+jXB[NW]+nx+1]
                    -operatorLB[SE]*u[sx+jXB[SE]+(nx+1)];
                for(Index i=7; i<operatorLB.size(); i++)
                {
                    rhs[sx+nx+1]-=
                               operatorLB[i]*u[sx+jXB[i]+(1+jYB[i])*(nx+1)];
                }
                ndiagL1[0]=operatorL[S];
                diagR[1]=operatorL[C];
                ndiagR1[1]=operatorL[N];                  
                ndiagR2[1]=operatorL[NE];
                rhs[sx+2*(nx+1)]=f[sx+2*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+2*(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+2*(nx+1)]
                    -operatorL[NW]*u[sx+jX[NW]+2*(nx+1)]
                    -operatorL[SE]*u[sx+jX[SE]+2*(nx+1)]
                    -operatorL[SW]*u[sx+(2+jY[SW])*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[sx+2*(nx+1)]-=
                                operatorL[i]*u[sx+jX[i]+(2+jY[i])*(nx+1)];
                }
                for(Index sy=3; sy<ny-2; sy++)
                {
                    ndiagL2[sy-3]=operatorL[SW];
                    ndiagL1[sy-2]=operatorL[S];
                    diagR[sy-1]=operatorL[C];
                    ndiagR1[sy-1]=operatorL[N];
                    ndiagR2[sy-1]=operatorL[NE];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
                        -operatorL[E]*u[sx+jX[E]+sy*(nx+1)]
                        -operatorL[NW]*u[sx+jX[NW]+sy*(nx+1)]
                        -operatorL[SE]*u[sx+jX[SE]+sy*(nx+1)];
                   for(Index i=9; i<operatorL.size(); i++)
                   {
                       rhs[sx+sy*(nx+1)]-=
                                operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                   }
                }
                ndiagL2[nx-5]=operatorL[SW];
                ndiagL1[nx-4]=operatorL[S];
                diagR[nx-3]=operatorL[C];
                ndiagR1[nx-3]=operatorL[N];
                rhs[sx+(ny-2)*(nx+1)]=f[sx+(ny-2)*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+(ny-2)*(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+(ny-2)*(nx+1)]
                    -operatorL[NW]*u[sx+jX[NW]+(ny-2)*(nx+1)]
                    -operatorL[SE]*u[sx+jX[SE]+(ny-2)*(nx+1)]
                    -operatorL[NE]*u[sx+(ny-2+jY[NE])*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[sx+(ny-2)*(nx+1)]-=
                               operatorL[i]*u[sx+jX[i]+(ny-2+jY[i])*(nx+1)];
                }
                operatorLB=stencil.getL(N,sx,ny-1,nx,ny);
                jXB=stencil.getJx(N,nx,ny);
                jYB=stencil.getJy(N,nx,ny);
                ndiagL2[nx-4]=operatorLB[SW];
                ndiagL1[nx-3]=operatorLB[S];
                diagR[nx-2]=operatorLB[C];
                rhs[sx+(ny-1)*(nx+1)]=f[sx+(ny-1)*(nx+1)]
                    -operatorLB[N]*u[sx+(ny-1+jYB[N])*(nx+1)] 
                    -operatorLB[W]*u[sx+jXB[W]+(ny-1)*(nx+1)]
                    -operatorLB[E]*u[sx+jXB[E]+(ny-1)*(nx+1)]
                    -operatorLB[NW]*u[sx+jXB[NW]+(ny-1)*(nx+1)]
                    -operatorLB[SE]*u[sx+jXB[SE]+(ny-1)*(nx+1)];
                for(Index i=9; i<operatorLB.size(); i++)
                {
                    rhs[sx+(ny-1)*(nx+1)]-=
                            operatorLB[i]*u[sx+jXB[i]+(ny-1+jYB[i])*(nx+1)];
                }
                yLRSolver(u,sx,nx,ny,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            }
////////////////last column//////////////////
            operatorLC=stencil.getL(SE,nx-1,1,nx,ny);
            jXC=stencil.getJx(SE,nx,ny);
            jYC=stencil.getJy(SE,nx,ny);
            rhs[nx-1+nx+1]=f[nx-1+nx+1]
                -operatorLC[S]*u[nx-1+(1+jYC[S])*(nx+1)]
                -operatorLC[W]*u[nx-2+nx+1]
                -operatorLC[E]*u[nx+nx+1]
                -operatorLC[NW]*u[nx-1+jXC[NW]+nx+1];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[nx-1+nx+1]-=
                            operatorLC[i]*u[nx-1+jXC[i]+(1+jYC[i])*(nx+1)];
            }
            operatorLB=stencil.getL(E,nx-1,2,nx,ny);
            jXB=stencil.getJx(E,nx,ny);
            jYB=stencil.getJy(E,nx,ny);
            ndiagL1[0]=operatorLB[S];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[N];                    
            ndiagR2[1]=operatorLB[NE];
            rhs[nx-1+2*(nx+1)]=f[nx-1+2*(nx+1)]
                -operatorLB[W]*u[nx-2+2*(nx+1)]
                -operatorLB[E]*u[nx+2*(nx+1)]
                -operatorLB[NW]*u[nx-1+jXB[NW]+2*(nx+1)]
                -operatorLB[SE]*u[nx-1+(2+jYB[SE])*(nx+1)];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[nx-1+2*(nx+1)]-=
                            operatorLB[i]*u[nx-1+jXB[i]+(2+jYB[i])*(nx+1)];
            }
            for(Index sy=3; sy<ny-2; sy++)  
            {
                ndiagL2[sy-3]=operatorLB[SE];
                ndiagL1[sy-2]=operatorLB[S];
                diagR[sy-1]=operatorLB[C];
                ndiagR1[sy-1]=operatorLB[N];                  
                ndiagR2[sy-1]=operatorLB[NE];
                rhs[nx-1+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorLB[W]*u[nx-1+jXB[W]+sy*(nx+1)]
                    -operatorLB[E]*u[nx-1+jXB[E]+sy*(nx+1)]
                    -operatorLB[NW]*u[nx-1+jXB[NW]+sy*(nx+1)];
                for(Index i=8; i<operatorLB.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                            operatorLB[i]*u[nx-1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
            }
            ndiagL2[ny-5]=operatorLB[SE];
            ndiagL1[ny-4]=operatorLB[S];
            diagR[ny-3]=operatorLB[C];
            ndiagR1[ny-3]=operatorLB[N];
            rhs[nx-1+(ny-2)*(nx+1)]=f[nx-1+(ny-2)*(nx+1)]
                -operatorLB[W]*u[nx-2+(ny-2)*(nx+1)]
                -operatorLB[E]*u[nx+(ny-2)*(nx+1)]
                -operatorLB[NW]*u[nx-1+jXB[NW]+(ny-2)*(nx+1)]
                -operatorLB[NE]*u[nx-1+(ny-2+jYB[NE])*(nx+1)];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[nx-1+(ny-2)*(nx+1)]-=
                        operatorLB[i]*u[nx-1+jXB[i]+(ny-2+jYB[i])*(nx+1)];
            }
            operatorLC=stencil.getL(NE,nx-1,ny-1,nx,ny);
            jXC=stencil.getJx(NE,nx,ny);
            jYC=stencil.getJy(NE,nx,ny);
            ndiagL2[ny-4]=operatorLC[NE];
            ndiagL1[ny-3]=operatorLC[S];
            diagR[ny-2]=operatorLC[C];
            rhs[nx-1+(ny-1)*(nx+1)]=f[nx-1+(ny-1)*(nx+1)]
                -operatorLC[W]*u[nx-2+(ny-1)*(nx+1)] 
                -operatorLC[E]*u[nx+(ny-1)*(nx+1)]
                -operatorLC[NW]*u[nx-1+jXC[NW]+(ny-1)*(nx+1)]
                -operatorLC[N]*u[nx-1+(ny-1+jYC[N])*(nx+1)];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[nx-1+(ny-1)*(nx+1)]-=
                        operatorLC[i]*u[nx-1+jXC[i]+(ny-1+jYC[i])*(nx+1)];
            }
            yLRSolver(u,nx-1,nx,ny,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            // process all even inner columns
            for(Index sx=2; sx < nx-1; sx+=2)
            {
                // set rhs
                operatorLB=stencil.getL(S,sx,1,nx,ny);
                jXB=stencil.getJx(S,nx,ny);
                jYB=stencil.getJy(S,nx,ny);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[N];
                ndiagR2[0]=operatorLB[NE];
                rhs[sx+nx+1]=f[sx+nx+1]
                    -operatorLB[S]*u[sx+(1+jYB[S])*(nx+1)]
                    -operatorLB[W]*u[sx+jXB[W]+nx+1] 
                    -operatorLB[E]*u[sx+jXB[E]+nx+1]
                    -operatorLB[NW]*u[sx+jXB[NW]+nx+1]
                    -operatorLB[SE]*u[sx+jXB[SE]+(nx+1)];
                for(Index i=7; i<operatorLB.size(); i++)
                {
                    rhs[sx+nx+1]-=operatorLB[i]*u[sx+jXB[i]+(1+jYB[i])*(nx+1)];
                }
                ndiagL1[0]=operatorL[S];
                diagR[1]=operatorL[C];
                ndiagR1[1]=operatorL[N];                  
                ndiagR2[1]=operatorL[NE];
                rhs[sx+2*(nx+1)]=f[sx+2*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+2*(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+2*(nx+1)]
                    -operatorL[NW]*u[sx+jX[NW]+2*(nx+1)]
                    -operatorL[SE]*u[sx+jX[SE]+2*(nx+1)]
                    -operatorL[SW]*u[sx+(2+jY[SW])*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[sx+2*(nx+1)]-=operatorL[i]*u[sx+jX[i]+(2+jY[i])*(nx+1)];
                }
                for(Index sy=3; sy<ny-2; sy++)
                {
                    ndiagL2[sy-3]=operatorL[SW];
                    ndiagL1[sy-2]=operatorL[S];
                    diagR[sy-1]=operatorL[C];
                    ndiagR1[sy-1]=operatorL[N];
                    ndiagR2[sy-1]=operatorL[NE];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
                        -operatorL[E]*u[sx+jX[E]+sy*(nx+1)]
                        -operatorL[NW]*u[sx+jX[NW]+sy*(nx+1)]
                        -operatorL[SE]*u[sx+jX[SE]+sy*(nx+1)];
                    for(Index i=9; i<operatorL.size(); i++)
                    {
                        rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                    }
                }
                ndiagL2[nx-5]=operatorL[SW];
                ndiagL1[nx-4]=operatorL[S];
                diagR[nx-3]=operatorL[C];
                ndiagR1[nx-3]=operatorL[N];
                rhs[sx+(ny-2)*(nx+1)]=f[sx+(ny-2)*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+(ny-2)*(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+(ny-2)*(nx+1)]
                    -operatorL[NW]*u[sx+jX[NW]+(ny-2)*(nx+1)]
                    -operatorL[SE]*u[sx+jX[SE]+(ny-2)*(nx+1)]
                    -operatorL[NE]*u[sx+(ny-2+jY[NE])*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[sx+(ny-2)*(nx+1)]-=
                                operatorL[i]*u[sx+jX[i]+(ny-2+jY[i])*(nx+1)];
                }
                operatorLB=stencil.getL(N,sx,ny-1,nx,ny);
                jXB=stencil.getJx(N,nx,ny);
                jYB=stencil.getJy(N,nx,ny);
                ndiagL2[nx-4]=operatorLB[SW];
                ndiagL1[nx-3]=operatorLB[S];
                diagR[nx-2]=operatorLB[C];
                rhs[sx+(ny-1)*(nx+1)]=f[sx+(ny-1)*(nx+1)]
                    -operatorLB[N]*u[sx+(ny-1+jYB[N])*(nx+1)] 
                    -operatorLB[W]*u[sx+jXB[W]+(ny-1)*(nx+1)]
                    -operatorLB[E]*u[sx+jXB[E]+(ny-1)*(nx+1)]
                    -operatorLB[NW]*u[sx+jXB[NW]+(ny-1)*(nx+1)]
                    -operatorLB[SE]*u[sx+jXB[SE]+(ny-1)*(nx+1)];
                for(Index i=9; i<operatorLB.size(); i++)
                {
                    rhs[sx+(ny-1)*(nx+1)]-=
                                operatorLB[i]*u[sx+jXB[i]+(ny-1+jYB[i])*(nx+1)];
                }
                yLRSolver(u,sx,nx,ny,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);    
            }
        }
        else // stencil not constant
        {
            NumericArray operatorL=stencil.getL(C,2,2,nx,ny);
            PositionArray jX=stencil.getJx(C,nx,ny);
            PositionArray jY=stencil.getJy(C,nx,ny);
            NumericArray operatorLB=stencil.getL(W,1,2,nx,ny);
            PositionArray jXB=stencil.getJx(W,nx,ny);
            PositionArray jYB=stencil.getJy(W,nx,ny);
            NumericArray operatorLC=stencil.getL(SW,1,1,nx,ny);
            PositionArray jXC=stencil.getJx(SW,nx,ny);
            PositionArray jYC=stencil.getJy(SW,nx,ny);
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[N];
            ndiagR2[0]=operatorLC[NW];
            rhs[1+nx+1]=f[1+nx+1]
                -operatorLC[S]*u[1+(1+jYC[S])*(nx+1)]
                -operatorLC[W]*u[nx+1]
                -operatorLC[E]*u[2+nx+1]
                -operatorLC[NE]*u[1+jXC[NE]+nx+1];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[1+nx+1]-=operatorLC[i]*u[1+jXC[i]+(1+jYC[i])*(nx+1)];
            }
            operatorLB=stencil.getL(W,1,2,nx,ny);
            jXB=stencil.getJx(W,nx,ny);
            jYB=stencil.getJy(W,nx,ny);
            ndiagL1[0]=operatorLB[S];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[N];                    
            ndiagR2[1]=operatorLB[NW];
            rhs[1+2*(nx+1)]=f[1+2*(nx+1)]
                -operatorLB[W]*u[2*(nx+1)]
                -operatorLB[E]*u[2+2*(nx+1)]
                -operatorLB[NE]*u[1+jXB[NE]+2*(nx+1)]
                -operatorLB[SE]*u[1+(2+jYB[SE])*(nx+1)];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[1+2*(nx+1)]-=operatorLB[i]*u[1+jXB[i]+(2+jYB[i])*(nx+1)];
            }
            for(Index sy=3; sy<ny-2; sy++)  
            {
                operatorLB=stencil.getL(W,1,sy,nx,ny);
                ndiagL2[sy-3]=operatorLB[SE];
                ndiagL1[sy-2]=operatorLB[S];
                diagR[sy-1]=operatorLB[C];
                ndiagR1[sy-1]=operatorLB[N];                  
                ndiagR2[sy-1]=operatorLB[NW];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorLB[W]*u[1+jXB[W]+sy*(nx+1)]
                    -operatorLB[E]*u[1+jXB[E]+sy*(nx+1)]
                    -operatorLB[NE]*u[1+jXB[NE]+sy*(nx+1)];
                for(Index i=8; i<operatorLB.size(); i++)
                {
                    rhs[1+sy*(nx+1)]-=
                                operatorLB[i]*u[1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
            }
            operatorLB=stencil.getL(W,1,ny-2,nx,ny);
            ndiagL2[ny-5]=operatorLB[SE];
            ndiagL1[ny-4]=operatorLB[S];
            diagR[ny-3]=operatorLB[C];
            ndiagR1[ny-3]=operatorLB[N];
            rhs[1+(ny-2)*(nx+1)]=f[1+(ny-2)*(nx+1)]
                -operatorLB[W]*u[(ny-2)*(nx+1)]
                -operatorLB[E]*u[2+(ny-2)*(nx+1)]
                -operatorLB[NE]*u[1+jXB[NE]+(ny-2)*(nx+1)]
                -operatorLB[NW]*u[1+(ny-2+jYB[NW])*(nx+1)];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[1+(ny-2)*(nx+1)]-=
                                operatorLB[i]*u[1+jXB[i]+(ny-2+jYB[i])*(nx+1)];
            }
            operatorLC=stencil.getL(NW,1,ny-1,nx,ny);
            jXC=stencil.getJx(NW,nx,ny);
            jYC=stencil.getJy(NW,nx,ny);
            ndiagL2[ny-4]=operatorLC[NE];
            ndiagL1[ny-3]=operatorLC[S];
            diagR[ny-2]=operatorLC[C];
            rhs[1+(ny-1)*(nx+1)]=f[1+(ny-1)*(nx+1)]
                -operatorLC[W]*u[(ny-1)*(nx+1)]
                -operatorLC[E]*u[2+(ny-1)*(nx+1)]
                -operatorLC[NW]*u[1+jXC[NW]+(ny-1)*(nx+1)]
                -operatorLC[N]*u[1+(ny-1+jYC[N])*(nx+1)];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[1+(ny-1)*(nx+1)]-=
                                operatorLC[i]*u[1+jXC[i]+(ny-1+jYC[i])*(nx+1)];
            }
            yLRSolver(u,1,nx,ny,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
////////////////////////////////////////////////////////////////////////////////
            // process odd inner columns
            for(Index sx=3; sx<nx-2; sx+=2)
            {
                //set rhs                  
                operatorLB=stencil.getL(S,sx,1,nx,ny);
                jXB=stencil.getJx(S,nx,ny);
                jYB=stencil.getJy(S,nx,ny);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[N];
                ndiagR2[0]=operatorLB[NE];
                rhs[sx+nx+1]=f[sx+nx+1]
                    -operatorLB[S]*u[sx+(1+jYB[S])*(nx+1)]
                    -operatorLB[W]*u[sx+jXB[W]+nx+1] 
                    -operatorLB[E]*u[sx+jXB[E]+nx+1]
                    -operatorLB[NW]*u[sx+jXB[NW]+nx+1]
                    -operatorLB[SE]*u[sx+jXB[SE]+(nx+1)];
                for(Index i=7; i<operatorLB.size(); i++)
                {
                    rhs[sx+nx+1]-=operatorLB[i]*u[sx+jXB[i]+(1+jYB[i])*(nx+1)];
                }
                operatorL=stencil.getL(C,sx,2,nx,ny);
                ndiagL1[0]=operatorL[S];
                diagR[1]=operatorL[C];
                ndiagR1[1]=operatorL[N];                  
                ndiagR2[1]=operatorL[NE];
                rhs[sx+2*(nx+1)]=f[sx+2*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+2*(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+2*(nx+1)]
                    -operatorL[NW]*u[sx+jX[NW]+2*(nx+1)]
                    -operatorL[SE]*u[sx+jX[SE]+2*(nx+1)]
                    -operatorL[SW]*u[sx+(2+jY[SW])*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[sx+2*(nx+1)]-=operatorL[i]*u[sx+jX[i]+(2+jY[i])*(nx+1)];
                }
                for(Index sy=3; sy<ny-2; sy++)
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    ndiagL2[sy-3]=operatorL[SW];
                    ndiagL1[sy-2]=operatorL[S];
                    diagR[sy-1]=operatorL[C];
                    ndiagR1[sy-1]=operatorL[N];
                    ndiagR2[sy-1]=operatorL[NE];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
                        -operatorL[E]*u[sx+jX[E]+sy*(nx+1)]
                        -operatorL[NW]*u[sx+jX[NW]+sy*(nx+1)]
                        -operatorL[SE]*u[sx+jX[SE]+sy*(nx+1)];
                    for(Index i=9; i<operatorL.size(); i++)
                    {
                        rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                    }
                }
                operatorL=stencil.getL(C,sx,ny-2,nx,ny);
                ndiagL2[nx-5]=operatorL[SW];
                ndiagL1[nx-4]=operatorL[S];
                diagR[nx-3]=operatorL[C];
                ndiagR1[nx-3]=operatorL[N];
                rhs[sx+(ny-2)*(nx+1)]=f[sx+(ny-2)*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+(ny-2)*(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+(ny-2)*(nx+1)]
                    -operatorL[NW]*u[sx+jX[NW]+(ny-2)*(nx+1)]
                    -operatorL[SE]*u[sx+jX[SE]+(ny-2)*(nx+1)]
                    -operatorL[NE]*u[sx+(ny-2+jY[NE])*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[sx+(ny-2)*(nx+1)]-=
                                   operatorL[i]*u[sx+jX[i]+(ny-2+jY[i])*(nx+1)];
                }
                operatorLB=stencil.getL(N,sx,ny-1,nx,ny);
                jXB=stencil.getJx(N,nx,ny);
                jYB=stencil.getJy(N,nx,ny);
                ndiagL2[nx-4]=operatorLB[SW];
                ndiagL1[nx-3]=operatorLB[S];
                diagR[nx-2]=operatorLB[C];
                rhs[sx+(ny-1)*(nx+1)]=f[sx+(ny-1)*(nx+1)]
                    -operatorLB[N]*u[sx+(ny-1+jYB[N])*(nx+1)] 
                    -operatorLB[W]*u[sx+jXB[W]+(ny-1)*(nx+1)]
                    -operatorLB[E]*u[sx+jXB[E]+(ny-1)*(nx+1)]
                    -operatorLB[NW]*u[sx+jXB[NW]+(ny-1)*(nx+1)]
                    -operatorLB[SE]*u[sx+jXB[SE]+(ny-1)*(nx+1)];
                for(Index i=9; i<operatorLB.size(); i++)
                {
                    rhs[sx+(ny-1)*(nx+1)]-=
                                operatorLB[i]*u[sx+jXB[i]+(ny-1+jYB[i])*(nx+1)];
                }
                yLRSolver(u,sx,nx,ny,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            }
////////////////last column//////////////////
            operatorLC=stencil.getL(SE,nx-1,1,nx,ny);
            jXC=stencil.getJx(SE,nx,ny);
            jYC=stencil.getJy(SE,nx,ny);
            rhs[nx-1+nx+1]=f[nx-1+nx+1]
                -operatorLC[S]*u[nx-1+(1+jYC[S])*(nx+1)]
                -operatorLC[W]*u[nx-2+nx+1]
                -operatorLC[E]*u[nx+nx+1]
                -operatorLC[NW]*u[nx-1+jXC[NW]+nx+1];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[nx-1+nx+1]-=operatorLC[i]*u[nx-1+jXC[i]+(1+jYC[i])*(nx+1)];
            }
            operatorLB=stencil.getL(E,nx-1,2,nx,ny);
            jXB=stencil.getJx(E,nx,ny);
            jYB=stencil.getJy(E,nx,ny);
            ndiagL1[0]=operatorLB[S];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[N];                    
            ndiagR2[1]=operatorLB[NE];
            rhs[nx-1+2*(nx+1)]=f[nx-1+2*(nx+1)]
                -operatorLB[W]*u[nx-2+2*(nx+1)]
                -operatorLB[E]*u[nx+2*(nx+1)]
                -operatorLB[NW]*u[nx-1+jXB[NW]+2*(nx+1)]
                -operatorLB[SE]*u[nx-1+(2+jYB[SE])*(nx+1)];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[nx-1+2*(nx+1)]-=
                                operatorLB[i]*u[nx-1+jXB[i]+(2+jYB[i])*(nx+1)];
            }
            for(Index sy=3; sy<ny-2; sy++)  
            {
                operatorLB=stencil.getL(E,nx-1,sy,nx,ny);
                ndiagL2[sy-3]=operatorLB[SE];
                ndiagL1[sy-2]=operatorLB[S];
                diagR[sy-1]=operatorLB[C];
                ndiagR1[sy-1]=operatorLB[N];                  
                ndiagR2[sy-1]=operatorLB[NE];
                rhs[nx-1+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorLB[W]*u[nx-1+jXB[W]+sy*(nx+1)]
                    -operatorLB[E]*u[nx-1+jXB[E]+sy*(nx+1)]
                    -operatorLB[NW]*u[nx-1+jXB[NW]+sy*(nx+1)];
                for(Index i=8; i<operatorLB.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                operatorLB[i]*u[nx-1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
            }
            operatorLB=stencil.getL(E,nx-1,ny-2,nx,ny);
            ndiagL2[ny-5]=operatorLB[SE];
            ndiagL1[ny-4]=operatorLB[S];
            diagR[ny-3]=operatorLB[C];
            ndiagR1[ny-3]=operatorLB[N];
            rhs[nx-1+(ny-2)*(nx+1)]=f[nx-1+(ny-2)*(nx+1)]
                -operatorLB[W]*u[nx-2+(ny-2)*(nx+1)]
                -operatorLB[E]*u[nx+(ny-2)*(nx+1)]
                -operatorLB[NW]*u[nx-1+jXB[NW]+(ny-2)*(nx+1)]
                -operatorLB[NE]*u[nx-1+(ny-2+jYB[NE])*(nx+1)];
            for(Index i=8; i<operatorLB.size(); i++)
            {
                rhs[nx-1+(ny-2)*(nx+1)]-=
                            operatorLB[i]*u[nx-1+jXB[i]+(ny-2+jYB[i])*(nx+1)];
            }
            operatorLC=stencil.getL(NE,nx-1,ny-1,nx,ny);
            jXC=stencil.getJx(NE,nx,ny);
            jYC=stencil.getJy(NE,nx,ny);
            ndiagL2[ny-4]=operatorLC[NE];
            ndiagL1[ny-3]=operatorLC[S];
            diagR[ny-2]=operatorLC[C];
            rhs[nx-1+(ny-1)*(nx+1)]=f[nx-1+(ny-1)*(nx+1)]
                -operatorLC[W]*u[nx-2+(ny-1)*(nx+1)] 
                -operatorLC[E]*u[nx+(ny-1)*(nx+1)]
                -operatorLC[NW]*u[nx-1+jXC[NW]+(ny-1)*(nx+1)]
                -operatorLC[N]*u[nx-1+(ny-1+jYC[N])*(nx+1)];
            for(Index i=7; i<operatorLC.size(); i++)
            {
                rhs[nx-1+(ny-1)*(nx+1)]-=
                            operatorLC[i]*u[nx-1+jXC[i]+(ny-1+jYC[i])*(nx+1)];
            }
            yLRSolver(u,nx-1,nx,ny,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            // process even inner columns
            for(Index sx=2; sx < nx-1; sx+=2)
            {
                // set rhs
                operatorLB=stencil.getL(S,sx,1,nx,ny);
                jXB=stencil.getJx(S,nx,ny);
                jYB=stencil.getJy(S,nx,ny);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[N];
                ndiagR2[0]=operatorLB[NE];
                rhs[sx+nx+1]=f[sx+nx+1]
                    -operatorLB[S]*u[sx+(1+jYB[S])*(nx+1)]
                    -operatorLB[W]*u[sx+jXB[W]+nx+1] 
                    -operatorLB[E]*u[sx+jXB[E]+nx+1]
                    -operatorLB[NW]*u[sx+jXB[NW]+nx+1]
                    -operatorLB[SE]*u[sx+jXB[SE]+(nx+1)];
                for(Index i=7; i<operatorLB.size(); i++)
                {
                    rhs[sx+nx+1]-=operatorLB[i]*u[sx+jXB[i]+(1+jYB[i])*(nx+1)];
                }
                operatorL=stencil.getL(C,sx,2,nx,ny);
                ndiagL1[0]=operatorL[S];
                diagR[1]=operatorL[C];
                ndiagR1[1]=operatorL[N];                  
                ndiagR2[1]=operatorL[NE];
                rhs[sx+2*(nx+1)]=f[sx+2*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+2*(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+2*(nx+1)]
                    -operatorL[NW]*u[sx+jX[NW]+2*(nx+1)]
                    -operatorL[SE]*u[sx+jX[SE]+2*(nx+1)]
                    -operatorL[SW]*u[sx+(2+jY[SW])*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[sx+2*(nx+1)]-=operatorL[i]*u[sx+jX[i]+(2+jY[i])*(nx+1)];
                }
                for(Index sy=3; sy<ny-2; sy++)
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    ndiagL2[sy-3]=operatorL[SW];
                    ndiagL1[sy-2]=operatorL[S];
                    diagR[sy-1]=operatorL[C];
                    ndiagR1[sy-1]=operatorL[N];
                    ndiagR2[sy-1]=operatorL[NE];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
                        -operatorL[E]*u[sx+jX[E]+sy*(nx+1)]
                        -operatorL[NW]*u[sx+jX[NW]+sy*(nx+1)]
                        -operatorL[SE]*u[sx+jX[SE]+sy*(nx+1)];
                    for(Index i=9; i<operatorL.size(); i++)
                    {
                        rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                    }
                }
                operatorL=stencil.getL(C,sx,ny-2,nx,ny);
                ndiagL2[nx-5]=operatorL[SW];
                ndiagL1[nx-4]=operatorL[S];
                diagR[nx-3]=operatorL[C];
                ndiagR1[nx-3]=operatorL[N];
                rhs[sx+(ny-2)*(nx+1)]=f[sx+(ny-2)*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+(ny-2)*(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+(ny-2)*(nx+1)]
                    -operatorL[NW]*u[sx+jX[NW]+(ny-2)*(nx+1)]
                    -operatorL[SE]*u[sx+jX[SE]+(ny-2)*(nx+1)]
                    -operatorL[NE]*u[sx+(ny-2+jY[NE])*(nx+1)];
                for(Index i=9; i<operatorL.size(); i++)
                {
                    rhs[sx+(ny-2)*(nx+1)]-=
                                operatorL[i]*u[sx+jX[i]+(ny-2+jY[i])*(nx+1)];
                }
                operatorLB=stencil.getL(N,sx,ny-1,nx,ny);
                jXB=stencil.getJx(N,nx,ny);
                jYB=stencil.getJy(N,nx,ny);
                ndiagL2[nx-4]=operatorLB[SW];
                ndiagL1[nx-3]=operatorLB[S];
                diagR[nx-2]=operatorLB[C];
                rhs[sx+(ny-1)*(nx+1)]=f[sx+(ny-1)*(nx+1)]
                    -operatorLB[N]*u[sx+(ny-1+jYB[N])*(nx+1)] 
                    -operatorLB[W]*u[sx+jXB[W]+(ny-1)*(nx+1)]
                    -operatorLB[E]*u[sx+jXB[E]+(ny-1)*(nx+1)]
                    -operatorLB[NW]*u[sx+jXB[NW]+(ny-1)*(nx+1)]
                    -operatorLB[SE]*u[sx+jXB[SE]+(ny-1)*(nx+1)];
                for(Index i=9; i<operatorLB.size(); i++)
                {
                    rhs[sx+(ny-1)*(nx+1)]-=
                                operatorLB[i]*u[sx+jXB[i]+(ny-1+jYB[i])*(nx+1)];
                }
                yLRSolver(u,sx,nx,ny,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);    
            }

        }
    }
    else //nx,ny to small
    {
        for(int i=0; i<2; i++)
        {
            gsRedBlack_.relax(u,f,stencil,nx,ny);
        }
    }
}
}
