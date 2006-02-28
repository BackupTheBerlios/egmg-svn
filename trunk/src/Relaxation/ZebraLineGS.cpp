/** \file ZebraLineGS.cpp
 * \author Andre Oeckerath
 * \brief ZebraLineGS.cpp contains the implementation of the class ZebraLineGS.
 * \see ZebraLineGS.h
 */
#include<iostream>
#include "ZebraLineGS.h"

namespace mg
{		
void ZebraLineGS::relax(
    std::valarray<Precision> &u,
    const std::valarray<Precision> &f, 
    const Stencil &stencil,
    const size_t nx,
    const size_t ny) const
{
    // valarrays needed for LR-decomposition of a tridiagonal matrix
    std::valarray<Precision> rhs(0.0,u.size());
    switch (stencil.size())
    {
        case 1:  // stencil of size 1
        {
            switch (direction_)
            {
                case ALTDIR:
                { 
                    ninepointxzebra(u, f, rhs, stencil, nx, ny);
                    ninepointyzebra(u, f, rhs, stencil, nx, ny);
                    break;
                }
                case XDIR:
                {
                    ninepointxzebra(u, f, rhs, stencil, nx, ny);
                    break;
                }
                case YDIR:
                {
                    ninepointyzebra(u, f, rhs, stencil, nx, ny);
                    break;
                }
                default:
                {
                    std::cerr << "Error in direction of the line relaxation!\n";
                    break;
                }
            }
            break;			
        }
	    case 2:  //stencil of size 2
        {
            switch (direction_)
            {
                case ALTDIR:
                { 
                    xzebra(u, f, rhs, stencil, nx, ny);
                    yzebra(u, f, rhs, stencil, nx, ny);
                    break;
                }
                case XDIR:
                {
                    xzebra(u, f, rhs, stencil, nx, ny);
                    break;
                }	
                case YDIR:
                {
                    yzebra(u, f, rhs, stencil, nx, ny);
                    break;
                }
                default:
                {
                    std::cerr << "Error in direction of the line relaxation!\n";
                    break;
                }
            }
            break;
        }
	    default:
        {
            std::cerr << "Stencil is too big (size>2)!" << std::endl;
            break;
        }
    }
}
/**
 * \todo clean up ninepointxzebra, ninepointyzebra, xzebra and yzebra by
 *       doint redudant things in seperate functions
 */
void ZebraLineGS::ninepointxzebra(
    std::valarray<Precision> &u,
    const std::valarray<Precision> &f, 
    std::valarray<Precision> &rhs,
    const Stencil &stencil,
    const size_t nx, 
    const size_t ny) const                      
{ 
    //valarrays needed for saving the tridiagonal matrix A of linear system
    //A u = rhs
    std::valarray<Precision> diagR(nx-1);
    std::valarray<Precision> ndiagR(nx-2);
    std::valarray<Precision> ndiagL(nx-2);
    if(stencil.isConstant())
    {
        // get const operator L
        const std::valarray<Precision> operatorL = stencil.getL(C,2,2,nx,ny);
        const std::valarray<int> jX = stencil.getJx(C);
        const std::valarray<int> jY = stencil.getJy(C);
            
        // for each line: correction of the rhs given by 
        // rhs = fv - [L[n]  0  L[s]]^t * u and elimination of the 
        // boundary condition in first and last inner point
        // odd lines
        for(size_t sy=1; sy<ny ; sy+=2) 
        {
            rhs[1+sy*(nx+1)] = f[1+sy*(nx+1)]
                    -operatorL[N]*u[1+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[1+(sy+jY[S])*(nx+1)]
                    -operatorL[W]*u[sy*(nx+1)];
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[1+sy*(nx+1)]-=
                                operatorL[i]*u[1+jX[i]+(sy+jY[i])*(nx+1)];
            }
            for(size_t sx=2; sx<nx-1; sx++)  
            {
                rhs[sx+sy*(nx+1)] = f[sx+sy*(nx+1)] 
                    -operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[sx+(sy+jY[S])*(nx+1)];
                for(size_t i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+sy*(nx+1)]-=
                                operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                }
            }
            rhs[(nx-1)+sy*(nx+1)] = f[nx-1+sy*(nx+1)]
                    -operatorL[N]*u[nx-1+(sy+jY[N])*(nx+1)] 
                    -operatorL[S]*u[nx-1+(sy+jY[S])*(nx+1)]
                    -operatorL[E]*u[nx+sy*(nx+1)];
            for(size_t i=5; i<operatorL.size(); i++)
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
            xLRSolver(u,sy,rhs,ndiagL,diagR,ndiagR);
        }
        // same for even lines 
        for(size_t sy=2; sy<ny ; sy+=2)
        {
            rhs[1+sy*(nx+1)] = f[1+sy*(nx+1)]
                -operatorL[N]*u[1+(sy+jY[N])*(nx+1)]
                -operatorL[S]*u[1+(sy+jY[S])*(nx+1)] 
                -operatorL[W]*u[sy*(nx+1)];
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[1+sy*(nx+1)]-=operatorL[i]*u[1+jX[i]+(sy+jY[i])*(nx+1)];
            }
            for(size_t sx=2; sx<nx-1; sx++) 
            {
                rhs[sx+sy*(nx+1)] = f[sx+sy*(nx+1)]
                    -operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[sx+(sy+jY[S])*(nx+1)];
                for(size_t i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                }
            }
            rhs[(nx-1)+sy*(nx+1)] = f[nx-1+sy*(nx+1)]
                -operatorL[N]*u[nx-1+(sy+jY[N])*(nx+1)]
                -operatorL[S]*u[nx-1+(sy+jY[S])*(nx+1)] 
                -operatorL[E]*u[nx+sy*(nx+1)];
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[nx-1+sy*(nx+1)]-=
                                   operatorL[i]*u[nx-1+jX[i]+(sy+jY[i])*(nx+1)];
            }
            diagR = operatorL[C];
            ndiagR = operatorL[E];
            ndiagL = operatorL[W];
            xLRSolver(u,sy,rhs,ndiagL,diagR,ndiagR);
        }
    }
    else
    {
        //Stencil ist not constant, so L needs to be evaluated in each grid
        //point no other change in the algorithm  
        std::valarray<Precision> operatorL = stencil.getL(C,2,2,nx,ny);
        std::valarray<int> jX = stencil.getJx(C);
        std::valarray<int> jY = stencil.getJy(C);
        if(nx > 2)
        {
            operatorL=stencil.getL(SW,1,1,nx,ny);
            diagR[0]=operatorL[C];
            ndiagR[0]=operatorL[E];
            rhs[1+nx+1]=f[1+nx+1]
                -operatorL[N]*u[1+(1+jY[N])*(nx+1)]
                -operatorL[S]*u[1+(1+jY[S])*(nx+1)] 
                -operatorL[W]*u[nx+1];
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[1+nx+1]-=operatorL[i]*u[1+jX[i]+(1+jY[i])*(nx+1)];
            }
            for(size_t sx=2; sx<nx-1; sx++)
            {
                operatorL=stencil.getL(S,sx,1,nx,ny);
                diagR[sx-1]=operatorL[C];
                ndiagR[sx-1]=operatorL[E];
                ndiagL[sx-2]=operatorL[W];
                rhs[sx+nx+1]=f[sx+nx+1]
                    -operatorL[N]*u[sx+(1+jY[N])*(nx+1)]
                    -operatorL[S]*u[sx+(1+jY[S])*(nx+1)];
                for(size_t i=5; i<operatorL.size(); i++)
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
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[nx-1+(nx+1)]-=operatorL[i]*u[nx-1+jX[i]+(1+jY[i])*(nx+1)];
            }
            xLRSolver(u,1,rhs,ndiagL,diagR,ndiagR);
            for(size_t sy=3; sy<ny-1 ; sy+=2)
            {
                operatorL=stencil.getL(W,1,sy,nx,ny);
                diagR[0]=operatorL[C];
                ndiagR[0]=operatorL[E];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorL[N]*u[1+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[1+(sy+jY[S])*(nx+1)] 
                    -operatorL[W]*u[sy*(nx+1)];
                for(size_t i=5; i<operatorL.size(); i++)
                {
                    rhs[1+sy*(nx+1)]-=operatorL[i]*u[1+jX[i]+(sy+jY[i])*(nx+1)];
                }
                for(size_t sx=2; sx<nx-1; sx++)
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    diagR[sx-1]=operatorL[C];
                    ndiagR[sx-1]=operatorL[E];
                    ndiagL[sx-2]=operatorL[W];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
                        -operatorL[S]*u[sx+(sy+jY[S])*(nx+1)];
                    for(size_t i=5; i<operatorL.size(); i++)
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
                for(size_t i=5; i<operatorL.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                   operatorL[i]*u[nx-1+jX[i]+(sy+jY[i])*(nx+1)];
                }
                xLRSolver(u,sy,rhs,ndiagL,diagR,ndiagR);
            }
            operatorL=stencil.getL(NW,1,ny-1,nx,ny);
            diagR[0]=operatorL[C];
            ndiagR[0]=operatorL[E];
            rhs[1+(ny-1)*(nx+1)] = f[1+(ny-1)*(nx+1)]
                -operatorL[N]*u[1+(ny-1+jY[N])*(nx+1)]
                -operatorL[S]*u[1+(ny-1+jY[S])*(nx+1)] 
                -operatorL[W]*u[(ny-1)*(nx+1)];
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[1+(ny-1)*(nx+1)]-=
                                    operatorL[i]*u[1+jX[i]+(ny-1+jY[i])*(nx+1)];
            }
            for(size_t sx=2; sx<nx-1; sx++)
            {
                operatorL=stencil.getL(N,sx,ny-1,nx,ny);
                diagR[sx-1]=operatorL[C];
                ndiagR[sx-1]=operatorL[E];
                ndiagL[sx-2]=operatorL[W];
                rhs[sx+(ny-1)*(nx+1)] = f[sx+(ny-1)*(nx+1)]
                    -operatorL[N]*u[sx+(ny-1+jY[N])*(nx+1)]
                    -operatorL[S]*u[sx+(ny-1+jY[S])*(nx+1)];
                for(size_t i=5; i<operatorL.size(); i++)
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
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[nx-1+(ny-1)*(nx+1)]-=
                                operatorL[i]*u[nx-1+jX[i]+(ny-1+jY[i])*(nx+1)];
            }
            xLRSolver(u,(ny-1),rhs,ndiagL,diagR,ndiagR);
            //even lines
            for(size_t sy=2; sy<ny ; sy+=2) 
            {                   
                operatorL=stencil.getL(W,1,sy,nx,ny);
                diagR[0]=operatorL[C];
                ndiagR[0]=operatorL[E];
                rhs[1+sy*(nx+1)] = f[1+sy*(nx+1)]
                    -operatorL[N]*u[1+(sy+jY[N])*(nx+1)]
                    -operatorL[S]*u[1+(sy+jY[S])*(nx+1)] 
                    -operatorL[W]*u[sy*(nx+1)];
                for(size_t i=5; i<operatorL.size(); i++)
                {
                    rhs[1+sy*(nx+1)]-=operatorL[i]*u[1+jX[i]+(sy+jY[i])*(nx+1)];
                }
                // L im Zentrum im Punkt (j/i)
                for(size_t sx=2; sx<nx-1; sx++)  
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    diagR[sx-1]=operatorL[C];
                    ndiagR[sx-1]=operatorL[E];
                    ndiagL[sx-2]=operatorL[W];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
                        -operatorL[S]*u[sx+(sy+jY[S])*(nx+1)];
                    for(size_t i=5; i<operatorL.size(); i++)
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
                for(size_t i=5; i<operatorL.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                   operatorL[i]*u[nx-1+jX[i]+(sy+jY[i])*(nx+1)];
                }
                xLRSolver(u,sy,rhs,ndiagL,diagR,ndiagR);
            }
        }       
        else // if nx and ny are too small do one GS_lex step
        {
            Precision temp=0;
            for(size_t sy=1; sy<ny; sy++)
            {
                operatorL=stencil.getL(C,1,sy,nx,ny);
                for(size_t sum=5; sum < operatorL.size(); sum++)
                {
                    temp -= operatorL[sum] * u[1+jX[sum]+(sy+jY[sum])*(nx+1)];
                }
                u[1+sy*(nx+1)] =1/operatorL[C]*(f[1+sy*(nx+1)]
                                            -operatorL[W]*u[1+jX[W]+sy*(nx+1)]
                                            -operatorL[E]*u[1+jX[E]+sy*(nx+1)] 
                                            -operatorL[N]*u[1+(sy+jY[N])*(nx+1)]
                                            -operatorL[S]*u[1+(sy+jY[S])*(nx+1)]
                                            -temp);
            }
        }
    }                       
}
void ZebraLineGS::ninepointyzebra(
    std::valarray<Precision> &u,
    const std::valarray<Precision> &f, 
    std::valarray<Precision> &rhs,
    const Stencil &stencil,
    const size_t nx, 
    const size_t ny) const
{ 
    //valarrays needed for saving the tridiagonal matrix A of linear system
    //A u = rhs
    std::valarray<Precision> diagR(ny-1);
    std::valarray<Precision> ndiagR(ny-2);
    std::valarray<Precision> ndiagL(ny-2);
    if(stencil.isConstant())
    {
        // get const operator L
        const std::valarray<Precision> operatorL=stencil.getL(C,2,2,nx,ny);
        const std::valarray<int> jX = stencil.getJx(C);
        const std::valarray<int> jY = stencil.getJy(C);
        // for each odd line correction of the rhs given by
        // rhs = fv + [1  0  1] * u and elimination of the 
        // boundary condition in first and last inner point  
        for(size_t sx=1; sx<nx ; sx+=2)
        {
            rhs[sx+(nx+1)]=f[sx+(nx+1)]
                -operatorL[W]*u[sx+jX[W]+(nx+1)]
                -operatorL[E]*u[sx+jX[E]+(nx+1)]
                -operatorL[S]*u[sx+(1+jY[S])*(nx+1)];
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[sx+nx+1]-=operatorL[i]*u[sx+jX[i]+(1+jY[i])*(nx+1)];
            }
            for(size_t sy=2; sy<ny-1; sy++)
            {
                rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+sy*(nx+1)];
                for(size_t i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                }
            }
            rhs[sx+(ny-1)*(nx+1)]=f[sx+(ny-1)*(nx+1)]
                -operatorL[W]*u[sx+jX[W]+(ny-1)*(nx+1)]
                -operatorL[E]*u[sx+jX[E]+(ny-1)*(nx+1)] 
                -operatorL[N]*u[sx+(ny-1+jY[N])*(nx+1)];
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[sx+(ny-1)*(nx+1)]-=
                                   operatorL[i]*u[sx+jX[i]+(ny-1+jY[i])*(nx+1)];
            }
            // set tridiagonalmatrix for solving A u = rhs
            // A[i][i] = L[c]; A[i-1][i] = L[s]; A[i+1][i] = L[n]
            diagR = operatorL[C];
            ndiagR = operatorL[N];
            ndiagL = operatorL[S];
            yLRSolver(u,sx,rhs,ndiagL,diagR,ndiagR);
        }
        // same for each even line
        for(size_t sx=2; sx<nx ; sx+=2) 
        {
            rhs[sx+(nx+1)]=f[sx+(nx+1)]
                -operatorL[W]*u[sx+jX[W]+(nx+1)]
                -operatorL[E]*u[sx+jX[E]+(nx+1)] 
                -operatorL[S]*u[sx+(1+jY[S])*(nx+1)]; 
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[sx+nx+1]-=operatorL[i]*u[sx+jX[i]+(1+jY[i])*(nx+1)];
            }
            for(size_t sy=2; sy<ny-1; sy++) 
            {
                rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+sy*(nx+1)];
                for(size_t i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+sy*(nx+1)]-=
                                    operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
                }
            }
            rhs[sx+(ny-1)*(nx+1)]=f[sx+(ny-1)*(nx+1)]
                -operatorL[W]*u[sx+jX[W]+(ny-1)*(nx+1)] 
                -operatorL[E]*u[sx+jX[E]+(ny-1)*(nx+1)]
                -operatorL[N]*u[sx+(ny-1+jY[N])*(nx+1)];
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[sx+(ny-1)*(nx+1)]-=
                                   operatorL[i]*u[sx+jX[i]+(ny-1+jY[i])*(nx+1)];
            }
            diagR=operatorL[C];
            ndiagR=operatorL[N];
            ndiagL=operatorL[S];
            yLRSolver(u,sx,rhs,ndiagL,diagR,ndiagR);
        }
    }
    else
    {
        //Stencil ist not constant, so L needs to be evaluated in each grid
        //point no other change in the algorithm  
        std::valarray<Precision> operatorL = stencil.getL(C,2,2,nx,ny);
        std::valarray<int> jX = stencil.getJx(C);
        std::valarray<int> jY = stencil.getJy(C);
        if(ny > 2)
        {
            operatorL=stencil.getL(SW,1,1,nx,ny);
            diagR[0]=operatorL[C];
            ndiagR[0]=operatorL[N];
            rhs[1+(nx+1)]=f[1+(nx+1)]
                -operatorL[W]*u[1+jX[W]+(nx+1)]
                -operatorL[E]*u[1+jX[E]+(nx+1)] 
                -operatorL[S]*u[1+(1+jY[S])*(nx+1)];
            for(size_t i=5; i<operatorL.size(); i++)
            {
                    rhs[1+nx+1]-=operatorL[i]*u[1+jX[i]+(1+jY[i])*(nx+1)];
            }
            for(size_t sy=2; sy<ny-1; sy++) 
            {
                operatorL = stencil.getL(W,1,sy,nx,ny);
                diagR[sy-1]=operatorL[C];
                ndiagR[sy-1]=operatorL[N];
                ndiagL[sy-2]=operatorL[S];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorL[W]*u[1+jX[W]+sy*(nx+1)]
                    -operatorL[E]*u[1+jX[E]+sy*(nx+1)];
                for(size_t i=5; i<operatorL.size(); i++)
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
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[1+(ny-1)*(nx+1)]-=
                                    operatorL[i]*u[1+jX[i]+(ny-1+jY[i])*(nx+1)];
            }
            yLRSolver(u,1,rhs,ndiagL,diagR,ndiagR);
            for(size_t sx=3; sx<nx-1 ; sx+=2)
            {
                operatorL=stencil.getL(S,sx,1,nx,ny);
                diagR[0]=operatorL[C];
                ndiagR[0]=operatorL[N];
                rhs[sx+(nx+1)]=f[sx+(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+(nx+1)] 
                    -operatorL[S]*u[sx+(1+jY[S])*(nx+1)];
                for(size_t i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+nx+1]-=operatorL[i]*u[sx+jX[i]+(1+jY[i])*(nx+1)];
                }
                for(size_t sy=2; sy<ny-1; sy++) 
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    diagR[sy-1]=operatorL[C];
                    ndiagR[sy-1]=operatorL[N];
                    ndiagL[sy-2]=operatorL[S];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
                        -operatorL[E]*u[sx+jX[E]+sy*(nx+1)];
                    for(size_t i=5; i<operatorL.size(); i++)
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
                for(size_t i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+(ny-1)*(nx+1)]-=
                                   operatorL[i]*u[sx+jX[i]+(ny-1+jY[i])*(nx+1)];
                }
                yLRSolver(u,sx,rhs,ndiagL,diagR,ndiagL);
            }
            operatorL=stencil.getL(SE,nx-1,1,nx,ny);
            diagR[0]=operatorL[C];
            ndiagR[0]=operatorL[N];
            rhs[nx-1+(nx+1)]=f[nx-1+(nx+1)]
                -operatorL[W]*u[nx-1+jX[W]+(nx+1)]
                -operatorL[E]*u[nx-1+jX[E]+(nx+1)] 
                -operatorL[S]*u[nx-1+(1+jY[S])*(nx+1)];
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[nx-1+nx+1]-=operatorL[i]*u[nx-1+jX[i]+(1+jY[i])*(nx+1)];
            }
            for(size_t sy=2; sy<ny-1; sy++) 
            {
                operatorL=stencil.getL(E,nx-1,sy,nx,ny);
                diagR[sy-1]=operatorL[C];
                ndiagR[sy-1]=operatorL[N];
                ndiagL[sy-2]=operatorL[S];
                rhs[nx-1+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorL[W]*u[nx-1+jX[W]+sy*(nx+1)]
                    -operatorL[E]*u[nx-1+jX[E]+sy*(nx+1)];
                for(size_t i=5; i<operatorL.size(); i++)
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
            for(size_t i=5; i<operatorL.size(); i++)
            {
                rhs[nx-1+(ny-1)*(nx+1)]-=
                                operatorL[i]*u[nx-1+jX[i]+(ny-1+jY[i])*(nx+1)];
            }
            yLRSolver(u,nx-1,rhs,ndiagL,diagR,ndiagR);
            for(size_t sx=2; sx<nx ; sx+=2)
            {
                operatorL=stencil.getL(S,sx,1,nx,ny);
                diagR[0]=operatorL[C];
                ndiagR[0]=operatorL[N];
                rhs[sx+(nx+1)]=f[sx+(nx+1)]
                    -operatorL[W]*u[sx+jX[W]+(nx+1)]
                    -operatorL[E]*u[sx+jX[E]+(nx+1)] 
                    -operatorL[S]*u[sx+(1+jY[S])*(nx+1)];
                for(size_t i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+nx+1]-=operatorL[i]*u[sx+jX[i]+(1+jY[i])*(nx+1)];
                }
                for(size_t sy=2; sy<ny-1; sy++) 
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    diagR[sy-1]=operatorL[C];
                    ndiagR[sy-1]=operatorL[N];                                     
                    ndiagL[sy-2]=operatorL[S];
                    rhs[sx+sy*(nx+1)]=f[sx+sy*(nx+1)]
                        -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
                        -operatorL[E]*u[sx+jX[E]+sy*(nx+1)]; 
                    for(size_t i=5; i<operatorL.size(); i++)
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
                for(size_t i=5; i<operatorL.size(); i++)
                {
                    rhs[sx+(ny-1)*(nx+1)]-=
                                   operatorL[i]*u[sx+jX[i]+(ny-1+jY[i])*(nx+1)];
                }
                yLRSolver(u,sx,rhs,ndiagL,diagR,ndiagR);
            }
        }
        else // if Nx and Ny are too small do one GS_lex step
        {
            for(size_t sx=1; sx<nx; sx++)
            {
                Precision temp=0;
                operatorL=stencil.getL(C,1,sx,nx,ny);
                for(size_t i=5; i<operatorL.size(); i++)
                {
                    temp-=operatorL[i]*u[1+jX[i]+(sx+jY[i])*(nx+1)];
                }
                u[1+sx*(nx+1)]=1/operatorL[C]*(f[1+sx*(nx+1)]
                                    -operatorL[W]*u[1+jX[W]+sx*(nx+1)]
                                    -operatorL[E]*u[1+jX[E]+sx*(nx+1)] 
                                    -operatorL[N]*u[1+(sx+jY[N])*(nx+1)]
                                    -operatorL[S]*u[1+(sx+jY[S])*(nx+1)]
                                    -temp);
            }
        }
    }                       
}
void ZebraLineGS::xzebra(
    std::valarray<Precision> &u,
    const std::valarray<Precision> &f, 
    std::valarray<Precision> &rhs,
    const Stencil &stencil,
    const size_t nx, 
    const size_t ny) const
{
    if((ny > 4) && (nx > 4))
    {   
        std::valarray<Precision> diagR(0.0,nx-1);
        std::valarray<Precision> ndiagR1(0.0,nx-2);
        std::valarray<Precision> ndiagL1(0.0,nx-2);
        std::valarray<Precision> ndiagR2(0.0,nx-3);
        std::valarray<Precision> ndiagL2(0.0,nx-3);
        if(stencil.isConstant())
        {
            // get const operator L
            const std::valarray<Precision> operatorL=stencil.getL(C,2,2,nx,ny);
            const std::valarray<int> jX=stencil.getJx(C);
            const std::valarray<int> jY=stencil.getJy(C);
            std::valarray<Precision> operatorLB=stencil.getL(S,2,1,nx,ny);
            std::valarray<int> jXB=stencil.getJx(S);
            std::valarray<int> jYB=stencil.getJy(S);        
            std::valarray<Precision> operatorLC=stencil.getL(SW,1,1,nx,ny);
            std::valarray<int> jXC=stencil.getJx(SW);
            std::valarray<int> jYC=stencil.getJy(SW);
            // set right hand side for line 1                   
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[E];
            ndiagR2[0]=operatorLC[NE];
            rhs[1+nx+1]=f[1+nx+1]
                -operatorLC[N]*u[1+(1+jYC[N])*(nx+1)]
                -operatorLC[S]*u[1+(1+jYC[S])*(nx+1)]
                -operatorLC[W]*u[nx+1]
                -operatorLC[NW]*u[1+(1+jYC[NW])*(nx+1)];  
            for(size_t i=7; i<operatorLC.size(); i++)
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
            for(size_t i=8; i<operatorLB.size(); i++)
            {
                rhs[2+nx+1]-=operatorLB[i]*u[2+jXB[i]+(1+jYB[i])*(nx+1)];
            }
            for(size_t sx=3; sx<nx-2; sx++)  
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
                for(size_t i=8; i<operatorLB.size(); i++)
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
            for(size_t i=8; i<operatorLB.size(); i++)
            {
                rhs[nx-2+nx+1]-=operatorLB[i]*u[nx-2+jXB[i]+(1+jYB[i])*(nx+1)];
            }
            operatorLC=stencil.get_L_se(nx-1,1,nx,ny);
            jXC=stencil.getJx(SE);
            jYC=stencil.getJy(SE);
            ndiagL2[nx-4]=operatorLC[NW];
            ndiagL1[nx-3]=operatorLC[W];
            diagR[nx-2]=operatorLC[C];
            rhs[(nx-1)+(nx+1)]=f[nx-1+nx+1]
                -operatorLC[N]*u[nx-1+(1+jYC[N])*(nx+1)] 
                -operatorLC[S]*u[nx-1+(1+jYC[S])*(nx+1)]
                -operatorLC[E]*u[nx+nx+1]
                -operatorLC[NE]*u[nx-1+(1+jYC[NE])*(nx+1)];
            for(size_t i=7; i<operatorLC.size(); i++)
            {
                rhs[nx-1+(nx+1)]-=
                                operatorLC[i]*u[nx-1+jXC[i]+(1+jYC[i])*(nx+1)];
            }
            xLRSolver(u,1,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            for(size_t sy=3; sy<ny-2; sy+=2)
            {                   
                operatorLB=stencil.getL(W,1,sy,nx,ny);
                jXB=stencil.getJx(W);
                jYB=stencil.getJy(W);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[E];
                ndiagR2[0]=operatorLB[NE];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorLB[N]*u[1+(sy+jYB[N])*(nx+1)]
                    -operatorLB[S]*u[1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[W]*u[sy*(nx+1)]
                    -operatorLB[NW]*u[1+(sy+jYB[NW])*(nx+1)]
                    -operatorLB[SE]*u[1+(sy+jYB[SE])*(nx+1)];
                for(size_t i=8; i<operatorLB.size(); i++)
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
                for(size_t i=9; i<operatorL.size(); i++)
                {
                    rhs[2+sy*(nx+1)]-=operatorL[i]*u[2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                for(size_t sx=3; sx<nx-2; sx++)  
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
                    for(size_t i=9; i<operatorL.size(); i++)
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
                for(size_t i=9; i<operatorL.size(); i++)
                {
                    rhs[nx-2+sy*(nx+1)]-=
                                   operatorL[i]*u[nx-2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                operatorLB=stencil.getL(E,nx-1,sy,nx,ny);
                jXB=stencil.getJx(E);
                jYB=stencil.getJy(E);
                ndiagL2[nx-4]=operatorLB[NW];
                ndiagL1[nx-3]=operatorLB[W];
                diagR[nx-2]=operatorLB[C];
                rhs[(nx-1)+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorLB[N]*u[nx-1+(sy+jYB[N])*(nx+1)] 
                    -operatorLB[S]*u[nx-1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[E]*u[nx+sy*(nx+1)]
                    -operatorLB[NE]*u[nx-1+(sy+jYB[NE])*(nx+1)]
                    -operatorLB[SE]*u[nx-1+(sy+jYB[SE])*(nx+1)];
                for(size_t i=9; i<operatorLB.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                operatorLB[i]*u[nx-1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
                xLRSolver(u,sy,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);    
            }
            //relax top line Zeile
            //set rhs in top line
            operatorLC=stencil.getL(NW,1,ny-1,nx,ny);
            jXC=stencil.getJx(NW);
            jYC=stencil.getJy(NW);
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[E];
            ndiagR2[0]=operatorLC[NW];
            rhs[1+(ny-1)*(nx+1)]=f[1+(ny-1)*(nx+1)]
                -operatorLC[N]*u[1+(ny-1+jYC[N])*(nx+1)]
                -operatorLC[S]*u[1+(ny-1+jYC[S])*(nx+1)]
                -operatorLC[W]*u[(ny-1)*(nx+1)]
                -operatorLC[NE]*u[1+(ny-1+jYC[NE])*(nx+1)];
            for(size_t i=7; i<operatorLC.size(); i++)
            {
                rhs[1+(ny-1)*(nx+1)]-=
                                operatorLC[i]*u[1+jXC[i]+(ny-1+jYC[i])*(nx+1)];
            }
            operatorLB=stencil.getL(N,2,ny-1,nx,ny);
            jXB=stencil.getJx(N);
            jYB=stencil.getJy(N);
            ndiagL1[0]=operatorLB[W];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[E];                    
            ndiagR2[1]=operatorLB[NE];
            rhs[2+(ny-1)*(nx+1)]=f[2+(ny-1)*(nx+1)]
                -operatorLB[N]*u[2+(ny-1+jYB[N])*(nx+1)]
                -operatorLB[S]*u[2+(ny-1+jYB[S])*(nx+1)]
                -operatorLB[SE]*u[2+(ny-1+jYB[SE])*(nx+1)]
                -operatorLB[NW]*u[2+jXB[NW]+(ny-1)*(nx+1)];
            for(size_t i=8; i<operatorLB.size(); i++)
            {
                rhs[2+(ny-1)*(nx+1)]-=
                                operatorLB[i]*u[2+jXB[i]+(ny-1+jYB[i])*(nx+1)];
            }
            for(size_t sx=3; sx<nx-2; sx++)  
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
                for(size_t i=8; i<operatorLB.size(); i++)
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
            for(size_t i=8; i<operatorLB.size(); i++)
            {
                rhs[nx-2+(ny-1)*(nx+1)]-=
                            operatorLB[i]*u[nx-2+jXB[i]+(ny-1+jYB[i])*(nx+1)];
            }
            operatorLC=stencil.getL(NE,nx-1,ny-1,nx,ny);
            jXC=stencil.getJx(NE);
            jYC=stencil.getJy(NE);
            ndiagL2[nx-4]=operatorLC[NW];
            ndiagL1[nx-3]=operatorLC[W];
            diagR[nx-2]=operatorLC[C];
            rhs[(nx-1)+(ny-1)*(nx+1)]=f[nx-1+(ny-1)*(nx+1)]
                -operatorLC[N]*u[nx-1+(ny-1+jYC[N])*(nx+1)]
                -operatorLC[S]*u[nx-1+(ny-1+jYC[S])*(nx+1)]
                -operatorLC[E]*u[nx+(ny-1)*(nx+1)]
                -operatorLC[NE]*u[nx-1+(ny-1+jYC[NE])*(nx+1)];
            for(size_t i=7; i<operatorLC.size(); i++)
            {
                rhs[nx-1+(ny-1)*(nx+1)]-=
                            operatorLC[i]*u[nx-1+jXC[i]+(ny-1+jYC[i])*(nx+1)];
            }
            xLRSolver(u,ny-1,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            // relax even inner lines
            for(size_t sy=2; sy < ny-1; sy+=2)
            {
                // set rhs           
                operatorLB=stencil.getL(W,1,sy,nx,ny);
                jXB=stencil.getJx(W);
                jYB=stencil.getJy(W);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[E];
                ndiagR2[0]=operatorLB[NE];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorLB[N]*u[1+(sy+jYB[N])*(nx+1)]
                    -operatorLB[S]*u[1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[W]*u[sy*(nx+1)]
                    -operatorLB[NW]*u[1+(sy+jYB[NW])*(nx+1)]
                    -operatorLB[SE]*u[1+(sy+jYB[SE])*(nx+1)];
                for(size_t i=8; i<operatorLB.size(); i++)
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
                for(size_t i=9; i<operatorL.size(); i++)
                {
                    rhs[2+sy*(nx+1)]-=operatorL[i]*u[2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                for(size_t sx=3; sx<nx-2; sx++)  
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
                    for(size_t i=9; i<operatorL.size(); i++)
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
                for(size_t i=9; i<operatorL.size(); i++)
                {
                    rhs[nx-2+sy*(nx+1)]-=
                                operatorL[i]*u[nx-2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                operatorLB=stencil.getL(E,nx-1,sy,nx,ny);
                jXB=stencil.getJx(E);
                jYB=stencil.getJy(E);
                ndiagL2[nx-4]=operatorLB[NW];
                ndiagL1[nx-3]=operatorLB[W];
                diagR[nx-2]=operatorLB[C];
                rhs[(nx-1)+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorLB[N]*u[nx-1+(sy+jYB[N])*(nx+1)] 
                    -operatorLB[S]*u[nx-1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[E]*u[nx+sy*(nx+1)]
                    -operatorLB[NE]*u[nx-1+(sy+jYB[NE])*(nx+1)]
                    -operatorLB[SE]*u[nx-1+(sy+jYB[SE])*(nx+1)];
                for(size_t i=9; i<operatorLB.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                operatorLB[i]*u[nx-1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
                xLRSolver(u,sy,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);    
            }
        }
        else // stencil not constant
        {
            std::valarray<Precision> operatorL=stencil.getL(C,2,2,nx,ny);
            std::valarray<int> jX=stencil.getJx(C);
            std::valarray<int> jY=stencil.getJy(C);
            std::valarray<Precision> operatorLB=stencil.getL(S,2,1,nx,ny);
            std::valarray<int> jXB=stencil.getJx(S);
            std::valarray<int> jYB=stencil.getJy(S);
            std::valarray<Precision> operatorLC=stencil.getL(SW,1,1,nx,ny);
            std::valarray<int> jXC=stencil.getJx(SW);
            std::valarray<int> jYC=stencil.getJy(SW);
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[E];
            ndiagR2[0]=operatorLC[NE];
            rhs[1+nx+1]=f[1+nx+1]
                -operatorLC[N]*u[1+(1+jYC[N])*(nx+1)]
                -operatorLC[S]*u[1+(1+jYC[S])*(nx+1)]
                -operatorLC[W]*u[nx+1]
                -operatorLC[NW]*u[1+(1+jYC[NW])*(nx+1)];
            for(size_t i=7; i<operatorLC.size(); i++)
            {
                rhs[1+nx+1]-=operatorLC[i]*u[1+jXC[i]+(1+jYC[i])*(nx+1)];
            }
            operatorLB=stencil.getL(S,2,1,nx,ny);
            jXB=stencil.getJx(S);
            jYB=stencil.getJy(S);
            ndiagL1[0]=operatorLB[W];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[E];                    
            ndiagR2[1]=operatorLB[SE];
            rhs[2+nx+1]=f[2+nx+1]
                -operatorLB[N]*u[2+(1+jYB[N])*(nx+1)]
                -operatorLB[S]*u[2+(1+jYB[S])*(nx+1)]
                -operatorLB[NE]*u[2+(1+jYB[NE])*(nx+1)]
                -operatorLB[NW]*u[2+jXB[NW]+nx+1];
            for(size_t i=8; i<operatorLB.size(); i++)
            {
                rhs[2+nx+1]-=operatorLB[i]*u[2+jXB[i]+(1+jYB[i])*(nx+1)];
            }
            for(size_t sx=3; sx<nx-2; sx++)  
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
                for(size_t i=8; i<operatorLB.size(); i++)
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
            for(size_t i=8; i<operatorLB.size(); i++)
            {
                rhs[nx-2+nx+1]-=operatorLB[i]*u[nx-2+jXB[i]+(1+jYB[i])*(nx+1)];
            }
            operatorLC=stencil.getL(SE,nx-1,1,nx,ny);
            jXC=stencil.getJx(SE);
            jYC=stencil.getJy(SE);
            ndiagL2[nx-4]=operatorLC[NW];
            ndiagL1[nx-3]=operatorLC[W];
            diagR[nx-2]=operatorLC[C];
            rhs[(nx-1)+(nx+1)]=f[nx-1+nx+1]
                -operatorLC[N]*u[nx-1+(1+jYC[N])*(nx+1)] 
                -operatorLC[S]*u[nx-1+(1+jYC[S])*(nx+1)]
                -operatorLC[E]*u[nx+nx+1]
                -operatorLC[NE]*u[nx-1+(1+jYC[NE])*(nx+1)];
            for(size_t i=7; i<operatorLC.size(); i++)
            {
                rhs[nx-1+(nx+1)]-=
                                operatorLC[i]*u[nx-1+jXC[i]+(1+jYC[i])*(nx+1)];
            }
            xLRSolver(u,1,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            // do odd inner lines
            for(size_t sy=3; sy < ny-2; sy+=2)
            {
                // set rhs                   
                operatorLB=stencil.getL(W,1,sy,nx,ny);
                jXB=stencil.getJx(W);
                jYB=stencil.getJy(W);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[E];
                ndiagR2[0]=operatorLB[NE];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorLB[N]*u[1+(sy+jYB[N])*(nx+1)]
                    -operatorLB[S]*u[1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[W]*u[sy*(nx+1)]
                    -operatorLB[NW]*u[1+(sy+jYB[NW])*(nx+1)]
                    -operatorLB[SE]*u[1+(sy+jYB[SE])*(nx+1)];
                for(size_t i=8; i<operatorLB.size(); i++)
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
                for(size_t i=9; i<operatorL.size(); i++)
                {
                    rhs[2+sy*(nx+1)]-=operatorL[i]*u[2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                for(size_t sx=3; sx<nx-2; sx++)  
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
                    for(size_t i=9; i<operatorL.size(); i++)
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
                for(size_t i=9; i<operatorL.size(); i++)
                {
                    rhs[nx-2+sy*(nx+1)]-=
                                operatorL[i]*u[nx-2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                operatorLB=stencil.getL(E,nx-1,sy,nx,ny);
                jXB=stencil.getJx(E);
                jYB=stencil.getJy(E);
                ndiagL2[nx-4]=operatorLB[NW];
                ndiagL1[nx-3]=operatorLB[W];
                diagR[nx-2]=operatorLB[C];
                rhs[(nx-1)+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorLB[N]*u[nx-1+(sy+jYB[N])*(nx+1)] 
                    -operatorLB[S]*u[nx-1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[E]*u[nx+sy*(nx+1)]
                    -operatorLB[NE]*u[nx-1+(sy+jYB[NE])*(nx+1)]
                    -operatorLB[SE]*u[nx-1+(sy+jYB[SE])*(nx+1)];
                for(size_t i=9; i<operatorLB.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                operatorLB[i]*u[nx-1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
                xLRSolver(u,sy,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            }
            //relax top line
            // set rhs in top line                
            operatorLC=stencil.getL(NW,1,ny-1,nx,ny);
            jXC=stencil.getJx(NW);
            jYC=stencil.getJy(NW);
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[E];
            ndiagR2[0]=operatorLC[NW];
            rhs[1+(ny-1)*(nx+1)]=f[1+(ny-1)*(nx+1)]
                -operatorLC[N]*u[1+(ny-1+jYC[N])*(nx+1)]
                -operatorLC[S]*u[1+(ny-1+jYC[S])*(nx+1)]
                -operatorLC[W]*u[(ny-1)*(nx+1)]
                -operatorLC[NE]*u[1+(ny-1+jYC[NE])*(nx+1)];
            for(size_t i=7; i<operatorLC.size(); i++)
            {
                rhs[1+(ny-1)*(nx+1)]-=
                            operatorLC[i]*u[1+jXC[i]+(ny-1+jYC[i])*(nx+1)];
            }
            operatorLB=stencil.getL(N,2,ny-1,nx,ny);
            jXB=stencil.getJx(N);
            jYB=stencil.getJy(N);
            ndiagL1[0]=operatorLB[W];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[E];                    
            ndiagR2[1]=operatorLB[NE];
            rhs[2+(ny-1)*(nx+1)]=f[2+(ny-1)*(nx+1)]
                -operatorLB[N]*u[2+(ny-1+jYB[N])*(nx+1)]
                -operatorLB[S]*u[2+(ny-1+jYB[S])*(nx+1)]
                -operatorLB[SE]*u[2+(ny-1+jYB[SE])*(nx+1)]
                -operatorLB[NW]*u[2+jXB[NW]+(ny-1)*(nx+1)];
            for(size_t i=8; i<operatorLB.size(); i++)
            {
                rhs[2+(ny-1)*(nx+1)]-=
                                operatorLB[i]*u[2+jXB[i]+(ny-1+jYB[i])*(nx+1)];
            }
            for(size_t sx=3; sx<nx-2; sx++)  
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
                for(size_t i=8; i<operatorLB.size(); i++)
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
            for(size_t i=8; i<operatorLB.size(); i++)
            {
                rhs[nx-2+(ny-1)*(nx+1)]-=
                            operatorLB[i]*u[nx-2+jXB[i]+(ny-1+jYB[i])*(nx+1)];
            }
            operatorLC=stencil.getL(NE,nx-1,ny-1,nx,ny);
            jXC=stencil.getJx(NE);
            jYC=stencil.getJy(NE);
            ndiagL2[nx-4]=operatorLC[NW];
            ndiagL1[nx-3]=operatorLC[W];
            diagR[nx-2]=operatorLC[C];
            rhs[(nx-1)+(ny-1)*(nx+1)]=f[nx-1+(ny-1)*(nx+1)]
                -operatorLC[N]*u[nx-1+(ny-1+jYC[N])*(nx+1)]
                -operatorLC[S]*u[nx-1+(ny-1+jYC[S])*(nx+1)]
                -operatorLC[E]*u[nx+(ny-1)*(nx+1)]
                -operatorLC[NE]*u[nx-1+(ny-1+jYC[NE])*(nx+1)];
            for(size_t i=7; i<operatorLC.size(); i++)
            {
                rhs[nx-1+(ny-1)*(nx+1)]-=
                            operatorLC[i]*u[nx-1+jXC[i]+(ny-1+jYC[i])*(nx+1)];
            }
            xLRSolver(u,ny-1,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            // relax even inner lines
            for(size_t sy=2; sy<ny-1; sy+=2)
            {
                // set rhs                   
                operatorLB=stencil.getL(W,1,sy,nx,ny);
                jXB=stencil.getJx(W);
                jYB=stencil.getJy(W);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[E];
                ndiagR2[0]=operatorLB[NE];
                rhs[1+sy*(nx+1)]=f[1+sy*(nx+1)]
                    -operatorLB[N]*u[1+(sy+jYB[N])*(nx+1)]
                    -operatorLB[S]*u[1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[W]*u[sy*(nx+1)]
                    -operatorLB[NW]*u[1+(sy+jYB[NW])*(nx+1)]
                    -operatorLB[SE]*u[1+(sy+jYB[SE])*(nx+1)];
                for(size_t i=8; i<operatorLB.size(); i++)
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
                for(size_t i=9; i<operatorL.size(); i++)
                {
                    rhs[2+sy*(nx+1)]-=operatorL[i]*u[2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                for(size_t sx=3; sx<nx-2; sx++)  
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
                    for(size_t i=9; i<operatorL.size(); i++)
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
                for(size_t i=9; i<operatorL.size(); i++)
                {
                    rhs[nx-2+sy*(nx+1)]-=
                                operatorL[i]*u[nx-2+jX[i]+(sy+jY[i])*(nx+1)];
                }
                operatorLB=stencil.getL(E,nx-1,sy,nx,ny);
                jXB=stencil.getJx(E);
                jYB=stencil.getJy(E);
                ndiagL2[nx-4]=operatorLB[NW];
                ndiagL1[nx-3]=operatorLB[W];
                diagR[nx-2]=operatorLB[C];
                rhs[(nx-1)+sy*(nx+1)]=f[nx-1+sy*(nx+1)]
                    -operatorLB[N]*u[nx-1+(sy+jYB[N])*(nx+1)] 
                    -operatorLB[S]*u[nx-1+(sy+jYB[S])*(nx+1)]
                    -operatorLB[E]*u[nx+sy*(nx+1)]
                    -operatorLB[NE]*u[nx-1+(sy+jYB[NE])*(nx+1)]
                    -operatorLB[SE]*u[nx-1+(sy+jYB[SE])*(nx+1)];
                for(size_t i=9; i<operatorLB.size(); i++)
                {
                    rhs[nx-1+sy*(nx+1)]-=
                                operatorLB[i]*u[nx-1+jXB[i]+(sy+jYB[i])*(nx+1)];
                }
                xLRSolver(u,sy,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);    
            }
        }
    }
    else 
    {
        for(int i=0; i<2; i++)
        {
            factor = 1.0/stencil.getCenter(SW,1,1,nx,ny);
            u[1*(nx+1)+1]+=factor*(f[1*(nx+1)+1]
                        -stencil.apply(u,SW,1,1,nx,ny));
            for (size_t sx=3;sx<(nx-1);sx+=2)
            {
                factor = 1.0/stencil.getCenter(S,sx,1,nx,ny);
                u[1*(nx+1)+sx]+=factor*(f[1*(nx+1)+sx]
                        -stencil.apply(u,S,sx,1,nx,ny));
            }
            factor = 1.0/stencil.getCenter(SE,(nx-1),1,nx,ny);
            u[1*(nx+1)+(nx-1)]+=factor*(f[1*(nx+1)+(nx-1)]
                        -stencil.apply(u,SE,(nx-1),1,nx,ny));
            for (size_t sy=3;sy<(ny-1);sy+=2)
            {
                factor = 1.0/stencil.getCenter(W,1,sy,nx,ny);
                u[sy*(nx+1)+1]+=factor*(f[sy*(nx+1)+1]
                        -stencil.apply(u,W,1,sy,nx,ny));
                for (size_t sx=3;sx<(nx-1);sx+=2)
                {
                    factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
                    u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
                            -stencil.apply(u,C,sx,sy,nx,ny));
                }
                factor = 1.0/stencil.getCenter(E,(nx-1),sy,nx,ny);
                u[sy*(nx+1)+(nx-1)]+=factor*(f[sy*(nx+1)+(nx-1)]
                        -stencil.apply(u,E,(nx-1),sy,nx,ny));    
            }
            //the missing red points in the center
            for (size_t sy=2;sy<(ny-1);sy+=2)
            {
                for (size_t sx=2;sx<(nx-1);sx+=2)
                {
                    factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
                    u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
                            -stencil.apply(u,C,sx,sy,nx,ny));
                }
            }
            factor = 1.0/stencil.getCenter(NW,1,(ny-1),nx,ny);
            u[(nx-1)*(nx+1)+1]+=factor*(f[(nx-1)*(nx+1)+1]
                        -stencil.apply(u,NW,1,(ny-1),nx,ny));
            for (size_t sx=3;sx<(nx-1);sx+=2)
            {
                factor = 1.0/stencil.getCenter(N,sx,(nx-1),nx,ny);
                u[(nx-1)*(nx+1)+sx]+=factor*(f[(nx-1)*(nx+1)+sx]
                        -stencil.apply(u,N,sx,(ny-1),nx,ny));
            }
            factor = 1.0/stencil.getCenter(NE,(nx-1),(nx-1),nx,ny);
            u[(nx-1)*(nx+1)+(nx-1)]+=factor*(f[(nx-1)*(nx+1)+(nx-1)]
                    -stencil.apply(u,NE,(nx-1),(ny-1),nx,ny));
            
            //do black points
            for (size_t sx=2;sx<(nx-1);sx+=2)
            {
                 factor = 1.0/stencil.getCenter(S,sx,1,nx,ny);
                    u[1*(nx+1)+sx]+=factor*(f[1*(nx+1)+sx]
                        -stencil.apply(u,S,sx,1,nx,ny));
            }
            for (size_t sy=2;sy<(ny-1);sy+=2)
            {
                factor = 1.0/stencil.getCenter(W,1,sy,nx,ny);
                u[sy*(nx+1)+1]+=factor*(f[sy*(nx+1)+1]
                        -stencil.apply(u,W,1,sy,nx,ny));    
                for (size_t sx=3;sx<(nx-1);sx+=2)
                {
                    factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
                    u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
                            -stencil.apply(u,C,sx,sy,nx,ny));
                }
                factor = 1.0/stencil.getCenter(E,(nx-1),sy,nx,ny);
                u[sy*(nx+1)+(nx-1)]+=factor*(f[sy*(nx+1)+(nx-1)]
                        -stencil.apply(u,E,(nx-1),sy,nx,ny));
            }
            for (size_t sy=3;sy<(ny-1);sy+=2)
            {
                for (size_t sx=2;sx<(nx-1);sx+=2)
                {
                    factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
                    u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
                            -stencil.apply(u,C,sx,sy,nx,ny));
                }
            }
            for (size_t sx=2;sx<(nx-1);sx+=2)
            {
                factor = 1.0/stencil.getCenter(N,sx,(ny-1),nx,ny);
                u[(nx-1)*(nx+1)+sx]+=factor*(f[(nx-1)*(nx+1)+sx]
                        -stencil.apply(u,N,sx,(ny-1),nx,ny));
            }
        }
    }
}
void ZebraLineGS::yzebra(
    std::valarray<Precision> &u,
    const std::valarray<Precision> &fv, 
    std::valarray<Precision> &rhs,
    const Stencil &stencil,
    const size_t Nx, 
    const size_t Ny) const
{
        if((Ny > 4) && (Nx > 4))
        {
            std::valarray<Precision> diagR(0.0,Ny-1);
            std::valarray<Precision> ndiagR1(0.0,Ny-2);
            std::valarray<Precision> ndiagL1(0.0,Ny-2);
            std::valarray<Precision> ndiagR2(0.0,Ny-3);
            std::valarray<Precision> ndiagL2(0.0,Ny-3);
                    
            if(stencil.isConstant() == true)
            {
                // get const operator L
                const std::valarray<Precision> L = stencil.get_L_c(2,2,Nx,Ny);
                const std::valarray<int> J_x = stencil.getJx(C);
                const std::valarray<int> J_y = stencil.getJy(C);
                
                std::valarray<Precision> L_b = stencil.get_L_w(1,2,Nx,Ny);
                std::valarray<int> J_b_x = stencil.getJx(W);
                std::valarray<int> J_b_y = stencil.getJy(W);
 
                std::valarray<Precision> L_c = stencil.get_L_sw(1,1,Nx,Ny);
                std::valarray<int> J_c_x = stencil.getJx(SW);
                std::valarray<int> J_c_y = stencil.getJy(SW);

                diagR[0] = L_c[C];
                ndiagR1[0] = L_c[N];
                ndiagR2[0] = L_c[NW];
                                
                rhs[1+Nx+1] = fv[1+Nx+1] - L_c[S] * u[1+(1+J_c_y[S])*(Nx+1)]
                    - L_c[W] * u[Nx+1] - L_c[E] * u[2+Nx+1] - L_c[NE] * u[1+J_c_x[NE]+Nx+1];
                    
                for(size_t sum=7; sum<L_c.size(); sum++)
                {
                    rhs[1+Nx+1] -= L_c[sum] * u[1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
                }

                L_b = stencil.get_L_w(1,2,Nx,Ny);
                J_b_x = stencil.getJx(W);
                J_b_y = stencil.getJy(W);

                ndiagL1[0] = L_b[S];
                diagR[1] = L_b[C];
                ndiagR1[1] = L_b[N];                    
                ndiagR2[1] = L_b[NW];

                rhs[1+2*(Nx+1)] = fv[1+2*(Nx+1)] - L_b[W] * u[2*(Nx+1)] - L_b[E] * u[2+2*(Nx+1)]
                  -L_b[NE] * u[1+J_b_x[NE]+2*(Nx+1)] - L_b[SE] * u[1+(2+J_b_y[SE])*(Nx+1)];

                for(size_t sum=8; sum<L_b.size(); sum++)
                {
                    rhs[1+2*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(2+J_b_y[sum])*(Nx+1)];
                }

                for(size_t j=3; j<Ny-2; j++)  
                {
                    ndiagL2[j-3] = L_b[SE];
                    ndiagL1[j-2] = L_b[S];
                    diagR[j-1] = L_b[C];
                    ndiagR1[j-1] = L_b[N];                  
                    ndiagR2[j-1] = L_b[NW];
                        
                    rhs[1+j*(Nx+1)] = fv[1+j*(Nx+1)] - L_b[W] * u[1+J_b_x[W]+j*(Nx+1)] - L_b[E] * u[1+J_b_x[E]+j*(Nx+1)]
                            -L_b[NE] * u[1+J_b_x[NE]+j*(Nx+1)];

                    for(size_t sum=8; sum<L_b.size(); sum++)
                    {
                      rhs[1+j*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(j+J_b_y[sum])*(Nx+1)];
                    }
                }
                    
                ndiagL2[Ny-5] = L_b[SE];
                ndiagL1[Ny-4] = L_b[S];
                diagR[Ny-3] = L_b[C];
                ndiagR1[Ny-3] = L_b[N];

                rhs[1+(Ny-2)*(Nx+1)] = fv[1+(Ny-2)*(Nx+1)] - L_b[W] * u[(Ny-2)*(Nx+1)] - L_b[E] * u[2+(Ny-2)*(Nx+1)]
                  -L_b[NE] * u[1+J_b_x[NE]+(Ny-2)*(Nx+1)] - L_b[NW] * u[1+(Ny-2+J_b_y[NW])*(Nx+1)];

                for(size_t sum=8; sum<L_b.size(); sum++)
                {
                  rhs[1+(Ny-2)*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(Ny-2+J_b_y[sum])*(Nx+1)];
                }
                    
                L_c = stencil.get_L_nw(1,Ny-1,Nx,Ny);
                J_c_x = stencil.getJx(NW);
                J_c_y = stencil.getJy(NW);

                ndiagL2[Ny-4] = L_c[NE];
                ndiagL1[Ny-3] = L_c[S];
                diagR[Ny-2] = L_c[C];

                rhs[1+(Ny-1)*(Nx+1)] = fv[1+(Ny-1)*(Nx+1)] - L_c[W] * u[(Ny-1)*(Nx+1)] - L_c[E] * u[2+(Ny-1)*(Nx+1)]
                        - L_c[NW] * u[1+J_c_x[NW]+(Ny-1)*(Nx+1)] - L_c[N] * u[1+(Ny-1+J_c_y[N])*(Nx+1)];

                for(size_t sum=7; sum<L_c.size(); sum++)
                {
                  rhs[1+(Ny-1)*(Nx+1)] -= L_c[sum] * u[1+J_c_x[sum]+(Ny-1+J_c_y[sum])*(Nx+1)];
                }

                // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<Ny-2; k++)  
                {
                    ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                    rhs[(k+1)*(Nx+1) + 1] = rhs[(k+1)*(Nx+1) + 1] - ndiagL1[k-1] * rhs[k*(Nx+1)+1]; 

                    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                    rhs[(k+2)*(Nx+1) + 1] = rhs[(k+2)*(Nx+1) + 1] - ndiagL2[k-1] * rhs[k*(Nx+1)+1];
                }
                ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
                diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
                rhs[(Ny-1)*(Nx+1) + 1] = rhs[(Ny-1)*(Nx+1) + 1] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ 1];

                // solve the linear system of equations R u = rhs
                u[1+(Nx+1)*(Ny-1)] = rhs[1+(Nx+1)*(Ny-1)] / diagR[Ny-2];
                for(size_t k=Ny-2; k>1; k--)
                {
                    u[1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[1+k*(Nx+1)] - ndiagR1[k-1] * u[1+(k+1)*(Nx+1)] );

                    rhs[(k-1)*(Nx+1)+1] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+1];
                }
                u[(Nx+1)+1] = 1/diagR[0] * ( rhs[(Nx+1)+1] - ndiagR1[0] * u[2*(Nx+1)+1] );



////////////////////////////////////////////////////////////////////////////////////////////////

                // durchlaufe alle inneren Spalten, Spaltenindex i              
                for(size_t i=3; i < Nx-2; i+=2)
                {
                    // setze rechte Seite                   
                    L_b = stencil.get_L_s(i,1,Nx,Ny);
                    J_b_x = stencil.getJx(S);
                    J_b_y = stencil.getJy(S);

                    diagR[0] = L_b[C];
                    ndiagR1[0] = L_b[N];
                    ndiagR2[0] = L_b[NE];
       
                    rhs[i+Nx+1] = fv[i+Nx+1] - L_b[S] * u[i+(1+J_b_y[S])*(Nx+1)] - L_b[W] * u[i+J_b_x[W]+Nx+1] 
                        - L_b[E] * u[i+J_b_x[E]+Nx+1] - L_b[NW] * u[i+J_b_x[NW]+Nx+1] - L_b[SE] * u[i+J_b_x[SE]+(Nx+1)];
                    
                    for(size_t sum=7; sum<L_b.size(); sum++)
                    {
                        rhs[i+Nx+1] -= L_b[sum] * u[i+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
                    }

                    ndiagL1[0] = L[S];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[N];                  
                    ndiagR2[1] = L[NE];

                    rhs[i+2*(Nx+1)] = fv[i+2*(Nx+1)] - L[W] * u[i+J_x[W]+2*(Nx+1)] - L[E] * u[i+J_x[E]+2*(Nx+1)]
                      -L[NW] * u[i+J_x[NW]+2*(Nx+1)] - L[SE] * u[i+J_x[SE]+2*(Nx+1)] - L[SW] * u[i+(2+J_y[SW])*(Nx+1)];

                    for(size_t sum=9; sum<L.size(); sum++)
                    {
                        rhs[i+2*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(2+J_y[sum])*(Nx+1)];
                    }
                        

                    for(size_t j=3; j<Ny-2; j++)
                    {
                       ndiagL2[j-3] = L[SW];
                       ndiagL1[j-2] = L[S];
                       diagR[j-1] = L[C];
                       ndiagR1[j-1] = L[N];
                       ndiagR2[j-1] = L[NE];
 
                       rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[W] * u[i+J_x[W]+j*(Nx+1)] - L[E] * u[i+J_x[E]+j*(Nx+1)]
                             -L[NW] * u[i+J_x[NW]+j*(Nx+1)] - L[SE] * u[i+J_x[SE]+j*(Nx+1)];
                       
                       for(size_t sum=9; sum<L.size(); sum++)
                       {
                           rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
                       }
                    }
                                        
                    ndiagL2[Nx-5] = L[SW];
                    ndiagL1[Nx-4] = L[S];
                    diagR[Nx-3] = L[C];
                    ndiagR1[Nx-3] = L[N];

                    rhs[i+(Ny-2)*(Nx+1)] = fv[i+(Ny-2)*(Nx+1)] - L[W] * u[i+J_x[W]+(Ny-2)*(Nx+1)] - L[E] * u[i+J_x[E]+(Ny-2)*(Nx+1)]
                         -L[NW] * u[i+J_x[NW]+(Ny-2)*(Nx+1)] - L[SE] * u[i+J_x[SE]+(Ny-2)*(Nx+1)] - L[NE] * u[i+(Ny-2+J_y[NE])*(Nx+1)];
                    
                    for(size_t sum=9; sum<L.size(); sum++)
                    {
                        rhs[i+(Ny-2)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-2+J_y[sum])*(Nx+1)];
                    }


                    L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
                    J_b_x = stencil.getJx(N);
                    J_b_y = stencil.getJy(N);

                    ndiagL2[Nx-4] = L_b[SW];
                    ndiagL1[Nx-3] = L_b[S];
                    diagR[Nx-2] = L_b[C];

                    rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L_b[N] * u[i+(Ny-1+J_b_y[N])*(Nx+1)] 
                        - L_b[W] * u[i+J_b_x[W]+(Ny-1)*(Nx+1)] - L_b[E] * u[i+J_b_x[E]+(Ny-1)*(Nx+1)] - L_b[NW] * u[i+J_b_x[NW]+(Ny-1)*(Nx+1)]
                        - L_b[SE] * u[i+J_b_x[SE]+(Ny-1)*(Nx+1)];

                    for(size_t sum=9; sum<L_b.size(); sum++)
                    {
                        rhs[i+(Ny-1)*(Nx+1)] -= L_b[sum] * u[i+J_b_x[sum]+(Ny-1+J_b_y[sum])*(Nx+1)];
                    }
                    
                    // LR-decomposition + transformation of the rhs
                    for(size_t k=1; k<Ny-2; k++)  
                    {
                        ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                        diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                        ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                        rhs[(k+1)*(Nx+1) + i] = rhs[(k+1)*(Nx+1) + i] - ndiagL1[k-1] * rhs[k*(Nx+1)+i]; 

                        ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                        ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                        diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                        rhs[(k+2)*(Nx+1) + i] = rhs[(k+2)*(Nx+1) + i] - ndiagL2[k-1] * rhs[k*(Nx+1)+i];
                    }
                    ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
                    diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
                    rhs[(Ny-1)*(Nx+1) + i] = rhs[(Ny-1)*(Nx+1) + i] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ i];

                    // solve the linear system of equations R u = rhs
                    u[i+(Nx+1)*(Ny-1)] = rhs[i+(Nx+1)*(Ny-1)] / diagR[Ny-2];

                    for(size_t k=Ny-2; k>1; k--)
                    {
                        u[i+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[i+k*(Nx+1)] - ndiagR1[k-1] * u[i+(k+1)*(Nx+1)] );

                        rhs[(k-1)*(Nx+1)+i] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+i];
                    }
                    u[(Nx+1)+i] = 1/diagR[0] * ( rhs[(Nx+1)+i] - ndiagR1[0] * u[2*(Nx+1)+i] );
                }

////////////////letzte Spalte//////////////////
                
                L_c = stencil.get_L_se(Nx-1,1,Nx,Ny);
                J_c_x = stencil.getJx(SE);
                J_c_y = stencil.getJy(SE);

                
                rhs[Nx-1+Nx+1] = fv[Nx-1+Nx+1] - L_c[S] * u[Nx-1+(1+J_c_y[S])*(Nx+1)]
                    - L_c[W] * u[Nx-2+Nx+1] - L_c[E] * u[Nx+Nx+1] - L_c[NW] * u[Nx-1+J_c_x[NW]+Nx+1];
                    
                for(size_t sum=7; sum<L_c.size(); sum++)
                {
                    rhs[Nx-1+Nx+1] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
                }

                L_b = stencil.get_L_e(Nx-1,2,Nx,Ny);
                J_b_x = stencil.getJx(E);
                J_b_y = stencil.getJy(E);

                ndiagL1[0] = L_b[S];
                diagR[1] = L_b[C];
                ndiagR1[1] = L_b[N];                    
                ndiagR2[1] = L_b[NE];

                rhs[Nx-1+2*(Nx+1)] = fv[Nx-1+2*(Nx+1)] - L_b[W] * u[Nx-2+2*(Nx+1)] - L_b[E] * u[Nx+2*(Nx+1)]
                  -L_b[NW] * u[Nx-1+J_b_x[NW]+2*(Nx+1)] - L_b[SE] * u[Nx-1+(2+J_b_y[SE])*(Nx+1)];

                for(size_t sum=8; sum<L_b.size(); sum++)
                {
                    rhs[Nx-1+2*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(2+J_b_y[sum])*(Nx+1)];
                }

                for(size_t j=3; j<Ny-2; j++)  
                {
                    ndiagL2[j-3] = L_b[SE];
                    ndiagL1[j-2] = L_b[S];
                    diagR[j-1] = L_b[C];
                    ndiagR1[j-1] = L_b[N];                  
                    ndiagR2[j-1] = L_b[NE];
                        
                    rhs[Nx-1+j*(Nx+1)] = fv[Nx-1+j*(Nx+1)] - L_b[W] * u[Nx-1+J_b_x[W]+j*(Nx+1)] - L_b[E] * u[Nx-1+J_b_x[E]+j*(Nx+1)]
                            -L_b[NW] * u[Nx-1+J_b_x[NW]+j*(Nx+1)];

                    for(size_t sum=8; sum<L_b.size(); sum++)
                    {
                        rhs[Nx-1+j*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(j+J_b_y[sum])*(Nx+1)];
                    }
                }
                    
                ndiagL2[Ny-5] = L_b[SE];
                ndiagL1[Ny-4] = L_b[S];
                diagR[Ny-3] = L_b[C];
                ndiagR1[Ny-3] = L_b[N];

                rhs[Nx-1+(Ny-2)*(Nx+1)] = fv[Nx-1+(Ny-2)*(Nx+1)] - L_b[W] * u[Nx-2+(Ny-2)*(Nx+1)] - L_b[E] * u[Nx+(Ny-2)*(Nx+1)]
                  - L_b[NW] * u[Nx-1+J_b_x[NW]+(Ny-2)*(Nx+1)] - L_b[NE] * u[Nx-1+(Ny-2+J_b_y[NE])*(Nx+1)];

                for(size_t sum=8; sum<L_b.size(); sum++)
                {
                  rhs[Nx-1+(Ny-2)*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(Ny-2+J_b_y[sum])*(Nx+1)];
                }
                    
                L_c = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
                J_c_x = stencil.getJx(NE);
                J_c_y = stencil.getJy(NE);

                ndiagL2[Ny-4] = L_c[NE];
                ndiagL1[Ny-3] = L_c[S];
                diagR[Ny-2] = L_c[C];

                rhs[Nx-1+(Ny-1)*(Nx+1)] = fv[Nx-1+(Ny-1)*(Nx+1)] - L_c[W] * u[Nx-2+(Ny-1)*(Nx+1)] 
                  - L_c[E] * u[Nx+(Ny-1)*(Nx+1)] - L_c[NW] * u[Nx-1+J_c_x[NW]+(Ny-1)*(Nx+1)] - L_c[N] * u[Nx-1+(Ny-1+J_c_y[N])*(Nx+1)];

                for(size_t sum=7; sum<L_c.size(); sum++)
                {
                    rhs[Nx-1+(Ny-1)*(Nx+1)] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(Ny-1+J_c_y[sum])*(Nx+1)];
                }

                // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<Ny-2; k++)  
                {
                    ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                    rhs[(k+1)*(Nx+1) + Nx-1] = rhs[(k+1)*(Nx+1) + Nx-1] - ndiagL1[k-1] * rhs[k*(Nx+1)+Nx-1]; 

                    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                    rhs[(k+2)*(Nx+1) + Nx-1] = rhs[(k+2)*(Nx+1) + Nx-1] - ndiagL2[k-1] * rhs[k*(Nx+1)+Nx-1];
                }
                ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
                diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
                rhs[(Ny-1)*(Nx+1) + Nx-1] = rhs[(Ny-1)*(Nx+1) + Nx-1] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ Nx-1];

                // solve the linear system of equations R u = rhs
                u[Nx-1+(Nx+1)*(Ny-1)] = rhs[Nx-1+(Nx+1)*(Ny-1)] / diagR[Ny-2];

                for(size_t k=Ny-2; k>1; k--)
                {
                    u[Nx-1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[Nx-1+k*(Nx+1)] - ndiagR1[k-1] * u[Nx-1+(k+1)*(Nx+1)] );

                    rhs[(k-1)*(Nx+1)+Nx-1] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+Nx-1];
                }
                u[(Nx+1)+Nx-1] = 1/diagR[0] * ( rhs[(Nx+1)+Nx-1] - ndiagR1[0] * u[2*(Nx+1)+Nx-1] );


                // durchlaufe alle inneren Spalten, Spaltenindex i              
                for(size_t i=2; i < Nx-1; i+=2)
                {
                    // setze rechte Seite                   
                    L_b = stencil.get_L_s(i,1,Nx,Ny);
                    J_b_x = stencil.getJx(S);
                    J_b_y = stencil.getJy(S);

                    diagR[0] = L_b[C];
                    ndiagR1[0] = L_b[N];
                    ndiagR2[0] = L_b[NE];
       
                    rhs[i+Nx+1] = fv[i+Nx+1] - L_b[S] * u[i+(1+J_b_y[S])*(Nx+1)] - L_b[W] * u[i+J_b_x[W]+Nx+1] 
                        - L_b[E] * u[i+J_b_x[E]+Nx+1] - L_b[NW] * u[i+J_b_x[NW]+Nx+1] - L_b[SE] * u[i+J_b_x[SE]+(Nx+1)];
                    
                    for(size_t sum=7; sum<L_b.size(); sum++)
                    {
                        rhs[i+Nx+1] -= L_b[sum] * u[i+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
                    }

                    ndiagL1[0] = L[S];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[N];                  
                    ndiagR2[1] = L[NE];

                    rhs[i+2*(Nx+1)] = fv[i+2*(Nx+1)] - L[W] * u[i+J_x[W]+2*(Nx+1)] - L[E] * u[i+J_x[E]+2*(Nx+1)]
                      -L[NW] * u[i+J_x[NW]+2*(Nx+1)] - L[SE] * u[i+J_x[SE]+2*(Nx+1)] - L[SW] * u[i+(2+J_y[SW])*(Nx+1)];

                    for(size_t sum=9; sum<L.size(); sum++)
                    {
                        rhs[i+2*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(2+J_y[sum])*(Nx+1)];
                    }
                        

                    for(size_t j=3; j<Ny-2; j++)
                    {
                       ndiagL2[j-3] = L[SW];
                       ndiagL1[j-2] = L[S];
                       diagR[j-1] = L[C];
                       ndiagR1[j-1] = L[N];
                       ndiagR2[j-1] = L[NE];
 
                       rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[W] * u[i+J_x[W]+j*(Nx+1)] - L[E] * u[i+J_x[E]+j*(Nx+1)]
                             -L[NW] * u[i+J_x[NW]+j*(Nx+1)] - L[SE] * u[i+J_x[SE]+j*(Nx+1)];
                       
                       for(size_t sum=9; sum<L.size(); sum++)
                       {
                           rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
                       }
                    }
                                        
                    ndiagL2[Nx-5] = L[SW];
                    ndiagL1[Nx-4] = L[S];
                    diagR[Nx-3] = L[C];
                    ndiagR1[Nx-3] = L[N];

                    rhs[i+(Ny-2)*(Nx+1)] = fv[i+(Ny-2)*(Nx+1)] - L[W] * u[i+J_x[W]+(Ny-2)*(Nx+1)] - L[E] * u[i+J_x[E]+(Ny-2)*(Nx+1)]
                         -L[NW] * u[i+J_x[NW]+(Ny-2)*(Nx+1)] - L[SE] * u[i+J_x[SE]+(Ny-2)*(Nx+1)] - L[NE] * u[i+(Ny-2+J_y[NE])*(Nx+1)];
                    
                    for(size_t sum=9; sum<L.size(); sum++)
                    {
                        rhs[i+(Ny-2)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-2+J_y[sum])*(Nx+1)];
                    }


                    L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
                    J_b_x = stencil.getJx(N);
                    J_b_y = stencil.getJy(N);

                    ndiagL2[Nx-4] = L_b[SW];
                    ndiagL1[Nx-3] = L_b[S];
                    diagR[Nx-2] = L_b[C];

                    rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L_b[N] * u[i+(Ny-1+J_b_y[N])*(Nx+1)] 
                        - L_b[W] * u[i+J_b_x[W]+(Ny-1)*(Nx+1)] - L_b[E] * u[i+J_b_x[E]+(Ny-1)*(Nx+1)] - L_b[NW] * u[i+J_b_x[NW]+(Ny-1)*(Nx+1)]
                        - L_b[SE] * u[i+J_b_x[SE]+(Ny-1)*(Nx+1)];

                    for(size_t sum=9; sum<L_b.size(); sum++)
                    {
                        rhs[i+(Ny-1)*(Nx+1)] -= L_b[sum] * u[i+J_b_x[sum]+(Ny-1+J_b_y[sum])*(Nx+1)];
                    }
                    
                    // LR-decomposition + transformation of the rhs
                    for(size_t k=1; k<Ny-2; k++)  
                    {
                        ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                        diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                        ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                        rhs[(k+1)*(Nx+1) + i] = rhs[(k+1)*(Nx+1) + i] - ndiagL1[k-1] * rhs[k*(Nx+1)+i]; 

                        ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                        ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                        diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                        rhs[(k+2)*(Nx+1) + i] = rhs[(k+2)*(Nx+1) + i] - ndiagL2[k-1] * rhs[k*(Nx+1)+i];
                    }
                    ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
                    diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
                    rhs[(Ny-1)*(Nx+1) + i] = rhs[(Ny-1)*(Nx+1) + i] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ i];

                    // solve the linear system of equations R u = rhs
                    u[i+(Nx+1)*(Ny-1)] = rhs[i+(Nx+1)*(Ny-1)] / diagR[Ny-2];

                    for(size_t k=Ny-2; k>1; k--)
                    {
                        u[i+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[i+k*(Nx+1)] - ndiagR1[k-1] * u[i+(k+1)*(Nx+1)] );

                        rhs[(k-1)*(Nx+1)+i] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+i];
                    }
                    u[(Nx+1)+i] = 1/diagR[0] * ( rhs[(Nx+1)+i] - ndiagR1[0] * u[2*(Nx+1)+i] );
                }

            }
            else // stencil not constant
            {
                std::valarray<Precision> L = stencil.get_L_c(2,2,Nx,Ny);
                std::valarray<int> J_x = stencil.getJx(C);
                std::valarray<int> J_y = stencil.getJy(C);
                
                std::valarray<Precision> L_b = stencil.get_L_w(1,2,Nx,Ny);
                std::valarray<int> J_b_x = stencil.getJx(W);
                std::valarray<int> J_b_y = stencil.getJy(W);
 
                std::valarray<Precision> L_c = stencil.get_L_sw(1,1,Nx,Ny);
                std::valarray<int> J_c_x = stencil.getJx(SW);
                std::valarray<int> J_c_y = stencil.getJy(SW);

                diagR[0] = L_c[C];
                ndiagR1[0] = L_c[N];
                ndiagR2[0] = L_c[NW];
                                
                rhs[1+Nx+1] = fv[1+Nx+1] - L_c[S] * u[1+(1+J_c_y[S])*(Nx+1)]
                    - L_c[W] * u[Nx+1] - L_c[E] * u[2+Nx+1] - L_c[NE] * u[1+J_c_x[NE]+Nx+1];
                    
                for(size_t sum=7; sum<L_c.size(); sum++)
                {
                    rhs[1+Nx+1] -= L_c[sum] * u[1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
                }

                L_b = stencil.get_L_w(1,2,Nx,Ny);
                J_b_x = stencil.getJx(W);
                J_b_y = stencil.getJy(W);

                ndiagL1[0] = L_b[S];
                diagR[1] = L_b[C];
                ndiagR1[1] = L_b[N];                    
                ndiagR2[1] = L_b[NW];

                rhs[1+2*(Nx+1)] = fv[1+2*(Nx+1)] - L_b[W] * u[2*(Nx+1)] - L_b[E] * u[2+2*(Nx+1)]
                  -L_b[NE] * u[1+J_b_x[NE]+2*(Nx+1)] - L_b[SE] * u[1+(2+J_b_y[SE])*(Nx+1)];

                for(size_t sum=8; sum<L_b.size(); sum++)
                {
                    rhs[1+2*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(2+J_b_y[sum])*(Nx+1)];
                }

                for(size_t j=3; j<Ny-2; j++)  
                {
                    L_b = stencil.get_L_w(1,j,Nx,Ny);
                    
                    ndiagL2[j-3] = L_b[SE];
                    ndiagL1[j-2] = L_b[S];
                    diagR[j-1] = L_b[C];
                    ndiagR1[j-1] = L_b[N];                  
                    ndiagR2[j-1] = L_b[NW];
                        
                    rhs[1+j*(Nx+1)] = fv[1+j*(Nx+1)] - L_b[W] * u[1+J_b_x[W]+j*(Nx+1)] - L_b[E] * u[1+J_b_x[E]+j*(Nx+1)]
                            -L_b[NE] * u[1+J_b_x[NE]+j*(Nx+1)];

                    for(size_t sum=8; sum<L_b.size(); sum++)
                    {
                      rhs[1+j*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(j+J_b_y[sum])*(Nx+1)];
                    }
                }

                L_b = stencil.get_L_w(1,Ny-2,Nx,Ny);
                    
                ndiagL2[Ny-5] = L_b[SE];
                ndiagL1[Ny-4] = L_b[S];
                diagR[Ny-3] = L_b[C];
                ndiagR1[Ny-3] = L_b[N];

                rhs[1+(Ny-2)*(Nx+1)] = fv[1+(Ny-2)*(Nx+1)] - L_b[W] * u[(Ny-2)*(Nx+1)] - L_b[E] * u[2+(Ny-2)*(Nx+1)]
                  -L_b[NE] * u[1+J_b_x[NE]+(Ny-2)*(Nx+1)] - L_b[NW] * u[1+(Ny-2+J_b_y[NW])*(Nx+1)];

                for(size_t sum=8; sum<L_b.size(); sum++)
                {
                  rhs[1+(Ny-2)*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(Ny-2+J_b_y[sum])*(Nx+1)];
                }
                    
                L_c = stencil.get_L_nw(1,Ny-1,Nx,Ny);
                J_c_x = stencil.getJx(NW);
                J_c_y = stencil.getJy(NW);

                ndiagL2[Ny-4] = L_c[NE];
                ndiagL1[Ny-3] = L_c[S];
                diagR[Ny-2] = L_c[C];

                rhs[1+(Ny-1)*(Nx+1)] = fv[1+(Ny-1)*(Nx+1)] - L_c[W] * u[(Ny-1)*(Nx+1)] - L_c[E] * u[2+(Ny-1)*(Nx+1)]
                        - L_c[NW] * u[1+J_c_x[NW]+(Ny-1)*(Nx+1)] - L_c[N] * u[1+(Ny-1+J_c_y[N])*(Nx+1)];

                for(size_t sum=7; sum<L_c.size(); sum++)
                {
                  rhs[1+(Ny-1)*(Nx+1)] -= L_c[sum] * u[1+J_c_x[sum]+(Ny-1+J_c_y[sum])*(Nx+1)];
                }

                // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<Ny-2; k++)  
                {
                    ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                    rhs[(k+1)*(Nx+1) + 1] = rhs[(k+1)*(Nx+1) + 1] - ndiagL1[k-1] * rhs[k*(Nx+1)+1]; 

                    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                    rhs[(k+2)*(Nx+1) + 1] = rhs[(k+2)*(Nx+1) + 1] - ndiagL2[k-1] * rhs[k*(Nx+1)+1];
                }
                ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
                diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
                rhs[(Ny-1)*(Nx+1) + 1] = rhs[(Ny-1)*(Nx+1) + 1] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ 1];

                // solve the linear system of equations R u = rhs
                u[1+(Nx+1)*(Ny-1)] = rhs[1+(Nx+1)*(Ny-1)] / diagR[Ny-2];
                for(size_t k=Ny-2; k>1; k--)
                {
                    u[1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[1+k*(Nx+1)] - ndiagR1[k-1] * u[1+(k+1)*(Nx+1)] );

                    rhs[(k-1)*(Nx+1)+1] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+1];
                }
                u[(Nx+1)+1] = 1/diagR[0] * ( rhs[(Nx+1)+1] - ndiagR1[0] * u[2*(Nx+1)+1] );



////////////////////////////////////////////////////////////////////////////////////////////////

                // durchlaufe alle inneren Spalten, Spaltenindex i              
                for(size_t i=3; i < Nx-2; i+=2)
                {
                    // setze rechte Seite                   
                    L_b = stencil.get_L_s(i,1,Nx,Ny);
                    J_b_x = stencil.getJx(S);
                    J_b_y = stencil.getJy(S);

                    diagR[0] = L_b[C];
                    ndiagR1[0] = L_b[N];
                    ndiagR2[0] = L_b[NE];
       
                    rhs[i+Nx+1] = fv[i+Nx+1] - L_b[S] * u[i+(1+J_b_y[S])*(Nx+1)] - L_b[W] * u[i+J_b_x[W]+Nx+1] 
                        - L_b[E] * u[i+J_b_x[E]+Nx+1] - L_b[NW] * u[i+J_b_x[NW]+Nx+1] - L_b[SE] * u[i+J_b_x[SE]+(Nx+1)];
                    
                    for(size_t sum=7; sum<L_b.size(); sum++)
                    {
                        rhs[i+Nx+1] -= L_b[sum] * u[i+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
                    }

                    L = stencil.get_L_c(i,2,Nx,Ny);

                    ndiagL1[0] = L[S];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[N];                  
                    ndiagR2[1] = L[NE];

                    rhs[i+2*(Nx+1)] = fv[i+2*(Nx+1)] - L[W] * u[i+J_x[W]+2*(Nx+1)] - L[E] * u[i+J_x[E]+2*(Nx+1)]
                      -L[NW] * u[i+J_x[NW]+2*(Nx+1)] - L[SE] * u[i+J_x[SE]+2*(Nx+1)] - L[SW] * u[i+(2+J_y[SW])*(Nx+1)];

                    for(size_t sum=9; sum<L.size(); sum++)
                    {
                        rhs[i+2*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(2+J_y[sum])*(Nx+1)];
                    }
                        

                    for(size_t j=3; j<Ny-2; j++)
                    {
                       L = stencil.get_L_c(i,j,Nx,Ny);
                       
                       ndiagL2[j-3] = L[SW];
                       ndiagL1[j-2] = L[S];
                       diagR[j-1] = L[C];
                       ndiagR1[j-1] = L[N];
                       ndiagR2[j-1] = L[NE];
 
                       rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[W] * u[i+J_x[W]+j*(Nx+1)] - L[E] * u[i+J_x[E]+j*(Nx+1)]
                             -L[NW] * u[i+J_x[NW]+j*(Nx+1)] - L[SE] * u[i+J_x[SE]+j*(Nx+1)];
                       
                       for(size_t sum=9; sum<L.size(); sum++)
                       {
                           rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
                       }
                    }

                    L = stencil.get_L_c(i,Ny-2,Nx,Ny);
                                        
                    ndiagL2[Nx-5] = L[SW];
                    ndiagL1[Nx-4] = L[S];
                    diagR[Nx-3] = L[C];
                    ndiagR1[Nx-3] = L[N];

                    rhs[i+(Ny-2)*(Nx+1)] = fv[i+(Ny-2)*(Nx+1)] - L[W] * u[i+J_x[W]+(Ny-2)*(Nx+1)] - L[E] * u[i+J_x[E]+(Ny-2)*(Nx+1)]
                         -L[NW] * u[i+J_x[NW]+(Ny-2)*(Nx+1)] - L[SE] * u[i+J_x[SE]+(Ny-2)*(Nx+1)] - L[NE] * u[i+(Ny-2+J_y[NE])*(Nx+1)];
                    
                    for(size_t sum=9; sum<L.size(); sum++)
                    {
                        rhs[i+(Ny-2)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-2+J_y[sum])*(Nx+1)];
                    }


                    L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
                    J_b_x = stencil.getJx(N);
                    J_b_y = stencil.getJy(N);

                    ndiagL2[Nx-4] = L_b[SW];
                    ndiagL1[Nx-3] = L_b[S];
                    diagR[Nx-2] = L_b[C];

                    rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L_b[N] * u[i+(Ny-1+J_b_y[N])*(Nx+1)] 
                        - L_b[W] * u[i+J_b_x[W]+(Ny-1)*(Nx+1)] - L_b[E] * u[i+J_b_x[E]+(Ny-1)*(Nx+1)] - L_b[NW] * u[i+J_b_x[NW]+(Ny-1)*(Nx+1)]
                        - L_b[SE] * u[i+J_b_x[SE]+(Ny-1)*(Nx+1)];

                    for(size_t sum=9; sum<L_b.size(); sum++)
                    {
                        rhs[i+(Ny-1)*(Nx+1)] -= L_b[sum] * u[i+J_b_x[sum]+(Ny-1+J_b_y[sum])*(Nx+1)];
                    }
                    
                    // LR-decomposition + transformation of the rhs
                    for(size_t k=1; k<Ny-2; k++)  
                    {
                        ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                        diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                        ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                        rhs[(k+1)*(Nx+1) + i] = rhs[(k+1)*(Nx+1) + i] - ndiagL1[k-1] * rhs[k*(Nx+1)+i]; 

                        ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                        ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                        diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                        rhs[(k+2)*(Nx+1) + i] = rhs[(k+2)*(Nx+1) + i] - ndiagL2[k-1] * rhs[k*(Nx+1)+i];
                    }
                    ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
                    diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
                    rhs[(Ny-1)*(Nx+1) + i] = rhs[(Ny-1)*(Nx+1) + i] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ i];

                    // solve the linear system of equations R u = rhs
                    u[i+(Nx+1)*(Ny-1)] = rhs[i+(Nx+1)*(Ny-1)] / diagR[Ny-2];

                    for(size_t k=Ny-2; k>1; k--)
                    {
                        u[i+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[i+k*(Nx+1)] - ndiagR1[k-1] * u[i+(k+1)*(Nx+1)] );

                        rhs[(k-1)*(Nx+1)+i] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+i];
                    }
                    u[(Nx+1)+i] = 1/diagR[0] * ( rhs[(Nx+1)+i] - ndiagR1[0] * u[2*(Nx+1)+i] );
                }

////////////////letzte Spalte//////////////////
                
                L_c = stencil.get_L_se(Nx-1,1,Nx,Ny);
                J_c_x = stencil.getJx(SE);
                J_c_y = stencil.getJy(SE);

                
                rhs[Nx-1+Nx+1] = fv[Nx-1+Nx+1] - L_c[S] * u[Nx-1+(1+J_c_y[S])*(Nx+1)]
                    - L_c[W] * u[Nx-2+Nx+1] - L_c[E] * u[Nx+Nx+1] - L_c[NW] * u[Nx-1+J_c_x[NW]+Nx+1];
                    
                for(size_t sum=7; sum<L_c.size(); sum++)
                {
                    rhs[Nx-1+Nx+1] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
                }

                L_b = stencil.get_L_e(Nx-1,2,Nx,Ny);
                J_b_x = stencil.getJx(E);
                J_b_y = stencil.getJy(E);

                ndiagL1[0] = L_b[S];
                diagR[1] = L_b[C];
                ndiagR1[1] = L_b[N];                    
                ndiagR2[1] = L_b[NE];

                rhs[Nx-1+2*(Nx+1)] = fv[Nx-1+2*(Nx+1)] - L_b[W] * u[Nx-2+2*(Nx+1)] - L_b[E] * u[Nx+2*(Nx+1)]
                  -L_b[NW] * u[Nx-1+J_b_x[NW]+2*(Nx+1)] - L_b[SE] * u[Nx-1+(2+J_b_y[SE])*(Nx+1)];

                for(size_t sum=8; sum<L_b.size(); sum++)
                {
                    rhs[Nx-1+2*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(2+J_b_y[sum])*(Nx+1)];
                }

                for(size_t j=3; j<Ny-2; j++)  
                {
                    L_b = stencil.get_L_e(Nx-1,j,Nx,Ny);

                    ndiagL2[j-3] = L_b[SE];
                    ndiagL1[j-2] = L_b[S];
                    diagR[j-1] = L_b[C];
                    ndiagR1[j-1] = L_b[N];                  
                    ndiagR2[j-1] = L_b[NE];
                        
                    rhs[Nx-1+j*(Nx+1)] = fv[Nx-1+j*(Nx+1)] - L_b[W] * u[Nx-1+J_b_x[W]+j*(Nx+1)] - L_b[E] * u[Nx-1+J_b_x[E]+j*(Nx+1)]
                            -L_b[NW] * u[Nx-1+J_b_x[NW]+j*(Nx+1)];

                    for(size_t sum=8; sum<L_b.size(); sum++)
                    {
                        rhs[Nx-1+j*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(j+J_b_y[sum])*(Nx+1)];
                    }
                }

                L_b = stencil.get_L_e(Nx-1,Ny-2,Nx,Ny);
                    
                ndiagL2[Ny-5] = L_b[SE];
                ndiagL1[Ny-4] = L_b[S];
                diagR[Ny-3] = L_b[C];
                ndiagR1[Ny-3] = L_b[N];

                rhs[Nx-1+(Ny-2)*(Nx+1)] = fv[Nx-1+(Ny-2)*(Nx+1)] - L_b[W] * u[Nx-2+(Ny-2)*(Nx+1)] - L_b[E] * u[Nx+(Ny-2)*(Nx+1)]
                  - L_b[NW] * u[Nx-1+J_b_x[NW]+(Ny-2)*(Nx+1)] - L_b[NE] * u[Nx-1+(Ny-2+J_b_y[NE])*(Nx+1)];

                for(size_t sum=8; sum<L_b.size(); sum++)
                {
                  rhs[Nx-1+(Ny-2)*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(Ny-2+J_b_y[sum])*(Nx+1)];
                }
                    
                L_c = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
                J_c_x = stencil.getJx(NE);
                J_c_y = stencil.getJy(NE);

                ndiagL2[Ny-4] = L_c[NE];
                ndiagL1[Ny-3] = L_c[S];
                diagR[Ny-2] = L_c[C];

                rhs[Nx-1+(Ny-1)*(Nx+1)] = fv[Nx-1+(Ny-1)*(Nx+1)] - L_c[W] * u[Nx-2+(Ny-1)*(Nx+1)] 
                  - L_c[E] * u[Nx+(Ny-1)*(Nx+1)] - L_c[NW] * u[Nx-1+J_c_x[NW]+(Ny-1)*(Nx+1)] - L_c[N] * u[Nx-1+(Ny-1+J_c_y[N])*(Nx+1)];

                for(size_t sum=7; sum<L_c.size(); sum++)
                {
                    rhs[Nx-1+(Ny-1)*(Nx+1)] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(Ny-1+J_c_y[sum])*(Nx+1)];
                }

                // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<Ny-2; k++)  
                {
                    ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                    rhs[(k+1)*(Nx+1) + Nx-1] = rhs[(k+1)*(Nx+1) + Nx-1] - ndiagL1[k-1] * rhs[k*(Nx+1)+Nx-1]; 

                    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                    rhs[(k+2)*(Nx+1) + Nx-1] = rhs[(k+2)*(Nx+1) + Nx-1] - ndiagL2[k-1] * rhs[k*(Nx+1)+Nx-1];
                }
                ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
                diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
                rhs[(Ny-1)*(Nx+1) + Nx-1] = rhs[(Ny-1)*(Nx+1) + Nx-1] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ Nx-1];

                // solve the linear system of equations R u = rhs
                u[Nx-1+(Nx+1)*(Ny-1)] = rhs[Nx-1+(Nx+1)*(Ny-1)] / diagR[Ny-2];

                for(size_t k=Ny-2; k>1; k--)
                {
                    u[Nx-1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[Nx-1+k*(Nx+1)] - ndiagR1[k-1] * u[Nx-1+(k+1)*(Nx+1)] );

                    rhs[(k-1)*(Nx+1)+Nx-1] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+Nx-1];
                }
                u[(Nx+1)+Nx-1] = 1/diagR[0] * ( rhs[(Nx+1)+Nx-1] - ndiagR1[0] * u[2*(Nx+1)+Nx-1] );


                // durchlaufe alle inneren Spalten, Spaltenindex i              
                for(size_t i=2; i < Nx-1; i+=2)
                {
                    // setze rechte Seite
                    L_b = stencil.get_L_s(i,1,Nx,Ny);
                    J_b_x = stencil.getJx(S);
                    J_b_y = stencil.getJy(S);

                    diagR[0] = L_b[C];
                    ndiagR1[0] = L_b[N];
                    ndiagR2[0] = L_b[NE];
       
                    rhs[i+Nx+1] = fv[i+Nx+1] - L_b[S] * u[i+(1+J_b_y[S])*(Nx+1)] - L_b[W] * u[i+J_b_x[W]+Nx+1] 
                        - L_b[E] * u[i+J_b_x[E]+Nx+1] - L_b[NW] * u[i+J_b_x[NW]+Nx+1] - L_b[SE] * u[i+J_b_x[SE]+(Nx+1)];
                    
                    for(size_t sum=7; sum<L_b.size(); sum++)
                    {
                        rhs[i+Nx+1] -= L_b[sum] * u[i+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
                    }

                    L = stencil.get_L_c(i,2,Nx,Ny);

                    ndiagL1[0] = L[S];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[N];                  
                    ndiagR2[1] = L[NE];

                    rhs[i+2*(Nx+1)] = fv[i+2*(Nx+1)] - L[W] * u[i+J_x[W]+2*(Nx+1)] - L[E] * u[i+J_x[E]+2*(Nx+1)]
                      -L[NW] * u[i+J_x[NW]+2*(Nx+1)] - L[SE] * u[i+J_x[SE]+2*(Nx+1)] - L[SW] * u[i+(2+J_y[SW])*(Nx+1)];

                    for(size_t sum=9; sum<L.size(); sum++)
                    {
                        rhs[i+2*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(2+J_y[sum])*(Nx+1)];
                    }
                        

                    for(size_t j=3; j<Ny-2; j++)
                    {
                       L = stencil.get_L_c(i,j,Nx,Ny);

                       ndiagL2[j-3] = L[SW];
                       ndiagL1[j-2] = L[S];
                       diagR[j-1] = L[C];
                       ndiagR1[j-1] = L[N];
                       ndiagR2[j-1] = L[NE];
 
                       rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[W] * u[i+J_x[W]+j*(Nx+1)] - L[E] * u[i+J_x[E]+j*(Nx+1)]
                             -L[NW] * u[i+J_x[NW]+j*(Nx+1)] - L[SE] * u[i+J_x[SE]+j*(Nx+1)];
                       
                       for(size_t sum=9; sum<L.size(); sum++)
                       {
                           rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
                       }
                    }

                    L = stencil.get_L_c(i,Ny-2,Nx,Ny);
                                        
                    ndiagL2[Nx-5] = L[SW];
                    ndiagL1[Nx-4] = L[S];
                    diagR[Nx-3] = L[C];
                    ndiagR1[Nx-3] = L[N];

                    rhs[i+(Ny-2)*(Nx+1)] = fv[i+(Ny-2)*(Nx+1)] - L[W] * u[i+J_x[W]+(Ny-2)*(Nx+1)] - L[E] * u[i+J_x[E]+(Ny-2)*(Nx+1)]
                         -L[NW] * u[i+J_x[NW]+(Ny-2)*(Nx+1)] - L[SE] * u[i+J_x[SE]+(Ny-2)*(Nx+1)] - L[NE] * u[i+(Ny-2+J_y[NE])*(Nx+1)];
                    
                    for(size_t sum=9; sum<L.size(); sum++)
                    {
                        rhs[i+(Ny-2)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-2+J_y[sum])*(Nx+1)];
                    }


                    L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
                    J_b_x = stencil.getJx(N);
                    J_b_y = stencil.getJy(N);

                    ndiagL2[Nx-4] = L_b[SW];
                    ndiagL1[Nx-3] = L_b[S];
                    diagR[Nx-2] = L_b[C];

                    rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L_b[N] * u[i+(Ny-1+J_b_y[N])*(Nx+1)] 
                        - L_b[W] * u[i+J_b_x[W]+(Ny-1)*(Nx+1)] - L_b[E] * u[i+J_b_x[E]+(Ny-1)*(Nx+1)] - L_b[NW] * u[i+J_b_x[NW]+(Ny-1)*(Nx+1)]
                        - L_b[SE] * u[i+J_b_x[SE]+(Ny-1)*(Nx+1)];

                    for(size_t sum=9; sum<L_b.size(); sum++)
                    {
                        rhs[i+(Ny-1)*(Nx+1)] -= L_b[sum] * u[i+J_b_x[sum]+(Ny-1+J_b_y[sum])*(Nx+1)];
                    }
                    
                    // LR-decomposition + transformation of the rhs
                    for(size_t k=1; k<Ny-2; k++)  
                    {
                        ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                        diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                        ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                        rhs[(k+1)*(Nx+1) + i] = rhs[(k+1)*(Nx+1) + i] - ndiagL1[k-1] * rhs[k*(Nx+1)+i]; 

                        ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                        ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                        diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                        rhs[(k+2)*(Nx+1) + i] = rhs[(k+2)*(Nx+1) + i] - ndiagL2[k-1] * rhs[k*(Nx+1)+i];
                    }
                    ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
                    diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
                    rhs[(Ny-1)*(Nx+1) + i] = rhs[(Ny-1)*(Nx+1) + i] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ i];

                    // solve the linear system of equations R u = rhs
                    u[i+(Nx+1)*(Ny-1)] = rhs[i+(Nx+1)*(Ny-1)] / diagR[Ny-2];

                    for(size_t k=Ny-2; k>1; k--)
                    {
                        u[i+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[i+k*(Nx+1)] - ndiagR1[k-1] * u[i+(k+1)*(Nx+1)] );

                        rhs[(k-1)*(Nx+1)+i] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+i];
                    }
                    u[(Nx+1)+i] = 1/diagR[0] * ( rhs[(Nx+1)+i] - ndiagR1[0] * u[2*(Nx+1)+i] );
                }

            }
        }
        else //parameter zu klein
        {
            
            for(int k=0; k<2; k++)
            {
                Precision factor = 1.0;

                /**
         * \todo in case of large stencils a four colour RB is needed
         */
        //first do the red points
        //south west corner
        factor = 1.0/stencil.get_center_sw(1,1,Nx,Ny);
        u[1*(Nx+1)+1]+=factor*(fv[1*(Nx+1)+1]
                    -stencil.apply_sw(u,1,1,Nx,Ny));
        //red points on south boarder
        for (size_t i=3;i<(Nx-1);i+=2)
        {
            factor = 1.0/stencil.get_center_s(i,1,Nx,Ny);
            u[1*(Nx+1)+i]+=factor*(fv[1*(Nx+1)+i]
                    -stencil.apply_s(u,i,1,Nx,Ny));
        }
        //south east corner
        factor = 1.0/stencil.get_center_se((Nx-1),1,Nx,Ny);
        u[1*(Nx+1)+(Nx-1)]+=factor*(fv[1*(Nx+1)+(Nx-1)]
                    -stencil.apply_se(u,(Nx-1),1,Nx,Ny));
        for (size_t j=3;j<(Ny-1);j+=2)
        {
            //west boarder point in j. line
            factor = 1.0/stencil.get_center_w(1,j,Nx,Ny);
            u[j*(Nx+1)+1]+=factor*(fv[j*(Nx+1)+1]
                    -stencil.apply_w(u,1,j,Nx,Ny));
            
            for (size_t i=3;i<(Nx-1);i+=2)
            {
                factor = 1.0/stencil.get_center_c(i,j,Nx,Ny);
                u[j*(Nx+1)+i]+=factor*(fv[j*(Nx+1)+i]
                        -stencil.apply_c(u,i,j,Nx,Ny));
            }
            
            //east boarder point in j. line
            factor = 1.0/stencil.get_center_e((Nx-1),j,Nx,Ny);
            u[j*(Nx+1)+(Nx-1)]+=factor*(fv[j*(Nx+1)+(Nx-1)]
                    -stencil.apply_e(u,(Nx-1),j,Nx,Ny));
            
        }
        
        //the missing red points in the center
        for (size_t j=2;j<(Ny-1);j+=2)
        {
            for (size_t i=2;i<(Nx-1);i+=2)
            {
                factor = 1.0/stencil.get_center_c(i,j,Nx,Ny);
                u[j*(Nx+1)+i]+=factor*(fv[j*(Nx+1)+i]
                        -stencil.apply_c(u,i,j,Nx,Ny));
            }
        }
        //north west corner
        
        factor = 1.0/stencil.get_center_nw(1,(Ny-1),Nx,Ny);
        u[(Nx-1)*(Nx+1)+1]+=factor*(fv[(Nx-1)*(Nx+1)+1]
                    -stencil.apply_nw(u,1,(Ny-1),Nx,Ny));
        //red points on north boarder
        for (size_t i=3;i<(Nx-1);i+=2)
        {
            factor = 1.0/stencil.get_center_n(i,(Nx-1),Nx,Ny);
            u[(Nx-1)*(Nx+1)+i]+=factor*(fv[(Nx-1)*(Nx+1)+i]
                    -stencil.apply_n(u,i,(Ny-1),Nx,Ny));
        }
        
        //north east corner
        factor = 1.0/stencil.get_center_ne((Nx-1),(Nx-1),Nx,Ny);
        u[(Nx-1)*(Nx+1)+(Nx-1)]+=factor*(fv[(Nx-1)*(Nx+1)+(Nx-1)]
                -stencil.apply_ne(u,(Nx-1),(Ny-1),Nx,Ny));
        
        //do black points
        //black points on south boarder
        
        for (size_t i=2;i<(Nx-1);i+=2)
        {
             factor = 1.0/stencil.get_center_s(i,1,Nx,Ny);
                u[1*(Nx+1)+i]+=factor*(fv[1*(Nx+1)+i]
                    -stencil.apply_s(u,i,1,Nx,Ny));
        }
        for (size_t j=2;j<(Ny-1);j+=2)
        {
            //west boarder points
            factor = 1.0/stencil.get_center_w(1,j,Nx,Ny);
            u[j*(Nx+1)+1]+=factor*(fv[j*(Nx+1)+1]
                    -stencil.apply_w(u,1,j,Nx,Ny));
            
            for (size_t i=3;i<(Nx-1);i+=2)
            {
                factor = 1.0/stencil.get_center_c(i,j,Nx,Ny);
                u[j*(Nx+1)+i]+=factor*(fv[j*(Nx+1)+i]
                        -stencil.apply_c(u,i,j,Nx,Ny));
            }
            
            //east boarder points
            factor = 1.0/stencil.get_center_e((Nx-1),j,Nx,Ny);
            u[j*(Nx+1)+(Nx-1)]+=factor*(fv[j*(Nx+1)+(Nx-1)]
                    -stencil.apply_e(u,(Nx-1),j,Nx,Ny));
        }
        for (size_t j=3;j<(Ny-1);j+=2)
        {
            for (size_t i=2;i<(Nx-1);i+=2)
            {
                factor = 1.0/stencil.get_center_c(i,j,Nx,Ny);
                u[j*(Nx+1)+i]+=factor*(fv[j*(Nx+1)+i]
                        -stencil.apply_c(u,i,j,Nx,Ny));
            }
        }
        //black points on north boarder
        for (size_t i=2;i<(Nx-1);i+=2)
        {
            factor = 1.0/stencil.get_center_n(i,(Ny-1),Nx,Ny);
            u[(Nx-1)*(Nx+1)+i]+=factor*(fv[(Nx-1)*(Nx+1)+i]
                    -stencil.apply_n(u,i,(Ny-1),Nx,Ny));
        }
            }


        }



}
void ZebraLineGS::xLRSolver(
    std::valarray<Precision>& u,
    const size_t sy,
    std::valarray<Precision>& rhs,
    std::valarray<Precision>& ndiagL,
    std::valarray<Precision>& diagR,
    const std::valarray<Precision>& ndiagR) const
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
void ZebraLineGS::yLRSolver(
    std::valarray<Precision>& u,
    const size_t sx,
    std::valarray<Precision>& rhs,
    std::valarray<Precision>& ndiagL,
    std::valarray<Precision>& diagR,
    const std::valarray<Precision>& ndiagR) const
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
void ZebraLineGS::xLRSolver(
    std::valarray<Precision>& u,
    const size_t sy,
    std::valarray<Precision>& rhs,
    std::valarray<Precision>& ndiagL1,
    std::valarray<Precision>& ndiagL2,
    std::valarray<Precision>& diagR,
    std::valarray<Precision>& ndiagR1,
    const std::valarray<Precision>& ndiagR2) const
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
}
