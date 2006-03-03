/** \file LineJAC.cpp
 * \author Andre Oeckerath
 * \brief LineJAC.cpp contains the implementation of the class LineJAC.
 * \see LineJAC.h
 * \todo check LRSovler calls for corretness
 * \todo clean up ninepointxline, ninepointyline, xline and yline by
 *       doing redudant things in seperate functions
 */
 
#include<iostream>
#include "LineJAC.h"
#include "../functions/residuum.h"

namespace mg
{
void LineJAC::relax(
    NumericArray &u,
    const NumericArray &f, 
    const Stencil &stencil,
    const size_t nx,
    const size_t ny) const
{
    // valarray needed for LR-decomposition of a tridiagonal matrix
    NumericArray resid=residuum(u,f,stencil,nx,ny);
    switch (stencil.size())
    {
        case 1:  // stencil of size 1
        {
            switch (direction_)
            {
                case ALTDIR:
                { 
                    ninepointxline(u,f,resid,stencil,nx,ny);
                    resid=residuum(u,f,stencil,nx,ny);
                    ninepointyline(u,f,resid,stencil,nx,ny);
                    break;
                }
                case XDIR:
                {
                    ninepointxline(u,f,resid,stencil,nx,ny);
                    break;
                }
                case YDIR:
                {
                    ninepointyline(u,f,resid,stencil,nx,ny);
                    break;
                }         
                default:
                {
                    std::cerr<<"Error in direction of the line relaxation!\n";
                    break;
                }
            }
            break;          
        }
        case 2:  // stencil of size 2
        {
            switch (direction_)
            {
                case ALTDIR:
                { 
                    xline(u,f,resid,stencil,nx,ny);
                    resid=residuum(u,f,stencil,nx,ny);
                    yline(u,f,resid,stencil,nx,ny);
                    break;
                }
                case XDIR:
                {
                    xline(u,f,resid,stencil,nx,ny);
                    break;
                }
                case YDIR:
                {
                    yline(u,f,resid,stencil,nx,ny);
                    break;
                }        
                default:
                {
                    std::cerr<<"Error in direction of the line relaxation!\n";
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
void LineJAC::ninepointxline(
    NumericArray &u,
    const NumericArray &f, 
    NumericArray resid,
    const Stencil &stencil,
    const size_t nx, 
    const size_t ny) const
                   
{ 
    NumericArray rhs(0.0,nx-1);
    NumericArray temp(0.0,(nx+1)*(ny+1));       
    //valarrays needed for saving the tridiagonal matrix A of linear system A u = rhs       
    NumericArray diagR(nx-1);
    NumericArray ndiagR(nx-2);
    NumericArray ndiagL(nx-2);
    if(stencil.isConstant())
    {
        // get const operator L
        const NumericArray operatorL=stencil.getL(C,2,2,nx,ny);
        const PositionArray jX=stencil.getJx(C);
        const PositionArray jY=stencil.getJy(C);
        for(size_t sy=1; sy<ny ; sy++) 
        {
            for(size_t sx=0; sx<nx-1; sx++)  
            {
                rhs[sx]=resid[sy*(nx+1)+sx+1];
            }                               
            // set tridiagonalmatrix for solving A u = rhs
            // A[i][i] = L[c]; A[i-1][i] = L[w]; A[i+1][i] = L[e]
            diagR=operatorL[C];
            ndiagR=operatorL[E];
            ndiagL=operatorL[W];
            xLRSolver(temp,sy,nx,rhs,ndiagL,diagR,ndiagR);
        }
        u+=omega_*temp;
    }   
    else
    {
        //Stencil ist not constant, so L needs to be evaluated in each grid point
        //no other change in the algorithm          
        NumericArray operatorL=stencil.getL(C,2,2,nx,ny);
        PositionArray jX=stencil.getJx(C);
        PositionArray jY=stencil.getJy(C);
        if(nx > 2)
        {
            operatorL=stencil.getL(SW,1,1,nx,ny);
            diagR[0]=operatorL[C];
            ndiagR[0]=operatorL[E];
            rhs[0]=resid[nx+1+1];
            for(size_t sx=2; sx<nx-1; sx++)
            {
                operatorL=stencil.getL(S,sx,1,nx,ny);
                diagR[sx-1]=operatorL[C];
                ndiagR[sx-1]=operatorL[E];
                ndiagL[sx-2]=operatorL[W];
                rhs[sx-1]=resid[nx+1+sx];
            }
            operatorL=stencil.getL(SE,nx-1,1,nx,ny);
            diagR[nx-2]=operatorL[C];
            ndiagL[nx-3]=operatorL[W];
            rhs[nx-2]=resid[nx+1+nx-1];
            xLRSolver(temp,1,nx,rhs,ndiagL,diagR,ndiagR);
            for(size_t sy=2; sy<ny-1 ; sy++)
            {
                operatorL=stencil.getL(W,1,sy,nx,ny);
                diagR[0]=operatorL[C];
                ndiagR[0]=operatorL[E];
                rhs[0]=resid[sy*(nx+1)+1];
                for(size_t sx=2; sx<nx-1; sx++)
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    diagR[sx-1]=operatorL[C];
                    ndiagR[sx-1]=operatorL[E];
                    ndiagL[sx-2]=operatorL[W];
                    rhs[sx-1]=resid[sy*(nx+1)+sx];
                }
                operatorL=stencil.getL(E,nx-1,sy,nx,ny);
                diagR[nx-2]=operatorL[C];
                ndiagL[nx-3]=operatorL[W];
                rhs[nx-2]=resid[sy*(nx+1)+nx-1];
                xLRSolver(temp,sy,nx,rhs,ndiagL,diagR,ndiagR);
            }
            operatorL=stencil.getL(NW,1,ny-1,nx,ny);
            diagR[0]=operatorL[C];
            ndiagR[0]=operatorL[E];
            rhs[0]=resid[(ny-1)*(nx+1)+1];
            for(size_t sx=2; sx<nx-1; sx++)
            {
                operatorL=stencil.getL(N,sx,ny-1,nx,ny);
                diagR[sx-1]=operatorL[C];
                ndiagR[sx-1]=operatorL[E];
                ndiagL[sx-2]=operatorL[W];
                rhs[sx-1]=resid[(ny-1)*(nx+1)+sx];
            }
            operatorL=stencil.getL(NE,nx-1,ny-1,nx,ny);
            diagR[nx-2]=operatorL[C];
            ndiagL[nx-3]=operatorL[W];
            rhs[nx-2]=resid[(ny-1)*(nx+1)+nx-1];
            xLRSolver(temp,ny-1,nx,rhs,ndiagL,diagR,ndiagR);
            u+=omega_*temp;
        }
        else  // if Nx and Ny are too small do one GS_lex step
        {
            jacobi_.relax(u,f,stencil,nx,ny);
        }
    }       
}
void LineJAC::ninepointyline(
    NumericArray &u,
    const NumericArray &f, 
    NumericArray resid,
    const Stencil &stencil,
    const size_t nx, 
    const size_t ny) const             
{ 
    NumericArray rhs(0.0,ny+1);
    NumericArray temp(0.0,(nx+1)*(ny+1));
    //valarrays needed for saving the tridiagonal matrix A of linear system A u = rhs
    NumericArray diagR(ny-1);
    NumericArray ndiagR(ny-2);
    NumericArray ndiagL(ny-2);
    if(stencil.isConstant())
    {
        // get const operator L
        const NumericArray operatorL=stencil.getL(C,2,2,nx,ny);
        const PositionArray jX=stencil.getJx(C);
        const PositionArray jY=stencil.getJy(C);
        // for each line: correction of the rhs given by 
        // rhs = fv - [L[w]  0  L[e]] * u and elimination of the 
        // boundary condition in first and last inner point
        for(size_t sx=1; sx<nx ; sx++) 
        {
            for(size_t sy=0; sy<ny-1; sy++)  
            {
                rhs[sy]=resid[(sy+1)*(nx+1)+sx];
            }
            // set tridiagonalmatrix for solving A u = rhs
            // A[i][i] = L[c]; A[i-1][i] = L[s]; A[i+1][i] = L[n]
            diagR=operatorL[C];
            ndiagR=operatorL[N];
            ndiagL=operatorL[S];
            yLRSolver(temp,sx,nx,ny,rhs,ndiagL,diagR,ndiagR);
        }
        u += omega_ * temp;
    }
    else
    {
        //Stencil ist not constant, so L needs to be evaluated in each grid point
        //no other change in the algorithm  
        NumericArray operatorL=stencil.getL(C,2,2,nx,ny);
        PositionArray jX=stencil.getJx(C);
        PositionArray jY=stencil.getJy(C);
        if(ny > 2)
        {
            operatorL=stencil.getL(SW,1,1,nx,ny);                
            diagR[0]=operatorL[C];
            ndiagR[0]=operatorL[N];           
            rhs[0]=resid[nx+1+1];
            for(size_t sy=2; sy<ny-1; sy++) 
            {
                operatorL=stencil.getL(W,1,sy,nx,ny);
                diagR[sy-1]=operatorL[C];
                ndiagR[sy-1]=operatorL[N];
                ndiagL[sy-2]=operatorL[S];
                rhs[sy-1]=resid[sy*(nx+1)+1];                   
            }
            operatorL=stencil.getL(NW,1,ny-1,nx,ny);
            diagR[ny-2]=operatorL[C];
            ndiagL[ny-3]=operatorL[S]; 
            rhs[ny-2]=resid[(ny-1)*(nx+1)+1];
            yLRSolver(temp,1,nx,ny,rhs,ndiagL,diagR,ndiagR);
            for(size_t sx=2; sx<nx-1 ; sx++)
            {
                operatorL=stencil.getL(S,sx,1,nx,ny);
                diagR[0]=operatorL[C];
                ndiagR[0]=operatorL[N];   
                rhs[0]=resid[nx+1+sx];
                for(size_t sy=2; sy<ny-1; sy++) 
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    diagR[sy-1]=operatorL[C];
                    ndiagR[sy-1]=operatorL[N];
                    ndiagL[sy-2]=operatorL[S];
                    rhs[sy-1]=resid[sy*(nx+1)+sx];
                }
                operatorL=stencil.getL(N,sx,ny-1,nx,ny);
                diagR[ny-2]=operatorL[C];
                ndiagL[ny-3]=operatorL[S];
                rhs[ny-2]=resid[(ny-1)*(nx+1)+sx];
                yLRSolver(temp,sx,nx,ny,rhs,ndiagL,diagR,ndiagR);
            }
            operatorL=stencil.getL(SE,nx-1,1,nx,ny);
            diagR[0]=operatorL[C];
            ndiagR[0]=operatorL[N];   
            rhs[0]=resid[nx+1+nx-1];
            for(size_t sy=2; sy<ny-1; sy++) 
            {
                operatorL=stencil.getL(E,nx-1,sy,nx,ny);
                diagR[sy-1]=operatorL[C];
                ndiagR[sy-1]=operatorL[N];
                ndiagL[sy-2]=operatorL[S];
                rhs[sy-1]=resid[sy*(nx+1)+nx-1];                    
            }
            operatorL=stencil.getL(NE,nx-1,ny-1,nx,ny);
            diagR[ny-2]=operatorL[C];
            ndiagL[ny-3]=operatorL[S];
            rhs[ny-2]=resid[(ny-1)*(nx+1)+nx-1];
            yLRSolver(temp,nx-1,nx,ny,rhs,ndiagL,diagR,ndiagR);
            u += omega_ * temp;
        }
        else // if Nx and Ny are too small do one GS_lex step
        {
            jacobi_.relax(u,f,stencil,nx,ny);
        }
    }       
}
void LineJAC::xline(
    NumericArray &u,
    const NumericArray &f, 
    NumericArray resid,
    const Stencil &stencil,
    const size_t nx, 
    const size_t ny) const

{
    if((ny > 4) && (nx > 4))
    {
        NumericArray rhs(0.0,nx-1);
        NumericArray temp(0.0,(nx+1)*(ny+1));
        NumericArray diagR(0.0,nx-1);
        NumericArray ndiagR1(0.0,nx-2);
        NumericArray ndiagL1(0.0,nx-2);
        NumericArray ndiagR2(0.0,nx-3);
        NumericArray ndiagL2(0.0,nx-3);
        if(stencil.isConstant())
        {
            // get const operator L
            const NumericArray operatorL=stencil.getL(C,2,2,nx,ny);
            const PositionArray jX=stencil.getJx(C);
            const PositionArray jY=stencil.getJy(C);
            NumericArray operatorLB=stencil.getL(S,2,1,nx,ny);
            PositionArray jXB=stencil.getJx(S);
            PositionArray jYB=stencil.getJy(S);
            NumericArray operatorLC=stencil.getL(SW,1,1,nx,ny);
            PositionArray jXC=stencil.getJx(SW);
            PositionArray jYC=stencil.getJy(SW);
            // set rhs for line 1                   
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[E];
            ndiagR2[0]=operatorLC[NE];         
            rhs[0]=resid[nx+1+1];
            ndiagL1[0]=operatorLB[W];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[E];                    
            ndiagR2[1]=operatorLB[SE];
            rhs[1]=resid[nx+1+2];
            for(size_t sx=3; sx<nx-2; sx++)  
            {
                ndiagL2[sx-3]=operatorLB[NW];
                ndiagL1[sx-2]=operatorLB[W];
                diagR[sx-1]=operatorLB[C];
                ndiagR1[sx-1]=operatorLB[E];                  
                ndiagR2[sx-1]=operatorLB[SE]; 
                rhs[sx-1]=resid[nx+1+sx];
            }
            ndiagL2[nx-5]=operatorLB[NW];
            ndiagL1[nx-4]=operatorLB[W];
            diagR[nx-3]=operatorLB[C];
            ndiagR1[nx-3]=operatorLB[E];
            rhs[nx-3]=resid[nx+1+nx-2];
            
            operatorLC=stencil.getL(SE,nx-1,1,nx,ny);
            jXC=stencil.getJx(SE);
            jYC=stencil.getJy(SE);
            ndiagL2[nx-4]=operatorLC[NW];
            ndiagL1[nx-3]=operatorLC[W];
            diagR[nx-2]=operatorLC[C];
            rhs[nx-2]=resid[nx+1+nx-1];
            xLRSolver(temp,1,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
////////////////////////////////////////////////////////////////////////////////
            // process all inner lines
            for(size_t sy=2; sy<ny-1; sy++)
            {
                // set rhs
                operatorLB=stencil.getL(W,1,sy,nx,ny);
                jXB=stencil.getJx(W);
                jYB=stencil.getJy(W);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[E];
                ndiagR2[0]=operatorLB[NE];
                rhs[0]=resid[sy*(nx+1)+1];
                ndiagL1[0]=operatorL[W];
                diagR[1]=operatorL[C];
                ndiagR1[1]=operatorL[E];                  
                ndiagR2[1]=operatorL[SE];
                rhs[1]=resid[sy*(nx+1)+2];
                for(size_t sx=3; sx<nx-2; sx++)  
                {
                    ndiagL2[sx-3]=operatorL[NW];
                    ndiagL1[sx-2]=operatorL[W];
                    diagR[sx-1]=operatorL[C];
                    ndiagR1[sx-1]=operatorL[E];                    
                    ndiagR2[sx-1]=operatorL[SE];
                    rhs[sx-1]=resid[sy*(nx+1)+sx];
                }
                ndiagL2[nx-5]=operatorL[NW];
                ndiagL1[nx-4]=operatorL[W];
                diagR[nx-3]=operatorL[C];
                ndiagR1[nx-3]=operatorL[E];
                rhs[nx-3]=resid[sy*(nx+1)+nx-2];
                operatorLB=stencil.getL(E,nx-1,sy,nx,ny);
                jXB=stencil.getJx(E);
                jYB=stencil.getJy(E);
                ndiagL2[nx-4]=operatorLB[NW];
                ndiagL1[nx-3]=operatorLB[W];
                diagR[nx-2]=operatorLB[C];
                rhs[nx-2]=resid[sy*(nx+1)+nx-1];
                xLRSolver(temp,sy,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            }
////////////////////////////////////////////////////////////////////////////////
            // set rhs in top line
            operatorLC=stencil.getL(NW,1,ny-1,nx,ny);
            jXC=stencil.getJx(NW);
            jYC=stencil.getJy(NW);
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[E];
            ndiagR2[0]=operatorLC[NW];         
            rhs[0]=resid[(ny-1)*(nx+1)+1];
            operatorLB=stencil.getL(N,2,ny-1,nx,ny);
            jXB=stencil.getJx(N);
            jYB=stencil.getJy(N);
            ndiagL1[0]=operatorLB[W];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[E];                    
            ndiagR2[1]=operatorLB[NE];
            rhs[1]=resid[(ny-1)*(nx+1)+2];
            for(size_t sx=3; sx<nx-2; sx++)  
            {
                ndiagL2[sx-3]=operatorLB[NW];
                ndiagL1[sx-2]=operatorLB[W];
                diagR[sx-1]=operatorLB[C];
                ndiagR1[sx-1]=operatorLB[E];                  
                ndiagR2[sx-1]=operatorLB[NE];
                rhs[sx-1]=resid[(ny-1)*(nx+1)+sx];
            }
            ndiagL2[nx-5]=operatorLB[NW];
            ndiagL1[nx-4]=operatorLB[W];
            diagR[nx-3]=operatorLB[C];
            ndiagR1[nx-3]=operatorLB[E];
            rhs[nx-3]=resid[(ny-1)*(nx+1)+nx-2];
            operatorLC=stencil.getL(NE,nx-1,ny-1,nx,ny);
            jXC=stencil.getJx(NE);
            jYC=stencil.getJy(NE);
            ndiagL2[nx-4]=operatorLC[NW];
            ndiagL1[nx-3]=operatorLC[W];
            diagR[nx-2]=operatorLC[C];
            rhs[nx-2]=resid[(ny-1)*(nx+1)+nx-1];
            xLRSolver(temp,ny-1,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            u += omega_ * temp;
        }
        else // not constant
        {
            NumericArray operatorL=stencil.getL(C,2,2,nx,ny);
            PositionArray jX=stencil.getJx(C);
            PositionArray jY=stencil.getJy(C);    
            NumericArray operatorLB=stencil.getL(S,2,1,nx,ny);
            PositionArray jXB=stencil.getJx(S);
            PositionArray jYB=stencil.getJy(S);
            NumericArray operatorLC=stencil.getL(SW,1,1,nx,ny);
            PositionArray jXC=stencil.getJx(SW);
            PositionArray jYC=stencil.getJy(SW);
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[E];
            ndiagR2[0]=operatorLC[NE];
            rhs[0]=resid[nx+1+1];
            ndiagL1[0]=operatorLB[W];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[E];                    
            ndiagR2[1]=operatorLB[SE];
            rhs[1]=resid[nx+1+2];
            for(size_t sx=3; sx<nx-2; sx++)  
            {
                operatorLB=stencil.getL(S,sx,1,nx,ny);   
                ndiagL2[sx-3]=operatorLB[NW];
                ndiagL1[sx-2]=operatorLB[W];
                diagR[sx-1]=operatorLB[C];
                ndiagR1[sx-1]=operatorLB[E];                  
                ndiagR2[sx-1]=operatorLB[SE];
                rhs[sx-1]=resid[nx+1+sx];
            }
            operatorLB=stencil.getL(S,nx-2,1,nx,ny);
            ndiagL2[nx-5]=operatorLB[NW];
            ndiagL1[nx-4]=operatorLB[W];
            diagR[nx-3]=operatorLB[C];
            ndiagR1[nx-3]=operatorLB[E];
            rhs[nx-3]=resid[nx+1+nx-2];
            operatorLC=stencil.getL(SE,nx-1,1,nx,ny);
            jXC=stencil.getJx(SE);
            jYC=stencil.getJy(SE);
            ndiagL2[nx-4]=operatorLC[NW];
            ndiagL1[nx-3]=operatorLC[W];
            diagR[nx-2]=operatorLC[C];
            rhs[nx-2]=resid[nx+1+nx-1];
            xLRSolver(temp,1,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
////////////////////////////////////////////////////////////////////////////////
            // process all inner lines                
            for(size_t sy=2; sy<ny-1; sy++)
            {
                // set rhs
                operatorLB=stencil.getL(W,1,sy,nx,ny);
                jXB=stencil.getJx(W);
                jYB=stencil.getJy(W);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[E];
                ndiagR2[0]=operatorLB[NE];     
                rhs[0]=resid[sy*(nx+1)+1];             
                operatorL=stencil.getL(C,2,sy,nx,ny);
                ndiagL1[0]=operatorL[W];
                diagR[1]=operatorL[C];
                ndiagR1[1]=operatorL[E];                  
                ndiagR2[1]=operatorL[SE];
                rhs[1]=resid[sy*(nx+1)+2];
                for(size_t sx=3; sx<nx-2; sx++)  
                {
                    operatorL=stencil.getL(C,sx,sy,nx,ny);
                    ndiagL2[sx-3]=operatorL[NW];
                    ndiagL1[sx-2]=operatorL[W];
                    diagR[sx-1]=operatorL[C];
                    ndiagR1[sx-1]=operatorL[E];                    
                    ndiagR2[sx-1]=operatorL[SE];
                    rhs[sx-1]=resid[sy*(nx+1)+sx];
                }
                operatorL=stencil.getL(C,nx-2,sy,nx,ny);
                ndiagL2[nx-5]=operatorL[NW];
                ndiagL1[nx-4]=operatorL[W];
                diagR[nx-3]=operatorL[C];
                ndiagR1[nx-3]=operatorL[E];
                rhs[nx-3]=resid[sy*(nx+1)+nx-2];
                operatorLB=stencil.getL(E,nx-1,sy,nx,ny);
                jXB=stencil.getJx(E);
                jYB=stencil.getJy(E);
                ndiagL2[nx-4]=operatorLB[NW];
                ndiagL1[nx-3]=operatorLB[W];
                diagR[nx-2]=operatorLB[C];
                rhs[nx-2]=resid[sy*(nx+1)+nx-1];
                xLRSolver(temp,sy,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            }
////////////////////////////////////////////////////////////////////////////////
            // set rhs in top line                
            operatorLC=stencil.getL(NW,1,ny-1,nx,ny);
            jXC=stencil.getJx(NW);
            jYC=stencil.getJy(NW);
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[E];
            ndiagR2[0]=operatorLC[NW];         
            rhs[0]=resid[(ny-1)*(nx+1)+1];
            operatorLB=stencil.getL(N,2,ny-1,nx,ny);
            jXB=stencil.getJx(N);
            jYB=stencil.getJy(N);
            ndiagL1[0]=operatorLB[W];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[E];                    
            ndiagR2[1]=operatorLB[NE];
            rhs[1]=resid[(ny-1)*(nx+1)+2];
            for(size_t sx=3; sx<nx-2; sx++)  
            {
                operatorLB=stencil.getL(N,sx,ny-1,nx,ny);
                ndiagL2[sx-3]=operatorLB[NW];
                ndiagL1[sx-2]=operatorLB[W];
                diagR[sx-1]=operatorLB[C];
                ndiagR1[sx-1]=operatorLB[E];                  
                ndiagR2[sx-1]=operatorLB[NE];
                rhs[sx-1]=resid[(ny-1)*(nx+1)+sx];
            }
            operatorLB=stencil.getL(N,nx-2,ny-1,nx,ny);
            ndiagL2[nx-5]=operatorLB[NW];
            ndiagL1[nx-4]=operatorLB[W];
            diagR[nx-3]=operatorLB[C];
            ndiagR1[nx-3]=operatorLB[E];
            rhs[nx-3]=resid[(ny-1)*(nx+1)+nx-2];
            operatorLC=stencil.getL(NE,nx-1,ny-1,nx,ny);
            jXC=stencil.getJx(NE);
            jYC=stencil.getJy(NE);
            ndiagL2[nx-4]=operatorLC[NW];
            ndiagL1[nx-3]=operatorLC[W];
            diagR[nx-2]=operatorLC[C];
            rhs[nx-2]=resid[(ny-1)*(nx+1)+nx-1];
            xLRSolver(temp,ny-1,nx,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            u += omega_ * temp;
        }
    }
    else 
    {
        for(int k=0; k<2; k++)
        {
            jacobi_.relax(u,f,stencil,nx,ny);
        }
    }
}
void LineJAC::yline(
    NumericArray &u,
    const NumericArray &f, 
    NumericArray resid,
    const Stencil &stencil,
    const size_t nx, 
    const size_t ny) const
{
    if((ny > 4) && (nx > 4))
    {
        NumericArray rhs(0.0,ny-1);
        NumericArray temp(0.0,(nx+1)*(ny+1));
        NumericArray diagR(0.0,ny-1);
        NumericArray ndiagR1(0.0,ny-2);
        NumericArray ndiagL1(0.0,ny-2);
        NumericArray ndiagR2(0.0,ny-3);
        NumericArray ndiagL2(0.0,ny-3);        
        if(stencil.isConstant())
        {
            // get const operator L
            const NumericArray operatorL=stencil.getL(C,2,2,nx,ny);
            const PositionArray jX=stencil.getJx(C);
            const PositionArray jY=stencil.getJy(C);
            NumericArray operatorLB=stencil.getL(W,1,2,nx,ny);
            PositionArray jXB=stencil.getJx(W);
            PositionArray jYB=stencil.getJy(W);
            NumericArray operatorLC=stencil.getL(SW,1,1,nx,ny);
            PositionArray jXC=stencil.getJx(SW);
            PositionArray jYC=stencil.getJy(SW);
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[N];
            ndiagR2[0]=operatorLC[NW];
            rhs[0]=resid[nx+1+1];
            ndiagL1[0]=operatorLB[S];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[N];                    
            ndiagR2[1]=operatorLB[NW];
            rhs[1]=resid[2*(nx+1)+1];
            for(size_t sy=3; sy<ny-2; sy++)  
            {
                ndiagL2[sy-3]=operatorLB[SE];
                ndiagL1[sy-2]=operatorLB[S];
                diagR[sy-1]=operatorLB[C];
                ndiagR1[sy-1]=operatorLB[N];                  
                ndiagR2[sy-1]=operatorLB[NW];
                rhs[sy-1]=resid[sy*(nx+1)+1];
            }
            ndiagL2[ny-5]=operatorLB[SE];
            ndiagL1[ny-4]=operatorLB[S];
            diagR[ny-3]=operatorLB[C];
            ndiagR1[ny-3]=operatorLB[N];
            rhs[ny-3]=resid[(ny-2)*(nx+1)+1];
            operatorLC=stencil.getL(NW,1,ny-1,nx,ny);
            jXC=stencil.getJx(NW);
            jYC=stencil.getJy(NW);
            ndiagL2[ny-4]=operatorLC[NE];
            ndiagL1[ny-3]=operatorLC[S];
            diagR[ny-2]=operatorLC[C];
            rhs[ny-2]=resid[(ny-1)*(nx+1)+1];
            yLRSolver(temp,1,nx,ny,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
////////////////////////////////////////////////////////////////////////////////
            // process all inner lines             
            for(size_t sx=2; sx<nx-1; sx++)
            {
                // set rhs
                operatorLB=stencil.getL(S,sx,1,nx,ny);
                jXB=stencil.getJx(S);
                jYB=stencil.getJy(S);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[N];
                ndiagR2[0]=operatorLB[NE];
                rhs[0]=resid[nx+1+sx];
                ndiagL1[0]=operatorL[S];
                diagR[1]=operatorL[C];
                ndiagR1[1]=operatorL[N];                  
                ndiagR2[1]=operatorL[NE];
                rhs[1]=resid[2*(nx+1)+sx];                     
                for(size_t sy=3; sy<ny-2; sy++)
                {
                   ndiagL2[sy-3]=operatorL[SW];
                   ndiagL1[sy-2]=operatorL[S];
                   diagR[sy-1]=operatorL[C];
                   ndiagR1[sy-1]=operatorL[N];
                   ndiagR2[sy-1]=operatorL[NE];
                   rhs[sy-1]=resid[sy*(nx+1)+sx];                       
                }
                ndiagL2[nx-5]=operatorL[SW];
                ndiagL1[nx-4]=operatorL[S];
                diagR[nx-3]=operatorL[C];
                ndiagR1[nx-3]=operatorL[N];
                rhs[ny-3]=resid[(ny-2)*(nx+1)+sx];
                operatorLB=stencil.getL(N,sx,ny-1,nx,ny);
                jXB=stencil.getJx(N);
                jYB=stencil.getJy(N);
                ndiagL2[nx-4]=operatorLB[SW];
                ndiagL1[nx-3]=operatorLB[S];
                diagR[nx-2]=operatorLB[C];
                rhs[ny-2]=resid[(ny-1)*(nx+1)+sx];
                yLRSolver(temp,sx,nx,ny,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
            }
////////////////last column//////////////////
            operatorLC=stencil.getL(SE,nx-1,1,nx,ny);
            jXC=stencil.getJx(SE);
            jYC=stencil.getJy(SE);
            rhs[0]=resid[nx+1+nx-1];
            operatorLB=stencil.getL(E,nx-1,2,nx,ny);
            jXB=stencil.getJx(E);
            jYB=stencil.getJy(E);
            ndiagL1[0]=operatorLB[S];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[N];                    
            ndiagR2[1]=operatorLB[NE];
            rhs[1]=resid[2*(nx+1)+nx-1];
            for(size_t sy=3; sy<ny-2; sy++)  
            {
                ndiagL2[sy-3]=operatorLB[SE];
                ndiagL1[sy-2]=operatorLB[S];
                diagR[sy-1]=operatorLB[C];
                ndiagR1[sy-1]=operatorLB[N];                  
                ndiagR2[sy-1]=operatorLB[NE];
                rhs[sy-1]=resid[sy*(nx+1)+nx-1];
            }
            ndiagL2[ny-5]=operatorLB[SE];
            ndiagL1[ny-4]=operatorLB[S];
            diagR[ny-3]=operatorLB[C];
            ndiagR1[ny-3]=operatorLB[N];
            rhs[ny-3]=resid[(ny-2)*(nx+1)+nx-1];
            operatorLC=stencil.getL(NE,nx-1,ny-1,nx,ny);
            jXC=stencil.getJx(NE);
            jYC=stencil.getJy(NE);
            ndiagL2[ny-4]=operatorLC[NE];
            ndiagL1[ny-3]=operatorLC[S];
            diagR[ny-2]=operatorLC[C];
            rhs[ny-2]=resid[(ny-1)*(nx+1)+nx-1];
            yLRSolver(
                temp,
                nx-1,nx,ny,
                rhs,
                ndiagL1,ndiagL2,
                diagR,ndiagR1,ndiagR2);
            u += omega_ * temp;
        }
        else // stencil not constant
        {
            NumericArray operatorL=stencil.getL(C,2,2,nx,ny);
            PositionArray jX=stencil.getJx(C);
            PositionArray jY=stencil.getJy(C);                
            NumericArray operatorLB=stencil.getL(W,1,2,nx,ny);
            PositionArray jXB=stencil.getJx(W);
            PositionArray jYB=stencil.getJy(W);
            NumericArray operatorLC=stencil.getL(SW,1,1,nx,ny);
            PositionArray jXC=stencil.getJx(SW);
            PositionArray jYC=stencil.getJy(SW);
            diagR[0]=operatorLC[C];
            ndiagR1[0]=operatorLC[N];
            ndiagR2[0]=operatorLC[NW];
            rhs[0]=resid[nx+1+1];
            ndiagL1[0]=operatorLB[S];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[N];                    
            ndiagR2[1]=operatorLB[NW];
            rhs[1]=resid[2*(nx+1)+1];
            for(size_t sy=3; sy<ny-2; sy++)  
            {
                operatorLB=stencil.getL(W,1,sy,nx,ny);
                ndiagL2[sy-3]=operatorLB[SE];
                ndiagL1[sy-2]=operatorLB[S];
                diagR[sy-1]=operatorLB[C];
                ndiagR1[sy-1]=operatorLB[N];                  
                ndiagR2[sy-1]=operatorLB[NW];
                rhs[sy-1]=resid[sy*(nx+1)+1];
            }
            operatorLB=stencil.getL(W,1,ny-2,nx,ny);
            ndiagL2[ny-5]=operatorLB[SE];
            ndiagL1[ny-4]=operatorLB[S];
            diagR[ny-3]=operatorLB[C];
            ndiagR1[ny-3]=operatorLB[N];
            rhs[ny-3]=resid[(ny-2)*(nx+1)+1];
            operatorLC=stencil.getL(NW,1,ny-1,nx,ny);
            jXC=stencil.getJx(NW);
            jYC=stencil.getJy(NW);
            ndiagL2[ny-4]=operatorLC[NE];
            ndiagL1[ny-3]=operatorLC[S];
            diagR[ny-2]=operatorLC[C];
            rhs[ny-2]=resid[(ny-1)*(nx+1)+1];
            yLRSolver(temp,1,nx,ny,rhs,ndiagL1,ndiagL2,diagR,ndiagR1,ndiagR2);
////////////////////////////////////////////////////////////////////////////////
            // process all inner columns  
            for(size_t sx=2; sx<nx-1; sx++)
            {
                // set rhs
                operatorLB=stencil.getL(S,sx,1,nx,ny);
                jXB=stencil.getJx(S);
                jYB=stencil.getJy(S);
                diagR[0]=operatorLB[C];
                ndiagR1[0]=operatorLB[N];
                ndiagR2[0]=operatorLB[NE];
                rhs[0]=resid[nx+1+sx];
                operatorL=stencil.getL(C,sx,2,nx,ny);
                ndiagL1[0]=operatorL[S];
                diagR[1]=operatorL[C];
                ndiagR1[1]=operatorL[N];                  
                ndiagR2[1]=operatorL[NE];
                rhs[1]=resid[2*(nx+1)+sx];  
                for(size_t sy=3; sy<ny-2; sy++)
                {
                   operatorL=stencil.getL(C,sx,sy,nx,ny);
                   ndiagL2[sy-3]=operatorL[SW];
                   ndiagL1[sy-2]=operatorL[S];
                   diagR[sy-1]=operatorL[C];
                   ndiagR1[sy-1]=operatorL[N];
                   ndiagR2[sy-1]=operatorL[NE];
                   rhs[sy-1]=resid[sy*(nx+1)+sx];
                }       
                operatorL=stencil.getL(C,sx,ny-2,nx,ny);
                ndiagL2[nx-5]=operatorL[SW];
                ndiagL1[nx-4]=operatorL[S];
                diagR[nx-3]=operatorL[C];
                ndiagR1[nx-3]=operatorL[N];
                rhs[ny-3]=resid[(ny-2)*(nx+1)+sx];
                operatorLB=stencil.getL(N,sx,ny-1,nx,ny);
                jXB=stencil.getJx(N);
                jYB=stencil.getJy(N);
                ndiagL2[nx-4]=operatorLB[SW];
                ndiagL1[nx-3]=operatorLB[S];
                diagR[nx-2]=operatorLB[C];
                rhs[ny-2]=resid[(ny-1)*(nx+1)+sx];
                yLRSolver(
                    temp,
                    sx,nx,ny,
                    rhs,
                    ndiagL1,ndiagL2,
                    diagR,ndiagR1,ndiagR2);
            }
////////////////last column//////////////////
            operatorLC=stencil.getL(SE,nx-1,1,nx,ny);
            jXC=stencil.getJx(SE);
            jYC=stencil.getJy(SE);
            rhs[0]=resid[nx+1+nx-1];
            operatorLB=stencil.getL(E,nx-1,2,nx,ny);
            jXB=stencil.getJx(E);
            jYB=stencil.getJy(E);
            ndiagL1[0]=operatorLB[S];
            diagR[1]=operatorLB[C];
            ndiagR1[1]=operatorLB[N];                    
            ndiagR2[1]=operatorLB[NE];
            rhs[1]=resid[2*(nx+1)+nx-1];
            for(size_t sy=3; sy<ny-2; sy++)  
            {
                operatorLB=stencil.getL(E,nx-1,sy,nx,ny);
                ndiagL2[sy-3]=operatorLB[SE];
                ndiagL1[sy-2]=operatorLB[S];
                diagR[sy-1]=operatorLB[C];
                ndiagR1[sy-1]=operatorLB[N];                  
                ndiagR2[sy-1]=operatorLB[NE];
                rhs[sy-1]=resid[sy*(nx+1)+nx-1];
            }
            operatorLB=stencil.getL(E,nx-1,ny-2,nx,ny);
            ndiagL2[ny-5]=operatorLB[SE];
            ndiagL1[ny-4]=operatorLB[S];
            diagR[ny-3]=operatorLB[C];
            ndiagR1[ny-3]=operatorLB[N];
            rhs[ny-3]=resid[(ny-2)*(nx+1)+nx-1];
            operatorLC=stencil.getL(NE,nx-1,ny-1,nx,ny);
            jXC=stencil.getJx(NE);
            jYC=stencil.getJy(NE);
            ndiagL2[ny-4]=operatorLC[NE];
            ndiagL1[ny-3]=operatorLC[S];
            diagR[ny-2]=operatorLC[C];
            rhs[ny-2]=resid[(ny-1)*(nx+1)+nx-1];
            yLRSolver(
                temp,
                nx-1,nx,ny,
                rhs,
                ndiagL1,ndiagL2,
                diagR,ndiagR1,ndiagR2);
            u += omega_ * temp;
        }
    }
    else //parameter zu klein
    {     
        for(int k=0; k<2; k++)
        {
            jacobi_.relax(u,f,stencil,nx,ny);
        }
    }
}
}
