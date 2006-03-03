/** \file ZebraLine.cpp
 * \author Andre Oeckerath
 * \brief ZebraLine.cpp contains the implementation of the class ZebraLine.
 * \see ZebraLine.h
 * \todo check LRSovler calls for corretness
 * \todo clean up ninepointxzebra, ninepointyzebra, xzebra and yzebra by
 *       doing redudant things in seperate functions
 */
#include "ZebraLine.h"
#include<iostream>
#include "../functions/residuum.h"


namespace mg
{			
void ZebraLine::relax(
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
                    ninepointxzebra(u,f,resid,stencil,nx,ny);
                    resid=residuum(u,f,stencil,nx,ny);
                    ninepointyzebra(u,f,resid,stencil,nx,ny);
                    break;
                }
    		    case XDIR:
                {
                    ninepointxzebra(u,f,resid,stencil,nx,ny);
                    break;
                }
                case YDIR:
                {
                    ninepointyzebra(u,f,resid,stencil,nx,ny);
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
                    xzebra(u,f,resid,stencil,nx,ny);
                    resid=residuum(u,f,stencil,nx,ny);
                    yzebra(u,f,resid,stencil,nx,ny);
                    break;
                }
                case XDIR:
                {
                    xzebra(u,f,resid,stencil,nx,ny);
                    break;
                }
                case YDIR:
                {
                    yzebra(u,f,resid,stencil,nx,ny);
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
void ZebraLine::ninepointxzebra(
    NumericArray &u,
    const NumericArray &f, 
    NumericArray resid,
    const Stencil &stencil,
    const size_t nx, 
    const size_t ny) const
{ 
    NumericArray rhs(0.0,nx+1);
    NumericArray temp(0.0,(nx+1)*(ny+1));
    
    //valarrays needed for saving the tridiagonal matrix A of linear system A u = rhs
    NumericArray diagR(nx-1);
    NumericArray ndiagR(nx-2);
    NumericArray ndiagL(nx-2);
    
    if(stencil.isConstant() == true)
    {
        // get const operator L
        const NumericArray L = stencil.get_L_c(2,2,nx,ny);
        const PositionArray J_x = stencil.getJx(C);
        const PositionArray J_y = stencil.getJy(C);
        
        // for each line: correction of the rhs given by rhs = fv - [L[n]  0  L[s]]^t * u and elimination of the 
        // boundary condition in first and last inner point
        for(size_t i=1; i<ny ; i+=2) 
        {
            for(size_t j=0; j<nx-1; j++)  
            {
                rhs[j] = resid[i*(nx+1)+j+1];
            }   
                                            
            // set tridiagonalmatrix for solving A u = rhs
            // A[i][i] = L[c]; A[i-1][i] = L[w]; A[i+1][i] = L[e]
            diagR = L[C];
            ndiagR = L[E];
            ndiagL = L[W];

            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<nx-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }

            // solve the linear system of equations R u = rhs
            temp[i*(nx+1)+(nx-1)] = rhs[nx-2] / diagR[nx-2];
            
            for(size_t j=nx-2; j>0; j--)
            {
                temp[i*(nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR[j-1] * temp[i*(nx+1)+j+1] );
            }
        }

        u += omega_ * temp;

        resid = residuum(u,f,stencil,nx,ny);

        temp = 0.0;

        // same for even lines 
        for(size_t i=2; i<ny ; i+=2)
        {
            for(size_t j=0; j<nx-1; j++)  
            {
                rhs[j] = resid[i*(nx+1)+j+1];
            }   

            diagR = L[C];
            ndiagR = L[E];
            ndiagL = L[W];

            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<nx-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }

            // solve the linear system of equations R u = rhs
            temp[i*(nx+1)+(nx-1)] = rhs[nx-2] / diagR[nx-2];
            
            for(size_t j=nx-2; j>0; j--)
            {
                temp[i*(nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR[j-1] * temp[i*(nx+1)+j+1] );
            }
        }

        u += omega_ * temp;


    }

    else
    {
        //Stencil ist not constant, so L needs to be evaluated in each grid point
        //no other change in the algorithm  
        NumericArray L = stencil.get_L_c(2,2,nx,ny);
        PositionArray J_x = stencil.getJx(C);
        PositionArray J_y = stencil.getJy(C);

        if(nx > 2)
        {
            L = stencil.get_L_sw(1,1,nx,ny);
                
            diagR[0] = L[C];
            ndiagR[0] = L[E];

            rhs[0] = resid[nx+1+1];

            for(size_t j=2; j<nx-1; j++)
            {
                L = stencil.get_L_s(j,1,nx,ny);
                diagR[j-1] = L[C];
                ndiagR[j-1] = L[E];
                ndiagL[j-2] = L[W];
                 
                rhs[j-1] = resid[nx+1+j];
            }

            L = stencil.get_L_se(nx-1,1,nx,ny);
            
            diagR[nx-2] = L[C];
            ndiagL[nx-3] = L[W];
                
            rhs[nx-2] = resid[nx+1+nx-1];

            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<nx-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }

            // solve the linear system of equations R u = rhs
            temp[nx+1+nx-1] = rhs[nx-2] / diagR[nx-2];
            
            for(size_t j=nx-2; j>0; j--)
            {
                temp[nx+1+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR[j-1] * temp[nx+1+j+1] );
            }               
            

            for(size_t i=3; i<ny-1 ; i+=2)
            {
                L = stencil.get_L_w(1,i,nx,ny);
                
                diagR[0] = L[C];
                ndiagR[0] = L[E];

                rhs[0] = resid[i*(nx+1)+1];

                for(size_t j=2; j<nx-1; j++)
                {
                    L = stencil.get_L_c(j,i,nx,ny);
                    diagR[j-1] = L[C];
                    ndiagR[j-1] = L[E];
                    ndiagL[j-2] = L[W];
                    
                    rhs[j-1] = resid[i*(nx+1)+j];
                }

                L = stencil.get_L_e(nx-1,i,nx,ny);
                
                diagR[nx-2] = L[C];
                ndiagL[nx-3] = L[W];
                
                rhs[nx-2] = resid[i*(nx+1)+nx-1];


                // LR-decomposition + transformation of the rhs vector
                for(size_t k=1; k<nx-1; k++)  
                {
                    ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                    diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                    rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
                }

                // solve the linear system of equations R u = rhs
                temp[i*(nx+1)+nx-1] = rhs[nx-2] / diagR[nx-2];
                
                for(size_t j=nx-2; j>0; j--)
                {
                    temp[i*(nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR[j-1] * temp[i*(nx+1)+j+1] );
                }
            }

            L = stencil.get_L_nw(1,ny-1,nx,ny);
                
            diagR[0] = L[C];
            ndiagR[0] = L[E];

            rhs[0] = resid[(ny-1)*(nx+1)+1];

            for(size_t j=2; j<nx-1; j++)
            {
                L = stencil.get_L_n(j,ny-1,nx,ny);
                diagR[j-1] = L[C];
                ndiagR[j-1] = L[E];
                ndiagL[j-2] = L[W];
                
                rhs[j-1] = resid[(ny-1)*(nx+1)+j];
            }

            L = stencil.get_L_ne(nx-1,ny-1,nx,ny);
                
            diagR[nx-2] = L[C];
            ndiagL[nx-3] = L[W];
            
            rhs[nx-2] = resid[(ny-1)*(nx+1)+nx-1];


            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<nx-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }

            // solve the linear system of equations R u = rhs
            temp[(ny-1)*(nx+1)+(nx-1)] = rhs[nx-2] / diagR[nx-2];
            
            for(size_t j=nx-2; j>0; j--)
            {
                temp[(ny-1)*(nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR[j-1] * temp[(ny-1)*(nx+1)+j+1] );
            }

            u += omega_ * temp;

            resid = residuum(u,f,stencil,nx,ny);

            temp = 0.0;

            //even lines
            for(size_t i=2; i<ny ; i+=2) 
            {                   
                L = stencil.get_L_w(1,i,nx,ny);
                
                diagR[0] = L[C];
                ndiagR[0] = L[E];

                rhs[0] = resid[i*(nx+1)+1];
                
                // L im Zentrum im Punkt (j/i)
                for(size_t j=2; j<nx-1; j++)  
                {
                    L = stencil.get_L_c(j,i,nx,ny);
                    diagR[j-1] = L[C];
                    ndiagR[j-1] = L[E];
                    ndiagL[j-2] = L[W];

                    rhs[j-1] = resid[i*(nx+1)+j];
                }
                
                L = stencil.get_L_e(nx-1,i,nx,ny);
                
                diagR[nx-2] = L[C];
                ndiagL[nx-3] = L[W];

                rhs[nx-2] = resid[i*(nx+1)+nx-1];

                // LR-decomposition + transformation of the rhs vector
                for(size_t k=1; k<nx-1; k++)  
                {
                    ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                    diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                    rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
                }

                // solve the linear system of equations R u = rhs
                temp[i*(nx+1)+nx-1] = rhs[nx-2] / diagR[nx-2];
                
                for(size_t j=nx-2; j>0; j--)
                {
                    temp[i*(nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR[j-1] * temp[i*(nx+1)+j+1] );
                }
            }

            u += omega_ * temp;
        }       
        else // if Nx and Ny are too small do one GS_lex step
        {
            gsRedBlack_.relax(u,f,stencil,nx,ny);
        }
    }                       
}
void ZebraLine::ninepointyzebra(
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
    
    if(stencil.isConstant() == true)
    {
        // get const operator L
        const NumericArray L = stencil.get_L_c(2,2,nx,ny);
        const PositionArray J_x = stencil.getJx(C);
        const PositionArray J_y = stencil.getJy(C);

        // for each line: correction of the rhs given by rhs = fv - [L[w]  0  L[e]] * u and elimination of the 
        // boundary condition in first and last inner point
        for(size_t i=1; i<nx ; i+=2) 
        {
            for(size_t j=0; j<ny-1; j++)  
            {
                rhs[j] = resid[(j+1)*(nx+1)+i];
            }

            // set tridiagonalmatrix for solving A u = rhs
            // A[i][i] = L[c]; A[i-1][i] = L[s]; A[i+1][i] = L[n]
            diagR = L[C];
            ndiagR = L[N];
            ndiagL = L[S];

            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<ny-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }
            
            // solve the linear system of equations R u = rhs
            temp[i+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];                

            for(size_t k=ny-2; k>0; k--)
            {
                temp[i+(nx+1)*k] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR[k-1] * temp[i+(nx+1)*(k+1)] );
            }
        }

        u += omega_ * temp;

        resid = residuum(u,f,stencil,nx,ny);

        temp = 0.0;


        // same for each even line
        for(size_t i=2; i<nx ; i+=2) 
        {
            for(size_t j=0; j<ny-1; j++)  
            {
                rhs[j] = resid[(j+1)*(nx+1)+i];
            }

            // set tridiagonalmatrix for solving A u = rhs
            // A[i][i] = L[c]; A[i-1][i] = L[s]; A[i+1][i] = L[n]
            diagR = L[C];
            ndiagR = L[N];
            ndiagL = L[S];

            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<ny-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }
            
            // solve the linear system of equations R u = rhs
            temp[i+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];                

            for(size_t k=ny-2; k>0; k--)
            {
                temp[i+(nx+1)*k] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR[k-1] * temp[i+(nx+1)*(k+1)] );
            }
        }

        u += omega_ * temp;
    }

    else
    {
        //Stencil ist not constant, so L needs to be evaluated in each grid point
        //no other change in the algorithm  
        NumericArray L = stencil.get_L_c(2,2,nx,ny);
        PositionArray J_x = stencil.getJx(C);
        PositionArray J_y = stencil.getJy(C);

        if(ny > 2)
        {
            L = stencil.get_L_sw(1,1,nx,ny);
                
            diagR[0] = L[C];
            ndiagR[0] = L[N];   
                              
            rhs[0] = resid[nx+1+1];
            
            for(size_t j=2; j<ny-1; j++) 
            {
                L = stencil.get_L_w(1,j,nx,ny);
                diagR[j-1] = L[C];
                ndiagR[j-1] = L[N];
                ndiagL[j-2] = L[S];

                rhs[j-1] = resid[j*(nx+1)+1];                   
            }
                
            L = stencil.get_L_nw(1,ny-1,nx,ny);
                
            diagR[ny-2] = L[C];
            ndiagL[ny-3] = L[S];
            
            rhs[ny-2] = resid[(ny-1)*(nx+1)+1];
            
            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<ny-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }

            // solve the linear system of equations R u = rhs
            temp[1+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];
            for(size_t k=ny-2; k>0; k--)
            {
                temp[1+(nx+1)*k] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR[k-1] * temp[1+(nx+1)*(k+1)] );
            }
            
            for(size_t i=3; i<nx-1 ; i+=2)
            {
                L = stencil.get_L_s(i,1,nx,ny);
                
                diagR[0] = L[C];
                ndiagR[0] = L[N];   
                
                rhs[0] = resid[nx+1+i];
                
                for(size_t j=2; j<ny-1; j++) 
                {
                    L = stencil.get_L_c(i,j,nx,ny);
                    diagR[j-1] = L[C];
                    ndiagR[j-1] = L[N];
                    ndiagL[j-2] = L[S];

                    rhs[j-1] = resid[j*(nx+1)+i];
                }
                
                L = stencil.get_L_n(i,ny-1,nx,ny);
                
                diagR[ny-2] = L[C];
                ndiagL[ny-3] = L[S];

                rhs[ny-2] = resid[(ny-1)*(nx+1)+i];
                
                // LR-decomposition + transformation of the rhs vector
                for(size_t k=1; k<ny-1; k++)  
                {
                    ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                    diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                    rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
                }
                
                // solve the linear system of equations R u = rhs
                temp[i+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];
                
                for(size_t k=ny-2; k>0; k--)
                {
                    temp[i+(nx+1)*k] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR[k-1] * temp[i+(nx+1)*(k+1)] );
                }
            }

            L = stencil.get_L_se(nx-1,1,nx,ny);
                
            diagR[0] = L[C];
            ndiagR[0] = L[N];   
                
            rhs[0] = resid[nx+1+nx-1];
                
            for(size_t j=2; j<ny-1; j++) 
            {
                L = stencil.get_L_e(nx-1,j,nx,ny);
                diagR[j-1] = L[C];
                ndiagR[j-1] = L[N];
                ndiagL[j-2] = L[S];

                rhs[j-1] = resid[j*(nx+1)+nx-1];                    
            }                   

            L = stencil.get_L_ne(nx-1,ny-1,nx,ny);
                
            diagR[ny-2] = L[C];
            ndiagL[ny-3] = L[S];

            rhs[ny-2] = resid[(ny-1)*(nx+1)+nx-1];

            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<ny-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }
                
            // solve the linear system of equations R u = rhs
            temp[nx-1+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];
                
            for(size_t k=ny-2; k>0; k--)
            {
                temp[nx-1+(nx+1)*k] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR[k-1] * temp[nx-1+(nx+1)*(k+1)] );
            }
            
            u += omega_ * temp;

            resid = residuum(u,f,stencil,nx,ny);

            temp = 0.0;

            for(size_t i=2; i<nx-1 ; i+=2)
            {
                L = stencil.get_L_s(i,1,nx,ny);
                
                diagR[0] = L[C];
                ndiagR[0] = L[N];   
                
                rhs[0] = resid[nx+1+i];
                
                for(size_t j=2; j<ny-1; j++) 
                {
                    L = stencil.get_L_c(i,j,nx,ny);
                    diagR[j-1] = L[C];
                    ndiagR[j-1] = L[N];
                    ndiagL[j-2] = L[S];

                    rhs[j-1] = resid[j*(nx+1)+i];
                }
                
                L = stencil.get_L_n(i,ny-1,nx,ny);
                
                diagR[ny-2] = L[C];
                ndiagL[ny-3] = L[S];

                rhs[ny-2] = resid[(ny-1)*(nx+1)+i];
                
                // LR-decomposition + transformation of the rhs vector
                for(size_t k=1; k<ny-1; k++)  
                {
                    ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                    diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                    rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
                }
                
                // solve the linear system of equations R u = rhs
                temp[i+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];
                
                for(size_t k=ny-2; k>0; k--)
                {
                    temp[i+(nx+1)*k] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR[k-1] * temp[i+(nx+1)*(k+1)] );
                }
            }

            u += omega_ * temp;
        }

        else // if Nx and Ny are too small do one GS_lex step
        {
            gsRedBlack_.relax(u,f,stencil,nx,ny);
        }
    }                       
}
void ZebraLine::xzebra(
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

        if(stencil.isConstant() == true)
        {
            // get const operator L
            const NumericArray L = stencil.get_L_c(2,2,nx,ny);
            const PositionArray J_x = stencil.getJx(C);
            const PositionArray J_y = stencil.getJy(C);
            
            NumericArray L_b = stencil.get_L_s(2,1,nx,ny);
            PositionArray J_b_x = stencil.getJx(S);
            PositionArray J_b_y = stencil.getJy(S);
                
            NumericArray L_c = stencil.get_L_sw(1,1,nx,ny);
            PositionArray J_c_x = stencil.getJx(SW);
            PositionArray J_c_y = stencil.getJy(SW);


            // setze rechte Seite für Zeile 1                   
            diagR[0] = L_c[C];
            ndiagR1[0] = L_c[E];
            ndiagR2[0] = L_c[NE];
                                
            rhs[0] = resid[nx+1+1];
            
            ndiagL1[0] = L_b[W];
            diagR[1] = L_b[C];
            ndiagR1[1] = L_b[E];                    
            ndiagR2[1] = L_b[SE];

            rhs[1] = resid[nx+1+2];

            for(size_t j=3; j<nx-2; j++)  
            {
                ndiagL2[j-3] = L_b[NW];
                ndiagL1[j-2] = L_b[W];
                diagR[j-1] = L_b[C];
                ndiagR1[j-1] = L_b[E];                  
                ndiagR2[j-1] = L_b[SE];
                    
                rhs[j-1] = resid[nx+1+j];
            }
                
            ndiagL2[nx-5] = L_b[NW];
            ndiagL1[nx-4] = L_b[W];
            diagR[nx-3] = L_b[C];
            ndiagR1[nx-3] = L_b[E];

            rhs[nx-3] = resid[nx+1+nx-2];
                
            L_c = stencil.get_L_se(nx-1,1,nx,ny);
            J_c_x = stencil.getJx(SE);
            J_c_y = stencil.getJy(SE);

            ndiagL2[nx-4] = L_c[NW];
            ndiagL1[nx-3] = L_c[W];
            diagR[nx-2] = L_c[C];

            rhs[nx-2] = resid[nx+1+nx-1];

            // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<nx-2; k++)  
            {
                ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1];
                
                ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];
            }

            ndiagL1[nx-2-1] = ndiagL1[nx-2-1]/diagR[nx-2-1];
            diagR[nx-2] -= ndiagL1[nx-2-1] * ndiagR1[nx-2-1];
            rhs[nx-2] = rhs[nx-2] - ndiagL1[nx-2-1] * rhs[nx-2-1];

            
            // solve the linear system of equations R u = rhs
            temp[nx+1+nx-1] = rhs[nx-2] / diagR[nx-2];
        
            for(size_t j=nx-2; j>1; j--)
            {
                temp[nx+1+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[nx+1+j+1] );
                
                rhs[j-2] -= ndiagR2[j-2] * temp[(nx+1)+j+1];
            }
            temp[nx+1+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[nx+1+1+1] );

   // durchlaufe ungerade innere Zeilen
            
            for(size_t i=3; i < ny-2; i+=2)
            {
                // setze rechte Seite                   
                L_b = stencil.get_L_w(1,i,nx,ny);
                J_b_x = stencil.getJx(W);
                J_b_y = stencil.getJy(W);

                diagR[0] = L_b[C];
                ndiagR1[0] = L_b[E];
                ndiagR2[0] = L_b[NE];
                                
                rhs[0] = resid[i*(nx+1)+1];
                                    
                ndiagL1[0] = L[W];
                diagR[1] = L[C];
                ndiagR1[1] = L[E];                  
                ndiagR2[1] = L[SE];

                rhs[1] = resid[i*(nx+1)+2];

                for(size_t j=3; j<nx-2; j++)  
                {
                    ndiagL2[j-3] = L[NW];
                    ndiagL1[j-2] = L[W];
                    diagR[j-1] = L[C];
                    ndiagR1[j-1] = L[E];                    
                    ndiagR2[j-1] = L[SE];
                    
                    rhs[j-1] = resid[i*(nx+1)+j];
                }
                
                ndiagL2[nx-5] = L[NW];
                ndiagL1[nx-4] = L[W];
                diagR[nx-3] = L[C];
                ndiagR1[nx-3] = L[E];

                rhs[nx-3] = resid[i*(nx+1)+nx-2];
                
                L_b = stencil.get_L_e(nx-1,i,nx,ny);
                J_b_x = stencil.getJx(E);
                J_b_y = stencil.getJy(E);

                ndiagL2[nx-4] = L_b[NW];
                ndiagL1[nx-3] = L_b[W];
                diagR[nx-2] = L_b[C];

                rhs[nx-2] = resid[i*(nx+1)+nx-1];
                
                // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<nx-2; k++)  
                {
                    ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                    rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1];
                    
                    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                    rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];                      
                }
                ndiagL1[nx-2-1] = ndiagL1[nx-2-1]/diagR[nx-2-1];
                diagR[nx-2] -= ndiagL1[nx-2-1] * ndiagR1[nx-2-1];
                rhs[nx-2] = rhs[nx-2] - ndiagL1[nx-3] * rhs[nx-3];
                
                // solve the linear system of equations R u = rhs
                temp[i*(nx+1)+nx-1] = rhs[nx-2] / diagR[nx-2];                                  
                for(size_t j=nx-2; j>1; j--)
                {
                    temp[i*(nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[i*(nx+1)+j+1] );                    
                    rhs[j-2] -= ndiagR2[j-2] * temp[i*(nx+1)+j+1];
                }
                temp[i*(nx+1)+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[i*(nx+1)+1+1] );
            }

    //relaxiere oberste Zeile
            
            // setze rechte Seite in oberster Zeile                 
            L_c = stencil.get_L_nw(1,ny-1,nx,ny);
            J_c_x = stencil.getJx(NW);
            J_c_y = stencil.getJy(NW);

            diagR[0] = L_c[C];
            ndiagR1[0] = L_c[E];
            ndiagR2[0] = L_c[NW];
                                
            rhs[0] = resid[(ny-1)*(nx+1)+1];
            
            L_b = stencil.get_L_n(2,ny-1,nx,ny);
            J_b_x = stencil.getJx(N);
            J_b_y = stencil.getJy(N);
            
            ndiagL1[0] = L_b[W];
            diagR[1] = L_b[C];
            ndiagR1[1] = L_b[E];                    
            ndiagR2[1] = L_b[NE];

            rhs[1] = resid[(ny-1)*(nx+1)+2];

            for(size_t j=3; j<nx-2; j++)  
            {
                ndiagL2[j-3] = L_b[NW];
                ndiagL1[j-2] = L_b[W];
                diagR[j-1] = L_b[C];
                ndiagR1[j-1] = L_b[E];                  
                ndiagR2[j-1] = L_b[NE];
                    
                rhs[j-1] = resid[(ny-1)*(nx+1)+j];
            }
                
            ndiagL2[nx-5] = L_b[NW];
            ndiagL1[nx-4] = L_b[W];
            diagR[nx-3] = L_b[C];
            ndiagR1[nx-3] = L_b[E];

            rhs[nx-3] = resid[(ny-1)*(nx+1)+nx-2];
                
            L_c = stencil.get_L_ne(nx-1,ny-1,nx,ny);
            J_c_x = stencil.getJx(NE);
            J_c_y = stencil.getJy(NE);

            ndiagL2[nx-4] = L_c[NW];
            ndiagL1[nx-3] = L_c[W];
            diagR[nx-2] = L_c[C];

            rhs[nx-2] = resid[(ny-1)*(nx+1)+nx-1];

            // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<nx-2; k++)  
            {
                ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1];
                    
                ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];                      
            }
            ndiagL1[nx-2-1] = ndiagL1[nx-2-1]/diagR[nx-2-1];
            diagR[nx-2] -= ndiagL1[nx-2-1] * ndiagR1[nx-2-1];
            rhs[nx-2] = rhs[nx-2] - ndiagL1[nx-3] * rhs[nx-3];

            // solve the linear system of equations R u = rhs
            temp[(ny-1)*(nx+1)+nx-1] = rhs[nx-2] / diagR[nx-2];
            for(size_t j=nx-2; j>1; j--)
            {
                temp[(ny-1)*(nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[(ny-1)*(nx+1)+j+1] );
                rhs[j-2] -= ndiagR2[j-2] * temp[(ny-1)*(nx+1)+j+1];
            }
            temp[(ny-1)*(nx+1)+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[(ny-1)*(nx+1)+1+1] );

            u += omega_ * temp;

            resid = residuum(u,f,stencil,nx,ny);

            temp = 0.0;

            // relaxiere gerade innere Zeilen

            for(size_t i=2; i < ny-1; i+=2)
            {
                // setze rechte Seite                   
                L_b = stencil.get_L_w(1,i,nx,ny);
                J_b_x = stencil.getJx(W);
                J_b_y = stencil.getJy(W);

                diagR[0] = L_b[C];
                ndiagR1[0] = L_b[E];
                ndiagR2[0] = L_b[NE];
                                
                rhs[0] = resid[i*(nx+1)+1];
                                    
                ndiagL1[0] = L[W];
                diagR[1] = L[C];
                ndiagR1[1] = L[E];                  
                ndiagR2[1] = L[SE];

                rhs[1] = resid[i*(nx+1)+2];

                for(size_t j=3; j<nx-2; j++)  
                {
                    ndiagL2[j-3] = L[NW];
                    ndiagL1[j-2] = L[W];
                    diagR[j-1] = L[C];
                    ndiagR1[j-1] = L[E];                    
                    ndiagR2[j-1] = L[SE];
                    
                    rhs[j-1] = resid[i*(nx+1)+j];
                }
                
                ndiagL2[nx-5] = L[NW];
                ndiagL1[nx-4] = L[W];
                diagR[nx-3] = L[C];
                ndiagR1[nx-3] = L[E];

                rhs[nx-3] = resid[i*(nx+1)+nx-2];
                
                L_b = stencil.get_L_e(nx-1,i,nx,ny);
                J_b_x = stencil.getJx(E);
                J_b_y = stencil.getJy(E);

                ndiagL2[nx-4] = L_b[NW];
                ndiagL1[nx-3] = L_b[W];
                diagR[nx-2] = L_b[C];

                rhs[nx-2] = resid[i*(nx+1)+nx-1];
                
                // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<nx-2; k++)  
                {
                    ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                    rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1];
                    
                    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                    rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];                      
                }
                ndiagL1[nx-2-1] = ndiagL1[nx-2-1]/diagR[nx-2-1];
                diagR[nx-2] -= ndiagL1[nx-2-1] * ndiagR1[nx-2-1];
                rhs[nx-2] = rhs[nx-2] - ndiagL1[nx-3] * rhs[nx-3];
                
                // solve the linear system of equations R u = rhs
                temp[i*(nx+1)+nx-1] = rhs[nx-2] / diagR[nx-2];                                  
                for(size_t j=nx-2; j>1; j--)
                {
                    temp[i*(nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[i*(nx+1)+j+1] );                    
                    rhs[j-2] -= ndiagR2[j-2] * temp[i*(nx+1)+j+1];
                }
                temp[i*(nx+1)+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[i*(nx+1)+1+1] );
            }

            u += omega_ * temp;
        }

                    
        else // stencil not constant
        {
            NumericArray L = stencil.get_L_c(2,2,nx,ny);
            PositionArray J_x = stencil.getJx(C);
            PositionArray J_y = stencil.getJy(C);
            
            NumericArray L_b = stencil.get_L_s(2,1,nx,ny);
            PositionArray J_b_x = stencil.getJx(S);
            PositionArray J_b_y = stencil.getJy(S);
                
            NumericArray L_c = stencil.get_L_sw(1,1,nx,ny);
            PositionArray J_c_x = stencil.getJx(SW);
            PositionArray J_c_y = stencil.getJy(SW);
            

            diagR[0] = L_c[C];
            ndiagR1[0] = L_c[E];
            ndiagR2[0] = L_c[NE];
                                
            rhs[0] = resid[nx+1+1];
            
            L_b = stencil.get_L_s(2,1,nx,ny);
            J_b_x = stencil.getJx(S);
            J_b_y = stencil.getJy(S);

            ndiagL1[0] = L_b[W];
            diagR[1] = L_b[C];
            ndiagR1[1] = L_b[E];                    
            ndiagR2[1] = L_b[SE];

            rhs[1] = resid[nx+1+2];

            for(size_t j=3; j<nx-2; j++)  
            {
                L_b = stencil.get_L_s(j,1,nx,ny);
                
                ndiagL2[j-3] = L_b[NW];
                ndiagL1[j-2] = L_b[W];
                diagR[j-1] = L_b[C];
                ndiagR1[j-1] = L_b[E];                  
                ndiagR2[j-1] = L_b[SE];
                    
                rhs[j-1] = resid[nx+1+j];
            }
            
            L_b = stencil.get_L_s(nx-2,1,nx,ny);

            ndiagL2[nx-5] = L_b[NW];
            ndiagL1[nx-4] = L_b[W];
            diagR[nx-3] = L_b[C];
            ndiagR1[nx-3] = L_b[E];

            rhs[nx-3] = resid[nx+1+nx-2];
                
            L_c = stencil.get_L_se(nx-1,1,nx,ny);
            J_c_x = stencil.getJx(SE);
            J_c_y = stencil.getJy(SE);

            ndiagL2[nx-4] = L_c[NW];
            ndiagL1[nx-3] = L_c[W];
            diagR[nx-2] = L_c[C];

            rhs[nx-2] = resid[nx+1+nx-1];

            // LR-decomposition + transformation of the rhs
            // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<nx-2; k++)  
            {
                ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1];
                
                ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];
            }

            ndiagL1[nx-2-1] = ndiagL1[nx-2-1]/diagR[nx-2-1];
            diagR[nx-2] -= ndiagL1[nx-2-1] * ndiagR1[nx-2-1];
            rhs[nx-2] = rhs[nx-2] - ndiagL1[nx-2-1] * rhs[nx-2-1];

            
            // solve the linear system of equations R u = rhs
            temp[nx+1+nx-1] = rhs[nx-2] / diagR[nx-2];
        
            for(size_t j=nx-2; j>1; j--)
            {
                temp[nx+1+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[nx+1+j+1] );
                
                rhs[j-2] -= ndiagR2[j-2] * temp[(nx+1)+j+1];
            }
            temp[nx+1+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[nx+1+1+1] );

   // durchlaufe ungerade innere Zeilen
            
            for(size_t i=3; i < ny-2; i+=2)
            {
                // setze rechte Seite                   
                L_b = stencil.get_L_w(1,i,nx,ny);
                J_b_x = stencil.getJx(W);
                J_b_y = stencil.getJy(W);

                diagR[0] = L_b[C];
                ndiagR1[0] = L_b[E];
                ndiagR2[0] = L_b[NE];
                                
                rhs[0] = resid[i*(nx+1)+1];
                
                L = stencil.get_L_c(2,i,nx,ny);

                ndiagL1[0] = L[W];
                diagR[1] = L[C];
                ndiagR1[1] = L[E];                  
                ndiagR2[1] = L[SE];

                rhs[1] = resid[i*(nx+1)+2];

                for(size_t j=3; j<nx-2; j++)  
                {
                    L = stencil.get_L_c(j,i,nx,ny);
                    
                    ndiagL2[j-3] = L[NW];
                    ndiagL1[j-2] = L[W];
                    diagR[j-1] = L[C];
                    ndiagR1[j-1] = L[E];                    
                    ndiagR2[j-1] = L[SE];
                    
                    rhs[j-1] = resid[i*(nx+1)+j];
                }
                
                L = stencil.get_L_c(nx-2,i,nx,ny);

                ndiagL2[nx-5] = L[NW];
                ndiagL1[nx-4] = L[W];
                diagR[nx-3] = L[C];
                ndiagR1[nx-3] = L[E];

                rhs[nx-3] = resid[i*(nx+1)+nx-2];
                
                L_b = stencil.get_L_e(nx-1,i,nx,ny);
                J_b_x = stencil.getJx(E);
                J_b_y = stencil.getJy(E);

                ndiagL2[nx-4] = L_b[NW];
                ndiagL1[nx-3] = L_b[W];
                diagR[nx-2] = L_b[C];

                rhs[nx-2] = resid[i*(nx+1)+nx-1];
                
                // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<nx-2; k++)  
                {
                    ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                    rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1];
                    
                    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                    rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];                      
                }
                ndiagL1[nx-2-1] = ndiagL1[nx-2-1]/diagR[nx-2-1];
                diagR[nx-2] -= ndiagL1[nx-2-1] * ndiagR1[nx-2-1];
                rhs[nx-2] = rhs[nx-2] - ndiagL1[nx-3] * rhs[nx-3];
                
                // solve the linear system of equations R u = rhs
                temp[i*(nx+1)+nx-1] = rhs[nx-2] / diagR[nx-2];                                  
                for(size_t j=nx-2; j>1; j--)
                {
                    temp[i*(nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[i*(nx+1)+j+1] );                    
                    rhs[j-2] -= ndiagR2[j-2] * temp[i*(nx+1)+j+1];
                }
                temp[i*(nx+1)+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[i*(nx+1)+1+1] );
            }

    //relaxiere oberste Zeile
            
            // setze rechte Seite in oberster Zeile                 
            L_c = stencil.get_L_nw(1,ny-1,nx,ny);
            J_c_x = stencil.getJx(NW);
            J_c_y = stencil.getJy(NW);

            diagR[0] = L_c[C];
            ndiagR1[0] = L_c[E];
            ndiagR2[0] = L_c[NW];
                                
            rhs[0] = resid[(ny-1)*(nx+1)+1];
            
            L_b = stencil.get_L_n(2,ny-1,nx,ny);

            J_b_x = stencil.getJx(N);
            J_b_y = stencil.getJy(N);
            
            ndiagL1[0] = L_b[W];
            diagR[1] = L_b[C];
            ndiagR1[1] = L_b[E];                    
            ndiagR2[1] = L_b[NE];

            rhs[1] = resid[(ny-1)*(nx+1)+2];

            for(size_t j=3; j<nx-2; j++)  
            {
                L_b = stencil.get_L_n(j,ny-1,nx,ny);
                
                ndiagL2[j-3] = L_b[NW];
                ndiagL1[j-2] = L_b[W];
                diagR[j-1] = L_b[C];
                ndiagR1[j-1] = L_b[E];                  
                ndiagR2[j-1] = L_b[NE];
                    
                rhs[j-1] = resid[(ny-1)*(nx+1)+j];
            }

            L_b = stencil.get_L_n(nx-2,ny-1,nx,ny);
                
            ndiagL2[nx-5] = L_b[NW];
            ndiagL1[nx-4] = L_b[W];
            diagR[nx-3] = L_b[C];
            ndiagR1[nx-3] = L_b[E];

            rhs[nx-3] = resid[(ny-1)*(nx+1)+nx-2];
                
            L_c = stencil.get_L_ne(nx-1,ny-1,nx,ny);
            J_c_x = stencil.getJx(NE);
            J_c_y = stencil.getJy(NE);

            ndiagL2[nx-4] = L_c[NW];
            ndiagL1[nx-3] = L_c[W];
            diagR[nx-2] = L_c[C];

            rhs[nx-2] = resid[(ny-1)*(nx+1)+nx-1];

            // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<nx-2; k++)  
            {
                ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1];
                    
                ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];                      
            }
            ndiagL1[nx-2-1] = ndiagL1[nx-2-1]/diagR[nx-2-1];
            diagR[nx-2] -= ndiagL1[nx-2-1] * ndiagR1[nx-2-1];
            rhs[nx-2] = rhs[nx-2] - ndiagL1[nx-3] * rhs[nx-3];

            // solve the linear system of equations R u = rhs
            temp[(ny-1)*(nx+1)+nx-1] = rhs[nx-2] / diagR[nx-2];
            for(size_t j=nx-2; j>1; j--)
            {
                temp[(ny-1)*(nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[(ny-1)*(nx+1)+j+1] );
                rhs[j-2] -= ndiagR2[j-2] * temp[(ny-1)*(nx+1)+j+1];
            }
            temp[(ny-1)*(nx+1)+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[(ny-1)*(nx+1)+1+1] );

            u += omega_ * temp;

            resid = residuum(u,f,stencil,nx,ny);

            temp = 0.0;

            // relaxiere gerade innere Zeilen

            for(size_t i=2; i < ny-1; i+=2)
            {
                // setze rechte Seite                   
                    L_b = stencil.get_L_w(1,i,nx,ny);
                    J_b_x = stencil.getJx(W);
                    J_b_y = stencil.getJy(W);

                    diagR[0] = L_b[C];
                    ndiagR1[0] = L_b[E];
                    ndiagR2[0] = L_b[NE];
                                    
                    rhs[0] = resid[i*(nx+1)+1];
                    
                    L = stencil.get_L_c(2,i,nx,ny);

                    ndiagL1[0] = L[W];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[E];                  
                    ndiagR2[1] = L[SE];

                    rhs[1] = resid[i*(nx+1)+2];

                    for(size_t j=3; j<nx-2; j++)  
                    { 
                        L = stencil.get_L_c(j,i,nx,ny);
 
                        ndiagL2[j-3] = L[NW];
                        ndiagL1[j-2] = L[W];
                        diagR[j-1] = L[C];
                        ndiagR1[j-1] = L[E];                    
                        ndiagR2[j-1] = L[SE];
                        
                        rhs[j-1] = resid[i*(nx+1)+j];
                    }
                    
                    L = stencil.get_L_c(nx-2,i,nx,ny);

                    ndiagL2[nx-5] = L[NW];
                    ndiagL1[nx-4] = L[W];
                    diagR[nx-3] = L[C];
                    ndiagR1[nx-3] = L[E];

                    rhs[nx-3] = resid[i*(nx+1)+nx-2];
                    
                    L_b = stencil.get_L_e(nx-1,i,nx,ny);
                    J_b_x = stencil.getJx(E);
                    J_b_y = stencil.getJy(E);

                    ndiagL2[nx-4] = L_b[NW];
                    ndiagL1[nx-3] = L_b[W];
                    diagR[nx-2] = L_b[C];

                    rhs[nx-2] = resid[i*(nx+1)+nx-1];
                    
                    // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<nx-2; k++)  
                {
                    ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
                    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                    rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1];
                    
                    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                    rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];                      
                }
                ndiagL1[nx-2-1] = ndiagL1[nx-2-1]/diagR[nx-2-1];
                diagR[nx-2] -= ndiagL1[nx-2-1] * ndiagR1[nx-2-1];
                rhs[nx-2] = rhs[nx-2] - ndiagL1[nx-3] * rhs[nx-3];
                
                // solve the linear system of equations R u = rhs
                temp[i*(nx+1)+nx-1] = rhs[nx-2] / diagR[nx-2];                                  
                for(size_t j=nx-2; j>1; j--)
                {
                    temp[i*(nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[i*(nx+1)+j+1] );                    
                    rhs[j-2] -= ndiagR2[j-2] * temp[i*(nx+1)+j+1];
                }
                temp[i*(nx+1)+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[i*(nx+1)+1+1] );
            }

            u += omega_ * temp;
        }

     }

     else 
     {
        for(int k=0; k<2; k++)
        {
            gsRedBlack_.relax(u,f,stencil,nx,ny);
        }
    }
}
void ZebraLine::yzebra(
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
                
        if(stencil.isConstant() == true)
        {
            // get const operator L
                const NumericArray L = stencil.get_L_c(2,2,nx,ny);
                const PositionArray J_x = stencil.getJx(C);
                const PositionArray J_y = stencil.getJy(C);
                
                NumericArray L_b = stencil.get_L_w(1,2,nx,ny);
                PositionArray J_b_x = stencil.getJx(W);
                PositionArray J_b_y = stencil.getJy(W);
 
                NumericArray L_c = stencil.get_L_sw(1,1,nx,ny);
                PositionArray J_c_x = stencil.getJx(SW);
                PositionArray J_c_y = stencil.getJy(SW);

                diagR[0] = L_c[C];
                ndiagR1[0] = L_c[N];
                ndiagR2[0] = L_c[NW];
                                
                rhs[0] = resid[nx+1+1];

                L_b = stencil.get_L_w(1,2,nx,ny);
                J_b_x = stencil.getJx(W);
                J_b_y = stencil.getJy(W);

                ndiagL1[0] = L_b[S];
                diagR[1] = L_b[C];
                ndiagR1[1] = L_b[N];                    
                ndiagR2[1] = L_b[NW];

                rhs[1] = resid[2*(nx+1)+1];

                for(size_t j=3; j<ny-2; j++)  
                {
                    ndiagL2[j-3] = L_b[SE];
                    ndiagL1[j-2] = L_b[S];
                    diagR[j-1] = L_b[C];
                    ndiagR1[j-1] = L_b[N];                  
                    ndiagR2[j-1] = L_b[NW];
                        
                    rhs[j-1] = resid[j*(nx+1)+1];
                }
                    
                ndiagL2[ny-5] = L_b[SE];
                ndiagL1[ny-4] = L_b[S];
                diagR[ny-3] = L_b[C];
                ndiagR1[ny-3] = L_b[N];

                rhs[ny-3] = resid[(ny-2)*(nx+1)+1];
                    
                L_c = stencil.get_L_nw(1,ny-1,nx,ny);
                J_c_x = stencil.getJx(NW);
                J_c_y = stencil.getJy(NW);

                ndiagL2[ny-4] = L_c[NE];
                ndiagL1[ny-3] = L_c[S];
                diagR[ny-2] = L_c[C];

                rhs[ny-2] = resid[(ny-1)*(nx+1)+1];

                // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<ny-2; k++)  
            {
                ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1]; 
                
                ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];
            }
            ndiagL1[ny-2-1] = ndiagL1[ny-2-1]/diagR[ny-2-1];
            diagR[ny-2] -= ndiagL1[ny-2-1] * ndiagR1[ny-2-1];
            rhs[ny-2] = rhs[ny-2] - ndiagL1[ny-2-1] * rhs[ny-3];

            // solve the linear system of equations R u = rhs
            temp[1+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];
            for(size_t k=ny-2; k>1; k--)
            {
                temp[1+k*(nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[1+(k+1)*(nx+1)] );

                rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(nx+1)+1];
            }
            temp[nx+1+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(nx+1)+1] );

////////////////////////////////////////////////////////////////////////////////////////////////

            // durchlaufe alle inneren Spalten, Spaltenindex i              
            for(size_t i=3; i < nx-2; i+=2)
            {
                // setze rechte Seite                   
                    L_b = stencil.get_L_s(i,1,nx,ny);
                    J_b_x = stencil.getJx(S);
                    J_b_y = stencil.getJy(S);

                    diagR[0] = L_b[C];
                    ndiagR1[0] = L_b[N];
                    ndiagR2[0] = L_b[NE];
       
                    rhs[0] = resid[nx+1+i];

                    ndiagL1[0] = L[S];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[N];                  
                    ndiagR2[1] = L[NE];

                    rhs[1] = resid[2*(nx+1)+i];                     

                    for(size_t j=3; j<ny-2; j++)
                    {
                       ndiagL2[j-3] = L[SW];
                       ndiagL1[j-2] = L[S];
                       diagR[j-1] = L[C];
                       ndiagR1[j-1] = L[N];
                       ndiagR2[j-1] = L[NE];
 
                       rhs[j-1] = resid[j*(nx+1)+i];                       
                    }
                                        
                    ndiagL2[nx-5] = L[SW];
                    ndiagL1[nx-4] = L[S];
                    diagR[nx-3] = L[C];
                    ndiagR1[nx-3] = L[N];

                    rhs[ny-3] = resid[(ny-2)*(nx+1)+i];

                    L_b = stencil.get_L_n(i,ny-1,nx,ny);
                    J_b_x = stencil.getJx(N);
                    J_b_y = stencil.getJy(N);

                    ndiagL2[nx-4] = L_b[SW];
                    ndiagL1[nx-3] = L_b[S];
                    diagR[nx-2] = L_b[C];

                    rhs[ny-2] = resid[(ny-1)*(nx+1)+i];
                    
                    // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<ny-2; k++)  
                {
                    ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];  
                    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                    rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1]; 
                    
                    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                    rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];
                }
                ndiagL1[ny-2-1] = ndiagL1[ny-2-1]/diagR[ny-2-1];
                diagR[ny-2] -= ndiagL1[ny-2-1] * ndiagR1[ny-2-1];
                rhs[ny-2] = rhs[ny-2] - ndiagL1[ny-2-1] * rhs[ny-3];

                // solve the linear system of equations R u = rhs
                temp[i+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];                    
                for(size_t k=ny-2; k>1; k--)
                {
                    temp[i+k*(nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[i+(k+1)*(nx+1)] );
                    rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(nx+1)+i];                        
                }
                temp[nx+1+i] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(nx+1)+i] );
            }

////////////////letzte Spalte//////////////////
            
            L_c = stencil.get_L_se(nx-1,1,nx,ny);
            J_c_x = stencil.getJx(SE);
            J_c_y = stencil.getJy(SE);

            rhs[0] = resid[nx+1+nx-1];

            L_b = stencil.get_L_e(nx-1,2,nx,ny);
            J_b_x = stencil.getJx(E);
            J_b_y = stencil.getJy(E);

            ndiagL1[0] = L_b[S];
            diagR[1] = L_b[C];
            ndiagR1[1] = L_b[N];                    
            ndiagR2[1] = L_b[NE];

            rhs[1] = resid[2*(nx+1)+nx-1];

            for(size_t j=3; j<ny-2; j++)  
            {
                ndiagL2[j-3] = L_b[SE];
                ndiagL1[j-2] = L_b[S];
                diagR[j-1] = L_b[C];
                ndiagR1[j-1] = L_b[N];                  
                ndiagR2[j-1] = L_b[NE];
                    
                rhs[j-1] = resid[j*(nx+1)+nx-1];
            }
                
            ndiagL2[ny-5] = L_b[SE];
            ndiagL1[ny-4] = L_b[S];
            diagR[ny-3] = L_b[C];
            ndiagR1[ny-3] = L_b[N];

            rhs[ny-3] = resid[(ny-2)*(nx+1)+nx-1];
                
            L_c = stencil.get_L_ne(nx-1,ny-1,nx,ny);
            J_c_x = stencil.getJx(NE);
            J_c_y = stencil.getJy(NE);

            ndiagL2[ny-4] = L_c[NE];
            ndiagL1[ny-3] = L_c[S];
            diagR[ny-2] = L_c[C];

            rhs[ny-2] = resid[(ny-1)*(nx+1)+nx-1];

            // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<ny-2; k++)  
            {
                ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1]; 
                    
                ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];
            }
            ndiagL1[ny-2-1] = ndiagL1[ny-2-1]/diagR[ny-2-1];
            diagR[ny-2] -= ndiagL1[ny-2-1] * ndiagR1[ny-2-1];
            rhs[ny-2] = rhs[ny-2] - ndiagL1[ny-2-1] * rhs[ny-3];

            // solve the linear system of equations R u = rhs
            temp[nx-1+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];
            for(size_t k=ny-2; k>1; k--)
            {
                temp[nx-1+k*(nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[nx-1+(k+1)*(nx+1)] );
                rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(nx+1)+nx-1];
            }
            temp[nx+1+nx-1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(nx+1)+nx-1] );
            
            u += omega_ * temp;

            resid = residuum(u,f,stencil,nx,ny);

            temp = 0.0;

            // durchlaufe alle inneren Spalten, Spaltenindex i              
            for(size_t i=2; i < nx-1; i+=2)
            {
                // setze rechte Seite                   
                    L_b = stencil.get_L_s(i,1,nx,ny);
                    J_b_x = stencil.getJx(S);
                    J_b_y = stencil.getJy(S);

                    diagR[0] = L_b[C];
                    ndiagR1[0] = L_b[N];
                    ndiagR2[0] = L_b[NE];
       
                    rhs[0] = resid[nx+1+i];

                    ndiagL1[0] = L[S];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[N];                  
                    ndiagR2[1] = L[NE];

                    rhs[1] = resid[2*(nx+1)+i];                     

                    for(size_t j=3; j<ny-2; j++)
                    {
                       ndiagL2[j-3] = L[SW];
                       ndiagL1[j-2] = L[S];
                       diagR[j-1] = L[C];
                       ndiagR1[j-1] = L[N];
                       ndiagR2[j-1] = L[NE];
 
                       rhs[j-1] = resid[j*(nx+1)+i];                       
                    }
                                        
                    ndiagL2[nx-5] = L[SW];
                    ndiagL1[nx-4] = L[S];
                    diagR[nx-3] = L[C];
                    ndiagR1[nx-3] = L[N];

                    rhs[ny-3] = resid[(ny-2)*(nx+1)+i];

                    L_b = stencil.get_L_n(i,ny-1,nx,ny);
                    J_b_x = stencil.getJx(N);
                    J_b_y = stencil.getJy(N);

                    ndiagL2[nx-4] = L_b[SW];
                    ndiagL1[nx-3] = L_b[S];
                    diagR[nx-2] = L_b[C];

                    rhs[ny-2] = resid[(ny-1)*(nx+1)+i];
                    
                    // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<ny-2; k++)  
                {
                    ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];  
                    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                    rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1]; 
                    
                    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                    rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];
                }
                ndiagL1[ny-2-1] = ndiagL1[ny-2-1]/diagR[ny-2-1];
                diagR[ny-2] -= ndiagL1[ny-2-1] * ndiagR1[ny-2-1];
                rhs[ny-2] = rhs[ny-2] - ndiagL1[ny-2-1] * rhs[ny-3];

                // solve the linear system of equations R u = rhs
                temp[i+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];                    
                for(size_t k=ny-2; k>1; k--)
                {
                    temp[i+k*(nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[i+(k+1)*(nx+1)] );
                    rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(nx+1)+i];                        
                }
                temp[nx+1+i] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(nx+1)+i] );
            }

            u += omega_ * temp;

        }
        else // stencil not constant
            {
                NumericArray L = stencil.get_L_c(2,2,nx,ny);
                PositionArray J_x = stencil.getJx(C);
                PositionArray J_y = stencil.getJy(C);
                
                NumericArray L_b = stencil.get_L_w(1,2,nx,ny);
                PositionArray J_b_x = stencil.getJx(W);
                PositionArray J_b_y = stencil.getJy(W);
 
                NumericArray L_c = stencil.get_L_sw(1,1,nx,ny);
                PositionArray J_c_x = stencil.getJx(SW);
                PositionArray J_c_y = stencil.getJy(SW);

                
                diagR[0] = L_c[C];
                ndiagR1[0] = L_c[N];
                ndiagR2[0] = L_c[NW];
                                
                rhs[0] = resid[nx+1+1];
                
                ndiagL1[0] = L_b[S];
                diagR[1] = L_b[C];
                ndiagR1[1] = L_b[N];                    
                ndiagR2[1] = L_b[NW];

                rhs[1] = resid[2*(nx+1)+1];

                for(size_t j=3; j<ny-2; j++)  
                {
                    L_b = stencil.get_L_w(1,j,nx,ny);

                    ndiagL2[j-3] = L_b[SE];
                    ndiagL1[j-2] = L_b[S];
                    diagR[j-1] = L_b[C];
                    ndiagR1[j-1] = L_b[N];                  
                    ndiagR2[j-1] = L_b[NW];
                        
                    rhs[j-1] = resid[j*(nx+1)+1];
                }

                L_b = stencil.get_L_w(1,ny-2,nx,ny);
                
                ndiagL2[ny-5] = L_b[SE];
                ndiagL1[ny-4] = L_b[S];
                diagR[ny-3] = L_b[C];
                ndiagR1[ny-3] = L_b[N];

                rhs[ny-3] = resid[(ny-2)*(nx+1)+1];
                    
                L_c = stencil.get_L_nw(1,ny-1,nx,ny);
                J_c_x = stencil.getJx(NW);
                J_c_y = stencil.getJy(NW);

                ndiagL2[ny-4] = L_c[NE];
                ndiagL1[ny-3] = L_c[S];
                diagR[ny-2] = L_c[C];

                rhs[ny-2] = resid[(ny-1)*(nx+1)+1];

                for(size_t k=1; k<ny-2; k++)  
                {
                    ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];  
                    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                    rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1]; 
                    
                    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                    rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];
                }
                ndiagL1[ny-2-1] = ndiagL1[ny-2-1]/diagR[ny-2-1];
                diagR[ny-2] -= ndiagL1[ny-2-1] * ndiagR1[ny-2-1];
                rhs[ny-2] = rhs[ny-2] - ndiagL1[ny-2-1] * rhs[ny-3];

                // solve the linear system of equations R u = rhs
            temp[1+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];
            for(size_t k=ny-2; k>1; k--)
            {
                temp[1+k*(nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[1+(k+1)*(nx+1)] );

                rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(nx+1)+1];
            }
            temp[nx+1+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(nx+1)+1] );

////////////////////////////////////////////////////////////////////////////////////////////////

            // durchlaufe alle inneren Spalten, Spaltenindex i              
            for(size_t i=3; i < nx-2; i+=2)
            {
                // setze rechte Seite                   
                    L_b = stencil.get_L_s(i,1,nx,ny);
                    J_b_x = stencil.getJx(S);
                    J_b_y = stencil.getJy(S);

                    diagR[0] = L_b[C];
                    ndiagR1[0] = L_b[N];
                    ndiagR2[0] = L_b[NE];
       
                    rhs[0] = resid[nx+1+i];

                    L = stencil.get_L_c(i,2,nx,ny);

                    ndiagL1[0] = L[S];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[N];                  
                    ndiagR2[1] = L[NE];

                    rhs[1] = resid[2*(nx+1)+i];  

                    for(size_t j=3; j<ny-2; j++)
                    {
                       L = stencil.get_L_c(i,j,nx,ny);
                       ndiagL2[j-3] = L[SW];
                       ndiagL1[j-2] = L[S];
                       diagR[j-1] = L[C];
                       ndiagR1[j-1] = L[N];
                       ndiagR2[j-1] = L[NE];
 
                       rhs[j-1] = resid[j*(nx+1)+i];
                    }
                                        
                    L = stencil.get_L_c(i,ny-2,nx,ny);
                    ndiagL2[nx-5] = L[SW];
                    ndiagL1[nx-4] = L[S];
                    diagR[nx-3] = L[C];
                    ndiagR1[nx-3] = L[N];

                    rhs[ny-3] = resid[(ny-2)*(nx+1)+i];

                    L_b = stencil.get_L_n(i,ny-1,nx,ny);
                    J_b_x = stencil.getJx(N);
                    J_b_y = stencil.getJy(N);

                    ndiagL2[nx-4] = L_b[SW];
                    ndiagL1[nx-3] = L_b[S];
                    diagR[nx-2] = L_b[C];

                    rhs[ny-2] = resid[(ny-1)*(nx+1)+i];
                    
                    // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<ny-2; k++)  
                {
                    ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];  
                    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                    rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1]; 
                    
                    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                    rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];
                }
                ndiagL1[ny-2-1] = ndiagL1[ny-2-1]/diagR[ny-2-1];
                diagR[ny-2] -= ndiagL1[ny-2-1] * ndiagR1[ny-2-1];
                rhs[ny-2] = rhs[ny-2] - ndiagL1[ny-2-1] * rhs[ny-3];

                // solve the linear system of equations R u = rhs
                temp[i+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];                    
                for(size_t k=ny-2; k>1; k--)
                {
                    temp[i+k*(nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[i+(k+1)*(nx+1)] );
                    rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(nx+1)+i];                        
                }
                temp[nx+1+i] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(nx+1)+i] );
            }

////////////////letzte Spalte//////////////////
            
            L_c = stencil.get_L_se(nx-1,1,nx,ny);
            J_c_x = stencil.getJx(SE);
            J_c_y = stencil.getJy(SE);

            rhs[0] = resid[nx+1+nx-1];

            L_b = stencil.get_L_e(nx-1,2,nx,ny);
            J_b_x = stencil.getJx(E);
            J_b_y = stencil.getJy(E);

            ndiagL1[0] = L_b[S];
            diagR[1] = L_b[C];
            ndiagR1[1] = L_b[N];                    
            ndiagR2[1] = L_b[NE];

            rhs[1] = resid[2*(nx+1)+nx-1];

            for(size_t j=3; j<ny-2; j++)  
            {
                L_b = stencil.get_L_e(nx-1,j,nx,ny);
                
                ndiagL2[j-3] = L_b[SE];
                ndiagL1[j-2] = L_b[S];
                diagR[j-1] = L_b[C];
                ndiagR1[j-1] = L_b[N];                  
                ndiagR2[j-1] = L_b[NE];
                    
                rhs[j-1] = resid[j*(nx+1)+nx-1];
            }
                
            L_b = stencil.get_L_e(nx-1,ny-2,nx,ny);

            ndiagL2[ny-5] = L_b[SE];
            ndiagL1[ny-4] = L_b[S];
            diagR[ny-3] = L_b[C];
            ndiagR1[ny-3] = L_b[N];

            rhs[ny-3] = resid[(ny-2)*(nx+1)+nx-1];
                
            L_c = stencil.get_L_ne(nx-1,ny-1,nx,ny);
            J_c_x = stencil.getJx(NE);
            J_c_y = stencil.getJy(NE);

            ndiagL2[ny-4] = L_c[NE];
            ndiagL1[ny-3] = L_c[S];
            diagR[ny-2] = L_c[C];

            rhs[ny-2] = resid[(ny-1)*(nx+1)+nx-1];

            // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<ny-2; k++)  
            {
                ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1]; 
                    
                ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];
            }
            ndiagL1[ny-2-1] = ndiagL1[ny-2-1]/diagR[ny-2-1];
            diagR[ny-2] -= ndiagL1[ny-2-1] * ndiagR1[ny-2-1];
            rhs[ny-2] = rhs[ny-2] - ndiagL1[ny-2-1] * rhs[ny-3];

            // solve the linear system of equations R u = rhs
            temp[nx-1+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];
            for(size_t k=ny-2; k>1; k--)
            {
                temp[nx-1+k*(nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[nx-1+(k+1)*(nx+1)] );
                rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(nx+1)+nx-1];
            }
            temp[nx+1+nx-1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(nx+1)+nx-1] );
            
            u += omega_ * temp;

            resid = residuum(u,f,stencil,nx,ny);

            temp = 0.0;

            // durchlaufe alle inneren Spalten, Spaltenindex i              
            for(size_t i=2; i < nx-1; i+=2)
            {
                // setze rechte Seite                   
                    L_b = stencil.get_L_s(i,1,nx,ny);
                    J_b_x = stencil.getJx(S);
                    J_b_y = stencil.getJy(S);

                    diagR[0] = L_b[C];
                    ndiagR1[0] = L_b[N];
                    ndiagR2[0] = L_b[NE];
       
                    rhs[0] = resid[nx+1+i];

                    L = stencil.get_L_c(i,2,nx,ny);

                    ndiagL1[0] = L[S];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[N];                  
                    ndiagR2[1] = L[NE];

                    rhs[1] = resid[2*(nx+1)+i];  

                    for(size_t j=3; j<ny-2; j++)
                    {
                       L = stencil.get_L_c(i,j,nx,ny);
                       ndiagL2[j-3] = L[SW];
                       ndiagL1[j-2] = L[S];
                       diagR[j-1] = L[C];
                       ndiagR1[j-1] = L[N];
                       ndiagR2[j-1] = L[NE];
 
                       rhs[j-1] = resid[j*(nx+1)+i];
                    }
                                        
                    L = stencil.get_L_c(i,ny-2,nx,ny);
                    ndiagL2[nx-5] = L[SW];
                    ndiagL1[nx-4] = L[S];
                    diagR[nx-3] = L[C];
                    ndiagR1[nx-3] = L[N];

                    rhs[ny-3] = resid[(ny-2)*(nx+1)+i];

                    L_b = stencil.get_L_n(i,ny-1,nx,ny);
                    J_b_x = stencil.getJx(N);
                    J_b_y = stencil.getJy(N);

                    ndiagL2[nx-4] = L_b[SW];
                    ndiagL1[nx-3] = L_b[S];
                    diagR[nx-2] = L_b[C];

                    rhs[ny-2] = resid[(ny-1)*(nx+1)+i];
                    
                    // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<ny-2; k++)  
                {
                    ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];  
                    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
                    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
                    rhs[k] = rhs[k] - ndiagL1[k-1] * rhs[k-1]; 
                    
                    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
                    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
                    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
                    rhs[k+1] = rhs[k+1] - ndiagL2[k-1] * rhs[k-1];
                }
                ndiagL1[ny-2-1] = ndiagL1[ny-2-1]/diagR[ny-2-1];
                diagR[ny-2] -= ndiagL1[ny-2-1] * ndiagR1[ny-2-1];
                rhs[ny-2] = rhs[ny-2] - ndiagL1[ny-2-1] * rhs[ny-3];

                // solve the linear system of equations R u = rhs
                temp[i+(nx+1)*(ny-1)] = rhs[ny-2] / diagR[ny-2];                    
                for(size_t k=ny-2; k>1; k--)
                {
                    temp[i+k*(nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[i+(k+1)*(nx+1)] );
                    rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(nx+1)+i];                        
                }
                temp[nx+1+i] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(nx+1)+i] );
            }

            u += omega_ * temp;
        }
    }
    else //parameter zu klein
    {
        for(int k=0; k<2; k++)
        {
            gsRedBlack_.relax(u,f,stencil,nx,ny);
        }
    }
}
}
