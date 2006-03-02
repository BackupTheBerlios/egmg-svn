/** \file ZebraLine.cpp
 * \author Andre Oeckerath
 * \brief ZebraLine.cpp contains the implementation of the class ZebraLine.
 * \see ZebraLine.h
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
void ZebraLine::ninepointxzebra(NumericArray &u, const NumericArray &fv, 
                        NumericArray resid, const Stencil &stencil, const size_t Nx, 
                        const size_t Ny) const
                   
{ 
    NumericArray rhs(0.0,Nx+1);
    NumericArray temp(0.0,(Nx+1)*(Ny+1));
    
    //valarrays needed for saving the tridiagonal matrix A of linear system A u = rhs
    NumericArray diagR(Nx-1);
    NumericArray ndiagR(Nx-2);
    NumericArray ndiagL(Nx-2);
    
    if(stencil.isConstant() == true)
    {
        // get const operator L
        const NumericArray L = stencil.get_L_c(2,2,Nx,Ny);
        const PositionArray J_x = stencil.getJx(C);
        const PositionArray J_y = stencil.getJy(C);
        
        // for each line: correction of the rhs given by rhs = fv - [L[n]  0  L[s]]^t * u and elimination of the 
        // boundary condition in first and last inner point
        for(size_t i=1; i<Ny ; i+=2) 
        {
            for(size_t j=0; j<Nx-1; j++)  
            {
                rhs[j] = resid[i*(Nx+1)+j+1];
            }   
                                            
            // set tridiagonalmatrix for solving A u = rhs
            // A[i][i] = L[c]; A[i-1][i] = L[w]; A[i+1][i] = L[e]
            diagR = L[C];
            ndiagR = L[E];
            ndiagL = L[W];

            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<Nx-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }

            // solve the linear system of equations R u = rhs
            temp[i*(Nx+1)+(Nx-1)] = rhs[Nx-2] / diagR[Nx-2];
            
            for(size_t j=Nx-2; j>0; j--)
            {
                temp[i*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR[j-1] * temp[i*(Nx+1)+j+1] );
            }
        }

        u += omega_ * temp;

        resid = residuum(u,fv,stencil,Nx,Ny);

        temp = 0.0;

        // same for even lines 
        for(size_t i=2; i<Ny ; i+=2)
        {
            for(size_t j=0; j<Nx-1; j++)  
            {
                rhs[j] = resid[i*(Nx+1)+j+1];
            }   

            diagR = L[C];
            ndiagR = L[E];
            ndiagL = L[W];

            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<Nx-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }

            // solve the linear system of equations R u = rhs
            temp[i*(Nx+1)+(Nx-1)] = rhs[Nx-2] / diagR[Nx-2];
            
            for(size_t j=Nx-2; j>0; j--)
            {
                temp[i*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR[j-1] * temp[i*(Nx+1)+j+1] );
            }
        }

        u += omega_ * temp;


    }

    else
    {
        //Stencil ist not constant, so L needs to be evaluated in each grid point
        //no other change in the algorithm  
        NumericArray L = stencil.get_L_c(2,2,Nx,Ny);
        PositionArray J_x = stencil.getJx(C);
        PositionArray J_y = stencil.getJy(C);

        if(Nx > 2)
        {
            L = stencil.get_L_sw(1,1,Nx,Ny);
                
            diagR[0] = L[C];
            ndiagR[0] = L[E];

            rhs[0] = resid[Nx+1+1];

            for(size_t j=2; j<Nx-1; j++)
            {
                L = stencil.get_L_s(j,1,Nx,Ny);
                diagR[j-1] = L[C];
                ndiagR[j-1] = L[E];
                ndiagL[j-2] = L[W];
                 
                rhs[j-1] = resid[Nx+1+j];
            }

            L = stencil.get_L_se(Nx-1,1,Nx,Ny);
            
            diagR[Nx-2] = L[C];
            ndiagL[Nx-3] = L[W];
                
            rhs[Nx-2] = resid[Nx+1+Nx-1];

            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<Nx-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }

            // solve the linear system of equations R u = rhs
            temp[Nx+1+Nx-1] = rhs[Nx-2] / diagR[Nx-2];
            
            for(size_t j=Nx-2; j>0; j--)
            {
                temp[Nx+1+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR[j-1] * temp[Nx+1+j+1] );
            }               
            

            for(size_t i=3; i<Ny-1 ; i+=2)
            {
                L = stencil.get_L_w(1,i,Nx,Ny);
                
                diagR[0] = L[C];
                ndiagR[0] = L[E];

                rhs[0] = resid[i*(Nx+1)+1];

                for(size_t j=2; j<Nx-1; j++)
                {
                    L = stencil.get_L_c(j,i,Nx,Ny);
                    diagR[j-1] = L[C];
                    ndiagR[j-1] = L[E];
                    ndiagL[j-2] = L[W];
                    
                    rhs[j-1] = resid[i*(Nx+1)+j];
                }

                L = stencil.get_L_e(Nx-1,i,Nx,Ny);
                
                diagR[Nx-2] = L[C];
                ndiagL[Nx-3] = L[W];
                
                rhs[Nx-2] = resid[i*(Nx+1)+Nx-1];


                // LR-decomposition + transformation of the rhs vector
                for(size_t k=1; k<Nx-1; k++)  
                {
                    ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                    diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                    rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
                }

                // solve the linear system of equations R u = rhs
                temp[i*(Nx+1)+Nx-1] = rhs[Nx-2] / diagR[Nx-2];
                
                for(size_t j=Nx-2; j>0; j--)
                {
                    temp[i*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR[j-1] * temp[i*(Nx+1)+j+1] );
                }
            }

            L = stencil.get_L_nw(1,Ny-1,Nx,Ny);
                
            diagR[0] = L[C];
            ndiagR[0] = L[E];

            rhs[0] = resid[(Ny-1)*(Nx+1)+1];

            for(size_t j=2; j<Nx-1; j++)
            {
                L = stencil.get_L_n(j,Ny-1,Nx,Ny);
                diagR[j-1] = L[C];
                ndiagR[j-1] = L[E];
                ndiagL[j-2] = L[W];
                
                rhs[j-1] = resid[(Ny-1)*(Nx+1)+j];
            }

            L = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
                
            diagR[Nx-2] = L[C];
            ndiagL[Nx-3] = L[W];
            
            rhs[Nx-2] = resid[(Ny-1)*(Nx+1)+Nx-1];


            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<Nx-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }

            // solve the linear system of equations R u = rhs
            temp[(Ny-1)*(Nx+1)+(Nx-1)] = rhs[Nx-2] / diagR[Nx-2];
            
            for(size_t j=Nx-2; j>0; j--)
            {
                temp[(Ny-1)*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR[j-1] * temp[(Ny-1)*(Nx+1)+j+1] );
            }

            u += omega_ * temp;

            resid = residuum(u,fv,stencil,Nx,Ny);

            temp = 0.0;

            //even lines
            for(size_t i=2; i<Ny ; i+=2) 
            {                   
                L = stencil.get_L_w(1,i,Nx,Ny);
                
                diagR[0] = L[C];
                ndiagR[0] = L[E];

                rhs[0] = resid[i*(Nx+1)+1];
                
                // L im Zentrum im Punkt (j/i)
                for(size_t j=2; j<Nx-1; j++)  
                {
                    L = stencil.get_L_c(j,i,Nx,Ny);
                    diagR[j-1] = L[C];
                    ndiagR[j-1] = L[E];
                    ndiagL[j-2] = L[W];

                    rhs[j-1] = resid[i*(Nx+1)+j];
                }
                
                L = stencil.get_L_e(Nx-1,i,Nx,Ny);
                
                diagR[Nx-2] = L[C];
                ndiagL[Nx-3] = L[W];

                rhs[Nx-2] = resid[i*(Nx+1)+Nx-1];

                // LR-decomposition + transformation of the rhs vector
                for(size_t k=1; k<Nx-1; k++)  
                {
                    ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                    diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                    rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
                }

                // solve the linear system of equations R u = rhs
                temp[i*(Nx+1)+Nx-1] = rhs[Nx-2] / diagR[Nx-2];
                
                for(size_t j=Nx-2; j>0; j--)
                {
                    temp[i*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR[j-1] * temp[i*(Nx+1)+j+1] );
                }
            }

            u += omega_ * temp;
        }       
        
        else // if Nx and Ny are too small do one GS_lex step
        {
            Precision temp=0;

            for(size_t k=1; k<Ny; k++)
            {
                L = stencil.get_L_c(1,k,Nx,Ny);

                for(size_t sum=5; sum < L.size(); sum++)
                {
                    temp -= L[sum] * u[1+J_x[sum]+(k+J_y[sum])*(Nx+1)];
                }
                u[1+k*(Nx+1)] = 1/L[C] * ( fv[1+k*(Nx+1)] - L[W] * u[1+J_x[W]+k*(Nx+1)] - L[E] * u[1+J_x[E]+k*(Nx+1)] 
                              - L[N] * u[1+(k+J_y[N])*(Nx+1)] - L[S] * u[1+(k+J_y[S])*(Nx+1)] - temp );
            }
        }
    }                       
}
void ZebraLine::ninepointyzebra(NumericArray &u, const NumericArray &fv, 
                        NumericArray resid, const Stencil &stencil, const size_t Nx, 
                        const size_t Ny) const
                   
{ 
    NumericArray rhs(0.0,Ny+1);
    NumericArray temp(0.0,(Nx+1)*(Ny+1));

    //valarrays needed for saving the tridiagonal matrix A of linear system A u = rhs
    NumericArray diagR(Ny-1);
    NumericArray ndiagR(Ny-2);
    NumericArray ndiagL(Ny-2);
    
    if(stencil.isConstant() == true)
    {
        // get const operator L
        const NumericArray L = stencil.get_L_c(2,2,Nx,Ny);
        const PositionArray J_x = stencil.getJx(C);
        const PositionArray J_y = stencil.getJy(C);

        // for each line: correction of the rhs given by rhs = fv - [L[w]  0  L[e]] * u and elimination of the 
        // boundary condition in first and last inner point
        for(size_t i=1; i<Nx ; i+=2) 
        {
            for(size_t j=0; j<Ny-1; j++)  
            {
                rhs[j] = resid[(j+1)*(Nx+1)+i];
            }

            // set tridiagonalmatrix for solving A u = rhs
            // A[i][i] = L[c]; A[i-1][i] = L[s]; A[i+1][i] = L[n]
            diagR = L[C];
            ndiagR = L[N];
            ndiagL = L[S];

            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<Ny-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }
            
            // solve the linear system of equations R u = rhs
            temp[i+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];                

            for(size_t k=Ny-2; k>0; k--)
            {
                temp[i+(Nx+1)*k] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR[k-1] * temp[i+(Nx+1)*(k+1)] );
            }
        }

        u += omega_ * temp;

        resid = residuum(u,fv,stencil,Nx,Ny);

        temp = 0.0;


        // same for each even line
        for(size_t i=2; i<Nx ; i+=2) 
        {
            for(size_t j=0; j<Ny-1; j++)  
            {
                rhs[j] = resid[(j+1)*(Nx+1)+i];
            }

            // set tridiagonalmatrix for solving A u = rhs
            // A[i][i] = L[c]; A[i-1][i] = L[s]; A[i+1][i] = L[n]
            diagR = L[C];
            ndiagR = L[N];
            ndiagL = L[S];

            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<Ny-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }
            
            // solve the linear system of equations R u = rhs
            temp[i+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];                

            for(size_t k=Ny-2; k>0; k--)
            {
                temp[i+(Nx+1)*k] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR[k-1] * temp[i+(Nx+1)*(k+1)] );
            }
        }

        u += omega_ * temp;
    }

    else
    {
        //Stencil ist not constant, so L needs to be evaluated in each grid point
        //no other change in the algorithm  
        NumericArray L = stencil.get_L_c(2,2,Nx,Ny);
        PositionArray J_x = stencil.getJx(C);
        PositionArray J_y = stencil.getJy(C);

        if(Ny > 2)
        {
            L = stencil.get_L_sw(1,1,Nx,Ny);
                
            diagR[0] = L[C];
            ndiagR[0] = L[N];   
                              
            rhs[0] = resid[Nx+1+1];
            
            for(size_t j=2; j<Ny-1; j++) 
            {
                L = stencil.get_L_w(1,j,Nx,Ny);
                diagR[j-1] = L[C];
                ndiagR[j-1] = L[N];
                ndiagL[j-2] = L[S];

                rhs[j-1] = resid[j*(Nx+1)+1];                   
            }
                
            L = stencil.get_L_nw(1,Ny-1,Nx,Ny);
                
            diagR[Ny-2] = L[C];
            ndiagL[Ny-3] = L[S];
            
            rhs[Ny-2] = resid[(Ny-1)*(Nx+1)+1];
            
            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<Ny-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }

            // solve the linear system of equations R u = rhs
            temp[1+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];
            for(size_t k=Ny-2; k>0; k--)
            {
                temp[1+(Nx+1)*k] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR[k-1] * temp[1+(Nx+1)*(k+1)] );
            }
            
            for(size_t i=3; i<Nx-1 ; i+=2)
            {
                L = stencil.get_L_s(i,1,Nx,Ny);
                
                diagR[0] = L[C];
                ndiagR[0] = L[N];   
                
                rhs[0] = resid[Nx+1+i];
                
                for(size_t j=2; j<Ny-1; j++) 
                {
                    L = stencil.get_L_c(i,j,Nx,Ny);
                    diagR[j-1] = L[C];
                    ndiagR[j-1] = L[N];
                    ndiagL[j-2] = L[S];

                    rhs[j-1] = resid[j*(Nx+1)+i];
                }
                
                L = stencil.get_L_n(i,Ny-1,Nx,Ny);
                
                diagR[Ny-2] = L[C];
                ndiagL[Ny-3] = L[S];

                rhs[Ny-2] = resid[(Ny-1)*(Nx+1)+i];
                
                // LR-decomposition + transformation of the rhs vector
                for(size_t k=1; k<Ny-1; k++)  
                {
                    ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                    diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                    rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
                }
                
                // solve the linear system of equations R u = rhs
                temp[i+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];
                
                for(size_t k=Ny-2; k>0; k--)
                {
                    temp[i+(Nx+1)*k] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR[k-1] * temp[i+(Nx+1)*(k+1)] );
                }
            }

            L = stencil.get_L_se(Nx-1,1,Nx,Ny);
                
            diagR[0] = L[C];
            ndiagR[0] = L[N];   
                
            rhs[0] = resid[Nx+1+Nx-1];
                
            for(size_t j=2; j<Ny-1; j++) 
            {
                L = stencil.get_L_e(Nx-1,j,Nx,Ny);
                diagR[j-1] = L[C];
                ndiagR[j-1] = L[N];
                ndiagL[j-2] = L[S];

                rhs[j-1] = resid[j*(Nx+1)+Nx-1];                    
            }                   

            L = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
                
            diagR[Ny-2] = L[C];
            ndiagL[Ny-3] = L[S];

            rhs[Ny-2] = resid[(Ny-1)*(Nx+1)+Nx-1];

            // LR-decomposition + transformation of the rhs vector
            for(size_t k=1; k<Ny-1; k++)  
            {
                ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
            }
                
            // solve the linear system of equations R u = rhs
            temp[Nx-1+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];
                
            for(size_t k=Ny-2; k>0; k--)
            {
                temp[Nx-1+(Nx+1)*k] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR[k-1] * temp[Nx-1+(Nx+1)*(k+1)] );
            }
            
            u += omega_ * temp;

            resid = residuum(u,fv,stencil,Nx,Ny);

            temp = 0.0;

            for(size_t i=2; i<Nx-1 ; i+=2)
            {
                L = stencil.get_L_s(i,1,Nx,Ny);
                
                diagR[0] = L[C];
                ndiagR[0] = L[N];   
                
                rhs[0] = resid[Nx+1+i];
                
                for(size_t j=2; j<Ny-1; j++) 
                {
                    L = stencil.get_L_c(i,j,Nx,Ny);
                    diagR[j-1] = L[C];
                    ndiagR[j-1] = L[N];
                    ndiagL[j-2] = L[S];

                    rhs[j-1] = resid[j*(Nx+1)+i];
                }
                
                L = stencil.get_L_n(i,Ny-1,Nx,Ny);
                
                diagR[Ny-2] = L[C];
                ndiagL[Ny-3] = L[S];

                rhs[Ny-2] = resid[(Ny-1)*(Nx+1)+i];
                
                // LR-decomposition + transformation of the rhs vector
                for(size_t k=1; k<Ny-1; k++)  
                {
                    ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
                    diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
                    rhs[k] = rhs[k] - ndiagL[k-1] * rhs[k-1];  
                }
                
                // solve the linear system of equations R u = rhs
                temp[i+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];
                
                for(size_t k=Ny-2; k>0; k--)
                {
                    temp[i+(Nx+1)*k] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR[k-1] * temp[i+(Nx+1)*(k+1)] );
                }
            }

            u += omega_ * temp;
        }

        else // if Nx and Ny are too small do one GS_lex step
        {
            for(size_t k=1; k<Nx; k++)
            {
                Precision temp=0;

                L = stencil.get_L_c(1,k,Nx,Ny);

                for(size_t sum=5; sum < L.size(); sum++)
                {
                    temp -= L[sum] * u[1+J_x[sum]+(k+J_y[sum])*(Nx+1)];
                }
                u[1+k*(Nx+1)] = 1/L[C] * ( fv[1+k*(Nx+1)] - L[W] * u[1+J_x[W]+k*(Nx+1)] - L[E] * u[1+J_x[E]+k*(Nx+1)] 
                              - L[N] * u[1+(k+J_y[N])*(Nx+1)] - L[S] * u[1+(k+J_y[S])*(Nx+1)] - temp );
            }
        }
    }                       
}
void ZebraLine::xzebra(NumericArray &u, const NumericArray &fv, 
                        NumericArray resid, const Stencil &stencil, const size_t Nx, 
                        const size_t Ny) const

{
    if((Ny > 4) && (Nx > 4))
    {   
        NumericArray rhs(0.0,Nx-1);
        NumericArray temp(0.0,(Nx+1)*(Ny+1));

        NumericArray diagR(0.0,Nx-1);
        NumericArray ndiagR1(0.0,Nx-2);
        NumericArray ndiagL1(0.0,Nx-2);
        NumericArray ndiagR2(0.0,Nx-3);
        NumericArray ndiagL2(0.0,Nx-3);

        if(stencil.isConstant() == true)
        {
            // get const operator L
            const NumericArray L = stencil.get_L_c(2,2,Nx,Ny);
            const PositionArray J_x = stencil.getJx(C);
            const PositionArray J_y = stencil.getJy(C);
            
            NumericArray L_b = stencil.get_L_s(2,1,Nx,Ny);
            PositionArray J_b_x = stencil.getJx(S);
            PositionArray J_b_y = stencil.getJy(S);
                
            NumericArray L_c = stencil.get_L_sw(1,1,Nx,Ny);
            PositionArray J_c_x = stencil.getJx(SW);
            PositionArray J_c_y = stencil.getJy(SW);


            // setze rechte Seite für Zeile 1                   
            diagR[0] = L_c[C];
            ndiagR1[0] = L_c[E];
            ndiagR2[0] = L_c[NE];
                                
            rhs[0] = resid[Nx+1+1];
            
            ndiagL1[0] = L_b[W];
            diagR[1] = L_b[C];
            ndiagR1[1] = L_b[E];                    
            ndiagR2[1] = L_b[SE];

            rhs[1] = resid[Nx+1+2];

            for(size_t j=3; j<Nx-2; j++)  
            {
                ndiagL2[j-3] = L_b[NW];
                ndiagL1[j-2] = L_b[W];
                diagR[j-1] = L_b[C];
                ndiagR1[j-1] = L_b[E];                  
                ndiagR2[j-1] = L_b[SE];
                    
                rhs[j-1] = resid[Nx+1+j];
            }
                
            ndiagL2[Nx-5] = L_b[NW];
            ndiagL1[Nx-4] = L_b[W];
            diagR[Nx-3] = L_b[C];
            ndiagR1[Nx-3] = L_b[E];

            rhs[Nx-3] = resid[Nx+1+Nx-2];
                
            L_c = stencil.get_L_se(Nx-1,1,Nx,Ny);
            J_c_x = stencil.getJx(SE);
            J_c_y = stencil.getJy(SE);

            ndiagL2[Nx-4] = L_c[NW];
            ndiagL1[Nx-3] = L_c[W];
            diagR[Nx-2] = L_c[C];

            rhs[Nx-2] = resid[Nx+1+Nx-1];

            // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<Nx-2; k++)  
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

            ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
            diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
            rhs[Nx-2] = rhs[Nx-2] - ndiagL1[Nx-2-1] * rhs[Nx-2-1];

            
            // solve the linear system of equations R u = rhs
            temp[Nx+1+Nx-1] = rhs[Nx-2] / diagR[Nx-2];
        
            for(size_t j=Nx-2; j>1; j--)
            {
                temp[Nx+1+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[Nx+1+j+1] );
                
                rhs[j-2] -= ndiagR2[j-2] * temp[(Nx+1)+j+1];
            }
            temp[Nx+1+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[Nx+1+1+1] );

   // durchlaufe ungerade innere Zeilen
            
            for(size_t i=3; i < Ny-2; i+=2)
            {
                // setze rechte Seite                   
                L_b = stencil.get_L_w(1,i,Nx,Ny);
                J_b_x = stencil.getJx(W);
                J_b_y = stencil.getJy(W);

                diagR[0] = L_b[C];
                ndiagR1[0] = L_b[E];
                ndiagR2[0] = L_b[NE];
                                
                rhs[0] = resid[i*(Nx+1)+1];
                                    
                ndiagL1[0] = L[W];
                diagR[1] = L[C];
                ndiagR1[1] = L[E];                  
                ndiagR2[1] = L[SE];

                rhs[1] = resid[i*(Nx+1)+2];

                for(size_t j=3; j<Nx-2; j++)  
                {
                    ndiagL2[j-3] = L[NW];
                    ndiagL1[j-2] = L[W];
                    diagR[j-1] = L[C];
                    ndiagR1[j-1] = L[E];                    
                    ndiagR2[j-1] = L[SE];
                    
                    rhs[j-1] = resid[i*(Nx+1)+j];
                }
                
                ndiagL2[Nx-5] = L[NW];
                ndiagL1[Nx-4] = L[W];
                diagR[Nx-3] = L[C];
                ndiagR1[Nx-3] = L[E];

                rhs[Nx-3] = resid[i*(Nx+1)+Nx-2];
                
                L_b = stencil.get_L_e(Nx-1,i,Nx,Ny);
                J_b_x = stencil.getJx(E);
                J_b_y = stencil.getJy(E);

                ndiagL2[Nx-4] = L_b[NW];
                ndiagL1[Nx-3] = L_b[W];
                diagR[Nx-2] = L_b[C];

                rhs[Nx-2] = resid[i*(Nx+1)+Nx-1];
                
                // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<Nx-2; k++)  
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
                ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
                diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
                rhs[Nx-2] = rhs[Nx-2] - ndiagL1[Nx-3] * rhs[Nx-3];
                
                // solve the linear system of equations R u = rhs
                temp[i*(Nx+1)+Nx-1] = rhs[Nx-2] / diagR[Nx-2];                                  
                for(size_t j=Nx-2; j>1; j--)
                {
                    temp[i*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[i*(Nx+1)+j+1] );                    
                    rhs[j-2] -= ndiagR2[j-2] * temp[i*(Nx+1)+j+1];
                }
                temp[i*(Nx+1)+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[i*(Nx+1)+1+1] );
            }

    //relaxiere oberste Zeile
            
            // setze rechte Seite in oberster Zeile                 
            L_c = stencil.get_L_nw(1,Ny-1,Nx,Ny);
            J_c_x = stencil.getJx(NW);
            J_c_y = stencil.getJy(NW);

            diagR[0] = L_c[C];
            ndiagR1[0] = L_c[E];
            ndiagR2[0] = L_c[NW];
                                
            rhs[0] = resid[(Ny-1)*(Nx+1)+1];
            
            L_b = stencil.get_L_n(2,Ny-1,Nx,Ny);
            J_b_x = stencil.getJx(N);
            J_b_y = stencil.getJy(N);
            
            ndiagL1[0] = L_b[W];
            diagR[1] = L_b[C];
            ndiagR1[1] = L_b[E];                    
            ndiagR2[1] = L_b[NE];

            rhs[1] = resid[(Ny-1)*(Nx+1)+2];

            for(size_t j=3; j<Nx-2; j++)  
            {
                ndiagL2[j-3] = L_b[NW];
                ndiagL1[j-2] = L_b[W];
                diagR[j-1] = L_b[C];
                ndiagR1[j-1] = L_b[E];                  
                ndiagR2[j-1] = L_b[NE];
                    
                rhs[j-1] = resid[(Ny-1)*(Nx+1)+j];
            }
                
            ndiagL2[Nx-5] = L_b[NW];
            ndiagL1[Nx-4] = L_b[W];
            diagR[Nx-3] = L_b[C];
            ndiagR1[Nx-3] = L_b[E];

            rhs[Nx-3] = resid[(Ny-1)*(Nx+1)+Nx-2];
                
            L_c = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
            J_c_x = stencil.getJx(NE);
            J_c_y = stencil.getJy(NE);

            ndiagL2[Nx-4] = L_c[NW];
            ndiagL1[Nx-3] = L_c[W];
            diagR[Nx-2] = L_c[C];

            rhs[Nx-2] = resid[(Ny-1)*(Nx+1)+Nx-1];

            // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<Nx-2; k++)  
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
            ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
            diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
            rhs[Nx-2] = rhs[Nx-2] - ndiagL1[Nx-3] * rhs[Nx-3];

            // solve the linear system of equations R u = rhs
            temp[(Ny-1)*(Nx+1)+Nx-1] = rhs[Nx-2] / diagR[Nx-2];
            for(size_t j=Nx-2; j>1; j--)
            {
                temp[(Ny-1)*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[(Ny-1)*(Nx+1)+j+1] );
                rhs[j-2] -= ndiagR2[j-2] * temp[(Ny-1)*(Nx+1)+j+1];
            }
            temp[(Ny-1)*(Nx+1)+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[(Ny-1)*(Nx+1)+1+1] );

            u += omega_ * temp;

            resid = residuum(u,fv,stencil,Nx,Ny);

            temp = 0.0;

            // relaxiere gerade innere Zeilen

            for(size_t i=2; i < Ny-1; i+=2)
            {
                // setze rechte Seite                   
                L_b = stencil.get_L_w(1,i,Nx,Ny);
                J_b_x = stencil.getJx(W);
                J_b_y = stencil.getJy(W);

                diagR[0] = L_b[C];
                ndiagR1[0] = L_b[E];
                ndiagR2[0] = L_b[NE];
                                
                rhs[0] = resid[i*(Nx+1)+1];
                                    
                ndiagL1[0] = L[W];
                diagR[1] = L[C];
                ndiagR1[1] = L[E];                  
                ndiagR2[1] = L[SE];

                rhs[1] = resid[i*(Nx+1)+2];

                for(size_t j=3; j<Nx-2; j++)  
                {
                    ndiagL2[j-3] = L[NW];
                    ndiagL1[j-2] = L[W];
                    diagR[j-1] = L[C];
                    ndiagR1[j-1] = L[E];                    
                    ndiagR2[j-1] = L[SE];
                    
                    rhs[j-1] = resid[i*(Nx+1)+j];
                }
                
                ndiagL2[Nx-5] = L[NW];
                ndiagL1[Nx-4] = L[W];
                diagR[Nx-3] = L[C];
                ndiagR1[Nx-3] = L[E];

                rhs[Nx-3] = resid[i*(Nx+1)+Nx-2];
                
                L_b = stencil.get_L_e(Nx-1,i,Nx,Ny);
                J_b_x = stencil.getJx(E);
                J_b_y = stencil.getJy(E);

                ndiagL2[Nx-4] = L_b[NW];
                ndiagL1[Nx-3] = L_b[W];
                diagR[Nx-2] = L_b[C];

                rhs[Nx-2] = resid[i*(Nx+1)+Nx-1];
                
                // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<Nx-2; k++)  
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
                ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
                diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
                rhs[Nx-2] = rhs[Nx-2] - ndiagL1[Nx-3] * rhs[Nx-3];
                
                // solve the linear system of equations R u = rhs
                temp[i*(Nx+1)+Nx-1] = rhs[Nx-2] / diagR[Nx-2];                                  
                for(size_t j=Nx-2; j>1; j--)
                {
                    temp[i*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[i*(Nx+1)+j+1] );                    
                    rhs[j-2] -= ndiagR2[j-2] * temp[i*(Nx+1)+j+1];
                }
                temp[i*(Nx+1)+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[i*(Nx+1)+1+1] );
            }

            u += omega_ * temp;
        }

                    
        else // stencil not constant
        {
            NumericArray L = stencil.get_L_c(2,2,Nx,Ny);
            PositionArray J_x = stencil.getJx(C);
            PositionArray J_y = stencil.getJy(C);
            
            NumericArray L_b = stencil.get_L_s(2,1,Nx,Ny);
            PositionArray J_b_x = stencil.getJx(S);
            PositionArray J_b_y = stencil.getJy(S);
                
            NumericArray L_c = stencil.get_L_sw(1,1,Nx,Ny);
            PositionArray J_c_x = stencil.getJx(SW);
            PositionArray J_c_y = stencil.getJy(SW);
            

            diagR[0] = L_c[C];
            ndiagR1[0] = L_c[E];
            ndiagR2[0] = L_c[NE];
                                
            rhs[0] = resid[Nx+1+1];
            
            L_b = stencil.get_L_s(2,1,Nx,Ny);
            J_b_x = stencil.getJx(S);
            J_b_y = stencil.getJy(S);

            ndiagL1[0] = L_b[W];
            diagR[1] = L_b[C];
            ndiagR1[1] = L_b[E];                    
            ndiagR2[1] = L_b[SE];

            rhs[1] = resid[Nx+1+2];

            for(size_t j=3; j<Nx-2; j++)  
            {
                L_b = stencil.get_L_s(j,1,Nx,Ny);
                
                ndiagL2[j-3] = L_b[NW];
                ndiagL1[j-2] = L_b[W];
                diagR[j-1] = L_b[C];
                ndiagR1[j-1] = L_b[E];                  
                ndiagR2[j-1] = L_b[SE];
                    
                rhs[j-1] = resid[Nx+1+j];
            }
            
            L_b = stencil.get_L_s(Nx-2,1,Nx,Ny);

            ndiagL2[Nx-5] = L_b[NW];
            ndiagL1[Nx-4] = L_b[W];
            diagR[Nx-3] = L_b[C];
            ndiagR1[Nx-3] = L_b[E];

            rhs[Nx-3] = resid[Nx+1+Nx-2];
                
            L_c = stencil.get_L_se(Nx-1,1,Nx,Ny);
            J_c_x = stencil.getJx(SE);
            J_c_y = stencil.getJy(SE);

            ndiagL2[Nx-4] = L_c[NW];
            ndiagL1[Nx-3] = L_c[W];
            diagR[Nx-2] = L_c[C];

            rhs[Nx-2] = resid[Nx+1+Nx-1];

            // LR-decomposition + transformation of the rhs
            // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<Nx-2; k++)  
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

            ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
            diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
            rhs[Nx-2] = rhs[Nx-2] - ndiagL1[Nx-2-1] * rhs[Nx-2-1];

            
            // solve the linear system of equations R u = rhs
            temp[Nx+1+Nx-1] = rhs[Nx-2] / diagR[Nx-2];
        
            for(size_t j=Nx-2; j>1; j--)
            {
                temp[Nx+1+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[Nx+1+j+1] );
                
                rhs[j-2] -= ndiagR2[j-2] * temp[(Nx+1)+j+1];
            }
            temp[Nx+1+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[Nx+1+1+1] );

   // durchlaufe ungerade innere Zeilen
            
            for(size_t i=3; i < Ny-2; i+=2)
            {
                // setze rechte Seite                   
                L_b = stencil.get_L_w(1,i,Nx,Ny);
                J_b_x = stencil.getJx(W);
                J_b_y = stencil.getJy(W);

                diagR[0] = L_b[C];
                ndiagR1[0] = L_b[E];
                ndiagR2[0] = L_b[NE];
                                
                rhs[0] = resid[i*(Nx+1)+1];
                
                L = stencil.get_L_c(2,i,Nx,Ny);

                ndiagL1[0] = L[W];
                diagR[1] = L[C];
                ndiagR1[1] = L[E];                  
                ndiagR2[1] = L[SE];

                rhs[1] = resid[i*(Nx+1)+2];

                for(size_t j=3; j<Nx-2; j++)  
                {
                    L = stencil.get_L_c(j,i,Nx,Ny);
                    
                    ndiagL2[j-3] = L[NW];
                    ndiagL1[j-2] = L[W];
                    diagR[j-1] = L[C];
                    ndiagR1[j-1] = L[E];                    
                    ndiagR2[j-1] = L[SE];
                    
                    rhs[j-1] = resid[i*(Nx+1)+j];
                }
                
                L = stencil.get_L_c(Nx-2,i,Nx,Ny);

                ndiagL2[Nx-5] = L[NW];
                ndiagL1[Nx-4] = L[W];
                diagR[Nx-3] = L[C];
                ndiagR1[Nx-3] = L[E];

                rhs[Nx-3] = resid[i*(Nx+1)+Nx-2];
                
                L_b = stencil.get_L_e(Nx-1,i,Nx,Ny);
                J_b_x = stencil.getJx(E);
                J_b_y = stencil.getJy(E);

                ndiagL2[Nx-4] = L_b[NW];
                ndiagL1[Nx-3] = L_b[W];
                diagR[Nx-2] = L_b[C];

                rhs[Nx-2] = resid[i*(Nx+1)+Nx-1];
                
                // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<Nx-2; k++)  
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
                ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
                diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
                rhs[Nx-2] = rhs[Nx-2] - ndiagL1[Nx-3] * rhs[Nx-3];
                
                // solve the linear system of equations R u = rhs
                temp[i*(Nx+1)+Nx-1] = rhs[Nx-2] / diagR[Nx-2];                                  
                for(size_t j=Nx-2; j>1; j--)
                {
                    temp[i*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[i*(Nx+1)+j+1] );                    
                    rhs[j-2] -= ndiagR2[j-2] * temp[i*(Nx+1)+j+1];
                }
                temp[i*(Nx+1)+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[i*(Nx+1)+1+1] );
            }

    //relaxiere oberste Zeile
            
            // setze rechte Seite in oberster Zeile                 
            L_c = stencil.get_L_nw(1,Ny-1,Nx,Ny);
            J_c_x = stencil.getJx(NW);
            J_c_y = stencil.getJy(NW);

            diagR[0] = L_c[C];
            ndiagR1[0] = L_c[E];
            ndiagR2[0] = L_c[NW];
                                
            rhs[0] = resid[(Ny-1)*(Nx+1)+1];
            
            L_b = stencil.get_L_n(2,Ny-1,Nx,Ny);

            J_b_x = stencil.getJx(N);
            J_b_y = stencil.getJy(N);
            
            ndiagL1[0] = L_b[W];
            diagR[1] = L_b[C];
            ndiagR1[1] = L_b[E];                    
            ndiagR2[1] = L_b[NE];

            rhs[1] = resid[(Ny-1)*(Nx+1)+2];

            for(size_t j=3; j<Nx-2; j++)  
            {
                L_b = stencil.get_L_n(j,Ny-1,Nx,Ny);
                
                ndiagL2[j-3] = L_b[NW];
                ndiagL1[j-2] = L_b[W];
                diagR[j-1] = L_b[C];
                ndiagR1[j-1] = L_b[E];                  
                ndiagR2[j-1] = L_b[NE];
                    
                rhs[j-1] = resid[(Ny-1)*(Nx+1)+j];
            }

            L_b = stencil.get_L_n(Nx-2,Ny-1,Nx,Ny);
                
            ndiagL2[Nx-5] = L_b[NW];
            ndiagL1[Nx-4] = L_b[W];
            diagR[Nx-3] = L_b[C];
            ndiagR1[Nx-3] = L_b[E];

            rhs[Nx-3] = resid[(Ny-1)*(Nx+1)+Nx-2];
                
            L_c = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
            J_c_x = stencil.getJx(NE);
            J_c_y = stencil.getJy(NE);

            ndiagL2[Nx-4] = L_c[NW];
            ndiagL1[Nx-3] = L_c[W];
            diagR[Nx-2] = L_c[C];

            rhs[Nx-2] = resid[(Ny-1)*(Nx+1)+Nx-1];

            // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<Nx-2; k++)  
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
            ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
            diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
            rhs[Nx-2] = rhs[Nx-2] - ndiagL1[Nx-3] * rhs[Nx-3];

            // solve the linear system of equations R u = rhs
            temp[(Ny-1)*(Nx+1)+Nx-1] = rhs[Nx-2] / diagR[Nx-2];
            for(size_t j=Nx-2; j>1; j--)
            {
                temp[(Ny-1)*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[(Ny-1)*(Nx+1)+j+1] );
                rhs[j-2] -= ndiagR2[j-2] * temp[(Ny-1)*(Nx+1)+j+1];
            }
            temp[(Ny-1)*(Nx+1)+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[(Ny-1)*(Nx+1)+1+1] );

            u += omega_ * temp;

            resid = residuum(u,fv,stencil,Nx,Ny);

            temp = 0.0;

            // relaxiere gerade innere Zeilen

            for(size_t i=2; i < Ny-1; i+=2)
            {
                // setze rechte Seite                   
                    L_b = stencil.get_L_w(1,i,Nx,Ny);
                    J_b_x = stencil.getJx(W);
                    J_b_y = stencil.getJy(W);

                    diagR[0] = L_b[C];
                    ndiagR1[0] = L_b[E];
                    ndiagR2[0] = L_b[NE];
                                    
                    rhs[0] = resid[i*(Nx+1)+1];
                    
                    L = stencil.get_L_c(2,i,Nx,Ny);

                    ndiagL1[0] = L[W];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[E];                  
                    ndiagR2[1] = L[SE];

                    rhs[1] = resid[i*(Nx+1)+2];

                    for(size_t j=3; j<Nx-2; j++)  
                    { 
                        L = stencil.get_L_c(j,i,Nx,Ny);
 
                        ndiagL2[j-3] = L[NW];
                        ndiagL1[j-2] = L[W];
                        diagR[j-1] = L[C];
                        ndiagR1[j-1] = L[E];                    
                        ndiagR2[j-1] = L[SE];
                        
                        rhs[j-1] = resid[i*(Nx+1)+j];
                    }
                    
                    L = stencil.get_L_c(Nx-2,i,Nx,Ny);

                    ndiagL2[Nx-5] = L[NW];
                    ndiagL1[Nx-4] = L[W];
                    diagR[Nx-3] = L[C];
                    ndiagR1[Nx-3] = L[E];

                    rhs[Nx-3] = resid[i*(Nx+1)+Nx-2];
                    
                    L_b = stencil.get_L_e(Nx-1,i,Nx,Ny);
                    J_b_x = stencil.getJx(E);
                    J_b_y = stencil.getJy(E);

                    ndiagL2[Nx-4] = L_b[NW];
                    ndiagL1[Nx-3] = L_b[W];
                    diagR[Nx-2] = L_b[C];

                    rhs[Nx-2] = resid[i*(Nx+1)+Nx-1];
                    
                    // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<Nx-2; k++)  
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
                ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
                diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
                rhs[Nx-2] = rhs[Nx-2] - ndiagL1[Nx-3] * rhs[Nx-3];
                
                // solve the linear system of equations R u = rhs
                temp[i*(Nx+1)+Nx-1] = rhs[Nx-2] / diagR[Nx-2];                                  
                for(size_t j=Nx-2; j>1; j--)
                {
                    temp[i*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[j-1] - ndiagR1[j-1] * temp[i*(Nx+1)+j+1] );                    
                    rhs[j-2] -= ndiagR2[j-2] * temp[i*(Nx+1)+j+1];
                }
                temp[i*(Nx+1)+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[i*(Nx+1)+1+1] );
            }

            u += omega_ * temp;
        }

     }

     else 
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
void ZebraLine::yzebra(NumericArray &u, const NumericArray &fv, 
                        NumericArray resid, const Stencil &stencil, const size_t Nx, 
                        const size_t Ny) const

{
    if((Ny > 4) && (Nx > 4))
    {

        NumericArray rhs(0.0,Ny-1);
        NumericArray temp(0.0,(Nx+1)*(Ny+1));


        NumericArray diagR(0.0,Ny-1);
        NumericArray ndiagR1(0.0,Ny-2);
        NumericArray ndiagL1(0.0,Ny-2);
        NumericArray ndiagR2(0.0,Ny-3);
        NumericArray ndiagL2(0.0,Ny-3);
                
        if(stencil.isConstant() == true)
        {
            // get const operator L
                const NumericArray L = stencil.get_L_c(2,2,Nx,Ny);
                const PositionArray J_x = stencil.getJx(C);
                const PositionArray J_y = stencil.getJy(C);
                
                NumericArray L_b = stencil.get_L_w(1,2,Nx,Ny);
                PositionArray J_b_x = stencil.getJx(W);
                PositionArray J_b_y = stencil.getJy(W);
 
                NumericArray L_c = stencil.get_L_sw(1,1,Nx,Ny);
                PositionArray J_c_x = stencil.getJx(SW);
                PositionArray J_c_y = stencil.getJy(SW);

                diagR[0] = L_c[C];
                ndiagR1[0] = L_c[N];
                ndiagR2[0] = L_c[NW];
                                
                rhs[0] = resid[Nx+1+1];

                L_b = stencil.get_L_w(1,2,Nx,Ny);
                J_b_x = stencil.getJx(W);
                J_b_y = stencil.getJy(W);

                ndiagL1[0] = L_b[S];
                diagR[1] = L_b[C];
                ndiagR1[1] = L_b[N];                    
                ndiagR2[1] = L_b[NW];

                rhs[1] = resid[2*(Nx+1)+1];

                for(size_t j=3; j<Ny-2; j++)  
                {
                    ndiagL2[j-3] = L_b[SE];
                    ndiagL1[j-2] = L_b[S];
                    diagR[j-1] = L_b[C];
                    ndiagR1[j-1] = L_b[N];                  
                    ndiagR2[j-1] = L_b[NW];
                        
                    rhs[j-1] = resid[j*(Nx+1)+1];
                }
                    
                ndiagL2[Ny-5] = L_b[SE];
                ndiagL1[Ny-4] = L_b[S];
                diagR[Ny-3] = L_b[C];
                ndiagR1[Ny-3] = L_b[N];

                rhs[Ny-3] = resid[(Ny-2)*(Nx+1)+1];
                    
                L_c = stencil.get_L_nw(1,Ny-1,Nx,Ny);
                J_c_x = stencil.getJx(NW);
                J_c_y = stencil.getJy(NW);

                ndiagL2[Ny-4] = L_c[NE];
                ndiagL1[Ny-3] = L_c[S];
                diagR[Ny-2] = L_c[C];

                rhs[Ny-2] = resid[(Ny-1)*(Nx+1)+1];

                // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<Ny-2; k++)  
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
            ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
            diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
            rhs[Ny-2] = rhs[Ny-2] - ndiagL1[Ny-2-1] * rhs[Ny-3];

            // solve the linear system of equations R u = rhs
            temp[1+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];
            for(size_t k=Ny-2; k>1; k--)
            {
                temp[1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[1+(k+1)*(Nx+1)] );

                rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(Nx+1)+1];
            }
            temp[Nx+1+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(Nx+1)+1] );

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
       
                    rhs[0] = resid[Nx+1+i];

                    ndiagL1[0] = L[S];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[N];                  
                    ndiagR2[1] = L[NE];

                    rhs[1] = resid[2*(Nx+1)+i];                     

                    for(size_t j=3; j<Ny-2; j++)
                    {
                       ndiagL2[j-3] = L[SW];
                       ndiagL1[j-2] = L[S];
                       diagR[j-1] = L[C];
                       ndiagR1[j-1] = L[N];
                       ndiagR2[j-1] = L[NE];
 
                       rhs[j-1] = resid[j*(Nx+1)+i];                       
                    }
                                        
                    ndiagL2[Nx-5] = L[SW];
                    ndiagL1[Nx-4] = L[S];
                    diagR[Nx-3] = L[C];
                    ndiagR1[Nx-3] = L[N];

                    rhs[Ny-3] = resid[(Ny-2)*(Nx+1)+i];

                    L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
                    J_b_x = stencil.getJx(N);
                    J_b_y = stencil.getJy(N);

                    ndiagL2[Nx-4] = L_b[SW];
                    ndiagL1[Nx-3] = L_b[S];
                    diagR[Nx-2] = L_b[C];

                    rhs[Ny-2] = resid[(Ny-1)*(Nx+1)+i];
                    
                    // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<Ny-2; k++)  
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
                ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
                diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
                rhs[Ny-2] = rhs[Ny-2] - ndiagL1[Ny-2-1] * rhs[Ny-3];

                // solve the linear system of equations R u = rhs
                temp[i+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];                    
                for(size_t k=Ny-2; k>1; k--)
                {
                    temp[i+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[i+(k+1)*(Nx+1)] );
                    rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(Nx+1)+i];                        
                }
                temp[Nx+1+i] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(Nx+1)+i] );
            }

////////////////letzte Spalte//////////////////
            
            L_c = stencil.get_L_se(Nx-1,1,Nx,Ny);
            J_c_x = stencil.getJx(SE);
            J_c_y = stencil.getJy(SE);

            rhs[0] = resid[Nx+1+Nx-1];

            L_b = stencil.get_L_e(Nx-1,2,Nx,Ny);
            J_b_x = stencil.getJx(E);
            J_b_y = stencil.getJy(E);

            ndiagL1[0] = L_b[S];
            diagR[1] = L_b[C];
            ndiagR1[1] = L_b[N];                    
            ndiagR2[1] = L_b[NE];

            rhs[1] = resid[2*(Nx+1)+Nx-1];

            for(size_t j=3; j<Ny-2; j++)  
            {
                ndiagL2[j-3] = L_b[SE];
                ndiagL1[j-2] = L_b[S];
                diagR[j-1] = L_b[C];
                ndiagR1[j-1] = L_b[N];                  
                ndiagR2[j-1] = L_b[NE];
                    
                rhs[j-1] = resid[j*(Nx+1)+Nx-1];
            }
                
            ndiagL2[Ny-5] = L_b[SE];
            ndiagL1[Ny-4] = L_b[S];
            diagR[Ny-3] = L_b[C];
            ndiagR1[Ny-3] = L_b[N];

            rhs[Ny-3] = resid[(Ny-2)*(Nx+1)+Nx-1];
                
            L_c = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
            J_c_x = stencil.getJx(NE);
            J_c_y = stencil.getJy(NE);

            ndiagL2[Ny-4] = L_c[NE];
            ndiagL1[Ny-3] = L_c[S];
            diagR[Ny-2] = L_c[C];

            rhs[Ny-2] = resid[(Ny-1)*(Nx+1)+Nx-1];

            // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<Ny-2; k++)  
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
            ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
            diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
            rhs[Ny-2] = rhs[Ny-2] - ndiagL1[Ny-2-1] * rhs[Ny-3];

            // solve the linear system of equations R u = rhs
            temp[Nx-1+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];
            for(size_t k=Ny-2; k>1; k--)
            {
                temp[Nx-1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[Nx-1+(k+1)*(Nx+1)] );
                rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(Nx+1)+Nx-1];
            }
            temp[Nx+1+Nx-1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(Nx+1)+Nx-1] );
            
            u += omega_ * temp;

            resid = residuum(u,fv,stencil,Nx,Ny);

            temp = 0.0;

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
       
                    rhs[0] = resid[Nx+1+i];

                    ndiagL1[0] = L[S];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[N];                  
                    ndiagR2[1] = L[NE];

                    rhs[1] = resid[2*(Nx+1)+i];                     

                    for(size_t j=3; j<Ny-2; j++)
                    {
                       ndiagL2[j-3] = L[SW];
                       ndiagL1[j-2] = L[S];
                       diagR[j-1] = L[C];
                       ndiagR1[j-1] = L[N];
                       ndiagR2[j-1] = L[NE];
 
                       rhs[j-1] = resid[j*(Nx+1)+i];                       
                    }
                                        
                    ndiagL2[Nx-5] = L[SW];
                    ndiagL1[Nx-4] = L[S];
                    diagR[Nx-3] = L[C];
                    ndiagR1[Nx-3] = L[N];

                    rhs[Ny-3] = resid[(Ny-2)*(Nx+1)+i];

                    L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
                    J_b_x = stencil.getJx(N);
                    J_b_y = stencil.getJy(N);

                    ndiagL2[Nx-4] = L_b[SW];
                    ndiagL1[Nx-3] = L_b[S];
                    diagR[Nx-2] = L_b[C];

                    rhs[Ny-2] = resid[(Ny-1)*(Nx+1)+i];
                    
                    // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<Ny-2; k++)  
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
                ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
                diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
                rhs[Ny-2] = rhs[Ny-2] - ndiagL1[Ny-2-1] * rhs[Ny-3];

                // solve the linear system of equations R u = rhs
                temp[i+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];                    
                for(size_t k=Ny-2; k>1; k--)
                {
                    temp[i+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[i+(k+1)*(Nx+1)] );
                    rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(Nx+1)+i];                        
                }
                temp[Nx+1+i] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(Nx+1)+i] );
            }

            u += omega_ * temp;

        }
        else // stencil not constant
            {
                NumericArray L = stencil.get_L_c(2,2,Nx,Ny);
                PositionArray J_x = stencil.getJx(C);
                PositionArray J_y = stencil.getJy(C);
                
                NumericArray L_b = stencil.get_L_w(1,2,Nx,Ny);
                PositionArray J_b_x = stencil.getJx(W);
                PositionArray J_b_y = stencil.getJy(W);
 
                NumericArray L_c = stencil.get_L_sw(1,1,Nx,Ny);
                PositionArray J_c_x = stencil.getJx(SW);
                PositionArray J_c_y = stencil.getJy(SW);

                
                diagR[0] = L_c[C];
                ndiagR1[0] = L_c[N];
                ndiagR2[0] = L_c[NW];
                                
                rhs[0] = resid[Nx+1+1];
                
                ndiagL1[0] = L_b[S];
                diagR[1] = L_b[C];
                ndiagR1[1] = L_b[N];                    
                ndiagR2[1] = L_b[NW];

                rhs[1] = resid[2*(Nx+1)+1];

                for(size_t j=3; j<Ny-2; j++)  
                {
                    L_b = stencil.get_L_w(1,j,Nx,Ny);

                    ndiagL2[j-3] = L_b[SE];
                    ndiagL1[j-2] = L_b[S];
                    diagR[j-1] = L_b[C];
                    ndiagR1[j-1] = L_b[N];                  
                    ndiagR2[j-1] = L_b[NW];
                        
                    rhs[j-1] = resid[j*(Nx+1)+1];
                }

                L_b = stencil.get_L_w(1,Ny-2,Nx,Ny);
                
                ndiagL2[Ny-5] = L_b[SE];
                ndiagL1[Ny-4] = L_b[S];
                diagR[Ny-3] = L_b[C];
                ndiagR1[Ny-3] = L_b[N];

                rhs[Ny-3] = resid[(Ny-2)*(Nx+1)+1];
                    
                L_c = stencil.get_L_nw(1,Ny-1,Nx,Ny);
                J_c_x = stencil.getJx(NW);
                J_c_y = stencil.getJy(NW);

                ndiagL2[Ny-4] = L_c[NE];
                ndiagL1[Ny-3] = L_c[S];
                diagR[Ny-2] = L_c[C];

                rhs[Ny-2] = resid[(Ny-1)*(Nx+1)+1];

                for(size_t k=1; k<Ny-2; k++)  
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
                ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
                diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
                rhs[Ny-2] = rhs[Ny-2] - ndiagL1[Ny-2-1] * rhs[Ny-3];

                // solve the linear system of equations R u = rhs
            temp[1+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];
            for(size_t k=Ny-2; k>1; k--)
            {
                temp[1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[1+(k+1)*(Nx+1)] );

                rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(Nx+1)+1];
            }
            temp[Nx+1+1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(Nx+1)+1] );

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
       
                    rhs[0] = resid[Nx+1+i];

                    L = stencil.get_L_c(i,2,Nx,Ny);

                    ndiagL1[0] = L[S];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[N];                  
                    ndiagR2[1] = L[NE];

                    rhs[1] = resid[2*(Nx+1)+i];  

                    for(size_t j=3; j<Ny-2; j++)
                    {
                       L = stencil.get_L_c(i,j,Nx,Ny);
                       ndiagL2[j-3] = L[SW];
                       ndiagL1[j-2] = L[S];
                       diagR[j-1] = L[C];
                       ndiagR1[j-1] = L[N];
                       ndiagR2[j-1] = L[NE];
 
                       rhs[j-1] = resid[j*(Nx+1)+i];
                    }
                                        
                    L = stencil.get_L_c(i,Ny-2,Nx,Ny);
                    ndiagL2[Nx-5] = L[SW];
                    ndiagL1[Nx-4] = L[S];
                    diagR[Nx-3] = L[C];
                    ndiagR1[Nx-3] = L[N];

                    rhs[Ny-3] = resid[(Ny-2)*(Nx+1)+i];

                    L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
                    J_b_x = stencil.getJx(N);
                    J_b_y = stencil.getJy(N);

                    ndiagL2[Nx-4] = L_b[SW];
                    ndiagL1[Nx-3] = L_b[S];
                    diagR[Nx-2] = L_b[C];

                    rhs[Ny-2] = resid[(Ny-1)*(Nx+1)+i];
                    
                    // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<Ny-2; k++)  
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
                ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
                diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
                rhs[Ny-2] = rhs[Ny-2] - ndiagL1[Ny-2-1] * rhs[Ny-3];

                // solve the linear system of equations R u = rhs
                temp[i+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];                    
                for(size_t k=Ny-2; k>1; k--)
                {
                    temp[i+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[i+(k+1)*(Nx+1)] );
                    rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(Nx+1)+i];                        
                }
                temp[Nx+1+i] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(Nx+1)+i] );
            }

////////////////letzte Spalte//////////////////
            
            L_c = stencil.get_L_se(Nx-1,1,Nx,Ny);
            J_c_x = stencil.getJx(SE);
            J_c_y = stencil.getJy(SE);

            rhs[0] = resid[Nx+1+Nx-1];

            L_b = stencil.get_L_e(Nx-1,2,Nx,Ny);
            J_b_x = stencil.getJx(E);
            J_b_y = stencil.getJy(E);

            ndiagL1[0] = L_b[S];
            diagR[1] = L_b[C];
            ndiagR1[1] = L_b[N];                    
            ndiagR2[1] = L_b[NE];

            rhs[1] = resid[2*(Nx+1)+Nx-1];

            for(size_t j=3; j<Ny-2; j++)  
            {
                L_b = stencil.get_L_e(Nx-1,j,Nx,Ny);
                
                ndiagL2[j-3] = L_b[SE];
                ndiagL1[j-2] = L_b[S];
                diagR[j-1] = L_b[C];
                ndiagR1[j-1] = L_b[N];                  
                ndiagR2[j-1] = L_b[NE];
                    
                rhs[j-1] = resid[j*(Nx+1)+Nx-1];
            }
                
            L_b = stencil.get_L_e(Nx-1,Ny-2,Nx,Ny);

            ndiagL2[Ny-5] = L_b[SE];
            ndiagL1[Ny-4] = L_b[S];
            diagR[Ny-3] = L_b[C];
            ndiagR1[Ny-3] = L_b[N];

            rhs[Ny-3] = resid[(Ny-2)*(Nx+1)+Nx-1];
                
            L_c = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
            J_c_x = stencil.getJx(NE);
            J_c_y = stencil.getJy(NE);

            ndiagL2[Ny-4] = L_c[NE];
            ndiagL1[Ny-3] = L_c[S];
            diagR[Ny-2] = L_c[C];

            rhs[Ny-2] = resid[(Ny-1)*(Nx+1)+Nx-1];

            // LR-decomposition + transformation of the rhs
            for(size_t k=1; k<Ny-2; k++)  
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
            ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
            diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
            rhs[Ny-2] = rhs[Ny-2] - ndiagL1[Ny-2-1] * rhs[Ny-3];

            // solve the linear system of equations R u = rhs
            temp[Nx-1+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];
            for(size_t k=Ny-2; k>1; k--)
            {
                temp[Nx-1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[Nx-1+(k+1)*(Nx+1)] );
                rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(Nx+1)+Nx-1];
            }
            temp[Nx+1+Nx-1] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(Nx+1)+Nx-1] );
            
            u += omega_ * temp;

            resid = residuum(u,fv,stencil,Nx,Ny);

            temp = 0.0;

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
       
                    rhs[0] = resid[Nx+1+i];

                    L = stencil.get_L_c(i,2,Nx,Ny);

                    ndiagL1[0] = L[S];
                    diagR[1] = L[C];
                    ndiagR1[1] = L[N];                  
                    ndiagR2[1] = L[NE];

                    rhs[1] = resid[2*(Nx+1)+i];  

                    for(size_t j=3; j<Ny-2; j++)
                    {
                       L = stencil.get_L_c(i,j,Nx,Ny);
                       ndiagL2[j-3] = L[SW];
                       ndiagL1[j-2] = L[S];
                       diagR[j-1] = L[C];
                       ndiagR1[j-1] = L[N];
                       ndiagR2[j-1] = L[NE];
 
                       rhs[j-1] = resid[j*(Nx+1)+i];
                    }
                                        
                    L = stencil.get_L_c(i,Ny-2,Nx,Ny);
                    ndiagL2[Nx-5] = L[SW];
                    ndiagL1[Nx-4] = L[S];
                    diagR[Nx-3] = L[C];
                    ndiagR1[Nx-3] = L[N];

                    rhs[Ny-3] = resid[(Ny-2)*(Nx+1)+i];

                    L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
                    J_b_x = stencil.getJx(N);
                    J_b_y = stencil.getJy(N);

                    ndiagL2[Nx-4] = L_b[SW];
                    ndiagL1[Nx-3] = L_b[S];
                    diagR[Nx-2] = L_b[C];

                    rhs[Ny-2] = resid[(Ny-1)*(Nx+1)+i];
                    
                    // LR-decomposition + transformation of the rhs
                for(size_t k=1; k<Ny-2; k++)  
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
                ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
                diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
                rhs[Ny-2] = rhs[Ny-2] - ndiagL1[Ny-2-1] * rhs[Ny-3];

                // solve the linear system of equations R u = rhs
                temp[i+(Nx+1)*(Ny-1)] = rhs[Ny-2] / diagR[Ny-2];                    
                for(size_t k=Ny-2; k>1; k--)
                {
                    temp[i+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[k-1] - ndiagR1[k-1] * temp[i+(k+1)*(Nx+1)] );
                    rhs[k-2] -= ndiagR2[k-2] * temp[(k+1)*(Nx+1)+i];                        
                }
                temp[Nx+1+i] = 1/diagR[0] * ( rhs[0] - ndiagR1[0] * temp[2*(Nx+1)+i] );
            }

            u += omega_ * temp;
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
}
