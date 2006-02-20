namespace mg
{
	void lineJAC::ninepointxline(std::valarray<precision> &u, const std::valarray<precision> &fv, 
		                    std::valarray<precision> resid, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const
					   
    { 
		std::valarray<precision> rhs(0.0,Nx-1);
		std::valarray<precision> temp(0.0,(Nx+1)*(Ny+1));		
		
		//valarrays needed for saving the tridiagonal matrix A of linear system A u = rhs		
		std::valarray<precision> diagR(Nx-1);
		std::valarray<precision> ndiagR(Nx-2);
		std::valarray<precision> ndiagL(Nx-2);

		if(stencil.is_constant() == true)
		{
			// get const operator L
			const std::valarray<precision> L = stencil.get_L_c(2,2,Nx,Ny);
			const std::valarray<int> J_x = stencil.get_J_x(c);
			const std::valarray<int> J_y = stencil.get_J_y(c);

			for(size_t i=1; i<Ny ; i++) 
			{
                for(size_t j=0; j<Nx-1; j++)  
				{
					rhs[j] = resid[i*(Nx+1)+j+1];
				}								
				// set tridiagonalmatrix for solving A u = rhs
				// A[i][i] = L[c]; A[i-1][i] = L[w]; A[i+1][i] = L[e]
				diagR = L[c];
				ndiagR = L[e];
				ndiagL = L[w];

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
			std::valarray<precision> L = stencil.get_L_c(2,2,Nx,Ny);
			std::valarray<int> J_x = stencil.get_J_x(c);
			std::valarray<int> J_y = stencil.get_J_y(c);
			
			if(Nx > 2)
			{
				L = stencil.get_L_sw(1,1,Nx,Ny);
				diagR[0] = L[c];
				ndiagR[0] = L[e];

				rhs[0] = resid[Nx+1+1];

				for(size_t j=2; j<Nx-1; j++)
				{
					L = stencil.get_L_s(j,1,Nx,Ny);
					diagR[j-1] = L[c];
					ndiagR[j-1] = L[e];
					ndiagL[j-2] = L[w];
					
					rhs[j-1] = resid[Nx+1+j];
				}

				L = stencil.get_L_se(Nx-1,1,Nx,Ny);
					
				diagR[Nx-2] = L[c];
				ndiagL[Nx-3] = L[w];
					
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
				
				for(size_t i=2; i<Ny-1 ; i++)
				{
					L = stencil.get_L_w(1,i,Nx,Ny);
					diagR[0] = L[c];
					ndiagR[0] = L[e];

					rhs[0] = resid[i*(Nx+1)+1];
					
					for(size_t j=2; j<Nx-1; j++)
					{
						L = stencil.get_L_c(j,i,Nx,Ny);
						diagR[j-1] = L[c];
						ndiagR[j-1] = L[e];
						ndiagL[j-2] = L[w];
						
						rhs[j-1] = resid[i*(Nx+1)+j];
					}

					L = stencil.get_L_e(Nx-1,i,Nx,Ny);
					
					diagR[Nx-2] = L[c];
					ndiagL[Nx-3] = L[w];
					
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
				diagR[0] = L[c];
				ndiagR[0] = L[e];

				rhs[0] = resid[(Ny-1)*(Nx+1)+1];
					
				for(size_t j=2; j<Nx-1; j++)
				{
					L = stencil.get_L_n(j,Ny-1,Nx,Ny);
					diagR[j-1] = L[c];
					ndiagR[j-1] = L[e];
					ndiagL[j-2] = L[w];
					 
					rhs[j-1] = resid[(Ny-1)*(Nx+1)+j];
				}

				L = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
					
				diagR[Nx-2] = L[c];
				ndiagL[Nx-3] = L[w];
					
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
			}
			
			
			else  // if Nx and Ny are too small do one GS_lex step
			{
				precision temp=0;

				for(size_t k=1; k<Ny; k++)
				{
					L = stencil.get_L_c(1,k,Nx,Ny);

					temp=0;

					for(size_t sum=5; sum < L.size(); sum++)
					{
						temp -= L[sum] * u[1+J_x[sum]+(k+J_y[sum])*(Nx+1)];
					}
					u[1+k*(Nx+1)] = 1/L[c] * ( fv[1+k*(Nx+1)] - L[w] * u[1+J_x[w]+k*(Nx+1)] - L[e] * u[1+J_x[e]+k*(Nx+1)] 
						          - L[n] * u[1+(k+J_y[n])*(Nx+1)] - L[s] * u[1+(k+J_y[s])*(Nx+1)] - temp );
				}
			}
		}		
    }
}

