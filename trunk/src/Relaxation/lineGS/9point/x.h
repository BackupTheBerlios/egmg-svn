namespace mg
{
	void ZebraLineGS::ninepointxline(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> &rhs, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const
					   
    { 
		//valarrays needed for saving the tridiagonal matrix A of linear system A u = rhs		
		std::valarray<Precision> diagR(Nx-1);
		std::valarray<Precision> ndiagR(Nx-2);
		std::valarray<Precision> ndiagL(Nx-2);

		if(stencil.isConstant() == true)
		{
			// get const operator L
			const std::valarray<Precision> L = stencil.get_L_c(2,2,Nx,Ny);
			const std::valarray<int> J_x = stencil.getJx(C);
			const std::valarray<int> J_y = stencil.getJy(C);
			
			// for each line: correction of the rhs given by rhs = fv - [L[n]  0  L[s]]^t * u and elimination of the 
			// boundary condition in first and last inner point
			for(size_t i=1; i<Ny ; i++) 
			{
				rhs[1+i*(Nx+1)] = fv[1+i*(Nx+1)] - L[N] * u[1+(i+J_y[N])*(Nx+1)] - L[S] * u[1+(i+J_y[S])*(Nx+1)]
					- L[W] * u[i*(Nx+1)];
				
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[1+i*(Nx+1)] -= L[sum] * u[1+J_x[sum]+(i+J_y[sum])*(Nx+1)];
				}
				

				for(size_t j=2; j<Nx-1; j++)  
				{
					rhs[j+i*(Nx+1)] = fv[j+i*(Nx+1)] - L[N] * u[j+(i+J_y[N])*(Nx+1)] - L[S] * u[j+(i+J_y[S])*(Nx+1)];

					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[j+i*(Nx+1)] -= L[sum] * u[j+J_x[sum]+(i+J_y[sum])*(Nx+1)];
					}

				}
				
				rhs[(Nx-1)+i*(Nx+1)] = fv[Nx-1+i*(Nx+1)] - L[N] * u[Nx-1+(i+J_y[N])*(Nx+1)] 
					- L[S] * u[Nx-1+(i+J_y[S])*(Nx+1)] - L[E] * u[Nx+i*(Nx+1)];

				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[Nx-1+i*(Nx+1)] -= L[sum] * u[Nx-1+J_x[sum]+(i+J_y[sum])*(Nx+1)];
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
					rhs[i*(Nx+1) + 1 + k] = rhs[i*(Nx+1) + 1 + k] - ndiagL[k-1] * rhs[i*(Nx+1)+k];  
				}

				// solve the linear system of equations R u = rhs
				u[i*(Nx+1)+(Nx-1)] = rhs[i*(Nx+1)+(Nx-1)] / diagR[Nx-2];
				
				for(size_t j=Nx-2; j>0; j--)
				{
					u[i*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[i*(Nx+1)+j] - ndiagR[j-1] * u[i*(Nx+1)+j+1] );
				}
			}
		}	

		else
		{
			//Stencil ist not constant, so L needs to be evaluated in each grid point
			//no other change in the algorithm			
			std::valarray<Precision> L = stencil.get_L_c(2,2,Nx,Ny);
			std::valarray<int> J_x = stencil.getJx(C);
			std::valarray<int> J_y = stencil.getJy(C);

			
			if(Nx > 2)
			{
				L = stencil.get_L_sw(1,1,Nx,Ny);
				diagR[0] = L[C];
				ndiagR[0] = L[E];

				rhs[1+Nx+1] = fv[1+Nx+1] - L[N] * u[1+(1+J_y[N])*(Nx+1)] - L[S] * u[1+(1+J_y[S])*(Nx+1)]
						- L[W] * u[Nx+1];
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[1+Nx+1] -= L[sum] * u[1+J_x[sum]+(1+J_y[sum])*(Nx+1)];
				}
					
				for(size_t j=2; j<Nx-1; j++)
				{
					L = stencil.get_L_s(j,1,Nx,Ny);
					diagR[j-1] = L[C];
					ndiagR[j-1] = L[E];
					ndiagL[j-2] = L[W];
					rhs[j+Nx+1] = fv[j+Nx+1] - L[N] * u[j+(1+J_y[N])*(Nx+1)] - L[S] * u[j+(1+J_y[S])*(Nx+1)];
						
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[j+Nx+1] -= L[sum] * u[j+J_x[sum]+(1+J_y[sum])*(Nx+1)];
					}
				}

				L = stencil.get_L_se(Nx-1,1,Nx,Ny);
					
				diagR[Nx-2] = L[C];
				ndiagL[Nx-3] = L[W];
					
				rhs[(Nx-1)+(Nx+1)] = fv[Nx-1+(Nx+1)] - L[N] * u[Nx-1+(1+J_y[N])*(Nx+1)] 
						- L[S] * u[Nx-1+(1+J_y[S])*(Nx+1)] - L[E] * u[Nx+(Nx+1)];
				
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[Nx-1+(Nx+1)] -= L[sum] * u[Nx-1+J_x[sum]+(1+J_y[sum])*(Nx+1)];
				}

				
				for(size_t k=1; k<Nx-1; k++)  
				{
					ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
					diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
					rhs[(Nx+1) + 1 + k] = rhs[(Nx+1) + 1 + k] - ndiagL[k-1] * rhs[(Nx+1)+k];  
				}
					
				
				u[(Nx+1)+(Nx-1)] = rhs[(Nx+1)+(Nx-1)] / diagR[Nx-2];

				for(size_t j=Nx-2; j>0; j--)
				{
					u[(Nx+1)+j] = 1/diagR[j-1] * ( rhs[(Nx+1)+j] - ndiagR[j-1] * u[(Nx+1)+j+1] );
				}
				
				
				for(size_t i=2; i<Ny-1 ; i++)
				{
					L = stencil.get_L_w(1,i,Nx,Ny);
					diagR[0] = L[C];
					ndiagR[0] = L[E];

					rhs[1+i*(Nx+1)] = fv[1+i*(Nx+1)] - L[N] * u[1+(i+J_y[N])*(Nx+1)] - L[S] * u[1+(i+J_y[S])*(Nx+1)]
						- L[W] * u[i*(Nx+1)];
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[1+i*(Nx+1)] -= L[sum] * u[1+J_x[sum]+(i+J_y[sum])*(Nx+1)];
					}
					
					for(size_t j=2; j<Nx-1; j++)
					{
						L = stencil.get_L_c(j,i,Nx,Ny);
						diagR[j-1] = L[C];
						ndiagR[j-1] = L[E];
						ndiagL[j-2] = L[W];
						rhs[j+i*(Nx+1)] = fv[j+i*(Nx+1)] - L[N] * u[j+(i+J_y[N])*(Nx+1)] - L[S] * u[j+(i+J_y[S])*(Nx+1)];
						
						for(size_t sum=5; sum<L.size(); sum++)
						{
							rhs[j+i*(Nx+1)] -= L[sum] * u[j+J_x[sum]+(i+J_y[sum])*(Nx+1)];
						}
					}

					L = stencil.get_L_e(Nx-1,i,Nx,Ny);
					
					diagR[Nx-2] = L[C];
					ndiagL[Nx-3] = L[W];
					
					rhs[(Nx-1)+i*(Nx+1)] = fv[Nx-1+i*(Nx+1)] - L[N] * u[Nx-1+(i+J_y[N])*(Nx+1)] 
						- L[S] * u[Nx-1+(i+J_y[S])*(Nx+1)] - L[E] * u[Nx+i*(Nx+1)];
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[Nx-1+i*(Nx+1)] -= L[sum] * u[Nx-1+J_x[sum]+(i+J_y[sum])*(Nx+1)];
					}

					
					for(size_t k=1; k<Nx-1; k++)  
					{
						ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
						diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
						rhs[i*(Nx+1) + 1 + k] = rhs[i*(Nx+1) + 1 + k] - ndiagL[k-1] * rhs[i*(Nx+1)+k];  
					}
					
					u[i*(Nx+1)+(Nx-1)] = rhs[i*(Nx+1)+(Nx-1)] / diagR[Nx-2];

					for(size_t j=Nx-2; j>0; j--)
					{
						u[i*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[i*(Nx+1)+j] - ndiagR[j-1] * u[i*(Nx+1)+j+1] );
					}
				}

				L = stencil.get_L_nw(1,Ny-1,Nx,Ny);
				diagR[0] = L[C];
				ndiagR[0] = L[E];

				rhs[1+(Ny-1)*(Nx+1)] = fv[1+(Ny-1)*(Nx+1)] - L[N] * u[1+(Ny-1+J_y[N])*(Nx+1)] - L[S] * u[1+(Ny-1+J_y[S])*(Nx+1)]
						- L[W] * u[(Ny-1)*(Nx+1)];
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[1+(Ny-1)*(Nx+1)] -= L[sum] * u[1+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
				}
					
				for(size_t j=2; j<Nx-1; j++)
				{
					L = stencil.get_L_n(j,Ny-1,Nx,Ny);
					diagR[j-1] = L[C];
					ndiagR[j-1] = L[E];
					ndiagL[j-2] = L[W];
					rhs[j+(Ny-1)*(Nx+1)] = fv[j+(Ny-1)*(Nx+1)] - L[N] * u[j+(Ny-1+J_y[N])*(Nx+1)] - L[S] * u[j+(Ny-1+J_y[S])*(Nx+1)];
						
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[j+(Ny-1)*(Nx+1)] -= L[sum] * u[j+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
					}
				}

				L = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
					
				diagR[Nx-2] = L[C];
				ndiagL[Nx-3] = L[W];
					
				rhs[(Nx-1)+(Ny-1)*(Nx+1)] = fv[Nx-1+(Ny-1)*(Nx+1)] - L[N] * u[Nx-1+(Ny-1+J_y[N])*(Nx+1)] 
						- L[S] * u[Nx-1+(Ny-1+J_y[S])*(Nx+1)] - L[E] * u[Nx+(Ny-1)*(Nx+1)];
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[Nx-1+(Ny-1)*(Nx+1)] -= L[sum] * u[Nx-1+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
				}

				
				for(size_t k=1; k<Nx-1; k++)  
				{
					ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
					diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
					rhs[(Ny-1)*(Nx+1) + 1 + k] = rhs[(Ny-1)*(Nx+1) + 1 + k] - ndiagL[k-1] * rhs[(Ny-1)*(Nx+1)+k];  
				}
					
				
				u[(Ny-1)*(Nx+1)+(Nx-1)] = rhs[(Ny-1)*(Nx+1)+(Nx-1)] / diagR[Nx-2];

				for(size_t j=Nx-2; j>0; j--)
				{
					u[(Ny-1)*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[(Ny-1)*(Nx+1)+j] - ndiagR[j-1] * u[(Ny-1)*(Nx+1)+j+1] );
				}
			}
			
			
			else  // if Nx and Ny are too small do one GS_lex step
			{
				Precision temp=0;

				for(size_t k=1; k<Ny; k++)
				{
					L = stencil.get_L_c(1,k,Nx,Ny);

					temp=0;

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
}

