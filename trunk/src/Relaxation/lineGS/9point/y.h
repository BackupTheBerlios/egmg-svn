namespace mg
{
	void ZebraLineGS::ninepointyline(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> &rhs, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const
					   
    { 
		//valarrays needed for saving the tridiagonal matrix A of linear system A u = rhs
		std::valarray<Precision> diagR(Ny-1);
		std::valarray<Precision> ndiagR(Ny-2);
		std::valarray<Precision> ndiagL(Ny-2);

		if(stencil.isConstant() == true)
		{
			// get const operator L
			const std::valarray<Precision> L = stencil.get_L_c(2,2,Nx,Ny);
			const std::valarray<int> J_x = stencil.getJx(C);
			const std::valarray<int> J_y = stencil.getJy(C);

			// for each line: correction of the rhs given by rhs = fv - [L[w]  0  L[e]] * u and elimination of the 
			// boundary condition in first and last inner point
			for(size_t i=1; i<Nx ; i++) 
			{
				rhs[i+(Nx+1)] = fv[i+(Nx+1)] - L[W] * u[i+J_x[W]+(Nx+1)] - L[E] * u[i+J_x[E]+(Nx+1)] 
					- L[S] * u[i+(1+J_y[S])*(Nx+1)];
				
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[i+Nx+1] -= L[sum] * u[i+J_x[sum]+(1+J_y[sum])*(Nx+1)];
				}
				
				for(size_t j=2; j<Ny-1; j++)  
				{
					rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[W] * u[i+J_x[W]+j*(Nx+1)] - L[E] * u[i+J_x[E]+j*(Nx+1)];

					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					}
				}
				rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L[W] * u[i+J_x[W]+(Ny-1)*(Nx+1)] 
					- L[E] * u[i+J_x[E]+(Ny-1)*(Nx+1)] - L[N] * u[i+(Ny-1+J_y[N])*(Nx+1)];
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[i+(Ny-1)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
				}

				// set tridiagonalmatrix for solving A u = rhs
				// A[i][i] = L[c]; A[i-1][i] = L[s]; A[i+1][i] = L[n]
				diagR = L[C];
				ndiagR = L[N];
				ndiagL = L[S];

				// LR-decomposition + transformation of the rhs vector
				for(size_t j=1; j<Ny-1; j++)  
				{
					ndiagL[j-1] = ndiagL[j-1]/diagR[j-1];  
					diagR[j] -= ndiagL[j-1] * ndiagR[j-1]; 
					rhs[(j+1)*(Nx+1) + i] = rhs[(j+1)*(Nx+1) + i] - ndiagL[j-1] * rhs[j*(Nx+1)+i];  
				}
				
				// solve the linear system of equations R u = rhs
				u[i+(Nx+1)*(Ny-1)] = rhs[i+(Nx+1)*(Ny-1)] / diagR[Ny-2];

				for(size_t k=Ny-2; k>0; k--)
				{
					u[i+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[i+k*(Nx+1)] - ndiagR[k-1] * u[i+(k+1)*(Nx+1)] );
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

			if(Ny > 2)
			{
				L = stencil.get_L_sw(1,1,Nx,Ny);
					
				diagR[0] = L[C];
				ndiagR[0] = L[N];	
				
				rhs[1+(Nx+1)] = fv[1+(Nx+1)] - L[W] * u[1+J_x[W]+(Nx+1)] - L[E] * u[1+J_x[E]+(Nx+1)] 
						- L[S] * u[1+(1+J_y[S])*(Nx+1)];
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[1+Nx+1] -= L[sum] * u[1+J_x[sum]+(1+J_y[sum])*(Nx+1)];
				}
				
				for(size_t j=2; j<Ny-1; j++) 
				{
					L = stencil.get_L_w(1,j,Nx,Ny);
					diagR[j-1] = L[C];
					ndiagR[j-1] = L[N];
					ndiagL[j-2] = L[S];

					rhs[1+j*(Nx+1)] = fv[1+j*(Nx+1)] - L[W] * u[1+J_x[W]+j*(Nx+1)] - L[E] * u[1+J_x[E]+j*(Nx+1)];
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[1+j*(Nx+1)] -= L[sum] * u[1+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					}
				}
					
				L = stencil.get_L_nw(1,Ny-1,Nx,Ny);
					
				diagR[Ny-2] = L[C];
				ndiagL[Ny-3] = L[S];

				rhs[1+(Ny-1)*(Nx+1)] = fv[1+(Ny-1)*(Nx+1)] - L[W] * u[1+J_x[W]+(Ny-1)*(Nx+1)] 
						- L[E] * u[1+J_x[E]+(Ny-1)*(Nx+1)] - L[N] * u[1+(Ny-1+J_y[N])*(Nx+1)];
					
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[1+(Ny-1)*(Nx+1)] -= L[sum] * u[1+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
				}
					
				for(size_t j=1; j<Ny-1; j++)  
				{
					ndiagL[j-1] = ndiagL[j-1]/diagR[j-1];  
					diagR[j] -= ndiagL[j-1] * ndiagR[j-1];
					rhs[(j+1)*(Nx+1) + 1] = rhs[(j+1)*(Nx+1) + 1] - ndiagL[j-1] * rhs[j*(Nx+1)+1];  
				}

				u[1+(Nx+1)*(Ny-1)] = rhs[1+(Nx+1)*(Ny-1)] / diagR[Ny-2];
				
				for(size_t k=Ny-2; k>0; k--)
				{
					u[1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[1+k*(Nx+1)] - ndiagR[k-1] * u[1+(k+1)*(Nx+1)] );
				}
				
				
				for(size_t i=2; i<Nx-1 ; i++)
				{
					L = stencil.get_L_s(i,1,Nx,Ny);
					
					diagR[0] = L[C];
					ndiagR[0] = L[N];	
					
					rhs[i+(Nx+1)] = fv[i+(Nx+1)] - L[W] * u[i+J_x[W]+(Nx+1)] - L[E] * u[i+J_x[E]+(Nx+1)] 
						- L[S] * u[i+(1+J_y[S])*(Nx+1)];
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[i+Nx+1] -= L[sum] * u[i+J_x[sum]+(1+J_y[sum])*(Nx+1)];
					}
					
					for(size_t j=2; j<Ny-1; j++) 
					{
						L = stencil.get_L_c(i,j,Nx,Ny);
						diagR[j-1] = L[C];
						ndiagR[j-1] = L[N];
						ndiagL[j-2] = L[S];

						rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[W] * u[i+J_x[W]+j*(Nx+1)] - L[E] * u[i+J_x[E]+j*(Nx+1)];
						for(size_t sum=5; sum<L.size(); sum++)
						{
							rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
						}
					}
					
					L = stencil.get_L_n(i,Ny-1,Nx,Ny);
					
					diagR[Ny-2] = L[C];
					ndiagL[Ny-3] = L[S];

					rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L[W] * u[i+J_x[W]+(Ny-1)*(Nx+1)] 
						- L[E] * u[i+J_x[E]+(Ny-1)*(Nx+1)] - L[N] * u[i+(Ny-1+J_y[N])*(Nx+1)];
					
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[i+(Ny-1)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
					}
					
					for(size_t j=1; j<Ny-1; j++)  
					{
						ndiagL[j-1] = ndiagL[j-1]/diagR[j-1];  
						diagR[j] -= ndiagL[j-1] * ndiagR[j-1];
						rhs[(j+1)*(Nx+1) + i] = rhs[(j+1)*(Nx+1) + i] - ndiagL[j-1] * rhs[j*(Nx+1)+i];  
					}

					u[i+(Nx+1)*(Ny-1)] = rhs[i+(Nx+1)*(Ny-1)] / diagR[Ny-2];
					
					for(size_t k=Ny-2; k>0; k--)
					{
						u[i+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[i+k*(Nx+1)] - ndiagR[k-1] * u[i+(k+1)*(Nx+1)] );
					}
				}

				L = stencil.get_L_se(Nx-1,1,Nx,Ny);
					
				diagR[0] = L[C];
				ndiagR[0] = L[N];	
					
				rhs[Nx-1+(Nx+1)] = fv[Nx-1+(Nx+1)] - L[W] * u[Nx-1+J_x[W]+(Nx+1)] - L[E] * u[Nx-1+J_x[E]+(Nx+1)] 
						- L[S] * u[Nx-1+(1+J_y[S])*(Nx+1)];
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[Nx-1+Nx+1] -= L[sum] * u[Nx-1+J_x[sum]+(1+J_y[sum])*(Nx+1)];
				}
					
				for(size_t j=2; j<Ny-1; j++) 
				{
					L = stencil.get_L_e(Nx-1,j,Nx,Ny);
					diagR[j-1] = L[C];
					ndiagR[j-1] = L[N];
					ndiagL[j-2] = L[S];

					rhs[Nx-1+j*(Nx+1)] = fv[Nx-1+j*(Nx+1)] - L[W] * u[Nx-1+J_x[W]+j*(Nx+1)] - L[E] * u[Nx-1+J_x[E]+j*(Nx+1)];
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[Nx-1+j*(Nx+1)] -= L[sum] * u[Nx-1+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					}
				}
					
				L = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
					
				diagR[Ny-2] = L[C];
				ndiagL[Ny-3] = L[S];

				rhs[Nx-1+(Ny-1)*(Nx+1)] = fv[Nx-1+(Ny-1)*(Nx+1)] - L[W] * u[Nx-1+J_x[W]+(Ny-1)*(Nx+1)] 
						- L[E] * u[Nx-1+J_x[E]+(Ny-1)*(Nx+1)] - L[N] * u[Nx-1+(Ny-1+J_y[N])*(Nx+1)];
					
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[Nx-1+(Ny-1)*(Nx+1)] -= L[sum] * u[Nx-1+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
				}
					
				for(size_t j=1; j<Ny-1; j++)  
				{
					ndiagL[j-1] = ndiagL[j-1]/diagR[j-1];  
					diagR[j] -= ndiagL[j-1] * ndiagR[j-1];
					rhs[(j+1)*(Nx+1) + Nx-1] = rhs[(j+1)*(Nx+1) + Nx-1] - ndiagL[j-1] * rhs[j*(Nx+1)+Nx-1];  
				}

				u[Nx-1+(Nx+1)*(Ny-1)] = rhs[Nx-1+(Nx+1)*(Ny-1)] / diagR[Ny-2];
					
				for(size_t k=Ny-2; k>0; k--)
				{
					u[Nx-1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[Nx-1+k*(Nx+1)] - ndiagR[k-1] * u[Nx-1+(k+1)*(Nx+1)] );
				}
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
}

