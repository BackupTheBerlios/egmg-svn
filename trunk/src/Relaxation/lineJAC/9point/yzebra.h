namespace mg
{
	void ZebraLineJAC::ninepointyzebra(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> resid, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const
					   
    { 
        std::valarray<Precision> rhs(0.0,Ny+1);
		std::valarray<Precision> temp(0.0,(Nx+1)*(Ny+1));

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
			std::valarray<Precision> L = stencil.get_L_c(2,2,Nx,Ny);
			std::valarray<int> J_x = stencil.getJx(C);
			std::valarray<int> J_y = stencil.getJy(C);

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
}

