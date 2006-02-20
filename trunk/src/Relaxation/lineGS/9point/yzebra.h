namespace mg
{
	void lineGS::ninepointyzebra(std::valarray<precision> &u, const std::valarray<precision> &fv, 
		                    std::valarray<precision> &rhs, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const
					   
    { 
		//valarrays needed for saving the tridiagonal matrix A of linear system A u = rhs
		std::valarray<precision> diagR(Ny-1);
		std::valarray<precision> ndiagR(Ny-2);
		std::valarray<precision> ndiagL(Ny-2);
		
		if(stencil.is_constant() == true)
		{
			// get const operator L
			const std::valarray<precision> L = stencil.get_L_c(2,2,Nx,Ny);
			const std::valarray<int> J_x = stencil.get_J_x(c);
			const std::valarray<int> J_y = stencil.get_J_y(c);
			
			// for each odd line correction of the rhs given by rhs = fv + [1  0  1] * u and elimination of the 
			// boundary condition in first and last inner point  
			for(size_t i=1; i<Nx ; i+=2) // i = Spalte
			{
				rhs[i+(Nx+1)] = fv[i+(Nx+1)] - L[w] * u[i+J_x[w]+(Nx+1)] - L[e] * u[i+J_x[e]+(Nx+1)] - L[s] * u[i+(1+J_y[s])*(Nx+1)];
				
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[i+Nx+1] -= L[sum] * u[i+J_x[sum]+(1+J_y[sum])*(Nx+1)];
				}

				for(size_t j=2; j<Ny-1; j++)  // j = Laufindex in einer Spalte
				{
					rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[w] * u[i+J_x[w]+j*(Nx+1)] - L[e] * u[i+J_x[e]+j*(Nx+1)];
					
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					}
				}

				rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L[w] * u[i+J_x[w]+(Ny-1)*(Nx+1)] - L[e] * u[i+J_x[e]+(Ny-1)*(Nx+1)] 
					- L[n] * u[i+(Ny-1+J_y[n])*(Nx+1)];
				
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[i+(Ny-1)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
				}

                // set tridiagonalmatrix for solving A u = rhs
				// A[i][i] = L[c]; A[i-1][i] = L[s]; A[i+1][i] = L[n]
				diagR = L[c];
				ndiagR = L[n];
				ndiagL = L[s];
				
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

			// same for each even line
			for(size_t i=2; i<Nx ; i+=2) 
			{
				rhs[i+(Nx+1)] = fv[i+(Nx+1)] - L[w] * u[i+J_x[w]+(Nx+1)] - L[e] * u[i+J_x[e]+(Nx+1)] 
					- L[s] * u[i+(1+J_y[s])*(Nx+1)]; 
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[i+Nx+1] -= L[sum] * u[i+J_x[sum]+(1+J_y[sum])*(Nx+1)];
				}

				for(size_t j=2; j<Ny-1; j++) 
				{
					rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[w] * u[i+J_x[w]+j*(Nx+1)] - L[e] * u[i+J_x[e]+j*(Nx+1)];

					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					}
				}
				
				rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L[w] * u[i+J_x[w]+(Ny-1)*(Nx+1)] 
					- L[e] * u[i+J_x[e]+(Ny-1)*(Nx+1)] - L[n] * u[i+(Ny-1+J_y[n])*(Nx+1)];

				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[i+(Ny-1)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
				}

				diagR = L[c];
				ndiagR = L[n];
				ndiagL = L[s];
				
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
		}

		else
		{
			//Stencil ist not constant, so L needs to be evaluated in each grid point
			//no other change in the algorithm	
			std::valarray<precision> L = stencil.get_L_c(2,2,Nx,Ny);
			std::valarray<int> J_x = stencil.get_J_x(c);
			std::valarray<int> J_y = stencil.get_J_y(c);

			if(Ny > 2)
			{
				L = stencil.get_L_sw(1,1,Nx,Ny);
					
				diagR[0] = L[c];
				ndiagR[0] = L[n];
					
				rhs[1+(Nx+1)] = fv[1+(Nx+1)] - L[w] * u[1+J_x[w]+(Nx+1)] - L[e] * u[1+J_x[e]+(Nx+1)] 
						- L[s] * u[1+(1+J_y[s])*(Nx+1)];
					
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[1+Nx+1] -= L[sum] * u[1+J_x[sum]+(1+J_y[sum])*(Nx+1)];
				}

				for(size_t j=2; j<Ny-1; j++) 
				{
					L = stencil.get_L_w(1,j,Nx,Ny);
						
					diagR[j-1] = L[c];
					ndiagR[j-1] = L[n];
					ndiagL[j-2] = L[s];
						
					rhs[1+j*(Nx+1)] = fv[1+j*(Nx+1)] - L[w] * u[1+J_x[w]+j*(Nx+1)] - L[e] * u[1+J_x[e]+j*(Nx+1)];

					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[1+j*(Nx+1)] -= L[sum] * u[1+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					}
						
				}
					
				L = stencil.get_L_nw(1,Ny-1,Nx,Ny);
					
				diagR[Ny-2] = L[c];
				ndiagL[Ny-3] = L[s];
					
				rhs[1+(Ny-1)*(Nx+1)] = fv[1+(Ny-1)*(Nx+1)] - L[w] * u[1+J_x[w]+(Ny-1)*(Nx+1)] 
						- L[e] * u[1+J_x[e]+(Ny-1)*(Nx+1)] - L[n] * u[1+(Ny-1+J_y[n])*(Nx+1)];

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
				
				for(size_t i=3; i<Nx-1 ; i+=2)
				{
					L = stencil.get_L_s(i,1,Nx,Ny);
					
					diagR[0] = L[c];
					ndiagR[0] = L[n];
					
					rhs[i+(Nx+1)] = fv[i+(Nx+1)] - L[w] * u[i+J_x[w]+(Nx+1)] - L[e] * u[i+J_x[e]+(Nx+1)] 
						- L[s] * u[i+(1+J_y[s])*(Nx+1)];
					
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[i+Nx+1] -= L[sum] * u[i+J_x[sum]+(1+J_y[sum])*(Nx+1)];
					}

					for(size_t j=2; j<Ny-1; j++) 
					{
						L = stencil.get_L_c(i,j,Nx,Ny);
						
						diagR[j-1] = L[c];
						ndiagR[j-1] = L[n];
						ndiagL[j-2] = L[s];
						
						rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[w] * u[i+J_x[w]+j*(Nx+1)] - L[e] * u[i+J_x[e]+j*(Nx+1)];

						for(size_t sum=5; sum<L.size(); sum++)
						{
							rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
						}
						
					}
					
					L = stencil.get_L_n(i,Ny-1,Nx,Ny);
					
					diagR[Ny-2] = L[c];
					ndiagL[Ny-3] = L[s];
					
					rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L[w] * u[i+J_x[w]+(Ny-1)*(Nx+1)] 
						- L[e] * u[i+J_x[e]+(Ny-1)*(Nx+1)] - L[n] * u[i+(Ny-1+J_y[n])*(Nx+1)];

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
					
				diagR[0] = L[c];
				ndiagR[0] = L[n];
					
				rhs[Nx-1+(Nx+1)] = fv[Nx-1+(Nx+1)] - L[w] * u[Nx-1+J_x[w]+(Nx+1)] - L[e] * u[Nx-1+J_x[e]+(Nx+1)] 
						- L[s] * u[Nx-1+(1+J_y[s])*(Nx+1)];
					
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[Nx-1+Nx+1] -= L[sum] * u[Nx-1+J_x[sum]+(1+J_y[sum])*(Nx+1)];
				}

				for(size_t j=2; j<Ny-1; j++) 
				{
					L = stencil.get_L_e(Nx-1,j,Nx,Ny);
					
					diagR[j-1] = L[c];
					ndiagR[j-1] = L[n];
					ndiagL[j-2] = L[s];
						
					rhs[Nx-1+j*(Nx+1)] = fv[Nx-1+j*(Nx+1)] - L[w] * u[Nx-1+J_x[w]+j*(Nx+1)] - L[e] * u[Nx-1+J_x[e]+j*(Nx+1)];

					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[Nx-1+j*(Nx+1)] -= L[sum] * u[Nx-1+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					}
						
				}
					
				L = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
					
				diagR[Ny-2] = L[c];
				ndiagL[Ny-3] = L[s];
					
				rhs[Nx-1+(Ny-1)*(Nx+1)] = fv[Nx-1+(Ny-1)*(Nx+1)] - L[w] * u[Nx-1+J_x[w]+(Ny-1)*(Nx+1)] 
						- L[e] * u[Nx-1+J_x[e]+(Ny-1)*(Nx+1)] - L[n] * u[Nx-1+(Ny-1+J_y[n])*(Nx+1)];

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

				
				for(size_t i=2; i<Nx ; i+=2)
				{
					L = stencil.get_L_s(i,1,Nx,Ny);
					
					diagR[0] = L[c];
					ndiagR[0] = L[n];
					
					rhs[i+(Nx+1)] = fv[i+(Nx+1)] - L[w] * u[i+J_x[w]+(Nx+1)] - L[e] * u[i+J_x[e]+(Nx+1)] 
						- L[s] * u[i+(1+J_y[s])*(Nx+1)];

					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[i+Nx+1] -= L[sum] * u[i+J_x[sum]+(1+J_y[sum])*(Nx+1)];
					}
					
					for(size_t j=2; j<Ny-1; j++) 
					{
						L = stencil.get_L_c(i,j,Nx,Ny);
						
						diagR[j-1] = L[c];
						ndiagR[j-1] = L[n];										
						ndiagL[j-2] = L[s];
						
						rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[w] * u[i+J_x[w]+j*(Nx+1)] - L[e] * u[i+J_x[e]+j*(Nx+1)]; 
						
						for(size_t sum=5; sum<L.size(); sum++)
						{
							rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
						}
					}

					L = stencil.get_L_n(i,Ny-1,Nx,Ny);
					
					diagR[Nx-2] = L[c];
					ndiagL[Nx-3] = L[s];

					rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L[w] * u[i+J_x[w]+(Ny-1)*(Nx+1)] 
						- L[e] * u[i+J_x[e]+(Ny-1)*(Nx+1)] - L[n] * u[i+(Ny-1+J_y[n])*(Nx+1)];
					
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
			}

			else // if Nx and Ny are too small do one GS_lex step
			{
				for(size_t k=1; k<Nx; k++)
				{
					precision temp=0;

					L = stencil.get_L_c(1,k,Nx,Ny);

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

