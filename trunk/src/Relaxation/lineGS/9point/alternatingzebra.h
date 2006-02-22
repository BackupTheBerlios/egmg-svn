/** \file alternatingzebra.h
 * \author Andre Oeckerath
 * \see lineGS.h
 */
namespace mg
{
	void lineGS::ninepointaltzebra(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> &rhs, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const
					   
    { 
				
		std::valarray<Precision> diagR(Nx-1);
		std::valarray<Precision> ndiagR(Nx-2);
		std::valarray<Precision> ndiagL(Nx-2);

		if(stencil.isConstant() == true )
		{
			const std::valarray<Precision> L = stencil.get_L_c(2,2,Nx,Ny);
			const std::valarray<int> J_x = stencil.getJx(c);
			const std::valarray<int> J_y = stencil.getJy(c);
			
			// begin linewise relaxation in x-direction
			// for each odd line correction of the rhs given by rhs = fv + [1  0  1]^t * u and elimination of the 
			// boundary condition in first and last inner point   
			for(size_t i=1; i<Ny ; i+=2)
			{
				rhs[1+i*(Nx+1)] = fv[1+i*(Nx+1)] - L[2] * u[1+(i+J_y[2])*(Nx+1)] - L[4] * u[1+(i+J_y[4])*(Nx+1)] 
					             - L[1] * u[i*(Nx+1)];
				
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[1+i*(Nx+1)] -= L[sum] * u[1+J_x[sum]+(i+J_y[sum])*(Nx+1)];
				}
				
				for(size_t j=2; j<Nx-1; j++)
				{
					rhs[j+i*(Nx+1)] = fv[j+i*(Nx+1)] - L[2] * u[j+(i+J_y[2])*(Nx+1)] - L[4] * u[j+(i+J_y[4])*(Nx+1)];
					
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[j+i*(Nx+1)] -= L[sum] * u[j+J_x[sum]+(i+J_y[sum])*(Nx+1)];
					}
				}
				rhs[(Nx-1)+i*(Nx+1)] = fv[Nx-1+i*(Nx+1)] - L[2] * u[Nx-1+(i+J_y[2])*(Nx+1)]	
					                  - L[4] * u[Nx-1+(i+J_y[4])*(Nx+1)] - L[3] * u[Nx+i*(Nx+1)];

				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[Nx-1+i*(Nx+1)] -= L[sum] * u[Nx-1+J_x[sum]+(i+J_y[sum])*(Nx+1)];
				}

				// set tridiagonalmatrix for solving A u = rhs
				diagR = L[0];
				ndiagR = L[3];
				ndiagL = L[1];

				
				// LR-decomposition + transformation of the rhs
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

			//even lines
			for(size_t i=2; i<Ny ; i+=2) 
			{
				rhs[1+i*(Nx+1)] = fv[1+i*(Nx+1)] - L[2] * u[1+(i+J_y[2])*(Nx+1)] - L[4] * u[1+(i+J_y[4])*(Nx+1)]
					             - L[1] * u[i*(Nx+1)];

				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[1+i*(Nx+1)] -= L[sum] * u[1+J_x[sum]+(i+J_y[sum])*(Nx+1)];
				}
				
				for(size_t j=2; j<Nx-1; j++)  
				{
					rhs[j+i*(Nx+1)] = fv[j+i*(Nx+1)] - L[2] * u[j+(i+J_y[2])*(Nx+1)] - L[4] * u[j+(i+J_y[4])*(Nx+1)]; 
					
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[j+i*(Nx+1)] -= L[sum] * u[j+J_x[sum]+(i+J_y[sum])*(Nx+1)];
					}
				}

				rhs[(Nx-1)+i*(Nx+1)] = fv[Nx-1+i*(Nx+1)] - L[2] * u[Nx-1+(i+J_y[2])*(Nx+1)] 
					- L[4] * u[Nx-1+(i+J_y[4])*(Nx+1)] - L[3] * u[Nx+i*(Nx+1)];

				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[Nx-1+i*(Nx+1)] -= L[sum] * u[Nx-1+J_x[sum]+(i+J_y[sum])*(Nx+1)];
				}

				diagR = L[0];
				ndiagR = L[3];
				ndiagL = L[1];

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

			
			diagR.resize(Ny-1);
			ndiagR.resize(Ny-2);
			ndiagL.resize(Ny-2);
			
			//begin linewise relaxation in y-direction
			for(size_t i=1; i<Nx ; i+=2) // i = Spalte
			{
				rhs[i+(Nx+1)] = fv[i+(Nx+1)] - L[1] * u[i+J_x[1]+(Nx+1)] - L[3] * u[i+J_x[3]+(Nx+1)] - L[4] * u[i+(1+J_y[4])*(Nx+1)];
				
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[i+Nx+1] -= L[sum] * u[i+J_x[sum]+(1+J_y[sum])*(Nx+1)];
				}

				for(size_t j=2; j<Ny-1; j++)  // j = Laufindex in einer Spalte
				{
					rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[1] * u[i+J_x[1]+j*(Nx+1)] - L[3] * u[i+J_x[3]+j*(Nx+1)];
					
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					}
				}

				rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L[1] * u[i+J_x[1]+(Ny-1)*(Nx+1)] - L[3] * u[i+J_x[3]+(Ny-1)*(Nx+1)] 
					- L[2] * u[i+(Ny-1+J_y[2])*(Nx+1)];
				
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[i+(Ny-1)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
				}


				// set tridiagonalmatrix for solving A u = rhs
				diagR = L[0];
				ndiagR = L[2];
				ndiagL = L[4];
				
				// LR-decomposition + transformation of the rhs
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
				rhs[i+(Nx+1)] = fv[i+(Nx+1)] - L[1] * u[i+J_x[1]+(Nx+1)] - L[3] * u[i+J_x[3]+(Nx+1)] 
					- L[4] * u[i+(1+J_y[4])*(Nx+1)]; 
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[i+Nx+1] -= L[sum] * u[i+J_x[sum]+(1+J_y[sum])*(Nx+1)];
				}

				for(size_t j=2; j<Ny-1; j++) 
				{
					rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[1] * u[i+J_x[1]+j*(Nx+1)] - L[3] * u[i+J_x[3]+j*(Nx+1)];

					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					}
				}
				
				rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L[1] * u[i+J_x[1]+(Ny-1)*(Nx+1)] 
					- L[3] * u[i+J_x[3]+(Ny-1)*(Nx+1)] - L[2] * u[i+(Ny-1+J_y[2])*(Nx+1)];

				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[i+(Ny-1)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
				}

				diagR = L[0];
				ndiagR = L[2];
				ndiagL = L[4];
				
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
			std::valarray<Precision> L = stencil.get_L_c(2,2,Nx,Ny);
			std::valarray<int> J_x = stencil.getJx(c);
			std::valarray<int> J_y = stencil.getJy(c);

			// begin linewise relaxation in x-direction
			// for each odd line correction of the rhs given by rhs = fv + [1  0  1]^t * u and elimination of the 
			// boundary condition in first and last inner point   
			
			diagR.resize(Nx-1);
			ndiagR.resize(Nx-2);
			ndiagL.resize(Nx-2);
			
			if( (Ny > 2) && (Nx > 2))
			{
				L = stencil.get_L_sw(1,1,Nx,Ny);
					
				diagR[0] = L[c];
				ndiagR[0] = L[e];

				rhs[1+Nx+1] = fv[1+Nx+1] - L[2] * u[1+(1+J_y[2])*(Nx+1)] - L[4] * u[1+(1+J_y[4])*(Nx+1)] 
						- L[1] * u[Nx+1];
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[1+Nx+1] -= L[sum] * u[1+J_x[sum]+(1+J_y[sum])*(Nx+1)];
				}

				for(size_t j=2; j<Nx-1; j++)
				{
					L = stencil.get_L_s(j,1,Nx,Ny);
					diagR[j-1] = L[c];
					ndiagR[j-1] = L[e];
					ndiagL[j-2] = L[w];
					rhs[j+Nx+1] = fv[j+Nx+1] - L[2] * u[j+(1+J_y[2])*(Nx+1)] - L[4] * u[j+(1+J_y[4])*(Nx+1)];
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[j+Nx+1] -= L[sum] * u[j+J_x[sum]+(1+J_y[sum])*(Nx+1)];
					}
				}

				L = stencil.get_L_se(Nx-1,1,Nx,Ny);
				
				diagR[Nx-2] = L[c];
				ndiagL[Nx-3] = L[w];
					
				rhs[(Nx-1)+(Nx+1)] = fv[Nx-1+(Nx+1)] - L[2] * u[Nx-1+(1+J_y[2])*(Nx+1)] - L[4] * u[Nx-1+(1+J_y[4])*(Nx+1)] 
						- L[3] * u[Nx+(Nx+1)];
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[Nx-1+(Nx+1)] -= L[sum] * u[Nx-1+J_x[sum]+(1+J_y[sum])*(Nx+1)];
				}

				// LR-decomposition + transformation of the rhs
				for(size_t k=1; k<Nx-1; k++)  
				{
					ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
					diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
					rhs[(Nx+1) + 1 + k] = rhs[(Nx+1) + 1 + k] - ndiagL[k-1] * rhs[(Nx+1)+k];  
				}

				// solve the linear system of equations R u = rhs
				u[(Nx+1)+(Nx-1)] = rhs[(Nx+1)+(Nx-1)] / diagR[Nx-2];

				for(size_t j=Nx-2; j>0; j--)
				{
					u[(Nx+1)+j] = 1/diagR[j-1] * ( rhs[(Nx+1)+j] - ndiagR[j-1] * u[(Nx+1)+j+1] );
				}
				
								
				for(size_t i=3; i<Ny-1 ; i+=2)
				{
					L = stencil.get_L_w(1,i,Nx,Ny);
					
					diagR[0] = L[0];
					ndiagR[0] = L[3];

					rhs[1+i*(Nx+1)] = fv[1+i*(Nx+1)] - L[2] * u[1+(i+J_y[2])*(Nx+1)] - L[4] * u[1+(i+J_y[4])*(Nx+1)]
						- L[1] * u[i*(Nx+1)];
					
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[1+i*(Nx+1)] -= L[sum] * u[1+J_x[sum]+(i+J_y[sum])*(Nx+1)];
					}

					// L im Zentrum im Punkt (j/i)
					
					for(size_t j=2; j<Nx-1; j++)
					{
						L = stencil.get_L_c(j,i,Nx,Ny);
						diagR[j-1] = L[0];
						ndiagR[j-1] = L[3];
						ndiagL[j-2] = L[1];
						rhs[j+i*(Nx+1)] = fv[j+i*(Nx+1)] - L[2] * u[j+(i+J_y[2])*(Nx+1)] - L[4] * u[j+(i+J_y[4])*(Nx+1)]; 
						
						for(size_t sum=5; sum<L.size(); sum++)
						{
							rhs[j+i*(Nx+1)] -= L[sum] * u[j+J_x[sum]+(i+J_y[sum])*(Nx+1)];
						}
					}

					// L am rechten Rand in Zeile i = Punkt (Nx-1/i)
					L = stencil.get_L_e(Nx-1,i,Nx,Ny);
					
					diagR[Nx-2] = L[0];
					ndiagL[Nx-3] = L[1];
					
					rhs[(Nx-1)+i*(Nx+1)] = fv[Nx-1+i*(Nx+1)] - L[2] * u[Nx-1+(i+J_y[2])*(Nx+1)] 
						- L[4] * u[Nx-1+(i+J_y[4])*(Nx+1)] - L[3] * u[Nx+i*(Nx+1)];
					
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[Nx-1+i*(Nx+1)] -= L[sum] * u[Nx-1+J_x[sum]+(i+J_y[sum])*(Nx+1)];
					}

					// LR-decomposition + transformation of the rhs
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

				L = stencil.get_L_nw(1,Ny-1,Nx,Ny);
					
				diagR[0] = L[0];
				ndiagR[0] = L[3];

				rhs[1+(Ny-1)*(Nx+1)] = fv[1+(Ny-1)*(Nx+1)] - L[2] * u[1+(Ny-1+J_y[2])*(Nx+1)] - L[4] * u[1+(Ny-1+J_y[4])*(Nx+1)] 
						- L[1] * u[(Ny-1)*(Nx+1)];
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[1+(Ny-1)*(Nx+1)] -= L[sum] * u[1+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
				}

				for(size_t j=2; j<Nx-1; j++)
				{
					L = stencil.get_L_n(j,Ny-1,Nx,Ny);
					diagR[j-1] = L[0];
					ndiagR[j-1] = L[3];
					ndiagL[j-2] = L[1];
					rhs[j+(Ny-1)*(Nx+1)] = fv[j+(Ny-1)*(Nx+1)] - L[2] * u[j+(Ny-1+J_y[2])*(Nx+1)] - L[4] * u[j+(Ny-1+J_y[4])*(Nx+1)];
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[j+(Ny-1)*(Nx+1)] -= L[sum] * u[j+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
					}
				}

				L = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
					
				diagR[Nx-2] = L[0];
				ndiagL[Nx-3] = L[1];
				
				rhs[(Nx-1)+(Ny-1)*(Nx+1)] = fv[Nx-1+(Ny-1)*(Nx+1)] - L[2] * u[Nx-1+(Ny-1+J_y[2])*(Nx+1)] - L[4] * u[Nx-1+(Ny-1+J_y[4])*(Nx+1)] 
				 	- L[3] * u[Nx+(Ny-1)*(Nx+1)];
				for(size_t sum=5; sum<L.size(); sum++)
				{
					rhs[Nx-1+(Ny-1)*(Nx+1)] -= L[sum] * u[Nx-1+J_x[sum]+(Ny-1+J_y[sum])*(Nx+1)];
				}


				// LR-decomposition + transformation of the rhs
				for(size_t k=1; k<Nx-1; k++)  
				{
					ndiagL[k-1] = ndiagL[k-1]/diagR[k-1];  
					diagR[k] -= ndiagL[k-1] * ndiagR[k-1]; 
					rhs[(Ny-1)*(Nx+1) + 1 + k] = rhs[(Ny-1)*(Nx+1) + 1 + k] - ndiagL[k-1] * rhs[(Ny-1)*(Nx+1)+k];  
				}

				// solve the linear system of equations R u = rhs
				u[(Ny-1)*(Nx+1)+(Nx-1)] = rhs[(Ny-1)*(Nx+1)+(Nx-1)] / diagR[Nx-2];

				for(size_t j=Nx-2; j>0; j--)
				{
					u[(Ny-1)*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[(Ny-1)*(Nx+1)+j] - ndiagR[j-1] * u[(Ny-1)*(Nx+1)+j+1] );
				}


				//even lines
				for(size_t i=2; i<Ny ; i+=2) 
				{					
					L = stencil.get_L_w(1,i,Nx,Ny);
					
					diagR[0] = L[0];
					ndiagR[0] = L[3];

					rhs[1+i*(Nx+1)] = fv[1+i*(Nx+1)] - L[2] * u[1+(i+J_y[2])*(Nx+1)] 
						- L[4] * u[1+(i+J_y[4])*(Nx+1)] - L[1] * u[i*(Nx+1)];
					
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[1+i*(Nx+1)] -= L[sum] * u[1+J_x[sum]+(i+J_y[sum])*(Nx+1)];
					}
					
					for(size_t j=2; j<Nx-1; j++)  
					{
						L = stencil.get_L_c(j,i,Nx,Ny);
						diagR[j-1] = L[0];
						ndiagR[j-1] = L[3];
						ndiagL[j-2] = L[1];
						rhs[j+i*(Nx+1)] = fv[j+i*(Nx+1)] - L[2] * u[j+(i+J_y[2])*(Nx+1)] - L[4] * u[j+(i+J_y[4])*(Nx+1)];
						
						for(size_t sum=5; sum<L.size(); sum++)
						{
							rhs[j+i*(Nx+1)] -= L[sum] * u[j+J_x[sum]+(i+J_y[sum])*(Nx+1)];
						}
					}

					L = stencil.get_L_e(Nx-1,i,Nx,Ny);
					
					diagR[Nx-2] = L[0];
					ndiagL[Nx-3] = L[1];

					rhs[(Nx-1)+i*(Nx+1)] = fv[Nx-1+i*(Nx+1)] - L[2] * u[Nx-1+(i+J_y[2])*(Nx+1)] 
						- L[4] * u[Nx-1+(i+J_y[4])*(Nx+1)] - L[3] * u[Nx+i*(Nx+1)];
					
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[Nx-1+i*(Nx+1)] -= L[sum] * u[Nx-1+J_x[sum]+(i+J_y[sum])*(Nx+1)];
					}

					// LR-decomposition + transformation of the rhs
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

				//begin linewise relaxation in y-direction
				diagR.resize(Ny-1);
			    ndiagR.resize(Ny-2);
			    ndiagL.resize(Ny-2);

				for(size_t i=1; i<Nx ; i+=2)
				{
					L = stencil.get_L_s(i,1,Nx,Ny);
					
					diagR[0] = L[0];
					ndiagR[0] = L[2];	
					
					rhs[i+(Nx+1)] = fv[i+(Nx+1)] - L[1] * u[i+J_x[1]+(Nx+1)] - L[3] * u[i+J_x[3]+(Nx+1)] 
						- L[4] * u[i+(1+J_y[4])*(Nx+1)];

					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[i+Nx+1] -= L[sum] * u[i+J_x[sum]+(1+J_y[sum])*(Nx+1)];
					}
					
					for(size_t j=2; j<Ny-1; j++) 
					{
						L = stencil.get_L_c(i,j,Nx,Ny);
						diagR[j-1] = L[0];
						ndiagR[j-1] = L[2];
						ndiagL[j-2] = L[4];
						rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[1] * u[i+J_x[1]+j*(Nx+1)] - L[3] * u[i+J_x[3]+j*(Nx+1)];
						
						for(size_t sum=5; sum<L.size(); sum++)
						{
							rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
						}
					}
					
					L = stencil.get_L_n(i,Ny-1,Nx,Ny);
					
					diagR[Nx-2] = L[0];
					ndiagL[Nx-3] = L[4];

					rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L[1] * u[i+J_x[1]+(Ny-1)*(Nx+1)] 
						- L[3] * u[i+J_x[3]+(Ny-1)*(Nx+1)] - L[2] * u[i+(Ny-1+J_y[2])*(Nx+1)];

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

				for(size_t i=2; i<Nx ; i+=2)
				{
					L = stencil.get_L_s(i,1,Nx,Ny);
					
					diagR[0] = L[0];
					ndiagR[0] = L[2];
					
					rhs[i+(Nx+1)] = fv[i+(Nx+1)] - L[1] * u[i+J_x[1]+(Nx+1)] - L[3] * u[i+J_x[3]+(Nx+1)] 
						- L[4] * u[i+(1+J_y[4])*(Nx+1)];
					
					for(size_t sum=5; sum<L.size(); sum++)
					{
						rhs[i+Nx+1] -= L[sum] * u[i+J_x[sum]+(1+J_y[sum])*(Nx+1)];
					}
					
					for(size_t j=2; j<Ny-1; j++) 
					{
						L = stencil.get_L_c(i,j,Nx,Ny);
						diagR[j-1] = L[0];
						ndiagR[j-1] = L[2];
						ndiagL[j-2] = L[4];
						rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[1] * u[i+J_x[1]+j*(Nx+1)] - L[3] * u[i+J_x[3]+j*(Nx+1)];
						
						for(size_t sum=5; sum<L.size(); sum++)
						{
							rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
						}
					}

					L = stencil.get_L_n(i,Ny-1,Nx,Ny);
					
					diagR[Nx-2] = L[0];
					ndiagL[Nx-3] = L[4];

					rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L[1] * u[i+J_x[1]+(Ny-1)*(Nx+1)] 
						- L[3] * u[i+J_x[3]+(Ny-1)*(Nx+1)] - L[2] * u[i+(Ny-1+J_y[2])*(Nx+1)];
					

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
			
			else
			{
				Precision temp=0;

				for(size_t k=1; k<Ny; k++)
				{
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

