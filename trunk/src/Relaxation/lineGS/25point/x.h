namespace mg
{
	void lineGS::xline(std::valarray<precision> &u, const std::valarray<precision> &fv, 
		                    std::valarray<precision> &rhs, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const

	{
		if((Ny > 4) && (Nx > 4))
		{
		  std::valarray<precision> diagR(0.0,Nx-1);
		  std::valarray<precision> ndiagR1(0.0,Nx-2);
		  std::valarray<precision> ndiagL1(0.0,Nx-2);
		  std::valarray<precision> ndiagR2(0.0,Nx-3);
		  std::valarray<precision> ndiagL2(0.0,Nx-3);

			if(stencil.is_constant() == true)
			{
				// get const operator L
				const std::valarray<precision> L = stencil.get_L_c(2,2,Nx,Ny);
				const std::valarray<int> J_x = stencil.get_J_x(c);
				const std::valarray<int> J_y = stencil.get_J_y(c);
				
				std::valarray<precision> L_b = stencil.get_L_s(2,1,Nx,Ny);
				std::valarray<int> J_b_x = stencil.get_J_x(s);
				std::valarray<int> J_b_y = stencil.get_J_y(s);
 
                std::valarray<precision> L_c = stencil.get_L_sw(1,1,Nx,Ny);
				std::valarray<int> J_c_x = stencil.get_J_x(sw);
				std::valarray<int> J_c_y = stencil.get_J_y(sw);

				// setze rechte Seite für Zeile 1					
				diagR[0] = L_c[c];
	            ndiagR1[0] = L_c[e];
	            ndiagR2[0] = L_c[ne];
		            				
				rhs[1+Nx+1] = fv[1+Nx+1] - L_c[n] * u[1+(1+J_c_y[n])*(Nx+1)] - L_c[s] * u[1+(1+J_c_y[s])*(Nx+1)]
     				- L_c[w] * u[Nx+1] - L_c[nw] * u[1+(1+J_c_y[nw])*(Nx+1)];
					
				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[1+Nx+1] -= L_c[sum] * u[1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
				}
				
				ndiagL1[0] = L_b[w];
				diagR[1] = L_b[c];
	            ndiagR1[1] = L_b[e];		            
		        ndiagR2[1] = L_b[se];

                rhs[2+Nx+1] = fv[2+Nx+1] - L_b[n] * u[2+(1+J_b_y[n])*(Nx+1)] - L_b[s] * u[2+(1+J_b_y[s])*(Nx+1)]
						    	-L_b[ne] * u[2+(1+J_b_y[ne])*(Nx+1)] - L_b[nw] * u[2+J_b_x[nw]+Nx+1];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
					rhs[2+Nx+1] -= L_b[sum] * u[2+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
				}

				for(size_t j=3; j<Nx-2; j++)  
				{
					ndiagL2[j-3] = L_b[nw];
					ndiagL1[j-2] = L_b[w];
				    diagR[j-1] = L_b[c];
		            ndiagR1[j-1] = L_b[e];		            
			        ndiagR2[j-1] = L_b[se];
						
					rhs[j+Nx+1] = fv[j+Nx+1] - L_b[n] * u[j+(1+J_b_y[n])*(Nx+1)] - L_b[s] * u[j+(1+J_b_y[s])*(Nx+1)]
							-L_b[ne] * u[j+(1+J_b_y[ne])*(Nx+1)];

					for(size_t sum=8; sum<L_b.size(); sum++)
					{
						rhs[j+Nx+1] -= L_b[sum] * u[j+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
					}
				}
					
				ndiagL2[Nx-5] = L_b[nw];
				ndiagL1[Nx-4] = L_b[w];
				diagR[Nx-3] = L_b[c];
		        ndiagR1[Nx-3] = L_b[e];

				rhs[Nx-2+Nx+1] = fv[Nx-2+Nx+1] - L_b[n] * u[Nx-2+(1+J_b_y[n])*(Nx+1)] - L_b[s] * u[Nx-2+(1+J_b_y[s])*(Nx+1)]
							-L_b[ne] * u[Nx-2+(1+J_b_y[ne])*(Nx+1)] - L_b[se] * u[Nx-2+J_b_x[se]+Nx+1];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
					rhs[Nx-2+Nx+1] -= L_b[sum] * u[Nx-2+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
				}
					
                L_c = stencil.get_L_se(Nx-1,1,Nx,Ny);
				J_c_x = stencil.get_J_x(se);
				J_c_y = stencil.get_J_y(se);

				ndiagL2[Nx-4] = L_c[nw];
				ndiagL1[Nx-3] = L_c[w];
				diagR[Nx-2] = L_c[c];

				rhs[(Nx-1)+(Nx+1)] = fv[Nx-1+Nx+1] - L_c[n] * u[Nx-1+(1+J_c_y[n])*(Nx+1)] 
						- L_c[s] * u[Nx-1+(1+J_c_y[s])*(Nx+1)] - L_c[e] * u[Nx+Nx+1] - L_c[ne] * u[Nx-1+(1+J_c_y[ne])*(Nx+1)];

				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[Nx-1+(Nx+1)] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
				}

				// LR-decomposition + transformation of the rhs
				for(size_t k=1; k<Nx-2; k++)  
				{
					ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
					diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
					ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
					rhs[(Nx+1) + 1 + k] = rhs[(Nx+1) + 1 + k] - ndiagL1[k-1] * rhs[(Nx+1)+k]; 

					ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
					ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
					diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
					rhs[(Nx+1) + 1 + k+1] = rhs[(Nx+1) + 1 + k+1] - ndiagL2[k-1] * rhs[(Nx+1)+k];
				}

				ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
				diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
				rhs[(Nx+1) + 1 + Nx-2] = rhs[(Nx+1) + 1 + Nx-2] - ndiagL1[Nx-2-1] * rhs[(Nx+1)+ Nx-2];

			    // solve the linear system of equations R u = rhs
			    u[(Nx+1)+(Nx-1)] = rhs[(Nx+1)+(Nx-1)] / diagR[Nx-2];
			
			    for(size_t j=Nx-2; j>1; j--)
				{
				   u[(Nx+1)+j] = 1/diagR[j-1] * ( rhs[(Nx+1)+j] - ndiagR1[j-1] * u[(Nx+1)+j+1] );

				   rhs[(Nx+1)+j-1] -= ndiagR2[j-2] * u[(Nx+1)+j+1];
				}
				u[(Nx+1)+1] = 1/diagR[0] * ( rhs[(Nx+1)+1] - ndiagR1[0] * u[(Nx+1)+1+1] );

/////////////////////////////////////////////////////////////////////////////////////////////////////////			
				// durchlaufe alle inneren Zeilen, Zeilenindex i				
				for(size_t i=2; i < Ny-1; i++)
				{
					// setze rechte Seite					
					L_b = stencil.get_L_w(1,i,Nx,Ny);
					J_b_x = stencil.get_J_x(w);
				    J_b_y = stencil.get_J_y(w);

					diagR[0] = L_b[c];
		            ndiagR1[0] = L_b[e];
		            ndiagR2[0] = L_b[ne];
		            				
					rhs[1+i*(Nx+1)] = fv[1+i*(Nx+1)] - L_b[n] * u[1+(i+J_b_y[n])*(Nx+1)] - L_b[s] * u[1+(i+J_b_y[s])*(Nx+1)]
					- L_b[w] * u[i*(Nx+1)] - L_b[nw] * u[1+(i+J_b_y[nw])*(Nx+1)] - L_b[se] * u[1+(i+J_b_y[se])*(Nx+1)];
					
					for(size_t sum=8; sum<L_b.size(); sum++)
					{
						rhs[1+i*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(i+J_b_y[sum])*(Nx+1)];
					}
                    					
					ndiagL1[0] = L[w];
					diagR[1] = L[c];
		            ndiagR1[1] = L[e];		            
			        ndiagR2[1] = L[se];

                    rhs[2+i*(Nx+1)] = fv[2+i*(Nx+1)] - L[n] * u[2+(i+J_y[n])*(Nx+1)] - L[s] * u[2+(i+J_y[s])*(Nx+1)]
						    	-L[ne] * u[2+(i+J_y[ne])*(Nx+1)] - L[sw] * u[2+(i+J_y[sw])*(Nx+1)] - L[nw] * u[2+J_x[nw]+i*(Nx+1)];

					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[2+i*(Nx+1)] -= L[sum] * u[2+J_x[sum]+(i+J_y[sum])*(Nx+1)];
					}

					for(size_t j=3; j<Nx-2; j++)  
					{
						ndiagL2[j-3] = L[nw];
						ndiagL1[j-2] = L[w];
					    diagR[j-1] = L[c];
		                ndiagR1[j-1] = L[e];		            
			            ndiagR2[j-1] = L[se];
						
						rhs[j+i*(Nx+1)] = fv[j+i*(Nx+1)] - L[n] * u[j+(i+J_y[n])*(Nx+1)] - L[s] * u[j+(i+J_y[s])*(Nx+1)]
							-L[ne] * u[j+(i+J_y[ne])*(Nx+1)] - L[sw] * u[j+(i+J_y[sw])*(Nx+1)];

						for(size_t sum=9; sum<L.size(); sum++)
						{
							rhs[j+i*(Nx+1)] -= L[sum] * u[j+J_x[sum]+(i+J_y[sum])*(Nx+1)];
						}
					}
					
					ndiagL2[Nx-5] = L[nw];
					ndiagL1[Nx-4] = L[w];
					diagR[Nx-3] = L[c];
		            ndiagR1[Nx-3] = L[e];

					rhs[Nx-2+i*(Nx+1)] = fv[Nx-2+i*(Nx+1)] - L[n] * u[Nx-2+(i+J_y[n])*(Nx+1)] - L[s] * u[Nx-2+(i+J_y[s])*(Nx+1)]
							-L[ne] * u[Nx-2+(i+J_y[ne])*(Nx+1)] - L[sw] * u[Nx-2+(i+J_y[sw])*(Nx+1)] - L[se] * u[Nx-2+J_x[se]+i*(Nx+1)];

					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[Nx-2+i*(Nx+1)] -= L[sum] * u[Nx-2+J_x[sum]+(i+J_y[sum])*(Nx+1)];
					}
					
                    L_b = stencil.get_L_e(Nx-1,i,Nx,Ny);
					J_b_x = stencil.get_J_x(e);
				    J_b_y = stencil.get_J_y(e);

					ndiagL2[Nx-4] = L_b[nw];
					ndiagL1[Nx-3] = L_b[w];
					diagR[Nx-2] = L_b[c];

					rhs[(Nx-1)+i*(Nx+1)] = fv[Nx-1+i*(Nx+1)] - L_b[n] * u[Nx-1+(i+J_b_y[n])*(Nx+1)] 
						- L_b[s] * u[Nx-1+(i+J_b_y[s])*(Nx+1)] - L_b[e] * u[Nx+i*(Nx+1)] - L_b[ne] * u[Nx-1+(i+J_b_y[ne])*(Nx+1)]
						- L_b[se] * u[Nx-1+(i+J_b_y[se])*(Nx+1)];

					for(size_t sum=9; sum<L_b.size(); sum++)
					{
						rhs[Nx-1+i*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(i+J_b_y[sum])*(Nx+1)];
					}
					
					// LR-decomposition + transformation of the rhs
				    for(size_t k=1; k<Nx-2; k++)  
					{
						ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
						diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
						ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
						rhs[i*(Nx+1) + 1 + k] = rhs[i*(Nx+1) + 1 + k] - ndiagL1[k-1] * rhs[i*(Nx+1)+k]; 

						ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
						ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
						diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
						rhs[i*(Nx+1) + 1 + k+1] = rhs[i*(Nx+1) + 1 + k+1] - ndiagL2[k-1] * rhs[i*(Nx+1)+k];
					}

					ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
					diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
					rhs[i*(Nx+1) + 1 + Nx-2] = rhs[i*(Nx+1) + 1 + Nx-2] - ndiagL1[Nx-2-1] * rhs[i*(Nx+1)+ Nx-2];

                    // solve the linear system of equations R u = rhs
				    u[i*(Nx+1)+(Nx-1)] = rhs[i*(Nx+1)+(Nx-1)] / diagR[Nx-2];
				
				    for(size_t j=Nx-2; j>1; j--)
					{
					   u[i*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[i*(Nx+1)+j] - ndiagR1[j-1] * u[i*(Nx+1)+j+1] );

					   rhs[i*(Nx+1)+j-1] -= ndiagR2[j-2] * u[i*(Nx+1)+j+1];
					}
					u[i*(Nx+1)+1] = 1/diagR[0] * ( rhs[i*(Nx+1)+1] - ndiagR1[0] * u[i*(Nx+1)+1+1] );
                }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// setze rechte Seite in oberster Zeile					
				L_c = stencil.get_L_nw(1,Ny-1,Nx,Ny);
				J_c_x = stencil.get_J_x(nw);
				J_c_y = stencil.get_J_y(nw);

				diagR[0] = L_c[c];
		        ndiagR1[0] = L_c[e];
		        ndiagR2[0] = L_c[nw];
		            				
				rhs[1+(Ny-1)*(Nx+1)] = fv[1+(Ny-1)*(Nx+1)] - L_c[n] * u[1+(Ny-1+J_c_y[n])*(Nx+1)] - L_c[s] * u[1+(Ny-1+J_c_y[s])*(Nx+1)]
				    - L_c[w] * u[(Ny-1)*(Nx+1)]  - L_c[ne] * u[1+(Ny-1+J_c_y[ne])*(Nx+1)];
					
				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[1+(Ny-1)*(Nx+1)] -= L_c[sum] * u[1+J_c_x[sum]+(Ny-1+J_c_y[sum])*(Nx+1)];
				}
				
				L_b = stencil.get_L_n(2,Ny-1,Nx,Ny);
				J_b_x = stencil.get_J_x(n);
				J_b_y = stencil.get_J_y(n);
				
                ndiagL1[0] = L_b[w];
				diagR[1] = L_b[c];
		        ndiagR1[1] = L_b[e];		            
			    ndiagR2[1] = L_b[ne];

                rhs[2+(Ny-1)*(Nx+1)] = fv[2+(Ny-1)*(Nx+1)] - L_b[n] * u[2+(Ny-1+J_b_y[n])*(Nx+1)] - L_b[s] * u[2+(Ny-1+J_b_y[s])*(Nx+1)]
				    	- L_b[se] * u[2+(Ny-1+J_b_y[se])*(Nx+1)] - L_b[nw] * u[2+J_b_x[nw]+(Ny-1)*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
					rhs[2+(Ny-1)*(Nx+1)] -= L_b[sum] * u[2+J_b_x[sum]+(Ny-1+J_b_y[sum])*(Nx+1)];
				}

				for(size_t j=3; j<Nx-2; j++)  
				{
					ndiagL2[j-3] = L_b[nw];
					ndiagL1[j-2] = L_b[w];
				    diagR[j-1] = L_b[c];
		            ndiagR1[j-1] = L_b[e];		            
			        ndiagR2[j-1] = L_b[ne];
						
					rhs[j+(Ny-1)*(Nx+1)] = fv[j+(Ny-1)*(Nx+1)] - L_b[n] * u[j+(Ny-1+J_b_y[n])*(Nx+1)] - L_b[s] * u[j+(Ny-1+J_b_y[s])*(Nx+1)]
				    		-L_b[se] * u[j+(Ny-1+J_b_y[se])*(Nx+1)];

					for(size_t sum=8; sum<L_b.size(); sum++)
					{
						rhs[j+(Ny-1)*(Nx+1)] -= L_b[sum] * u[j+J_b_x[sum]+(Ny-1+J_b_y[sum])*(Nx+1)];
					}
				}
					
				ndiagL2[Nx-5] = L_b[nw];
				ndiagL1[Nx-4] = L_b[w];
				diagR[Nx-3] = L_b[c];
		        ndiagR1[Nx-3] = L_b[e];

				rhs[Nx-2+(Ny-1)*(Nx+1)] = fv[Nx-2+(Ny-1)*(Nx+1)] - L_b[n] * u[Nx-2+(Ny-1+J_b_y[n])*(Nx+1)] - L_b[s] * u[Nx-2+(Ny-1+J_b_y[s])*(Nx+1)]
						-L_b[se] * u[Nx-2+(Ny-1+J_b_y[se])*(Nx+1)] - L_b[ne] * u[Nx-2+J_b_x[ne]+(Ny-1)*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
					rhs[Nx-2+(Ny-1)*(Nx+1)] -= L_b[sum] * u[Nx-2+J_b_x[sum]+(Ny-1+J_b_y[sum])*(Nx+1)];
				}
					
                L_c = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
				J_c_x = stencil.get_J_x(ne);
				J_c_y = stencil.get_J_y(ne);

				ndiagL2[Nx-4] = L_c[nw];
				ndiagL1[Nx-3] = L_c[w];
				diagR[Nx-2] = L_c[c];

				rhs[(Nx-1)+(Ny-1)*(Nx+1)] = fv[Nx-1+(Ny-1)*(Nx+1)] - L_c[n] * u[Nx-1+(Ny-1+J_c_y[n])*(Nx+1)]
  						- L_c[s] * u[Nx-1+(Ny-1+J_c_y[s])*(Nx+1)] - L_c[e] * u[Nx+(Ny-1)*(Nx+1)] - L_c[ne] * u[Nx-1+(Ny-1+J_c_y[ne])*(Nx+1)];

				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[Nx-1+(Ny-1)*(Nx+1)] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(Ny-1+J_c_y[sum])*(Nx+1)];
				}

				// LR-decomposition + transformation of the rhs
				for(size_t k=1; k<Nx-2; k++)  
				{
					ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
					diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
					ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
					rhs[(Ny-1)*(Nx+1) + 1 + k] = rhs[(Ny-1)*(Nx+1) + 1 + k] - ndiagL1[k-1] * rhs[(Ny-1)*(Nx+1)+k]; 

					ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
					ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
					diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
					rhs[(Ny-1)*(Nx+1) + 1 + k+1] = rhs[(Ny-1)*(Nx+1) + 1 + k+1] - ndiagL2[k-1] * rhs[(Ny-1)*(Nx+1)+k];
				}

				ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
				diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
				rhs[(Ny-1)*(Nx+1) + 1 + Nx-2] = rhs[(Ny-1)*(Nx+1) + 1 + Nx-2] - ndiagL1[Nx-2-1] * rhs[(Ny-1)*(Nx+1)+ Nx-2];

                // solve the linear system of equations R u = rhs
				u[(Ny-1)*(Nx+1)+(Nx-1)] = rhs[(Ny-1)*(Nx+1)+(Nx-1)] / diagR[Nx-2];
				
				for(size_t j=Nx-2; j>1; j--)
				{
				   u[(Ny-1)*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[(Ny-1)*(Nx+1)+j] - ndiagR1[j-1] * u[(Ny-1)*(Nx+1)+j+1] );

				   rhs[(Ny-1)*(Nx+1)+j-1] -= ndiagR2[j-2] * u[(Ny-1)*(Nx+1)+j+1];
				}
				u[(Ny-1)*(Nx+1)+1] = 1/diagR[0] * ( rhs[(Ny-1)*(Nx+1)+1] - ndiagR1[0] * u[(Ny-1)*(Nx+1)+1+1] );
			}

			else // not constant
			{
				std::valarray<precision> L = stencil.get_L_c(2,2,Nx,Ny);
				std::valarray<int> J_x = stencil.get_J_x(c);
				std::valarray<int> J_y = stencil.get_J_y(c);
				
				std::valarray<precision> L_b = stencil.get_L_s(2,1,Nx,Ny);
				std::valarray<int> J_b_x = stencil.get_J_x(s);
				std::valarray<int> J_b_y = stencil.get_J_y(s);
 
                std::valarray<precision> L_c = stencil.get_L_sw(1,1,Nx,Ny);
				std::valarray<int> J_c_x = stencil.get_J_x(sw);
				std::valarray<int> J_c_y = stencil.get_J_y(sw);

				
				diagR[0] = L_c[c];
	            ndiagR1[0] = L_c[e];
	            ndiagR2[0] = L_c[ne];
		            				
				rhs[1+Nx+1] = fv[1+Nx+1] - L_c[n] * u[1+(1+J_c_y[n])*(Nx+1)] - L_c[s] * u[1+(1+J_c_y[s])*(Nx+1)]
     				- L_c[w] * u[Nx+1] - L_c[nw] * u[1+(1+J_c_y[nw])*(Nx+1)];
					
				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[1+Nx+1] -= L_c[sum] * u[1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
				}
				
				ndiagL1[0] = L_b[w];
				diagR[1] = L_b[c];
	            ndiagR1[1] = L_b[e];		            
		        ndiagR2[1] = L_b[se];

                rhs[2+Nx+1] = fv[2+Nx+1] - L_b[n] * u[2+(1+J_b_y[n])*(Nx+1)] - L_b[s] * u[2+(1+J_b_y[s])*(Nx+1)]
						    	-L_b[ne] * u[2+(1+J_b_y[ne])*(Nx+1)] - L_b[nw] * u[2+J_b_x[nw]+Nx+1];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
					rhs[2+Nx+1] -= L_b[sum] * u[2+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
				}

				for(size_t j=3; j<Nx-2; j++)  
				{
					
					L_b = stencil.get_L_s(j,1,Nx,Ny);
					
					ndiagL2[j-3] = L_b[nw];
					ndiagL1[j-2] = L_b[w];
				    diagR[j-1] = L_b[c];
		            ndiagR1[j-1] = L_b[e];		            
			        ndiagR2[j-1] = L_b[se];
						
					rhs[j+Nx+1] = fv[j+Nx+1] - L_b[n] * u[j+(1+J_b_y[n])*(Nx+1)] - L_b[s] * u[j+(1+J_b_y[s])*(Nx+1)]
							-L_b[ne] * u[j+(1+J_b_y[ne])*(Nx+1)];

					for(size_t sum=8; sum<L_b.size(); sum++)
					{
						rhs[j+Nx+1] -= L_b[sum] * u[j+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
					}
				}
					
				L_b = stencil.get_L_s(Nx-2,1,Nx,Ny);

				ndiagL2[Nx-5] = L_b[nw];
				ndiagL1[Nx-4] = L_b[w];
				diagR[Nx-3] = L_b[c];
		        ndiagR1[Nx-3] = L_b[e];

				rhs[Nx-2+Nx+1] = fv[Nx-2+Nx+1] - L_b[n] * u[Nx-2+(1+J_b_y[n])*(Nx+1)] - L_b[s] * u[Nx-2+(1+J_b_y[s])*(Nx+1)]
							-L_b[ne] * u[Nx-2+(1+J_b_y[ne])*(Nx+1)] - L_b[se] * u[Nx-2+J_b_x[se]+Nx+1];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
					rhs[Nx-2+Nx+1] -= L_b[sum] * u[Nx-2+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
				}
					
                L_c = stencil.get_L_se(Nx-1,1,Nx,Ny);
				J_c_x = stencil.get_J_x(se);
				J_c_y = stencil.get_J_y(se);

				ndiagL2[Nx-4] = L_c[nw];
				ndiagL1[Nx-3] = L_c[w];
				diagR[Nx-2] = L_c[c];

				rhs[(Nx-1)+(Nx+1)] = fv[Nx-1+Nx+1] - L_c[n] * u[Nx-1+(1+J_c_y[n])*(Nx+1)] 
						- L_c[s] * u[Nx-1+(1+J_c_y[s])*(Nx+1)] - L_c[e] * u[Nx+Nx+1] - L_c[ne] * u[Nx-1+(1+J_c_y[ne])*(Nx+1)];

				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[Nx-1+(Nx+1)] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
				}

				// LR-decomposition + transformation of the rhs
				for(size_t k=1; k<Nx-2; k++)  
				{
					ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
					diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
					ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
					rhs[(Nx+1) + 1 + k] = rhs[(Nx+1) + 1 + k] - ndiagL1[k-1] * rhs[(Nx+1)+k]; 

					ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
					ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
					diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
					rhs[(Nx+1) + 1 + k+1] = rhs[(Nx+1) + 1 + k+1] - ndiagL2[k-1] * rhs[(Nx+1)+k];
				}

				ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
				diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
				rhs[(Nx+1) + 1 + Nx-2] = rhs[(Nx+1) + 1 + Nx-2] - ndiagL1[Nx-2-1] * rhs[(Nx+1)+ Nx-2];

			    // solve the linear system of equations R u = rhs
			    u[(Nx+1)+(Nx-1)] = rhs[(Nx+1)+(Nx-1)] / diagR[Nx-2];
			
			    for(size_t j=Nx-2; j>1; j--)
				{
				   u[(Nx+1)+j] = 1/diagR[j-1] * ( rhs[(Nx+1)+j] - ndiagR1[j-1] * u[(Nx+1)+j+1] );

				   rhs[(Nx+1)+j-1] -= ndiagR2[j-2] * u[(Nx+1)+j+1];
				}
				u[(Nx+1)+1] = 1/diagR[0] * ( rhs[(Nx+1)+1] - ndiagR1[0] * u[(Nx+1)+1+1] );

/////////////////////////////////////////////////////////////////////////////////////////////////////////			
				// durchlaufe alle inneren Zeilen, Zeilenindex i				
				for(size_t i=2; i < Ny-1; i++)
				{
					// setze rechte Seite					
					L_b = stencil.get_L_w(1,i,Nx,Ny);
					J_b_x = stencil.get_J_x(w);
				    J_b_y = stencil.get_J_y(w);

					diagR[0] = L_b[c];
		            ndiagR1[0] = L_b[e];
		            ndiagR2[0] = L_b[ne];
		            				
					rhs[1+i*(Nx+1)] = fv[1+i*(Nx+1)] - L_b[n] * u[1+(i+J_b_y[n])*(Nx+1)] - L_b[s] * u[1+(i+J_b_y[s])*(Nx+1)]
					- L_b[w] * u[i*(Nx+1)] - L_b[nw] * u[1+(i+J_b_y[nw])*(Nx+1)] - L_b[se] * u[1+(i+J_b_y[se])*(Nx+1)];
					
					for(size_t sum=8; sum<L_b.size(); sum++)
					{
						rhs[1+i*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(i+J_b_y[sum])*(Nx+1)];
					}
                    					
					L = stencil.get_L_c(2,i,Nx,Ny);
					
					ndiagL1[0] = L[w];
					diagR[1] = L[c];
		            ndiagR1[1] = L[e];		            
			        ndiagR2[1] = L[se];

                    rhs[2+i*(Nx+1)] = fv[2+i*(Nx+1)] - L[n] * u[2+(i+J_y[n])*(Nx+1)] - L[s] * u[2+(i+J_y[s])*(Nx+1)]
						    	-L[ne] * u[2+(i+J_y[ne])*(Nx+1)] - L[sw] * u[2+(i+J_y[sw])*(Nx+1)] - L[nw] * u[2+J_x[nw]+i*(Nx+1)];

					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[2+i*(Nx+1)] -= L[sum] * u[2+J_x[sum]+(i+J_y[sum])*(Nx+1)];
					}

					for(size_t j=3; j<Nx-2; j++)  
					{
						L = stencil.get_L_c(j,i,Nx,Ny);
						
						ndiagL2[j-3] = L[nw];
						ndiagL1[j-2] = L[w];
					    diagR[j-1] = L[c];
		                ndiagR1[j-1] = L[e];		            
			            ndiagR2[j-1] = L[se];
						
						rhs[j+i*(Nx+1)] = fv[j+i*(Nx+1)] - L[n] * u[j+(i+J_y[n])*(Nx+1)] - L[s] * u[j+(i+J_y[s])*(Nx+1)]
							-L[ne] * u[j+(i+J_y[ne])*(Nx+1)] - L[sw] * u[j+(i+J_y[sw])*(Nx+1)];

						for(size_t sum=9; sum<L.size(); sum++)
						{
							rhs[j+i*(Nx+1)] -= L[sum] * u[j+J_x[sum]+(i+J_y[sum])*(Nx+1)];
						}
					}
					
                    L = stencil.get_L_c(Nx-2,i,Nx,Ny);

					ndiagL2[Nx-5] = L[nw];
					ndiagL1[Nx-4] = L[w];
					diagR[Nx-3] = L[c];
		            ndiagR1[Nx-3] = L[e];

					rhs[Nx-2+i*(Nx+1)] = fv[Nx-2+i*(Nx+1)] - L[n] * u[Nx-2+(i+J_y[n])*(Nx+1)] - L[s] * u[Nx-2+(i+J_y[s])*(Nx+1)]
							-L[ne] * u[Nx-2+(i+J_y[ne])*(Nx+1)] - L[sw] * u[Nx-2+(i+J_y[sw])*(Nx+1)] - L[se] * u[Nx-2+J_x[se]+i*(Nx+1)];

					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[Nx-2+i*(Nx+1)] -= L[sum] * u[Nx-2+J_x[sum]+(i+J_y[sum])*(Nx+1)];
					}
					
                    L_b = stencil.get_L_e(Nx-1,i,Nx,Ny);
					J_b_x = stencil.get_J_x(e);
				    J_b_y = stencil.get_J_y(e);

					ndiagL2[Nx-4] = L_b[nw];
					ndiagL1[Nx-3] = L_b[w];
					diagR[Nx-2] = L_b[c];

					rhs[(Nx-1)+i*(Nx+1)] = fv[Nx-1+i*(Nx+1)] - L_b[n] * u[Nx-1+(i+J_b_y[n])*(Nx+1)] 
						- L_b[s] * u[Nx-1+(i+J_b_y[s])*(Nx+1)] - L_b[e] * u[Nx+i*(Nx+1)] - L_b[ne] * u[Nx-1+(i+J_b_y[ne])*(Nx+1)]
						- L_b[se] * u[Nx-1+(i+J_b_y[se])*(Nx+1)];

					for(size_t sum=9; sum<L_b.size(); sum++)
					{
						rhs[Nx-1+i*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(i+J_b_y[sum])*(Nx+1)];
					}
					
					// LR-decomposition + transformation of the rhs
				    for(size_t k=1; k<Nx-2; k++)  
					{
						ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
						diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
						ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
						rhs[i*(Nx+1) + 1 + k] = rhs[i*(Nx+1) + 1 + k] - ndiagL1[k-1] * rhs[i*(Nx+1)+k]; 

						ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
						ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
						diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
						rhs[i*(Nx+1) + 1 + k+1] = rhs[i*(Nx+1) + 1 + k+1] - ndiagL2[k-1] * rhs[i*(Nx+1)+k];
					}

					ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
					diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
					rhs[i*(Nx+1) + 1 + Nx-2] = rhs[i*(Nx+1) + 1 + Nx-2] - ndiagL1[Nx-2-1] * rhs[i*(Nx+1)+ Nx-2];

                    // solve the linear system of equations R u = rhs
				    u[i*(Nx+1)+(Nx-1)] = rhs[i*(Nx+1)+(Nx-1)] / diagR[Nx-2];
				
				    for(size_t j=Nx-2; j>1; j--)
					{
					   u[i*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[i*(Nx+1)+j] - ndiagR1[j-1] * u[i*(Nx+1)+j+1] );

					   rhs[i*(Nx+1)+j-1] -= ndiagR2[j-2] * u[i*(Nx+1)+j+1];
					}
					u[i*(Nx+1)+1] = 1/diagR[0] * ( rhs[i*(Nx+1)+1] - ndiagR1[0] * u[i*(Nx+1)+1+1] );
                }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// setze rechte Seite in oberster Zeile					
				L_c = stencil.get_L_nw(1,Ny-1,Nx,Ny);
				J_c_x = stencil.get_J_x(nw);
				J_c_y = stencil.get_J_y(nw);

				diagR[0] = L_c[c];
		        ndiagR1[0] = L_c[e];
		        ndiagR2[0] = L_c[nw];
		            				
				rhs[1+(Ny-1)*(Nx+1)] = fv[1+(Ny-1)*(Nx+1)] - L_c[n] * u[1+(Ny-1+J_c_y[n])*(Nx+1)] - L_c[s] * u[1+(Ny-1+J_c_y[s])*(Nx+1)]
				    - L_c[w] * u[(Ny-1)*(Nx+1)]  - L_c[ne] * u[1+(Ny-1+J_c_y[ne])*(Nx+1)];
					
				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[1+(Ny-1)*(Nx+1)] -= L_c[sum] * u[1+J_c_x[sum]+(Ny-1+J_c_y[sum])*(Nx+1)];
				}
				
				L_b = stencil.get_L_n(2,Ny-1,Nx,Ny);
				J_b_x = stencil.get_J_x(n);
				J_b_y = stencil.get_J_y(n);
				
                ndiagL1[0] = L_b[w];
				diagR[1] = L_b[c];
		        ndiagR1[1] = L_b[e];		            
			    ndiagR2[1] = L_b[ne];

                rhs[2+(Ny-1)*(Nx+1)] = fv[2+(Ny-1)*(Nx+1)] - L_b[n] * u[2+(Ny-1+J_b_y[n])*(Nx+1)] - L_b[s] * u[2+(Ny-1+J_b_y[s])*(Nx+1)]
				    	- L_b[se] * u[2+(Ny-1+J_b_y[se])*(Nx+1)] - L_b[nw] * u[2+J_b_x[nw]+(Ny-1)*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
					rhs[2+(Ny-1)*(Nx+1)] -= L_b[sum] * u[2+J_b_x[sum]+(Ny-1+J_b_y[sum])*(Nx+1)];
				}

				for(size_t j=3; j<Nx-2; j++)  
				{
					L_b = stencil.get_L_n(j,Ny-1,Nx,Ny);
					
					ndiagL2[j-3] = L_b[nw];
					ndiagL1[j-2] = L_b[w];
				    diagR[j-1] = L_b[c];
		            ndiagR1[j-1] = L_b[e];		            
			        ndiagR2[j-1] = L_b[ne];
						
					rhs[j+(Ny-1)*(Nx+1)] = fv[j+(Ny-1)*(Nx+1)] - L_b[n] * u[j+(Ny-1+J_b_y[n])*(Nx+1)] - L_b[s] * u[j+(Ny-1+J_b_y[s])*(Nx+1)]
				    		-L_b[se] * u[j+(Ny-1+J_b_y[se])*(Nx+1)];

					for(size_t sum=8; sum<L_b.size(); sum++)
					{
						rhs[j+(Ny-1)*(Nx+1)] -= L_b[sum] * u[j+J_b_x[sum]+(Ny-1+J_b_y[sum])*(Nx+1)];
					}
				}
					
				L_b = stencil.get_L_n(Nx-2,Ny-1,Nx,Ny);

				ndiagL2[Nx-5] = L_b[nw];
				ndiagL1[Nx-4] = L_b[w];
				diagR[Nx-3] = L_b[c];
		        ndiagR1[Nx-3] = L_b[e];

				rhs[Nx-2+(Ny-1)*(Nx+1)] = fv[Nx-2+(Ny-1)*(Nx+1)] - L_b[n] * u[Nx-2+(Ny-1+J_b_y[n])*(Nx+1)] - L_b[s] * u[Nx-2+(Ny-1+J_b_y[s])*(Nx+1)]
						-L_b[se] * u[Nx-2+(Ny-1+J_b_y[se])*(Nx+1)] - L_b[ne] * u[Nx-2+J_b_x[ne]+(Ny-1)*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
					rhs[Nx-2+(Ny-1)*(Nx+1)] -= L_b[sum] * u[Nx-2+J_b_x[sum]+(Ny-1+J_b_y[sum])*(Nx+1)];
				}
					
                L_c = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
				J_c_x = stencil.get_J_x(ne);
				J_c_y = stencil.get_J_y(ne);

				ndiagL2[Nx-4] = L_c[nw];
				ndiagL1[Nx-3] = L_c[w];
				diagR[Nx-2] = L_c[c];

				rhs[(Nx-1)+(Ny-1)*(Nx+1)] = fv[Nx-1+(Ny-1)*(Nx+1)] - L_c[n] * u[Nx-1+(Ny-1+J_c_y[n])*(Nx+1)]
  						- L_c[s] * u[Nx-1+(Ny-1+J_c_y[s])*(Nx+1)] - L_c[e] * u[Nx+(Ny-1)*(Nx+1)] - L_c[ne] * u[Nx-1+(Ny-1+J_c_y[ne])*(Nx+1)];

				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[Nx-1+(Ny-1)*(Nx+1)] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(Ny-1+J_c_y[sum])*(Nx+1)];
				}

				// LR-decomposition + transformation of the rhs
				for(size_t k=1; k<Nx-2; k++)  
				{
					ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
					diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
					ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
					rhs[(Ny-1)*(Nx+1) + 1 + k] = rhs[(Ny-1)*(Nx+1) + 1 + k] - ndiagL1[k-1] * rhs[(Ny-1)*(Nx+1)+k]; 

					ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
					ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
					diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
					rhs[(Ny-1)*(Nx+1) + 1 + k+1] = rhs[(Ny-1)*(Nx+1) + 1 + k+1] - ndiagL2[k-1] * rhs[(Ny-1)*(Nx+1)+k];
				}

				ndiagL1[Nx-2-1] = ndiagL1[Nx-2-1]/diagR[Nx-2-1];
				diagR[Nx-2] -= ndiagL1[Nx-2-1] * ndiagR1[Nx-2-1];
				rhs[(Ny-1)*(Nx+1) + 1 + Nx-2] = rhs[(Ny-1)*(Nx+1) + 1 + Nx-2] - ndiagL1[Nx-2-1] * rhs[(Ny-1)*(Nx+1)+ Nx-2];

                // solve the linear system of equations R u = rhs
				u[(Ny-1)*(Nx+1)+(Nx-1)] = rhs[(Ny-1)*(Nx+1)+(Nx-1)] / diagR[Nx-2];
				
				for(size_t j=Nx-2; j>1; j--)
				{
				   u[(Ny-1)*(Nx+1)+j] = 1/diagR[j-1] * ( rhs[(Ny-1)*(Nx+1)+j] - ndiagR1[j-1] * u[(Ny-1)*(Nx+1)+j+1] );

				   rhs[(Ny-1)*(Nx+1)+j-1] -= ndiagR2[j-2] * u[(Ny-1)*(Nx+1)+j+1];
				}
				u[(Ny-1)*(Nx+1)+1] = 1/diagR[0] * ( rhs[(Ny-1)*(Nx+1)+1] - ndiagR1[0] * u[(Ny-1)*(Nx+1)+1+1] );
				
			}

		}

		else 
		{
			
			for(int k=0; k<2; k++)
			{
				precision factor = 1.0;

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
