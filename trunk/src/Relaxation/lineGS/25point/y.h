namespace mg
{
	void ZebraLineGS::yline(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> &rhs, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const

	{
	    if((Ny > 4) && (Nx > 4))
		{
			std::valarray<Precision> diagR(0.0,Ny-1);
		    std::valarray<Precision> ndiagR1(0.0,Ny-2);
		    std::valarray<Precision> ndiagL1(0.0,Ny-2);
			std::valarray<Precision> ndiagR2(0.0,Ny-3);
		    std::valarray<Precision> ndiagL2(0.0,Ny-3);
	                
            if(stencil.isConstant() == true)
			{
				// get const operator L
				const std::valarray<Precision> L = stencil.get_L_c(2,2,Nx,Ny);
				const std::valarray<int> J_x = stencil.getJx(C);
				const std::valarray<int> J_y = stencil.getJy(C);
				
				std::valarray<Precision> L_b = stencil.get_L_w(1,2,Nx,Ny);
				std::valarray<int> J_b_x = stencil.getJx(W);
				std::valarray<int> J_b_y = stencil.getJy(W);
 
                std::valarray<Precision> L_c = stencil.get_L_sw(1,1,Nx,Ny);
				std::valarray<int> J_c_x = stencil.getJx(SW);
				std::valarray<int> J_c_y = stencil.getJy(SW);

                diagR[0] = L_c[C];
	            ndiagR1[0] = L_c[N];
	            ndiagR2[0] = L_c[NW];
                                
                rhs[1+Nx+1] = fv[1+Nx+1] - L_c[S] * u[1+(1+J_c_y[S])*(Nx+1)]
     				- L_c[W] * u[Nx+1] - L_c[E] * u[2+Nx+1] - L_c[NE] * u[1+J_c_x[NE]+Nx+1];
					
				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[1+Nx+1] -= L_c[sum] * u[1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
				}

                ndiagL1[0] = L_b[S];
				diagR[1] = L_b[C];
	            ndiagR1[1] = L_b[N];		            
		        ndiagR2[1] = L_b[NW];

                rhs[1+2*(Nx+1)] = fv[1+2*(Nx+1)] - L_b[W] * u[2*(Nx+1)] - L_b[E] * u[2+2*(Nx+1)]
				  -L_b[NE] * u[1+J_b_x[NE]+2*(Nx+1)] - L_b[SE] * u[1+(2+J_b_y[SE])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				    rhs[1+2*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(2+J_b_y[sum])*(Nx+1)];
				}

				for(size_t j=3; j<Ny-2; j++)  
				{
					ndiagL2[j-3] = L_b[SE];
					ndiagL1[j-2] = L_b[S];
			  	    diagR[j-1] = L_b[C];
  		            ndiagR1[j-1] = L_b[N];		            
			        ndiagR2[j-1] = L_b[NW];
						
					rhs[1+j*(Nx+1)] = fv[1+j*(Nx+1)] - L_b[W] * u[1+J_b_x[W]+j*(Nx+1)] - L_b[E] * u[1+J_b_x[E]+j*(Nx+1)]
							-L_b[NE] * u[1+J_b_x[NE]+j*(Nx+1)];

					for(size_t sum=8; sum<L_b.size(); sum++)
					{
					  rhs[1+j*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(j+J_b_y[sum])*(Nx+1)];
					}
				}
					
				ndiagL2[Ny-5] = L_b[SE];
				ndiagL1[Ny-4] = L_b[S];
				diagR[Ny-3] = L_b[C];
		        ndiagR1[Ny-3] = L_b[N];

				rhs[1+(Ny-2)*(Nx+1)] = fv[1+(Ny-2)*(Nx+1)] - L_b[W] * u[(Ny-2)*(Nx+1)] - L_b[E] * u[2+(Ny-2)*(Nx+1)]
				  -L_b[NE] * u[1+J_b_x[NE]+(Ny-2)*(Nx+1)] - L_b[NW] * u[1+(Ny-2+J_b_y[NW])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				  rhs[1+(Ny-2)*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(Ny-2+J_b_y[sum])*(Nx+1)];
				}
					
                L_c = stencil.get_L_nw(1,Ny-1,Nx,Ny);
				J_c_x = stencil.getJx(NW);
				J_c_y = stencil.getJy(NW);

				ndiagL2[Ny-4] = L_c[NE];
				ndiagL1[Ny-3] = L_c[S];
				diagR[Ny-2] = L_c[C];

				rhs[1+(Ny-1)*(Nx+1)] = fv[1+(Ny-1)*(Nx+1)] - L_c[W] * u[(Ny-1)*(Nx+1)] - L_c[E] * u[2+(Ny-1)*(Nx+1)]
					    - L_c[NW] * u[1+J_c_x[NW]+(Ny-1)*(Nx+1)] - L_c[N] * u[1+(Ny-1+J_c_y[N])*(Nx+1)];

				for(size_t sum=7; sum<L_c.size(); sum++)
				{
				  rhs[1+(Ny-1)*(Nx+1)] -= L_c[sum] * u[1+J_c_x[sum]+(Ny-1+J_c_y[sum])*(Nx+1)];
				}

				// LR-decomposition + transformation of the rhs
				for(size_t k=1; k<Ny-2; k++)  
				{
					ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
					diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
					ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
					rhs[(k+1)*(Nx+1) + 1] = rhs[(k+1)*(Nx+1) + 1] - ndiagL1[k-1] * rhs[k*(Nx+1)+1]; 

					ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
					ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
					diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
					rhs[(k+2)*(Nx+1) + 1] = rhs[(k+2)*(Nx+1) + 1] - ndiagL2[k-1] * rhs[k*(Nx+1)+1];
				}
				ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
				diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
				rhs[(Ny-1)*(Nx+1) + 1] = rhs[(Ny-1)*(Nx+1) + 1] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ 1];

			    // solve the linear system of equations R u = rhs
			    u[1+(Nx+1)*(Ny-1)] = rhs[1+(Nx+1)*(Ny-1)] / diagR[Ny-2];
				for(size_t k=Ny-2; k>1; k--)
				{
					u[1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[1+k*(Nx+1)] - ndiagR1[k-1] * u[1+(k+1)*(Nx+1)] );

				    rhs[(k-1)*(Nx+1)+1] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+1];
				}
				u[(Nx+1)+1] = 1/diagR[0] * ( rhs[(Nx+1)+1] - ndiagR1[0] * u[2*(Nx+1)+1] );



////////////////////////////////////////////////////////////////////////////////////////////////

                // durchlaufe alle inneren Spalten, Spaltenindex i				
				for(size_t i=2; i < Nx-1; i++)
				{
                    // setze rechte Seite					
					L_b = stencil.get_L_s(i,1,Nx,Ny);
					J_b_x = stencil.getJx(S);
				    J_b_y = stencil.getJy(S);

					diagR[0] = L_b[C];
        		    ndiagR1[0] = L_b[N];
		            ndiagR2[0] = L_b[NE];
       
                    rhs[i+Nx+1] = fv[i+Nx+1] - L_b[S] * u[i+(1+J_b_y[S])*(Nx+1)] - L_b[W] * u[i+J_b_x[W]+Nx+1] 
						- L_b[E] * u[i+J_b_x[E]+Nx+1] - L_b[NW] * u[i+J_b_x[NW]+Nx+1] - L_b[SE] * u[i+J_b_x[SE]+(Nx+1)];
					
		        	for(size_t sum=7; sum<L_b.size(); sum++)
			        {
				       	rhs[i+Nx+1] -= L_b[sum] * u[i+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
				    }

                    ndiagL1[0] = L[S];
					diagR[1] = L[C];
		            ndiagR1[1] = L[N];		            
 			        ndiagR2[1] = L[NE];

                    rhs[i+2*(Nx+1)] = fv[i+2*(Nx+1)] - L[W] * u[i+J_x[W]+2*(Nx+1)] - L[E] * u[i+J_x[E]+2*(Nx+1)]
					  -L[NW] * u[i+J_x[NW]+2*(Nx+1)] - L[SE] * u[i+J_x[SE]+2*(Nx+1)] - L[SW] * u[i+(2+J_y[SW])*(Nx+1)];

					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[i+2*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(2+J_y[sum])*(Nx+1)];
					}
    				    

					for(size_t j=3; j<Ny-2; j++)
					{
					   ndiagL2[j-3] = L[SW];
					   ndiagL1[j-2] = L[S];
					   diagR[j-1] = L[C];
					   ndiagR1[j-1] = L[N];
					   ndiagR2[j-1] = L[NE];
 
                       rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[W] * u[i+J_x[W]+j*(Nx+1)] - L[E] * u[i+J_x[E]+j*(Nx+1)]
					         -L[NW] * u[i+J_x[NW]+j*(Nx+1)] - L[SE] * u[i+J_x[SE]+j*(Nx+1)];
					   
					   for(size_t sum=9; sum<L.size(); sum++)
					   {
						   rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					   }
					}
                                        
                    ndiagL2[Nx-5] = L[SW];
					ndiagL1[Nx-4] = L[S];
					diagR[Nx-3] = L[C];
                    ndiagR1[Nx-3] = L[N];

					rhs[i+(Ny-2)*(Nx+1)] = fv[i+(Ny-2)*(Nx+1)] - L[W] * u[i+J_x[W]+(Ny-2)*(Nx+1)] - L[E] * u[i+J_x[E]+(Ny-2)*(Nx+1)]
					     -L[NW] * u[i+J_x[NW]+(Ny-2)*(Nx+1)] - L[SE] * u[i+J_x[SE]+(Ny-2)*(Nx+1)] - L[NE] * u[i+(Ny-2+J_y[NE])*(Nx+1)];
                    
					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[i+(Ny-2)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-2+J_y[sum])*(Nx+1)];
					}


					L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
					J_b_x = stencil.getJx(N);
				    J_b_y = stencil.getJy(N);

					ndiagL2[Nx-4] = L_b[SW];
					ndiagL1[Nx-3] = L_b[S];
					diagR[Nx-2] = L_b[C];

					rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L_b[N] * u[i+(Ny-1+J_b_y[N])*(Nx+1)] 
						- L_b[W] * u[i+J_b_x[W]+(Ny-1)*(Nx+1)] - L_b[E] * u[i+J_b_x[E]+(Ny-1)*(Nx+1)] - L_b[NW] * u[i+J_b_x[NW]+(Ny-1)*(Nx+1)]
						- L_b[SE] * u[i+J_b_x[SE]+(Ny-1)*(Nx+1)];

					for(size_t sum=9; sum<L_b.size(); sum++)
					{
						rhs[i+(Ny-1)*(Nx+1)] -= L_b[sum] * u[i+J_b_x[sum]+(Ny-1+J_b_y[sum])*(Nx+1)];
					}
					
					// LR-decomposition + transformation of the rhs
					for(size_t k=1; k<Ny-2; k++)  
					{
						ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
					    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
					    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
					    rhs[(k+1)*(Nx+1) + i] = rhs[(k+1)*(Nx+1) + i] - ndiagL1[k-1] * rhs[k*(Nx+1)+i]; 

					    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
					    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
					    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
					    rhs[(k+2)*(Nx+1) + i] = rhs[(k+2)*(Nx+1) + i] - ndiagL2[k-1] * rhs[k*(Nx+1)+i];
					}
					ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
				    diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
				    rhs[(Ny-1)*(Nx+1) + i] = rhs[(Ny-1)*(Nx+1) + i] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ i];

					// solve the linear system of equations R u = rhs
					u[i+(Nx+1)*(Ny-1)] = rhs[i+(Nx+1)*(Ny-1)] / diagR[Ny-2];

				    for(size_t k=Ny-2; k>1; k--)
					{
						u[i+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[i+k*(Nx+1)] - ndiagR1[k-1] * u[i+(k+1)*(Nx+1)] );

						rhs[(k-1)*(Nx+1)+i] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+i];
					}
					u[(Nx+1)+i] = 1/diagR[0] * ( rhs[(Nx+1)+i] - ndiagR1[0] * u[2*(Nx+1)+i] );
				}

////////////////letzte Spalte//////////////////
                
                L_c = stencil.get_L_se(Nx-1,1,Nx,Ny);
				J_c_x = stencil.getJx(SE);
				J_c_y = stencil.getJy(SE);

				
                rhs[Nx-1+Nx+1] = fv[Nx-1+Nx+1] - L_c[S] * u[Nx-1+(1+J_c_y[S])*(Nx+1)]
     				- L_c[W] * u[Nx-2+Nx+1] - L_c[E] * u[Nx+Nx+1] - L_c[NW] * u[Nx-1+J_c_x[NW]+Nx+1];
					
				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[Nx-1+Nx+1] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
				}

                L_b = stencil.get_L_e(Nx-1,2,Nx,Ny);
				J_b_x = stencil.getJx(E);
        		J_b_y = stencil.getJy(E);

                ndiagL1[0] = L_b[S];
				diagR[1] = L_b[C];
	            ndiagR1[1] = L_b[N];		            
		        ndiagR2[1] = L_b[NE];

                rhs[Nx-1+2*(Nx+1)] = fv[Nx-1+2*(Nx+1)] - L_b[W] * u[Nx-2+2*(Nx+1)] - L_b[E] * u[Nx+2*(Nx+1)]
				  -L_b[NW] * u[Nx-1+J_b_x[NW]+2*(Nx+1)] - L_b[SE] * u[Nx-1+(2+J_b_y[SE])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				    rhs[Nx-1+2*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(2+J_b_y[sum])*(Nx+1)];
				}

				for(size_t j=3; j<Ny-2; j++)  
				{
					ndiagL2[j-3] = L_b[SE];
					ndiagL1[j-2] = L_b[S];
			  	    diagR[j-1] = L_b[C];
  		            ndiagR1[j-1] = L_b[N];		            
			        ndiagR2[j-1] = L_b[NE];
						
					rhs[Nx-1+j*(Nx+1)] = fv[Nx-1+j*(Nx+1)] - L_b[W] * u[Nx-1+J_b_x[W]+j*(Nx+1)] - L_b[E] * u[Nx-1+J_b_x[E]+j*(Nx+1)]
							-L_b[NW] * u[Nx-1+J_b_x[NW]+j*(Nx+1)];

					for(size_t sum=8; sum<L_b.size(); sum++)
					{
						rhs[Nx-1+j*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(j+J_b_y[sum])*(Nx+1)];
					}
				}
					
				ndiagL2[Ny-5] = L_b[SE];
				ndiagL1[Ny-4] = L_b[S];
				diagR[Ny-3] = L_b[C];
		        ndiagR1[Ny-3] = L_b[N];

				rhs[Nx-1+(Ny-2)*(Nx+1)] = fv[Nx-1+(Ny-2)*(Nx+1)] - L_b[W] * u[Nx-2+(Ny-2)*(Nx+1)] - L_b[E] * u[Nx+(Ny-2)*(Nx+1)]
				  - L_b[NW] * u[Nx-1+J_b_x[NW]+(Ny-2)*(Nx+1)] - L_b[NE] * u[Nx-1+(Ny-2+J_b_y[NE])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				  rhs[Nx-1+(Ny-2)*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(Ny-2+J_b_y[sum])*(Nx+1)];
				}
					
                L_c = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
				J_c_x = stencil.getJx(NE);
				J_c_y = stencil.getJy(NE);

				ndiagL2[Ny-4] = L_c[NE];
				ndiagL1[Ny-3] = L_c[S];
				diagR[Ny-2] = L_c[C];

				rhs[Nx-1+(Ny-1)*(Nx+1)] = fv[Nx-1+(Ny-1)*(Nx+1)] - L_c[W] * u[Nx-2+(Ny-1)*(Nx+1)] 
				  - L_c[E] * u[Nx+(Ny-1)*(Nx+1)] - L_c[NW] * u[Nx-1+J_c_x[NW]+(Ny-1)*(Nx+1)] - L_c[N] * u[Nx-1+(Ny-1+J_c_y[N])*(Nx+1)];

				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[Nx-1+(Ny-1)*(Nx+1)] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(Ny-1+J_c_y[sum])*(Nx+1)];
				}

				// LR-decomposition + transformation of the rhs
				for(size_t k=1; k<Ny-2; k++)  
				{
					ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
					diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
					ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
					rhs[(k+1)*(Nx+1) + Nx-1] = rhs[(k+1)*(Nx+1) + Nx-1] - ndiagL1[k-1] * rhs[k*(Nx+1)+Nx-1]; 

					ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
					ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
					diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
					rhs[(k+2)*(Nx+1) + Nx-1] = rhs[(k+2)*(Nx+1) + Nx-1] - ndiagL2[k-1] * rhs[k*(Nx+1)+Nx-1];
				}
				ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
				diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
				rhs[(Ny-1)*(Nx+1) + Nx-1] = rhs[(Ny-1)*(Nx+1) + Nx-1] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ Nx-1];

			    // solve the linear system of equations R u = rhs
			    u[Nx-1+(Nx+1)*(Ny-1)] = rhs[Nx-1+(Nx+1)*(Ny-1)] / diagR[Ny-2];

				for(size_t k=Ny-2; k>1; k--)
				{
					u[Nx-1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[Nx-1+k*(Nx+1)] - ndiagR1[k-1] * u[Nx-1+(k+1)*(Nx+1)] );

				    rhs[(k-1)*(Nx+1)+Nx-1] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+Nx-1];
				}
				u[(Nx+1)+Nx-1] = 1/diagR[0] * ( rhs[(Nx+1)+Nx-1] - ndiagR1[0] * u[2*(Nx+1)+Nx-1] );

			}
			else // stencil not constant
			{
				std::valarray<Precision> L = stencil.get_L_c(2,2,Nx,Ny);
				std::valarray<int> J_x = stencil.getJx(C);
				std::valarray<int> J_y = stencil.getJy(C);
				
				std::valarray<Precision> L_b = stencil.get_L_w(1,2,Nx,Ny);
				std::valarray<int> J_b_x = stencil.getJx(W);
				std::valarray<int> J_b_y = stencil.getJy(W);
 
                std::valarray<Precision> L_c = stencil.get_L_sw(1,1,Nx,Ny);
				std::valarray<int> J_c_x = stencil.getJx(SW);
				std::valarray<int> J_c_y = stencil.getJy(SW);

                diagR[0] = L_c[C];
	            ndiagR1[0] = L_c[N];
	            ndiagR2[0] = L_c[NW];
                                
                rhs[1+Nx+1] = fv[1+Nx+1] - L_c[S] * u[1+(1+J_c_y[S])*(Nx+1)]
     				- L_c[W] * u[Nx+1] - L_c[E] * u[2+Nx+1] - L_c[NE] * u[1+J_c_x[NE]+Nx+1];
					
				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[1+Nx+1] -= L_c[sum] * u[1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
				}

                
                ndiagL1[0] = L_b[S];
				diagR[1] = L_b[C];
	            ndiagR1[1] = L_b[N];		            
		        ndiagR2[1] = L_b[NW];

                rhs[1+2*(Nx+1)] = fv[1+2*(Nx+1)] - L_b[W] * u[2*(Nx+1)] - L_b[E] * u[2+2*(Nx+1)]
				  -L_b[NE] * u[1+J_b_x[NE]+2*(Nx+1)] - L_b[SE] * u[1+(2+J_b_y[SE])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				    rhs[1+2*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(2+J_b_y[sum])*(Nx+1)];
				}

				for(size_t j=3; j<Ny-2; j++)  
				{
					L_b = stencil.get_L_w(1,j,Nx,Ny);

					ndiagL2[j-3] = L_b[SE];
					ndiagL1[j-2] = L_b[S];
			  	    diagR[j-1] = L_b[C];
  		            ndiagR1[j-1] = L_b[N];		            
			        ndiagR2[j-1] = L_b[NW];
						
					rhs[1+j*(Nx+1)] = fv[1+j*(Nx+1)] - L_b[W] * u[1+J_b_x[W]+j*(Nx+1)] - L_b[E] * u[1+J_b_x[E]+j*(Nx+1)]
							-L_b[NE] * u[1+J_b_x[NE]+j*(Nx+1)];

					for(size_t sum=8; sum<L_b.size(); sum++)
					{
					  rhs[1+j*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(j+J_b_y[sum])*(Nx+1)];
					}
				}

				L_b = stencil.get_L_w(1,Ny-2,Nx,Ny);
				
				ndiagL2[Ny-5] = L_b[SE];
				ndiagL1[Ny-4] = L_b[S];
				diagR[Ny-3] = L_b[C];
		        ndiagR1[Ny-3] = L_b[N];

				rhs[1+(Ny-2)*(Nx+1)] = fv[1+(Ny-2)*(Nx+1)] - L_b[W] * u[(Ny-2)*(Nx+1)] - L_b[E] * u[2+(Ny-2)*(Nx+1)]
				  -L_b[NE] * u[1+J_b_x[NE]+(Ny-2)*(Nx+1)] - L_b[NW] * u[1+(Ny-2+J_b_y[NW])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				  rhs[1+(Ny-2)*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(Ny-2+J_b_y[sum])*(Nx+1)];
				}
					
                L_c = stencil.get_L_nw(1,Ny-1,Nx,Ny);
				J_c_x = stencil.getJx(NW);
				J_c_y = stencil.getJy(NW);

				ndiagL2[Ny-4] = L_c[NE];
				ndiagL1[Ny-3] = L_c[S];
				diagR[Ny-2] = L_c[C];

				rhs[1+(Ny-1)*(Nx+1)] = fv[1+(Ny-1)*(Nx+1)] - L_c[W] * u[(Ny-1)*(Nx+1)] - L_c[E] * u[2+(Ny-1)*(Nx+1)]
					    - L_c[NW] * u[1+J_c_x[NW]+(Ny-1)*(Nx+1)] - L_c[N] * u[1+(Ny-1+J_c_y[N])*(Nx+1)];

				for(size_t sum=7; sum<L_c.size(); sum++)
				{
				  rhs[1+(Ny-1)*(Nx+1)] -= L_c[sum] * u[1+J_c_x[sum]+(Ny-1+J_c_y[sum])*(Nx+1)];
				}

				// LR-decomposition + transformation of the rhs
				for(size_t k=1; k<Ny-2; k++)  
				{
					ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
					diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
					ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
					rhs[(k+1)*(Nx+1) + 1] = rhs[(k+1)*(Nx+1) + 1] - ndiagL1[k-1] * rhs[k*(Nx+1)+1]; 

					ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
					ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
					diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
					rhs[(k+2)*(Nx+1) + 1] = rhs[(k+2)*(Nx+1) + 1] - ndiagL2[k-1] * rhs[k*(Nx+1)+1];
				}
				ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
				diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
				rhs[(Ny-1)*(Nx+1) + 1] = rhs[(Ny-1)*(Nx+1) + 1] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ 1];

			    // solve the linear system of equations R u = rhs
			    u[1+(Nx+1)*(Ny-1)] = rhs[1+(Nx+1)*(Ny-1)] / diagR[Ny-2];
				for(size_t k=Ny-2; k>1; k--)
				{
					u[1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[1+k*(Nx+1)] - ndiagR1[k-1] * u[1+(k+1)*(Nx+1)] );

				    rhs[(k-1)*(Nx+1)+1] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+1];
				}
				u[(Nx+1)+1] = 1/diagR[0] * ( rhs[(Nx+1)+1] - ndiagR1[0] * u[2*(Nx+1)+1] );



////////////////////////////////////////////////////////////////////////////////////////////////

                // durchlaufe alle inneren Spalten, Spaltenindex i				
				for(size_t i=2; i < Nx-1; i++)
				{
                    // setze rechte Seite					
					L_b = stencil.get_L_s(i,1,Nx,Ny);
					J_b_x = stencil.getJx(S);
				    J_b_y = stencil.getJy(S);

					diagR[0] = L_b[C];
        		    ndiagR1[0] = L_b[N];
		            ndiagR2[0] = L_b[NE];
       
                    rhs[i+Nx+1] = fv[i+Nx+1] - L_b[S] * u[i+(1+J_b_y[S])*(Nx+1)] - L_b[W] * u[i+J_b_x[W]+Nx+1] 
						- L_b[E] * u[i+J_b_x[E]+Nx+1] - L_b[NW] * u[i+J_b_x[NW]+Nx+1] - L_b[SE] * u[i+J_b_x[SE]+(Nx+1)];
					
		        	for(size_t sum=7; sum<L_b.size(); sum++)
			        {
				       	rhs[i+Nx+1] -= L_b[sum] * u[i+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
				    }

                    L = stencil.get_L_c(i,2,Nx,Ny);

                    ndiagL1[0] = L[S];
					diagR[1] = L[C];
		            ndiagR1[1] = L[N];		            
 			        ndiagR2[1] = L[NE];

                    rhs[i+2*(Nx+1)] = fv[i+2*(Nx+1)] - L[W] * u[i+J_x[W]+2*(Nx+1)] - L[E] * u[i+J_x[E]+2*(Nx+1)]
					  -L[NW] * u[i+J_x[NW]+2*(Nx+1)] - L[SE] * u[i+J_x[SE]+2*(Nx+1)] - L[SW] * u[i+(2+J_y[SW])*(Nx+1)];

					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[i+2*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(2+J_y[sum])*(Nx+1)];
					}
    				    

					for(size_t j=3; j<Ny-2; j++)
					{
					   L = stencil.get_L_c(i,j,Nx,Ny);
					   ndiagL2[j-3] = L[SW];
					   ndiagL1[j-2] = L[S];
					   diagR[j-1] = L[C];
					   ndiagR1[j-1] = L[N];
					   ndiagR2[j-1] = L[NE];
 
                       rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[W] * u[i+J_x[W]+j*(Nx+1)] - L[E] * u[i+J_x[E]+j*(Nx+1)]
					         -L[NW] * u[i+J_x[NW]+j*(Nx+1)] - L[SE] * u[i+J_x[SE]+j*(Nx+1)];
					   
					   for(size_t sum=9; sum<L.size(); sum++)
					   {
						   rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					   }
					}
                                        
					L = stencil.get_L_c(i,Ny-2,Nx,Ny);
                    ndiagL2[Nx-5] = L[SW];
					ndiagL1[Nx-4] = L[S];
					diagR[Nx-3] = L[C];
                    ndiagR1[Nx-3] = L[N];

					rhs[i+(Ny-2)*(Nx+1)] = fv[i+(Ny-2)*(Nx+1)] - L[W] * u[i+J_x[W]+(Ny-2)*(Nx+1)] - L[E] * u[i+J_x[E]+(Ny-2)*(Nx+1)]
					     -L[NW] * u[i+J_x[NW]+(Ny-2)*(Nx+1)] - L[SE] * u[i+J_x[SE]+(Ny-2)*(Nx+1)] - L[NE] * u[i+(Ny-2+J_y[NE])*(Nx+1)];
                    
					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[i+(Ny-2)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-2+J_y[sum])*(Nx+1)];
					}


					L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
					J_b_x = stencil.getJx(N);
				    J_b_y = stencil.getJy(N);

					ndiagL2[Nx-4] = L_b[SW];
					ndiagL1[Nx-3] = L_b[S];
					diagR[Nx-2] = L_b[C];

					rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L_b[N] * u[i+(Ny-1+J_b_y[N])*(Nx+1)] 
						- L_b[W] * u[i+J_b_x[W]+(Ny-1)*(Nx+1)] - L_b[E] * u[i+J_b_x[E]+(Ny-1)*(Nx+1)] - L_b[NW] * u[i+J_b_x[NW]+(Ny-1)*(Nx+1)]
						- L_b[SE] * u[i+J_b_x[SE]+(Ny-1)*(Nx+1)];

					for(size_t sum=9; sum<L_b.size(); sum++)
					{
						rhs[i+(Ny-1)*(Nx+1)] -= L_b[sum] * u[i+J_b_x[sum]+(Ny-1+J_b_y[sum])*(Nx+1)];
					}
					
					// LR-decomposition + transformation of the rhs
					for(size_t k=1; k<Ny-2; k++)  
					{
						ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
					    diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
					    ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
					    rhs[(k+1)*(Nx+1) + i] = rhs[(k+1)*(Nx+1) + i] - ndiagL1[k-1] * rhs[k*(Nx+1)+i]; 

					    ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
					    ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
					    diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
					    rhs[(k+2)*(Nx+1) + i] = rhs[(k+2)*(Nx+1) + i] - ndiagL2[k-1] * rhs[k*(Nx+1)+i];
					}
					ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
				    diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
				    rhs[(Ny-1)*(Nx+1) + i] = rhs[(Ny-1)*(Nx+1) + i] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ i];

					// solve the linear system of equations R u = rhs
					u[i+(Nx+1)*(Ny-1)] = rhs[i+(Nx+1)*(Ny-1)] / diagR[Ny-2];

				    for(size_t k=Ny-2; k>1; k--)
					{
						u[i+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[i+k*(Nx+1)] - ndiagR1[k-1] * u[i+(k+1)*(Nx+1)] );

						rhs[(k-1)*(Nx+1)+i] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+i];
					}
					u[(Nx+1)+i] = 1/diagR[0] * ( rhs[(Nx+1)+i] - ndiagR1[0] * u[2*(Nx+1)+i] );
				}

////////////////letzte Spalte//////////////////
                
                L_c = stencil.get_L_se(Nx-1,1,Nx,Ny);
				J_c_x = stencil.getJx(SE);
				J_c_y = stencil.getJy(SE);

				
                rhs[Nx-1+Nx+1] = fv[Nx-1+Nx+1] - L_c[S] * u[Nx-1+(1+J_c_y[S])*(Nx+1)]
     				- L_c[W] * u[Nx-2+Nx+1] - L_c[E] * u[Nx+Nx+1] - L_c[NW] * u[Nx-1+J_c_x[NW]+Nx+1];
					
				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[Nx-1+Nx+1] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
				}

                L_b = stencil.get_L_e(Nx-1,2,Nx,Ny);
				J_b_x = stencil.getJx(E);
        		J_b_y = stencil.getJy(E);

                ndiagL1[0] = L_b[S];
				diagR[1] = L_b[C];
	            ndiagR1[1] = L_b[N];		            
		        ndiagR2[1] = L_b[NE];

                rhs[Nx-1+2*(Nx+1)] = fv[Nx-1+2*(Nx+1)] - L_b[W] * u[Nx-2+2*(Nx+1)] - L_b[E] * u[Nx+2*(Nx+1)]
				  -L_b[NW] * u[Nx-1+J_b_x[NW]+2*(Nx+1)] - L_b[SE] * u[Nx-1+(2+J_b_y[SE])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				    rhs[Nx-1+2*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(2+J_b_y[sum])*(Nx+1)];
				}

				for(size_t j=3; j<Ny-2; j++)  
				{
					L_b = stencil.get_L_e(Nx-1,j,Nx,Ny);
					
					ndiagL2[j-3] = L_b[SE];
					ndiagL1[j-2] = L_b[S];
			  	    diagR[j-1] = L_b[C];
  		            ndiagR1[j-1] = L_b[N];		            
			        ndiagR2[j-1] = L_b[NE];
						
					rhs[Nx-1+j*(Nx+1)] = fv[Nx-1+j*(Nx+1)] - L_b[W] * u[Nx-1+J_b_x[W]+j*(Nx+1)] - L_b[E] * u[Nx-1+J_b_x[E]+j*(Nx+1)]
							-L_b[NW] * u[Nx-1+J_b_x[NW]+j*(Nx+1)];

					for(size_t sum=8; sum<L_b.size(); sum++)
					{
						rhs[Nx-1+j*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(j+J_b_y[sum])*(Nx+1)];
					}
				}
					
				L_b = stencil.get_L_e(Nx-1,Ny-2,Nx,Ny);

				ndiagL2[Ny-5] = L_b[SE];
				ndiagL1[Ny-4] = L_b[S];
				diagR[Ny-3] = L_b[C];
		        ndiagR1[Ny-3] = L_b[N];

				rhs[Nx-1+(Ny-2)*(Nx+1)] = fv[Nx-1+(Ny-2)*(Nx+1)] - L_b[W] * u[Nx-2+(Ny-2)*(Nx+1)] - L_b[E] * u[Nx+(Ny-2)*(Nx+1)]
				  - L_b[NW] * u[Nx-1+J_b_x[NW]+(Ny-2)*(Nx+1)] - L_b[NE] * u[Nx-1+(Ny-2+J_b_y[NE])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				  rhs[Nx-1+(Ny-2)*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(Ny-2+J_b_y[sum])*(Nx+1)];
				}
					
                L_c = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
				J_c_x = stencil.getJx(NE);
				J_c_y = stencil.getJy(NE);

				ndiagL2[Ny-4] = L_c[NE];
				ndiagL1[Ny-3] = L_c[S];
				diagR[Ny-2] = L_c[C];

				rhs[Nx-1+(Ny-1)*(Nx+1)] = fv[Nx-1+(Ny-1)*(Nx+1)] - L_c[W] * u[Nx-2+(Ny-1)*(Nx+1)] 
				  - L_c[E] * u[Nx+(Ny-1)*(Nx+1)] - L_c[NW] * u[Nx-1+J_c_x[NW]+(Ny-1)*(Nx+1)] - L_c[N] * u[Nx-1+(Ny-1+J_c_y[N])*(Nx+1)];

				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[Nx-1+(Ny-1)*(Nx+1)] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(Ny-1+J_c_y[sum])*(Nx+1)];
				}

				// LR-decomposition + transformation of the rhs
				for(size_t k=1; k<Ny-2; k++)  
				{
					ndiagL1[k-1] = ndiagL1[k-1]/diagR[k-1];
					diagR[k] -= ndiagL1[k-1] * ndiagR1[k-1];
					ndiagR1[k] -= ndiagL1[k-1] * ndiagR2[k-1];
					rhs[(k+1)*(Nx+1) + Nx-1] = rhs[(k+1)*(Nx+1) + Nx-1] - ndiagL1[k-1] * rhs[k*(Nx+1)+Nx-1]; 

					ndiagL2[k-1] = ndiagL2[k-1]/diagR[k-1];
					ndiagL1[k] -= ndiagL2[k-1] * ndiagR1[k-1];
					diagR[k+1] -= ndiagL2[k-1] * ndiagR2[k-1];
					rhs[(k+2)*(Nx+1) + Nx-1] = rhs[(k+2)*(Nx+1) + Nx-1] - ndiagL2[k-1] * rhs[k*(Nx+1)+Nx-1];
				}
				ndiagL1[Ny-2-1] = ndiagL1[Ny-2-1]/diagR[Ny-2-1];
				diagR[Ny-2] -= ndiagL1[Ny-2-1] * ndiagR1[Ny-2-1];
				rhs[(Ny-1)*(Nx+1) + Nx-1] = rhs[(Ny-1)*(Nx+1) + Nx-1] - ndiagL1[Ny-2-1] * rhs[(Ny-2)*(Nx+1)+ Nx-1];

			    // solve the linear system of equations R u = rhs
			    u[Nx-1+(Nx+1)*(Ny-1)] = rhs[Nx-1+(Nx+1)*(Ny-1)] / diagR[Ny-2];

				for(size_t k=Ny-2; k>1; k--)
				{
					u[Nx-1+k*(Nx+1)] = 1/diagR[k-1] * ( rhs[Nx-1+k*(Nx+1)] - ndiagR1[k-1] * u[Nx-1+(k+1)*(Nx+1)] );

				    rhs[(k-1)*(Nx+1)+Nx-1] -= ndiagR2[k-2] * u[(k+1)*(Nx+1)+Nx-1];
				}
				u[(Nx+1)+Nx-1] = 1/diagR[0] * ( rhs[(Nx+1)+Nx-1] - ndiagR1[0] * u[2*(Nx+1)+Nx-1] );

     		}
		}
		else //parameter zu klein
		{
			  
			for(int k=0; k<2; k++)
			{
				Precision factor = 1.0;
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
