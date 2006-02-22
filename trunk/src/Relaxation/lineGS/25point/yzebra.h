namespace mg
{
	void lineGS::yzebra(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
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
				const std::valarray<int> J_x = stencil.getJx(c);
				const std::valarray<int> J_y = stencil.getJy(c);
				
				std::valarray<Precision> L_b = stencil.get_L_w(1,2,Nx,Ny);
				std::valarray<int> J_b_x = stencil.getJx(w);
				std::valarray<int> J_b_y = stencil.getJy(w);
 
                std::valarray<Precision> L_c = stencil.get_L_sw(1,1,Nx,Ny);
				std::valarray<int> J_c_x = stencil.getJx(sw);
				std::valarray<int> J_c_y = stencil.getJy(sw);

                diagR[0] = L_c[c];
	            ndiagR1[0] = L_c[n];
	            ndiagR2[0] = L_c[nw];
                                
                rhs[1+Nx+1] = fv[1+Nx+1] - L_c[s] * u[1+(1+J_c_y[s])*(Nx+1)]
     				- L_c[w] * u[Nx+1] - L_c[e] * u[2+Nx+1] - L_c[ne] * u[1+J_c_x[ne]+Nx+1];
					
				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[1+Nx+1] -= L_c[sum] * u[1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
				}

                L_b = stencil.get_L_w(1,2,Nx,Ny);
				J_b_x = stencil.getJx(w);
        		J_b_y = stencil.getJy(w);

                ndiagL1[0] = L_b[s];
				diagR[1] = L_b[c];
	            ndiagR1[1] = L_b[n];		            
		        ndiagR2[1] = L_b[nw];

                rhs[1+2*(Nx+1)] = fv[1+2*(Nx+1)] - L_b[w] * u[2*(Nx+1)] - L_b[e] * u[2+2*(Nx+1)]
				  -L_b[ne] * u[1+J_b_x[ne]+2*(Nx+1)] - L_b[se] * u[1+(2+J_b_y[se])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				    rhs[1+2*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(2+J_b_y[sum])*(Nx+1)];
				}

				for(size_t j=3; j<Ny-2; j++)  
				{
					ndiagL2[j-3] = L_b[se];
					ndiagL1[j-2] = L_b[s];
			  	    diagR[j-1] = L_b[c];
  		            ndiagR1[j-1] = L_b[n];		            
			        ndiagR2[j-1] = L_b[nw];
						
					rhs[1+j*(Nx+1)] = fv[1+j*(Nx+1)] - L_b[w] * u[1+J_b_x[w]+j*(Nx+1)] - L_b[e] * u[1+J_b_x[e]+j*(Nx+1)]
							-L_b[ne] * u[1+J_b_x[ne]+j*(Nx+1)];

					for(size_t sum=8; sum<L_b.size(); sum++)
					{
					  rhs[1+j*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(j+J_b_y[sum])*(Nx+1)];
					}
				}
					
				ndiagL2[Ny-5] = L_b[se];
				ndiagL1[Ny-4] = L_b[s];
				diagR[Ny-3] = L_b[c];
		        ndiagR1[Ny-3] = L_b[n];

				rhs[1+(Ny-2)*(Nx+1)] = fv[1+(Ny-2)*(Nx+1)] - L_b[w] * u[(Ny-2)*(Nx+1)] - L_b[e] * u[2+(Ny-2)*(Nx+1)]
				  -L_b[ne] * u[1+J_b_x[ne]+(Ny-2)*(Nx+1)] - L_b[nw] * u[1+(Ny-2+J_b_y[nw])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				  rhs[1+(Ny-2)*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(Ny-2+J_b_y[sum])*(Nx+1)];
				}
					
                L_c = stencil.get_L_nw(1,Ny-1,Nx,Ny);
				J_c_x = stencil.getJx(nw);
				J_c_y = stencil.getJy(nw);

				ndiagL2[Ny-4] = L_c[ne];
				ndiagL1[Ny-3] = L_c[s];
				diagR[Ny-2] = L_c[c];

				rhs[1+(Ny-1)*(Nx+1)] = fv[1+(Ny-1)*(Nx+1)] - L_c[w] * u[(Ny-1)*(Nx+1)] - L_c[e] * u[2+(Ny-1)*(Nx+1)]
					    - L_c[nw] * u[1+J_c_x[nw]+(Ny-1)*(Nx+1)] - L_c[n] * u[1+(Ny-1+J_c_y[n])*(Nx+1)];

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
				for(size_t i=3; i < Nx-2; i+=2)
				{
                    // setze rechte Seite					
					L_b = stencil.get_L_s(i,1,Nx,Ny);
					J_b_x = stencil.getJx(s);
				    J_b_y = stencil.getJy(s);

					diagR[0] = L_b[c];
        		    ndiagR1[0] = L_b[n];
		            ndiagR2[0] = L_b[ne];
       
                    rhs[i+Nx+1] = fv[i+Nx+1] - L_b[s] * u[i+(1+J_b_y[s])*(Nx+1)] - L_b[w] * u[i+J_b_x[w]+Nx+1] 
						- L_b[e] * u[i+J_b_x[e]+Nx+1] - L_b[nw] * u[i+J_b_x[nw]+Nx+1] - L_b[se] * u[i+J_b_x[se]+(Nx+1)];
					
		        	for(size_t sum=7; sum<L_b.size(); sum++)
			        {
				       	rhs[i+Nx+1] -= L_b[sum] * u[i+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
				    }

                    ndiagL1[0] = L[s];
					diagR[1] = L[c];
		            ndiagR1[1] = L[n];		            
 			        ndiagR2[1] = L[ne];

                    rhs[i+2*(Nx+1)] = fv[i+2*(Nx+1)] - L[w] * u[i+J_x[w]+2*(Nx+1)] - L[e] * u[i+J_x[e]+2*(Nx+1)]
					  -L[nw] * u[i+J_x[nw]+2*(Nx+1)] - L[se] * u[i+J_x[se]+2*(Nx+1)] - L[sw] * u[i+(2+J_y[sw])*(Nx+1)];

					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[i+2*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(2+J_y[sum])*(Nx+1)];
					}
    				    

					for(size_t j=3; j<Ny-2; j++)
					{
					   ndiagL2[j-3] = L[sw];
					   ndiagL1[j-2] = L[s];
					   diagR[j-1] = L[c];
					   ndiagR1[j-1] = L[n];
					   ndiagR2[j-1] = L[ne];
 
                       rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[w] * u[i+J_x[w]+j*(Nx+1)] - L[e] * u[i+J_x[e]+j*(Nx+1)]
					         -L[nw] * u[i+J_x[nw]+j*(Nx+1)] - L[se] * u[i+J_x[se]+j*(Nx+1)];
					   
					   for(size_t sum=9; sum<L.size(); sum++)
					   {
						   rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					   }
					}
                                        
                    ndiagL2[Nx-5] = L[sw];
					ndiagL1[Nx-4] = L[s];
					diagR[Nx-3] = L[c];
                    ndiagR1[Nx-3] = L[n];

					rhs[i+(Ny-2)*(Nx+1)] = fv[i+(Ny-2)*(Nx+1)] - L[w] * u[i+J_x[w]+(Ny-2)*(Nx+1)] - L[e] * u[i+J_x[e]+(Ny-2)*(Nx+1)]
					     -L[nw] * u[i+J_x[nw]+(Ny-2)*(Nx+1)] - L[se] * u[i+J_x[se]+(Ny-2)*(Nx+1)] - L[ne] * u[i+(Ny-2+J_y[ne])*(Nx+1)];
                    
					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[i+(Ny-2)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-2+J_y[sum])*(Nx+1)];
					}


					L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
					J_b_x = stencil.getJx(n);
				    J_b_y = stencil.getJy(n);

					ndiagL2[Nx-4] = L_b[sw];
					ndiagL1[Nx-3] = L_b[s];
					diagR[Nx-2] = L_b[c];

					rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L_b[n] * u[i+(Ny-1+J_b_y[n])*(Nx+1)] 
						- L_b[w] * u[i+J_b_x[w]+(Ny-1)*(Nx+1)] - L_b[e] * u[i+J_b_x[e]+(Ny-1)*(Nx+1)] - L_b[nw] * u[i+J_b_x[nw]+(Ny-1)*(Nx+1)]
						- L_b[se] * u[i+J_b_x[se]+(Ny-1)*(Nx+1)];

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
				J_c_x = stencil.getJx(se);
				J_c_y = stencil.getJy(se);

				
                rhs[Nx-1+Nx+1] = fv[Nx-1+Nx+1] - L_c[s] * u[Nx-1+(1+J_c_y[s])*(Nx+1)]
     				- L_c[w] * u[Nx-2+Nx+1] - L_c[e] * u[Nx+Nx+1] - L_c[nw] * u[Nx-1+J_c_x[nw]+Nx+1];
					
				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[Nx-1+Nx+1] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
				}

                L_b = stencil.get_L_e(Nx-1,2,Nx,Ny);
				J_b_x = stencil.getJx(e);
        		J_b_y = stencil.getJy(e);

                ndiagL1[0] = L_b[s];
				diagR[1] = L_b[c];
	            ndiagR1[1] = L_b[n];		            
		        ndiagR2[1] = L_b[ne];

                rhs[Nx-1+2*(Nx+1)] = fv[Nx-1+2*(Nx+1)] - L_b[w] * u[Nx-2+2*(Nx+1)] - L_b[e] * u[Nx+2*(Nx+1)]
				  -L_b[nw] * u[Nx-1+J_b_x[nw]+2*(Nx+1)] - L_b[se] * u[Nx-1+(2+J_b_y[se])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				    rhs[Nx-1+2*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(2+J_b_y[sum])*(Nx+1)];
				}

				for(size_t j=3; j<Ny-2; j++)  
				{
					ndiagL2[j-3] = L_b[se];
					ndiagL1[j-2] = L_b[s];
			  	    diagR[j-1] = L_b[c];
  		            ndiagR1[j-1] = L_b[n];		            
			        ndiagR2[j-1] = L_b[ne];
						
					rhs[Nx-1+j*(Nx+1)] = fv[Nx-1+j*(Nx+1)] - L_b[w] * u[Nx-1+J_b_x[w]+j*(Nx+1)] - L_b[e] * u[Nx-1+J_b_x[e]+j*(Nx+1)]
							-L_b[nw] * u[Nx-1+J_b_x[nw]+j*(Nx+1)];

					for(size_t sum=8; sum<L_b.size(); sum++)
					{
						rhs[Nx-1+j*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(j+J_b_y[sum])*(Nx+1)];
					}
				}
					
				ndiagL2[Ny-5] = L_b[se];
				ndiagL1[Ny-4] = L_b[s];
				diagR[Ny-3] = L_b[c];
		        ndiagR1[Ny-3] = L_b[n];

				rhs[Nx-1+(Ny-2)*(Nx+1)] = fv[Nx-1+(Ny-2)*(Nx+1)] - L_b[w] * u[Nx-2+(Ny-2)*(Nx+1)] - L_b[e] * u[Nx+(Ny-2)*(Nx+1)]
				  - L_b[nw] * u[Nx-1+J_b_x[nw]+(Ny-2)*(Nx+1)] - L_b[ne] * u[Nx-1+(Ny-2+J_b_y[ne])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				  rhs[Nx-1+(Ny-2)*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(Ny-2+J_b_y[sum])*(Nx+1)];
				}
					
                L_c = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
				J_c_x = stencil.getJx(ne);
				J_c_y = stencil.getJy(ne);

				ndiagL2[Ny-4] = L_c[ne];
				ndiagL1[Ny-3] = L_c[s];
				diagR[Ny-2] = L_c[c];

				rhs[Nx-1+(Ny-1)*(Nx+1)] = fv[Nx-1+(Ny-1)*(Nx+1)] - L_c[w] * u[Nx-2+(Ny-1)*(Nx+1)] 
				  - L_c[e] * u[Nx+(Ny-1)*(Nx+1)] - L_c[nw] * u[Nx-1+J_c_x[nw]+(Ny-1)*(Nx+1)] - L_c[n] * u[Nx-1+(Ny-1+J_c_y[n])*(Nx+1)];

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


				// durchlaufe alle inneren Spalten, Spaltenindex i				
				for(size_t i=2; i < Nx-1; i+=2)
				{
                    // setze rechte Seite					
					L_b = stencil.get_L_s(i,1,Nx,Ny);
					J_b_x = stencil.getJx(s);
				    J_b_y = stencil.getJy(s);

					diagR[0] = L_b[c];
        		    ndiagR1[0] = L_b[n];
		            ndiagR2[0] = L_b[ne];
       
                    rhs[i+Nx+1] = fv[i+Nx+1] - L_b[s] * u[i+(1+J_b_y[s])*(Nx+1)] - L_b[w] * u[i+J_b_x[w]+Nx+1] 
						- L_b[e] * u[i+J_b_x[e]+Nx+1] - L_b[nw] * u[i+J_b_x[nw]+Nx+1] - L_b[se] * u[i+J_b_x[se]+(Nx+1)];
					
		        	for(size_t sum=7; sum<L_b.size(); sum++)
			        {
				       	rhs[i+Nx+1] -= L_b[sum] * u[i+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
				    }

                    ndiagL1[0] = L[s];
					diagR[1] = L[c];
		            ndiagR1[1] = L[n];		            
 			        ndiagR2[1] = L[ne];

                    rhs[i+2*(Nx+1)] = fv[i+2*(Nx+1)] - L[w] * u[i+J_x[w]+2*(Nx+1)] - L[e] * u[i+J_x[e]+2*(Nx+1)]
					  -L[nw] * u[i+J_x[nw]+2*(Nx+1)] - L[se] * u[i+J_x[se]+2*(Nx+1)] - L[sw] * u[i+(2+J_y[sw])*(Nx+1)];

					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[i+2*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(2+J_y[sum])*(Nx+1)];
					}
    				    

					for(size_t j=3; j<Ny-2; j++)
					{
					   ndiagL2[j-3] = L[sw];
					   ndiagL1[j-2] = L[s];
					   diagR[j-1] = L[c];
					   ndiagR1[j-1] = L[n];
					   ndiagR2[j-1] = L[ne];
 
                       rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[w] * u[i+J_x[w]+j*(Nx+1)] - L[e] * u[i+J_x[e]+j*(Nx+1)]
					         -L[nw] * u[i+J_x[nw]+j*(Nx+1)] - L[se] * u[i+J_x[se]+j*(Nx+1)];
					   
					   for(size_t sum=9; sum<L.size(); sum++)
					   {
						   rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					   }
					}
                                        
                    ndiagL2[Nx-5] = L[sw];
					ndiagL1[Nx-4] = L[s];
					diagR[Nx-3] = L[c];
                    ndiagR1[Nx-3] = L[n];

					rhs[i+(Ny-2)*(Nx+1)] = fv[i+(Ny-2)*(Nx+1)] - L[w] * u[i+J_x[w]+(Ny-2)*(Nx+1)] - L[e] * u[i+J_x[e]+(Ny-2)*(Nx+1)]
					     -L[nw] * u[i+J_x[nw]+(Ny-2)*(Nx+1)] - L[se] * u[i+J_x[se]+(Ny-2)*(Nx+1)] - L[ne] * u[i+(Ny-2+J_y[ne])*(Nx+1)];
                    
					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[i+(Ny-2)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-2+J_y[sum])*(Nx+1)];
					}


					L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
					J_b_x = stencil.getJx(n);
				    J_b_y = stencil.getJy(n);

					ndiagL2[Nx-4] = L_b[sw];
					ndiagL1[Nx-3] = L_b[s];
					diagR[Nx-2] = L_b[c];

					rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L_b[n] * u[i+(Ny-1+J_b_y[n])*(Nx+1)] 
						- L_b[w] * u[i+J_b_x[w]+(Ny-1)*(Nx+1)] - L_b[e] * u[i+J_b_x[e]+(Ny-1)*(Nx+1)] - L_b[nw] * u[i+J_b_x[nw]+(Ny-1)*(Nx+1)]
						- L_b[se] * u[i+J_b_x[se]+(Ny-1)*(Nx+1)];

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

			}
			else // stencil not constant
			{
				std::valarray<Precision> L = stencil.get_L_c(2,2,Nx,Ny);
				std::valarray<int> J_x = stencil.getJx(c);
				std::valarray<int> J_y = stencil.getJy(c);
				
				std::valarray<Precision> L_b = stencil.get_L_w(1,2,Nx,Ny);
				std::valarray<int> J_b_x = stencil.getJx(w);
				std::valarray<int> J_b_y = stencil.getJy(w);
 
                std::valarray<Precision> L_c = stencil.get_L_sw(1,1,Nx,Ny);
				std::valarray<int> J_c_x = stencil.getJx(sw);
				std::valarray<int> J_c_y = stencil.getJy(sw);

                diagR[0] = L_c[c];
	            ndiagR1[0] = L_c[n];
	            ndiagR2[0] = L_c[nw];
                                
                rhs[1+Nx+1] = fv[1+Nx+1] - L_c[s] * u[1+(1+J_c_y[s])*(Nx+1)]
     				- L_c[w] * u[Nx+1] - L_c[e] * u[2+Nx+1] - L_c[ne] * u[1+J_c_x[ne]+Nx+1];
					
				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[1+Nx+1] -= L_c[sum] * u[1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
				}

                L_b = stencil.get_L_w(1,2,Nx,Ny);
				J_b_x = stencil.getJx(w);
        		J_b_y = stencil.getJy(w);

                ndiagL1[0] = L_b[s];
				diagR[1] = L_b[c];
	            ndiagR1[1] = L_b[n];		            
		        ndiagR2[1] = L_b[nw];

                rhs[1+2*(Nx+1)] = fv[1+2*(Nx+1)] - L_b[w] * u[2*(Nx+1)] - L_b[e] * u[2+2*(Nx+1)]
				  -L_b[ne] * u[1+J_b_x[ne]+2*(Nx+1)] - L_b[se] * u[1+(2+J_b_y[se])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				    rhs[1+2*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(2+J_b_y[sum])*(Nx+1)];
				}

				for(size_t j=3; j<Ny-2; j++)  
				{
					L_b = stencil.get_L_w(1,j,Nx,Ny);
					
					ndiagL2[j-3] = L_b[se];
					ndiagL1[j-2] = L_b[s];
			  	    diagR[j-1] = L_b[c];
  		            ndiagR1[j-1] = L_b[n];		            
			        ndiagR2[j-1] = L_b[nw];
						
					rhs[1+j*(Nx+1)] = fv[1+j*(Nx+1)] - L_b[w] * u[1+J_b_x[w]+j*(Nx+1)] - L_b[e] * u[1+J_b_x[e]+j*(Nx+1)]
							-L_b[ne] * u[1+J_b_x[ne]+j*(Nx+1)];

					for(size_t sum=8; sum<L_b.size(); sum++)
					{
					  rhs[1+j*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(j+J_b_y[sum])*(Nx+1)];
					}
				}

				L_b = stencil.get_L_w(1,Ny-2,Nx,Ny);
					
				ndiagL2[Ny-5] = L_b[se];
				ndiagL1[Ny-4] = L_b[s];
				diagR[Ny-3] = L_b[c];
		        ndiagR1[Ny-3] = L_b[n];

				rhs[1+(Ny-2)*(Nx+1)] = fv[1+(Ny-2)*(Nx+1)] - L_b[w] * u[(Ny-2)*(Nx+1)] - L_b[e] * u[2+(Ny-2)*(Nx+1)]
				  -L_b[ne] * u[1+J_b_x[ne]+(Ny-2)*(Nx+1)] - L_b[nw] * u[1+(Ny-2+J_b_y[nw])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				  rhs[1+(Ny-2)*(Nx+1)] -= L_b[sum] * u[1+J_b_x[sum]+(Ny-2+J_b_y[sum])*(Nx+1)];
				}
					
                L_c = stencil.get_L_nw(1,Ny-1,Nx,Ny);
				J_c_x = stencil.getJx(nw);
				J_c_y = stencil.getJy(nw);

				ndiagL2[Ny-4] = L_c[ne];
				ndiagL1[Ny-3] = L_c[s];
				diagR[Ny-2] = L_c[c];

				rhs[1+(Ny-1)*(Nx+1)] = fv[1+(Ny-1)*(Nx+1)] - L_c[w] * u[(Ny-1)*(Nx+1)] - L_c[e] * u[2+(Ny-1)*(Nx+1)]
					    - L_c[nw] * u[1+J_c_x[nw]+(Ny-1)*(Nx+1)] - L_c[n] * u[1+(Ny-1+J_c_y[n])*(Nx+1)];

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
				for(size_t i=3; i < Nx-2; i+=2)
				{
                    // setze rechte Seite					
					L_b = stencil.get_L_s(i,1,Nx,Ny);
					J_b_x = stencil.getJx(s);
				    J_b_y = stencil.getJy(s);

					diagR[0] = L_b[c];
        		    ndiagR1[0] = L_b[n];
		            ndiagR2[0] = L_b[ne];
       
                    rhs[i+Nx+1] = fv[i+Nx+1] - L_b[s] * u[i+(1+J_b_y[s])*(Nx+1)] - L_b[w] * u[i+J_b_x[w]+Nx+1] 
						- L_b[e] * u[i+J_b_x[e]+Nx+1] - L_b[nw] * u[i+J_b_x[nw]+Nx+1] - L_b[se] * u[i+J_b_x[se]+(Nx+1)];
					
		        	for(size_t sum=7; sum<L_b.size(); sum++)
			        {
				       	rhs[i+Nx+1] -= L_b[sum] * u[i+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
				    }

					L = stencil.get_L_c(i,2,Nx,Ny);

                    ndiagL1[0] = L[s];
					diagR[1] = L[c];
		            ndiagR1[1] = L[n];		            
 			        ndiagR2[1] = L[ne];

                    rhs[i+2*(Nx+1)] = fv[i+2*(Nx+1)] - L[w] * u[i+J_x[w]+2*(Nx+1)] - L[e] * u[i+J_x[e]+2*(Nx+1)]
					  -L[nw] * u[i+J_x[nw]+2*(Nx+1)] - L[se] * u[i+J_x[se]+2*(Nx+1)] - L[sw] * u[i+(2+J_y[sw])*(Nx+1)];

					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[i+2*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(2+J_y[sum])*(Nx+1)];
					}
    				    

					for(size_t j=3; j<Ny-2; j++)
					{
					   L = stencil.get_L_c(i,j,Nx,Ny);
					   
					   ndiagL2[j-3] = L[sw];
					   ndiagL1[j-2] = L[s];
					   diagR[j-1] = L[c];
					   ndiagR1[j-1] = L[n];
					   ndiagR2[j-1] = L[ne];
 
                       rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[w] * u[i+J_x[w]+j*(Nx+1)] - L[e] * u[i+J_x[e]+j*(Nx+1)]
					         -L[nw] * u[i+J_x[nw]+j*(Nx+1)] - L[se] * u[i+J_x[se]+j*(Nx+1)];
					   
					   for(size_t sum=9; sum<L.size(); sum++)
					   {
						   rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					   }
					}

					L = stencil.get_L_c(i,Ny-2,Nx,Ny);
                                        
                    ndiagL2[Nx-5] = L[sw];
					ndiagL1[Nx-4] = L[s];
					diagR[Nx-3] = L[c];
                    ndiagR1[Nx-3] = L[n];

					rhs[i+(Ny-2)*(Nx+1)] = fv[i+(Ny-2)*(Nx+1)] - L[w] * u[i+J_x[w]+(Ny-2)*(Nx+1)] - L[e] * u[i+J_x[e]+(Ny-2)*(Nx+1)]
					     -L[nw] * u[i+J_x[nw]+(Ny-2)*(Nx+1)] - L[se] * u[i+J_x[se]+(Ny-2)*(Nx+1)] - L[ne] * u[i+(Ny-2+J_y[ne])*(Nx+1)];
                    
					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[i+(Ny-2)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-2+J_y[sum])*(Nx+1)];
					}


					L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
					J_b_x = stencil.getJx(n);
				    J_b_y = stencil.getJy(n);

					ndiagL2[Nx-4] = L_b[sw];
					ndiagL1[Nx-3] = L_b[s];
					diagR[Nx-2] = L_b[c];

					rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L_b[n] * u[i+(Ny-1+J_b_y[n])*(Nx+1)] 
						- L_b[w] * u[i+J_b_x[w]+(Ny-1)*(Nx+1)] - L_b[e] * u[i+J_b_x[e]+(Ny-1)*(Nx+1)] - L_b[nw] * u[i+J_b_x[nw]+(Ny-1)*(Nx+1)]
						- L_b[se] * u[i+J_b_x[se]+(Ny-1)*(Nx+1)];

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
				J_c_x = stencil.getJx(se);
				J_c_y = stencil.getJy(se);

				
                rhs[Nx-1+Nx+1] = fv[Nx-1+Nx+1] - L_c[s] * u[Nx-1+(1+J_c_y[s])*(Nx+1)]
     				- L_c[w] * u[Nx-2+Nx+1] - L_c[e] * u[Nx+Nx+1] - L_c[nw] * u[Nx-1+J_c_x[nw]+Nx+1];
					
				for(size_t sum=7; sum<L_c.size(); sum++)
				{
					rhs[Nx-1+Nx+1] -= L_c[sum] * u[Nx-1+J_c_x[sum]+(1+J_c_y[sum])*(Nx+1)];
				}

                L_b = stencil.get_L_e(Nx-1,2,Nx,Ny);
				J_b_x = stencil.getJx(e);
        		J_b_y = stencil.getJy(e);

                ndiagL1[0] = L_b[s];
				diagR[1] = L_b[c];
	            ndiagR1[1] = L_b[n];		            
		        ndiagR2[1] = L_b[ne];

                rhs[Nx-1+2*(Nx+1)] = fv[Nx-1+2*(Nx+1)] - L_b[w] * u[Nx-2+2*(Nx+1)] - L_b[e] * u[Nx+2*(Nx+1)]
				  -L_b[nw] * u[Nx-1+J_b_x[nw]+2*(Nx+1)] - L_b[se] * u[Nx-1+(2+J_b_y[se])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				    rhs[Nx-1+2*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(2+J_b_y[sum])*(Nx+1)];
				}

				for(size_t j=3; j<Ny-2; j++)  
				{
					L_b = stencil.get_L_e(Nx-1,j,Nx,Ny);

					ndiagL2[j-3] = L_b[se];
					ndiagL1[j-2] = L_b[s];
			  	    diagR[j-1] = L_b[c];
  		            ndiagR1[j-1] = L_b[n];		            
			        ndiagR2[j-1] = L_b[ne];
						
					rhs[Nx-1+j*(Nx+1)] = fv[Nx-1+j*(Nx+1)] - L_b[w] * u[Nx-1+J_b_x[w]+j*(Nx+1)] - L_b[e] * u[Nx-1+J_b_x[e]+j*(Nx+1)]
							-L_b[nw] * u[Nx-1+J_b_x[nw]+j*(Nx+1)];

					for(size_t sum=8; sum<L_b.size(); sum++)
					{
						rhs[Nx-1+j*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(j+J_b_y[sum])*(Nx+1)];
					}
				}

				L_b = stencil.get_L_e(Nx-1,Ny-2,Nx,Ny);
					
				ndiagL2[Ny-5] = L_b[se];
				ndiagL1[Ny-4] = L_b[s];
				diagR[Ny-3] = L_b[c];
		        ndiagR1[Ny-3] = L_b[n];

				rhs[Nx-1+(Ny-2)*(Nx+1)] = fv[Nx-1+(Ny-2)*(Nx+1)] - L_b[w] * u[Nx-2+(Ny-2)*(Nx+1)] - L_b[e] * u[Nx+(Ny-2)*(Nx+1)]
				  - L_b[nw] * u[Nx-1+J_b_x[nw]+(Ny-2)*(Nx+1)] - L_b[ne] * u[Nx-1+(Ny-2+J_b_y[ne])*(Nx+1)];

				for(size_t sum=8; sum<L_b.size(); sum++)
				{
				  rhs[Nx-1+(Ny-2)*(Nx+1)] -= L_b[sum] * u[Nx-1+J_b_x[sum]+(Ny-2+J_b_y[sum])*(Nx+1)];
				}
					
                L_c = stencil.get_L_ne(Nx-1,Ny-1,Nx,Ny);
				J_c_x = stencil.getJx(ne);
				J_c_y = stencil.getJy(ne);

				ndiagL2[Ny-4] = L_c[ne];
				ndiagL1[Ny-3] = L_c[s];
				diagR[Ny-2] = L_c[c];

				rhs[Nx-1+(Ny-1)*(Nx+1)] = fv[Nx-1+(Ny-1)*(Nx+1)] - L_c[w] * u[Nx-2+(Ny-1)*(Nx+1)] 
				  - L_c[e] * u[Nx+(Ny-1)*(Nx+1)] - L_c[nw] * u[Nx-1+J_c_x[nw]+(Ny-1)*(Nx+1)] - L_c[n] * u[Nx-1+(Ny-1+J_c_y[n])*(Nx+1)];

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


				// durchlaufe alle inneren Spalten, Spaltenindex i				
				for(size_t i=2; i < Nx-1; i+=2)
				{
                    // setze rechte Seite
					L_b = stencil.get_L_s(i,1,Nx,Ny);
					J_b_x = stencil.getJx(s);
				    J_b_y = stencil.getJy(s);

					diagR[0] = L_b[c];
        		    ndiagR1[0] = L_b[n];
		            ndiagR2[0] = L_b[ne];
       
                    rhs[i+Nx+1] = fv[i+Nx+1] - L_b[s] * u[i+(1+J_b_y[s])*(Nx+1)] - L_b[w] * u[i+J_b_x[w]+Nx+1] 
						- L_b[e] * u[i+J_b_x[e]+Nx+1] - L_b[nw] * u[i+J_b_x[nw]+Nx+1] - L_b[se] * u[i+J_b_x[se]+(Nx+1)];
					
		        	for(size_t sum=7; sum<L_b.size(); sum++)
			        {
				       	rhs[i+Nx+1] -= L_b[sum] * u[i+J_b_x[sum]+(1+J_b_y[sum])*(Nx+1)];
				    }

					L = stencil.get_L_c(i,2,Nx,Ny);

                    ndiagL1[0] = L[s];
					diagR[1] = L[c];
		            ndiagR1[1] = L[n];		            
 			        ndiagR2[1] = L[ne];

                    rhs[i+2*(Nx+1)] = fv[i+2*(Nx+1)] - L[w] * u[i+J_x[w]+2*(Nx+1)] - L[e] * u[i+J_x[e]+2*(Nx+1)]
					  -L[nw] * u[i+J_x[nw]+2*(Nx+1)] - L[se] * u[i+J_x[se]+2*(Nx+1)] - L[sw] * u[i+(2+J_y[sw])*(Nx+1)];

					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[i+2*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(2+J_y[sum])*(Nx+1)];
					}
    				    

					for(size_t j=3; j<Ny-2; j++)
					{
					   L = stencil.get_L_c(i,j,Nx,Ny);

					   ndiagL2[j-3] = L[sw];
					   ndiagL1[j-2] = L[s];
					   diagR[j-1] = L[c];
					   ndiagR1[j-1] = L[n];
					   ndiagR2[j-1] = L[ne];
 
                       rhs[i+j*(Nx+1)] = fv[i+j*(Nx+1)] - L[w] * u[i+J_x[w]+j*(Nx+1)] - L[e] * u[i+J_x[e]+j*(Nx+1)]
					         -L[nw] * u[i+J_x[nw]+j*(Nx+1)] - L[se] * u[i+J_x[se]+j*(Nx+1)];
					   
					   for(size_t sum=9; sum<L.size(); sum++)
					   {
						   rhs[i+j*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(j+J_y[sum])*(Nx+1)];
					   }
					}

					L = stencil.get_L_c(i,Ny-2,Nx,Ny);
                                        
                    ndiagL2[Nx-5] = L[sw];
					ndiagL1[Nx-4] = L[s];
					diagR[Nx-3] = L[c];
                    ndiagR1[Nx-3] = L[n];

					rhs[i+(Ny-2)*(Nx+1)] = fv[i+(Ny-2)*(Nx+1)] - L[w] * u[i+J_x[w]+(Ny-2)*(Nx+1)] - L[e] * u[i+J_x[e]+(Ny-2)*(Nx+1)]
					     -L[nw] * u[i+J_x[nw]+(Ny-2)*(Nx+1)] - L[se] * u[i+J_x[se]+(Ny-2)*(Nx+1)] - L[ne] * u[i+(Ny-2+J_y[ne])*(Nx+1)];
                    
					for(size_t sum=9; sum<L.size(); sum++)
					{
						rhs[i+(Ny-2)*(Nx+1)] -= L[sum] * u[i+J_x[sum]+(Ny-2+J_y[sum])*(Nx+1)];
					}


					L_b = stencil.get_L_n(i,Ny-1,Nx,Ny);
					J_b_x = stencil.getJx(n);
				    J_b_y = stencil.getJy(n);

					ndiagL2[Nx-4] = L_b[sw];
					ndiagL1[Nx-3] = L_b[s];
					diagR[Nx-2] = L_b[c];

					rhs[i+(Ny-1)*(Nx+1)] = fv[i+(Ny-1)*(Nx+1)] - L_b[n] * u[i+(Ny-1+J_b_y[n])*(Nx+1)] 
						- L_b[w] * u[i+J_b_x[w]+(Ny-1)*(Nx+1)] - L_b[e] * u[i+J_b_x[e]+(Ny-1)*(Nx+1)] - L_b[nw] * u[i+J_b_x[nw]+(Ny-1)*(Nx+1)]
						- L_b[se] * u[i+J_b_x[se]+(Ny-1)*(Nx+1)];

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
