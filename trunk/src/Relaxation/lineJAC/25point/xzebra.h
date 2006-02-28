namespace mg
{
	void ZebraLineJAC::xzebra(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> resid, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const

	{
        if((Ny > 4) && (Nx > 4))
		{	
			std::valarray<Precision> rhs(0.0,Nx-1);
		    std::valarray<Precision> temp(0.0,(Nx+1)*(Ny+1));

			std::valarray<Precision> diagR(0.0,Nx-1);
		    std::valarray<Precision> ndiagR1(0.0,Nx-2);
		    std::valarray<Precision> ndiagL1(0.0,Nx-2);
		    std::valarray<Precision> ndiagR2(0.0,Nx-3);
		    std::valarray<Precision> ndiagL2(0.0,Nx-3);

			if(stencil.isConstant() == true)
			{
				// get const operator L
				const std::valarray<Precision> L = stencil.get_L_c(2,2,Nx,Ny);
				const std::valarray<int> J_x = stencil.getJx(C);
				const std::valarray<int> J_y = stencil.getJy(C);
				
				std::valarray<Precision> L_b = stencil.get_L_s(2,1,Nx,Ny);
				std::valarray<int> J_b_x = stencil.getJx(S);
				std::valarray<int> J_b_y = stencil.getJy(S);
                    
                std::valarray<Precision> L_c = stencil.get_L_sw(1,1,Nx,Ny);
				std::valarray<int> J_c_x = stencil.getJx(SW);
				std::valarray<int> J_c_y = stencil.getJy(SW);


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
				std::valarray<Precision> L = stencil.get_L_c(2,2,Nx,Ny);
				std::valarray<int> J_x = stencil.getJx(C);
				std::valarray<int> J_y = stencil.getJy(C);
				
				std::valarray<Precision> L_b = stencil.get_L_s(2,1,Nx,Ny);
				std::valarray<int> J_b_x = stencil.getJx(S);
				std::valarray<int> J_b_y = stencil.getJy(S);
                    
                std::valarray<Precision> L_c = stencil.get_L_sw(1,1,Nx,Ny);
				std::valarray<int> J_c_x = stencil.getJx(SW);
				std::valarray<int> J_c_y = stencil.getJy(SW);
				

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
}

