/** \file lineJAC.cpp
 * \author Andre Oeckerath
 * \see lineJAC.h
 */
#include "lineJAC.h"
#include<iostream>
#include "../functions/residuum.h"
#include "./lineJAC/9point/x.h"
#include "./lineJAC/9point/y.h"
#include "./lineJAC/25point/x.h"
#include "./lineJAC/25point/y.h"
#include "./lineJAC/9point/xzebra.h"
#include "./lineJAC/9point/yzebra.h"
#include "./lineJAC/25point/xzebra.h"
#include "./lineJAC/25point/yzebra.h"


namespace mg
{
			
	void lineJAC::relax(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		               const Stencil &stencil, const size_t Nx, const size_t Ny) const
	{
		// valarray needed for LR-decomposition of a tridiagonal matrix
		

		std::valarray<Precision> resid(0.0,(Nx+1)*(Ny+1));
		resid = residuum(u,fv,stencil,Nx,Ny);

		switch (stencil.size())
		{
		case 1:  // stern der größe 1
			{
				switch (direction_)
				{
				case 0:
					{ 
					    ninepointxline(u, fv, resid, stencil, Nx, Ny);
                        resid = residuum(u,fv,stencil,Nx,Ny);
                        ninepointyline(u, fv, resid, stencil, Nx, Ny);
						break;
					}

			    case 1:
					{
						ninepointxline(u, fv, resid, stencil, Nx, Ny);
						break;
					}
					
                case 2:
					{
						ninepointyline(u, fv, resid, stencil, Nx, Ny);
						break;
					}

                case 3:
					{ 
					    ninepointxzebra(u, fv, resid, stencil, Nx, Ny);
                        resid = residuum(u,fv,stencil,Nx,Ny);
                        ninepointyzebra(u, fv, resid, stencil, Nx, Ny);
						break;
					}

			    case 4:
					{
						ninepointxzebra(u, fv, resid, stencil, Nx, Ny);
						break;
					}
					
                case 5:
					{
						ninepointyzebra(u, fv, resid, stencil, Nx, Ny);
						break;
					}
				               
                default:
					{
						std::cout << "Error in direction of the line relaxation!\n";
						break;
					}
				}
		        break;			
			}

	    case 2:  // stern der größe 2
			{
				switch (direction_)
				{
				case 0:
					{ 
						xline(u, fv, resid, stencil, Nx, Ny);
                        resid = residuum(u,fv,stencil,Nx,Ny);
						yline(u, fv, resid, stencil, Nx, Ny);
						break;
					}

			    case 1:
					{
						xline(u, fv, resid, stencil, Nx, Ny);
						break;
					}
					
                case 2:
					{
						yline(u, fv, resid, stencil, Nx, Ny);
						break;
					}
					
				case 3:
					{ 
					    xzebra(u, fv, resid, stencil, Nx, Ny);
                        resid = residuum(u,fv,stencil,Nx,Ny);
                        yzebra(u, fv, resid, stencil, Nx, Ny);
						break;
					}

			    case 4:
					{
						xzebra(u, fv, resid, stencil, Nx, Ny);
						break;
					}
					
                case 5:
					{
						yzebra(u, fv, resid, stencil, Nx, Ny);
						break;
					}

                default:
					{
						std::cout << "Error in direction of the line relaxation!\n";
						break;
					}
				}
				break;
			}

	    default:
			{
				std::cout << "Sterngröße nicht behandelbar!" << std::endl;
				break;
			}


		}

	}


}
