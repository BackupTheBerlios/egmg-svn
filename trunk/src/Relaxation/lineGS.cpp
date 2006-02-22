/** \file lineGS.cpp
 * \author Andre Oeckerath
 * \see lineGS.h
 */
#include "lineGS.h"
#include<iostream>
#include "./lineGS/9point/x.h"
#include "./lineGS/9point/y.h"
#include "./lineGS/9point/xzebra.h"
#include "./lineGS/9point/yzebra.h"
#include "./lineGS/25point/x.h"
#include "./lineGS/25point/y.h"
#include "./lineGS/25point/xzebra.h"
#include "./lineGS/25point/yzebra.h"

namespace mg
{
			
	void lineGS::relax(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		               const Stencil &stencil, const size_t Nx, const size_t Ny) const
	{
		// valarrays needed for LR-decomposition of a tridiagonal matrix
		std::valarray<Precision> rhs(u);
		rhs=0;
		
		
		switch (stencil.size())
		{
		case 1:  // stern der größe 1
			{
				switch (direction_)
				{
				case 0:
					{ 
						ninepointxline(u, fv, rhs, stencil, Nx, Ny);
						ninepointyline(u, fv, rhs, stencil, Nx, Ny);
						break;
					}

			    case 1:
					{
						ninepointxline(u, fv, rhs, stencil, Nx, Ny);
						break;
					}
					
                case 2:
					{
						ninepointyline(u, fv, rhs, stencil, Nx, Ny);
						break;
					}

				case 3:
					{ 
						ninepointxzebra(u, fv, rhs, stencil, Nx, Ny);
                        ninepointyzebra(u, fv, rhs, stencil, Nx, Ny);
						break;
					}

				case 4:
					{
						ninepointxzebra(u, fv, rhs, stencil, Nx, Ny);
			            break;
					}
		        
                case 5:
					{
						ninepointyzebra(u, fv, rhs, stencil, Nx, Ny);
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
						xline(u, fv, rhs, stencil, Nx, Ny);
						yline(u, fv, rhs, stencil, Nx, Ny);
						break;
					}

			    case 1:
					{
						xline(u, fv, rhs, stencil, Nx, Ny);
						break;
					}
					
                case 2:
					{
						yline(u, fv, rhs, stencil, Nx, Ny);
						break;
					}

				case 3:
					{ 
						xzebra(u, fv, rhs, stencil, Nx, Ny);
						yzebra(u, fv, rhs, stencil, Nx, Ny);
						break;
					}

				case 4:
					{
						xzebra(u, fv, rhs, stencil, Nx, Ny);
			            break;
					}
		        
                case 5:
					{
						yzebra(u, fv, rhs, stencil, Nx, Ny);
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
