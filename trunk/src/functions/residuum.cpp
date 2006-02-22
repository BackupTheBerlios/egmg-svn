/** \file residuum.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementaion of the function residuum.
 * \see residuum.h
 */
#include "residuum.h"

namespace mg
{
	std::valarray<Precision> residuum(const std::valarray<Precision>& u,
								   const std::valarray<Precision>& fv,
								   const Stencil& stencil,
								   const size_t Nx,const size_t Ny)
	{
		std::valarray<Precision> result(0.0,u.size());
		if (stencil.size() < 2)
		{
			for (size_t j=1; j<Ny;j++)
				for(size_t i=1; i<Nx;i++)
				{
					result[j*(Nx+1)+i]=fv[j*(Nx+1)+i]
							-stencil.apply_c(u,i,j,Nx,Ny);
				}
		}
		else
		{
			//south west corner
			result[1*(Nx+1)+1]=fv[1*(Nx+1)+1]
					-stencil.apply_sw(u,1,1,Nx,Ny);
			//south east corner
			result[1*(Nx+1)+(Nx-1)]=fv[1*(Nx+1)+(Nx-1)]
					-stencil.apply_se(u,(Nx-1),1,Nx,Ny);
			//north west corner
			result[(Nx-1)*(Nx+1)+1]=fv[(Nx-1)*(Nx+1)+1]
					-stencil.apply_nw(u,1,(Ny-1),Nx,Ny);
			//north east corner
			result[(Nx-1)*(Nx+1)+(Nx-1)]=fv[(Nx-1)*(Nx+1)+(Nx-1)]
				-stencil.apply_ne(u,(Nx-1),(Ny-1),Nx,Ny);
			//south boarder
			for (size_t i=2;i<(Nx-1);i++)
			{
				result[1*(Nx+1)+i]=fv[1*(Nx+1)+i]
					-stencil.apply_s(u,i,1,Nx,Ny);
			}
			//north boarder
			for (size_t i=2;i<(Nx-1);i++)
			{
				result[(Nx-1)*(Nx+1)+i]=fv[(Nx-1)*(Nx+1)+i]
					-stencil.apply_n(u,i,(Ny-1),Nx,Ny);
			}
			//west boarder
			for (size_t j=2;j<(Ny-1);j++)
			{
				result[j*(Nx+1)+1]=fv[j*(Nx+1)+1]
					-stencil.apply_w(u,1,j,Nx,Ny);
			}
			//east boarder;
			for (size_t j=2;j<(Ny-1);j++)
			{
				result[j*(Nx+1)+(Nx-1)]=fv[j*(Nx+1)+(Nx-1)]
					-stencil.apply_e(u,(Nx-1),j,Nx,Ny);
			}
			//the center
			for (size_t j=2;j<(Ny-1);j++)
			{
				for (size_t i=2;i<(Nx-1);i++)
				{
					result[j*(Nx+1)+i]=fv[j*(Nx+1)+i]
							-stencil.apply_c(u,i,j,Nx,Ny);
				}
			}
		}
		return result;
	}	
}
