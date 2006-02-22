/** \file max_residuum.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementaion of the function max_residuum
 * \see max_residuum.h
 */
#include "max_residuum.h"

namespace mg
{
	Precision max_residuum(const std::valarray<Precision>& u,
						const std::valarray<Precision>& fv,
						const Stencil& stencil,
						const size_t Nx,const size_t Ny)
	{
		Precision result=0;
		if (stencil.size() < 2)
		{
			for (size_t j=1; j<Ny;j++)
				for(size_t i=1; i<Nx;i++)
				{
					Precision temp_res=fv[j*(Nx+1)+i]
							-stencil.apply_c(u,i,j,Nx,Ny);
					result = std::fabs(temp_res) > result
									? std::fabs(temp_res) : result;
				}
		}
		else
		{
			//south west corner
			Precision temp_res=fv[1*(Nx+1)+1]
					-stencil.apply_sw(u,1,1,Nx,Ny);
			result = std::fabs(temp_res) > result ? std::fabs(temp_res) : result;
			//south east corner
			temp_res=fv[1*(Nx+1)+(Nx-1)]
					-stencil.apply_se(u,(Nx-1),1,Nx,Ny);
			result = std::fabs(temp_res) > result ? std::fabs(temp_res) : result;
			//north west corner
			temp_res=fv[(Nx-1)*(Nx+1)+1]
					-stencil.apply_nw(u,1,(Ny-1),Nx,Ny);
			result = std::fabs(temp_res) > result ? std::fabs(temp_res) : result;
			//north east corner
			temp_res=fv[(Nx-1)*(Nx+1)+(Nx-1)]
				-stencil.apply_ne(u,(Nx-1),(Ny-1),Nx,Ny);
			result = std::fabs(temp_res) > result ? std::fabs(temp_res) : result;
			//south boarder
			for (size_t i=2;i<(Nx-1);i++)
			{
				temp_res=fv[1*(Nx+1)+i]
					-stencil.apply_s(u,i,1,Nx,Ny);
				result = std::fabs(temp_res) > result
									? std::fabs(temp_res) : result;
			}
			//north boarder
			for (size_t i=2;i<(Nx-1);i++)
			{
				temp_res=fv[(Nx-1)*(Nx+1)+i]
					-stencil.apply_n(u,i,(Ny-1),Nx,Ny);
				result = std::fabs(temp_res) > result
									? std::fabs(temp_res) : result;
			}
			//west boarder
			for (size_t j=2;j<(Ny-1);j++)
			{
				temp_res=fv[j*(Nx+1)+1]
					-stencil.apply_w(u,1,j,Nx,Ny);
				result = std::fabs(temp_res) > result
									? std::fabs(temp_res) : result;
			}
			//east boarder
			for (size_t j=2;j<(Ny-1);j++)
			{
				temp_res=fv[j*(Nx+1)+(Nx-1)]
					-stencil.apply_e(u,(Nx-1),j,Nx,Ny);
				result = std::fabs(temp_res) > result
									? std::fabs(temp_res) : result;
			}
			//the center
			for (size_t j=2;j<(Ny-1);j++)
			{
				for (size_t i=2;i<(Nx-1);i++)
				{
					temp_res=fv[j*(Nx+1)+i]
							-stencil.apply_c(u,i,j,Nx,Ny);
					result = std::fabs(temp_res) > result
									? std::fabs(temp_res) : result;
				}
			}
		}
		return result;
	}
	
}
