/** \file twonorm_residuum.cpp
 * \author Matthias Rettenmeier
 * \brief contains the implimentaion of the function twonorm_residuum
 * \see twonorm_residuum.h
 */
#include "twonorm_residuum.h"

namespace mg
{
	Precision twonorm_residuum(const std::valarray<Precision>& u,
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
					temp_res=temp_res*temp_res;
                                        result += temp_res;
				}
		}
		else
		{
			//south west corner
			Precision temp_res=fv[1*(Nx+1)+1]
					-stencil.apply_sw(u,1,1,Nx,Ny);
			temp_res=temp_res*temp_res;
			result += temp_res;
			//south east corner
			temp_res=fv[1*(Nx+1)+(Nx-1)]
					-stencil.apply_se(u,(Nx-1),1,Nx,Ny);
			temp_res=temp_res*temp_res;
			result += temp_res;
			//north west corner
			temp_res=fv[(Nx-1)*(Nx+1)+1]
					-stencil.apply_nw(u,1,(Ny-1),Nx,Ny);
			temp_res=temp_res*temp_res;
			result += temp_res;
			//north east corner
			temp_res=fv[(Nx-1)*(Nx+1)+(Nx-1)]
				-stencil.apply_ne(u,(Nx-1),(Ny-1),Nx,Ny);
			temp_res=temp_res*temp_res;
			result += temp_res;
			//south boarder
			for (size_t i=2;i<(Nx-1);i++)
			{
				temp_res=fv[1*(Nx+1)+i]
					-stencil.apply_s(u,i,1,Nx,Ny);
				temp_res=temp_res*temp_res;
                                result += temp_res;
			}
			//north boarder
			for (size_t i=2;i<(Nx-1);i++)
			{
				temp_res=fv[(Nx-1)*(Nx+1)+i]
					-stencil.apply_n(u,i,(Ny-1),Nx,Ny);
				temp_res=temp_res*temp_res;
                                result += temp_res;
			}
			//west boarder
			for (size_t j=2;j<(Ny-1);j++)
			{
				temp_res=fv[j*(Nx+1)+1]
					-stencil.apply_w(u,1,j,Nx,Ny);
				temp_res=temp_res*temp_res;
                                result += temp_res;
			}
			//east boarder
			for (size_t j=2;j<(Ny-1);j++)
			{
				temp_res=fv[j*(Nx+1)+(Nx-1)]
					-stencil.apply_e(u,(Nx-1),j,Nx,Ny);
				temp_res=temp_res*temp_res;
                                result += temp_res;
			}
			//the center
			for (size_t j=2;j<(Ny-1);j++)
			{
				for (size_t i=2;i<(Nx-1);i++)
				{
					temp_res=fv[j*(Nx+1)+i]
							-stencil.apply_c(u,i,j,Nx,Ny);
					temp_res=temp_res*temp_res;
                                        result += temp_res;
				}
			}
		}
                result = std::sqrt(result);
		return result;
	}

}
