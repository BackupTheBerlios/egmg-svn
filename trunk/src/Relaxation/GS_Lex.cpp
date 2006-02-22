/** \file GS_Lex.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementaion of the class GS_Lex.
 * \see GS_Lex.h
 */
#include "GS_Lex.h"

namespace mg
{
	void GS_Lex::relax(std::valarray<Precision>& u,
					const std::valarray<Precision>& fv,
					const Stencil& stencil,
					const size_t Nx,const size_t Ny) const
	{
		Precision factor = 1.0;
		if (stencil.size() < 2)
		{
			for (size_t j=1;j<Ny;j++)
				for (size_t i=1;i<Nx;i++)
				{
					factor = 1.0/stencil.get_center_c(i,j,Nx,Ny);
					u[j*(Nx+1)+i]+=factor*(fv[j*(Nx+1)+i]
							-stencil.apply_c(u,i,j,Nx,Ny));
				}
		}
		else
		{
			//south west corner
			factor = 1.0/stencil.get_center_sw(1,1,Nx,Ny);
			u[1*(Nx+1)+1]+=factor*(fv[1*(Nx+1)+1]
					-stencil.apply_sw(u,1,1,Nx,Ny));
			//south boarder
			for (size_t i=2;i<(Nx-1);i++)
			{
				factor = 1.0/stencil.get_center_s(i,1,Nx,Ny);
				u[1*(Nx+1)+i]+=factor*(fv[1*(Nx+1)+i]
					-stencil.apply_s(u,i,1,Nx,Ny));
			}
			//south east corner
			factor = 1.0/stencil.get_center_se((Nx-1),1,Nx,Ny);
			u[1*(Nx+1)+(Nx-1)]+=factor*(fv[1*(Nx+1)+(Nx-1)]
					-stencil.apply_se(u,(Nx-1),1,Nx,Ny));
			//everything up to north west corner
			for (size_t j=2;j<(Ny-1);j++)
			{
				//west boarder point in j. line
				factor = 1.0/stencil.get_center_w(1,j,Nx,Ny);
				u[j*(Nx+1)+1]+=factor*(fv[j*(Nx+1)+1]
					-stencil.apply_w(u,1,j,Nx,Ny));
				//center points
				for (size_t i=2;i<(Nx-1);i++)
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
			//north west corner
			factor = 1.0/stencil.get_center_nw(1,(Nx-1),Nx,Ny);
			u[(Nx-1)*(Nx+1)+1]+=factor*(fv[(Nx-1)*(Nx+1)+1]
					-stencil.apply_nw(u,1,(Ny-1),Nx,Ny));
			//north boarder
			for (size_t i=2;i<(Nx-1);i++)
			{
				factor = 1.0/stencil.get_center_n(i,(Nx-1),Nx,Ny);
				u[(Nx-1)*(Nx+1)+i]+=factor*(fv[(Nx-1)*(Nx+1)+i]
					-stencil.apply_n(u,i,(Ny-1),Nx,Ny));
			}
			//north east corner
			factor = 1.0/stencil.get_center_ne((Nx-1),(Nx-1),Nx,Ny);
			u[(Nx-1)*(Nx+1)+(Nx-1)]+=factor*(fv[(Nx-1)*(Nx+1)+(Nx-1)]
				-stencil.apply_ne(u,(Nx-1),(Ny-1),Nx,Ny));
		}	
	}
}
