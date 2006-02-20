/** \file GS_RB.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementaion of the class GS_RB.
 * \see GS_RB.h
 */
#include "GS_RB.h"

namespace mg
{

void GS_RB::relax(std::valarray<precision>& u,
				const std::valarray<precision>& fv,
				const Stencil& stencil,
				const size_t Nx,const size_t Ny) const
{
    precision factor = 1.0;
	if (stencil.size() < 2)
	{
		//first do the red points
		for (size_t j=1;j<Ny;j+=2)
			for (size_t i=1;i<Nx;i+=2)
			{
				factor = 1.0/stencil.get_center_c(i,j,Nx,Ny);
				u[j*(Nx+1)+i]+=factor*(fv[j*(Nx+1)+i]
							-stencil.apply_c(u,i,j,Nx,Ny));
			}
		for (size_t j=2;j<(Ny-1);j+=2)
			for (size_t i=2;i<(Nx-1);i+=2)
			{
				factor = 1.0/stencil.get_center_c(i,j,Nx,Ny);
				u[j*(Nx+1)+i]+=factor*(fv[j*(Nx+1)+i]
							-stencil.apply_c(u,i,j,Nx,Ny));
			}
		//do black points
		for (size_t j=1;j<Ny;j+=2)
			for (size_t i=2;i<(Nx-1);i+=2)
			{
				factor = 1.0/stencil.get_center_c(i,j,Nx,Ny);
				u[j*(Nx+1)+i]+=factor*(fv[j*(Nx+1)+i]
							-stencil.apply_c(u,i,j,Nx,Ny));
			}
		for (size_t j=2;j<(Ny-1);j+=2)
			for (size_t i=1;i<Nx;i+=2)
			{
				factor = 1.0/stencil.get_center_c(i,j,Nx,Ny);
				u[j*(Nx+1)+i]+=factor*(fv[j*(Nx+1)+i]
							-stencil.apply_c(u,i,j,Nx,Ny));
			}
	}
	else
	{
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

