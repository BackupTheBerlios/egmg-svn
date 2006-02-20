/** \file Dendy_interpolation.cpp
 * \author Benedikt Engbroks
 * \brief Contains the implementation of the class Dendy_interpolation
 * \see Dendy_interpolation.h
 */
#include "Dendy_interpolation.h"

namespace mg
{
// This function only works for max. compact 9-point stencils
std::valarray<precision> Dendy_interpolation::prolong(
		const std::valarray<precision>& u, const Stencil& stencil, const size_t Nx,const size_t Ny) const
{
	register const size_t Nx_new = 2*Nx;
	register const size_t Ny_new = 2*Ny;
	std::valarray<precision> result((Nx_new+1)*(Ny_new+1));
	std::valarray<int> J_x = stencil.get_J_x(c);
	std::valarray<int> J_y = stencil.get_J_y(c);
	std::valarray<precision> L = stencil.get_L_c(0,0,Nx,Ny);
	std::valarray<size_t> posi(9);
	// posi[0] = No. of the sw-element of the stencil
	// posi[1] = No. of the s-element of the stencil
	// posi[2] = No. of the se-element of the stencil
	// posi[3] = No. of the w-element of the stencil
	// posi[4] = No. of the c-element of the stencil (i.e. 0)
	// posi[5] = No. of the e-element of the stencil
	// posi[6] = No. of the nw-element of the stencil
	// posi[7] = No. of the n-element of the stencil
	// posi[8] = No. of the ne-element of the stencil
	// (see also the definitions of the constants in the .h-file)
	for (size_t j=0;j<J_x.size();j++)
	{
		posi[(J_x[j]+1)+3*(J_y[j]+1)] = j;
	}

	precision scale = 0;
	precision weight1 = 0;
	precision weight2 = 0;
	precision erg = 0;
	
	//"interpolation" of coarse grid points
	for (size_t j=0;j<=Ny;j++)
		for (size_t i=0;i<=Nx;i++)
			result[2*j*(Nx_new+1)+2*i]=u[j*(Nx+1)+i];
	
	//interpolation of fine grid points on coarse grid lines
	for (size_t j=0;j<=Ny_new;j+=2)
		for (size_t i=1;i<=Nx_new;i+=2)
		{
			L = stencil.get_L_c(i,j,Nx,Ny);
			scale = -L[0];
			weight1 = 0;
			weight2 = 0;
			if (posi[S] != 0) scale -= L[posi[S]];
			if (posi[N] != 0) scale -= L[posi[N]];
			if (posi[W] != 0) weight1 += L[posi[W]];
			if (posi[SW] != 0) weight1 += L[posi[SW]];
			if (posi[NW] != 0) weight1 += L[posi[NW]];
			if (posi[E] != 0) weight2 += L[posi[E]];
			if (posi[SE] != 0) weight2 += L[posi[SE]];
			if (posi[NE] != 0) weight2 += L[posi[NE]];
			result[j*(Nx_new+1)+i]=(weight1*result[j*(Nx_new+1)+i-1] + weight2*result[j*(Nx_new+1)+i+1]) / scale;
		}
	
	//interpolation of fine grid points on fine grid lines and coarse grid columns
	for (size_t j=1;j<=Ny_new;j+=2)
		for (size_t i=0;i<=Nx_new;i+=2)
		{
			L = stencil.get_L_c(i,j,Nx,Ny);
			scale = -L[0];
			weight1 = 0;
			weight2 = 0;
			if (posi[W] != 0) scale -= L[posi[W]];
			if (posi[E] != 0) scale -= L[posi[E]];
			if (posi[S] != 0) weight1 += L[posi[S]];
			if (posi[SW] != 0) weight1 += L[posi[SW]];
			if (posi[SE] != 0) weight1 += L[posi[SE]];
			if (posi[N] != 0) weight2 += L[posi[N]];
			if (posi[NW] != 0) weight2 += L[posi[NW]];
			if (posi[NE] != 0) weight2 += L[posi[NE]];
			result[j*(Nx_new+1)+i]=(weight1*result[(j-1)*(Nx_new+1)+i] + weight2*result[(j+1)*(Nx_new+1)+i]) / scale;
		}
			
	//interpolation of fine grid points on fine grid lines and fine grid columns
	for (size_t j=1;j<=Ny_new;j+=2)
		for (size_t i=1;i<=Nx_new;i+=2)
		{
			L = stencil.get_L_c(i,j,Nx,Ny);
			erg = 0;
			scale = -L[0];
			if (posi[W] != 0) erg += L[posi[W]]*result[j*(Nx_new+1)+i-1];
			if (posi[E] != 0) erg += L[posi[E]]*result[j*(Nx_new+1)+i+1];
			if (posi[S] != 0) erg += L[posi[S]]*result[(j-1)*(Nx_new+1)+i];
			if (posi[SW] != 0) erg += L[posi[SW]]*result[(j-1)*(Nx_new+1)+i-1];
			if (posi[SE] != 0) erg += L[posi[SE]]*result[(j-1)*(Nx_new+1)+i+1];
			if (posi[N] != 0) erg += L[posi[N]]*result[(j+1)*(Nx_new+1)+i];
			if (posi[NW] != 0) erg += L[posi[NW]]*result[(j+1)*(Nx_new+1)+i-1];
			if (posi[NE] != 0) erg += L[posi[NE]]*result[(j+1)*(Nx_new+1)+i+1];
			result[j*(Nx_new+1)+i]=erg / scale;
		}
	return result;
}

}
