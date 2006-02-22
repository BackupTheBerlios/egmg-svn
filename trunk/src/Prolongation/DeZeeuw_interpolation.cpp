/** \file DeZeeuw_interpolation.cpp
 * \author Benedikt Engbroks
 * \brief Contains the implementation of the class DeZeeuw_interpolation
 * \see DeZeeuw_interpolation.h
 */
#include "DeZeeuw_interpolation.h"
#include "../Stencil/Stencil.h"
#include <algorithm>

namespace mg
{
// This function only works for max. compact 9-point stencils
std::valarray<Precision> DeZeeuw_interpolation::prolongate(
		const std::valarray<Precision>& u, const Stencil& stencil, const size_t Nx,const size_t Ny) const
{
	register const size_t Nx_new = 2*Nx;
	register const size_t Ny_new = 2*Ny;
	std::valarray<Precision> result((Nx_new+1)*(Ny_new+1));
	std::valarray<int> J_x = stencil.getJx(c);
	std::valarray<int> J_y = stencil.getJy(c);
	std::valarray<Precision> L = stencil.get_L_c(0,0,Nx,Ny);
	std::valarray<size_t> posi(9);
	std::valarray<Precision> mS(9);
	std::valarray<Precision> mT(9);
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

	Precision scale = 0;
	Precision weight1 = 0;
	Precision weight2 = 0;
	Precision erg = 0;
	Precision symsum;
	Precision d_w;
	Precision d_e;
	Precision d_n;
	Precision d_s;
	Precision sigma1;
	Precision c_1;
	Precision w_w;
	Precision w_e;
	Precision sigma2;
	Precision c_2;
	Precision w_s;
	Precision w_n;

	//"interpolation" of coarse grid points
	for (size_t j=0;j<=Ny;j++)
		for (size_t i=0;i<=Nx;i++)
			result[2*j*(Nx_new+1)+2*i]=u[j*(Nx+1)+i];
	
	//interpolation of fine grid points on coarse grid lines
	for (size_t j=0;j<=Ny_new;j+=2)
		for (size_t i=1;i<=Nx_new;i+=2)
		{
			L = stencil.get_L_c(i,j,Nx,Ny);
			symsum = 0;
			
			// Divide the stencil defined by L und posi into a symmetric and an antisymmetric part.
			for (size_t k=0; k<9; k++)
			{
				mS[k] = 0.5*(L[posi[k]] + L[posi[8-k]]);
				symsum += mS[k];
				mT[k] = 0.5*(L[posi[k]] - L[posi[8-k]]);
			}
			d_w = std::max(std::fabs(mS[SW] + mS[W] + mS[NW]), std::max(std::fabs(mS[SW]), std::fabs(mS[NW])));
			d_e = std::max(std::fabs(mS[SE] + mS[E] + mS[NE]), std::max(std::fabs(mS[SE]), std::fabs(mS[NE])));
			d_n = std::max(std::fabs(mS[NW] + mS[N] + mS[NE]), std::max(std::fabs(mS[NW]), std::fabs(mS[NE])));
			d_s = std::max(std::fabs(mS[SW] + mS[S] + mS[SE]), std::max(std::fabs(mS[SW]), std::fabs(mS[SE])));
			sigma1 = 0.5 * std::min(1.0, std::fabs(1 - symsum/L[0]));
			c_1 = mT[SE] + mT[E] + mT[NE] - mT[SW] - mT[W] - mT[NW];
			w_w = sigma1 * (1 + (d_w-d_e)/(d_w+d_e) + c_1 / (d_w+d_e+d_n+d_s));
			w_e = 2*sigma1-w_w;
			weight1 = std::min(2*sigma1, std::max(w_w, 0.0));
			weight2 = std::min(2*sigma1, std::max(w_e, 0.0));
			result[j*(Nx_new+1)+i]=(weight1*result[j*(Nx_new+1)+i-1] + weight2*result[j*(Nx_new+1)+i+1]);
		}
	
	//interpolation of fine grid points on fine grid lines and coarse grid columns
	for (size_t j=1;j<=Ny_new;j+=2)
		for (size_t i=0;i<=Nx_new;i+=2)
		{
			L = stencil.get_L_c(i,j,Nx,Ny);
			symsum = 0;
			
			// Divide the stencil defined by L und posi into a symmetric and an antisymmetric part.
			for (int k=0; k<9; k++)
			{
				mS[k] = 0.5*(L[posi[k]] + L[posi[8-k]]);
				symsum += mS[k];
				mT[k] = 0.5*(L[posi[k]] - L[posi[8-k]]);
			}
			d_w = std::max(std::fabs(mS[SW] + mS[W] + mS[NW]), std::max(std::fabs(mS[SW]), std::fabs(mS[NW])));
			d_e = std::max(std::fabs(mS[SE] + mS[E] + mS[NE]), std::max(std::fabs(mS[SE]), std::fabs(mS[NE])));
			d_n = std::max(std::fabs(mS[NW] + mS[N] + mS[NE]), std::max(std::fabs(mS[NW]), std::fabs(mS[NE])));
			d_s = std::max(std::fabs(mS[SW] + mS[S] + mS[SE]), std::max(std::fabs(mS[SW]), std::fabs(mS[SE])));
			sigma2 = 0.5 * std::min(1.0, std::fabs(1 - symsum/L[0]));
			c_2 = mT[NW] + mT[N] + mT[NE] - mT[SW] - mT[S] - mT[SE];
			w_n = sigma2 * (1 + (d_s-d_n)/(d_s+d_n) + c_2 / (d_w+d_e+d_n+d_s));
			w_s = 2*sigma2-w_n;
			weight1 = std::min(2*sigma2, std::max(w_s, 0.0));
			weight2 = std::min(2*sigma2, std::max(w_n, 0.0));
			result[j*(Nx_new+1)+i]=(weight1*result[(j-1)*(Nx_new+1)+i] + weight2*result[(j+1)*(Nx_new+1)+i]);
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
