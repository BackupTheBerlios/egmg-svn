/** \file Bicubic_interpolation.cpp
 * \author Benedikt Engbroks
 * \brief Contains the implementaion of the class Bicubic_interpolation
 * \see Bicubic_interpolation.h
 */
#include "Bicubic_interpolation.h"

namespace mg
{
	
std::valarray<Precision> Bicubic_interpolation::prolongate(
		const std::valarray<Precision>& u, const Stencil&,
									const size_t Nx,const size_t Ny) const
{
	const size_t Nx_new = 2*Nx;
	const size_t Ny_new = 2*Ny;
	std::valarray<Precision> result((Nx_new+1)*(Ny_new+1));
	
	//"interpolation" of coarse grid points
	for (size_t j=0;j<=Ny;j++)
		for (size_t i=0;i<=Nx;i++)
			result[2*j*(Nx_new+1)+2*i]=u[j*(Nx+1)+i];
		
	if (Nx_new > 3 && Ny_new > 3)
	{
		//interpolation of fine grid points on coarse grid lines
		for (size_t j=0;j<=Ny_new;j+=2)
		{
			result[j*(Nx_new+1)+1]=1./8 * (3*result[j*(Nx_new+1)] + 6*result[j*(Nx_new+1)+2] - result[j*(Nx_new+1)+4]);
			for (size_t i=3;i<=Nx_new-2;i+=2)
				result[j*(Nx_new+1)+i]= 1./16 * (-result[j*(Nx_new+1)+i-3] + 9*result[j*(Nx_new+1)+i-1] + 9*result[j*(Nx_new+1)+i+1] - result[j*(Nx_new+1)+i+3]);
			result[j*(Nx_new+1)+Nx_new-1]= 1./8 * (-result[j*(Nx_new+1)+Nx_new-4] + 6*result[j*(Nx_new+1)+Nx_new-2] + 3*result[j*(Nx_new+1)+Nx_new]);
		}
		
		//interpolation of fine grid points on fine grid lines
		for (size_t i=0;i<=Nx_new;i++)
		{
			result[Nx_new+1+i] = 1./8 * (3*result[i] + 6*result[2*(Nx_new+1)+i] - result[4*(Nx_new+1)+i]); 
			for (size_t j=3;j<=Ny_new-2;j+=2)
				result[j*(Nx_new+1)+i]=1./16*(-result[(j-3)*(Nx_new+1)+i] + 9*result[(j-1)*(Nx_new+1)+i] + 9*result[(j+1)*(Nx_new+1)+i] - result[(j+3)*(Nx_new+1)+i]);
			result[(Ny_new-1)*(Nx_new+1)+i] = 1./8 * (-result[(Ny_new-4)*(Nx_new+1)+i] + 6*result[(Ny_new-2)*(Nx_new+1)+i] + 3*result[Ny_new*(Nx_new+1)+i]); 		
		}
	}
	else
	{
		//bilinear interpolation on a grid with less than 4 fine grid points in each direction

		//interpolation of fine grid points on coarse grid lines
		for (size_t j=0;j<=Ny_new;j+=2)
			for (size_t i=1;i<=Nx_new;i+=2)
				result[j*(Nx_new+1)+i]=1./2*(result[j*(Nx_new+1)+i-1]+result[j*(Nx_new+1)+i+1]);
		
		//interpolation of fine grid points on fine grid lines
		for (size_t j=1;j<=Ny_new;j+=2)
			for (size_t i=0;i<=Nx_new;i++)
				result[j*(Nx_new+1)+i]=1./2*(result[(j-1)*(Nx_new+1)+i]+result[(j+1)*(Nx_new+1)+i]);
	}
	return result;
}

}
