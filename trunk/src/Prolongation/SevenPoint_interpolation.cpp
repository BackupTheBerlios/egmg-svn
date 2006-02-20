/** \file SevenPoint_interpolation.cpp
 * \author Benedikt Engbroks
 * \brief Contains the implementaion of the class SevenPoint_interpolation
 * \see SevenPoint_interpolation.h
 */
#include "SevenPoint_interpolation.h"

namespace mg
{
	
std::valarray<precision> SevenPoint_interpolation::prolong(
		const std::valarray<precision>& u, const Stencil&,
									const size_t Nx,const size_t Ny) const
{
	const size_t Nx_new = 2*Nx;
	const size_t Ny_new = 2*Ny;
	std::valarray<precision> result((Nx_new+1)*(Ny_new+1));
	//"interpolation" of coarse grid points
	for (size_t j=0;j<=Ny;j++)
		for (size_t i=0;i<=Nx;i++)
			result[2*j*(Nx_new+1)+2*i]=u[j*(Nx+1)+i];
	//interpolation of fine grid points on coarse grid lines
	for (size_t j=0;j<=Ny_new;j+=2)
		for (size_t i=1;i<=Nx_new;i+=2)
			result[j*(Nx_new+1)+i]=1./2*(result[j*(Nx_new+1)+i-1]
			+result[j*(Nx_new+1)+i+1]);
	//interpolation of fine grid points on fine grid lines and coarse
	//grid columns
	for (size_t j=1;j<=Ny_new;j+=2)
		for (size_t i=0;i<=Nx_new;i+=2)
			result[j*(Nx_new+1)+i]=1./2*(result[(j-1)*(Nx_new+1)+i]
			+result[(j+1)*(Nx_new+1)+i]);
	//interpolation of fine grid points on fine grid lines and fine grid columns
	for (size_t j=1;j<=Ny_new;j+=2)
		for (size_t i=1;i<=Nx_new;i+=2)
			result[j*(Nx_new+1)+i]=1./2*(result[(j-1)*(Nx_new+1)+i-1]
			+result[(j+1)*(Nx_new+1)+i+1]);
	return result;
}

}
