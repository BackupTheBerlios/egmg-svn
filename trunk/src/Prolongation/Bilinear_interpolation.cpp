/** \file Bilinear_interpolation.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class Bilinear_interpolation
 */
#include "Bilinear_interpolation.h"

namespace mg
{
	
std::valarray<precision> Bilinear_interpolation::prolong(
		const std::valarray<precision>& u, const Stencil&,
									const size_t Nx,const size_t Ny) const
{
	const size_t Nx_new = 2*Nx;
	const size_t Ny_new = 2*Ny;
	std::valarray<precision> result(0.0,(Nx_new+1)*(Nx_new+1));
	//"interpolation" of coarse grid points
	for (size_t j=0;j<=Ny;j++)
		for (size_t i=0;i<=Nx;i++)
			result[2*j*(Nx_new+1)+2*i]=u[j*(Nx+1)+i];
	//interpolation of fine grid points on coarse grid lines
	for (size_t j=0;j<=Ny_new;j+=2)
		for (size_t i=1;i<=Nx_new;i+=2)
			result[j*(Nx_new+1)+i]=1./2*(result[j*(Nx_new+1)+i-1]
											+result[j*(Nx_new+1)+i+1]);
	//interpolation of fine grid points on fine grid lines
	for (size_t j=1;j<=Ny_new;j+=2)
		for (size_t i=0;i<=Nx_new;i++)
			result[j*(Nx_new+1)+i]=1./2*(result[(j-1)*(Nx_new+1)+i]
										+result[(j+1)*(Nx_new+1)+i]);
	return result;
}

}
