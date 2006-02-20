/** \file Injection.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class Injection
 */
#include "Injection.h"
#include <stdexcept>

namespace mg
{
std::valarray<precision> Injection::restriction(
		const std::valarray<precision>& u, const Stencil&,const Prolongation&,
		const size_t Nx, const size_t Ny) const
{
	//if it is not possible to do standart coarsening throw an exeption
	if ((Nx%2 != 0) || (Ny%2 != 0))
		throw std::domain_error("u");
	const size_t Nx_new = Nx/2;
	const size_t Ny_new = Ny/2;
	std::valarray<precision> result(0.0,(Nx_new+1)*(Ny_new+1));
	for (size_t j=0;j<=Ny_new;j++)
	{
		for(size_t i=0;i<=Nx_new;i++)
			result[j*(Nx_new+1)+i]=weight_*u[2*j*(Nx+1)+2*i];
	}
	return result;
}
}
