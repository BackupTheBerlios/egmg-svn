/** \file Half_Weighting.cpp
 * \author Mareike Riepl, Michaela Anghel, Daniela Steffes-Lai
 * \brief Contains the implementation of the class Half Weighting
 */
#include "Half_Weighting.h"
#include <stdexcept>

namespace mg
{
NumericArray Half_Weighting::restriction(
		const NumericArray& u, const Stencil&,const Prolongation&,
		const size_t Nx, const size_t Ny) const
{
	//if it is not possible to do standart coarsening throw an exeption
	if ((Nx%2 != 0) || (Ny%2 != 0))
		throw std::domain_error("u");
	const size_t Nx_new = Nx/2;
	const size_t Ny_new = Ny/2;
	NumericArray result(0.0,(Nx_new+1)*(Ny_new+1));
	//do injection on the boarders
	for (size_t j=0;j<=Ny_new;j++)
	{
		result[j*(Nx_new+1)]=u[2*j*(Nx+1)];				// north boarder
		result[j*(Nx_new+1)+Nx_new]=u[2*j*(Nx+1)+Nx];  // west boarder
		result[j]=u[2*j];							   // south boarder
	    result[Ny_new*(Nx_new+1)+j]=u[Ny*(Nx+1)+2*j]; // east boarder
	}
	//do half-weighting on the centered points
	for (size_t j=1;j<Ny_new;j++)
	{
		for(size_t i=1;i<Nx_new;i++)
			result[j*(Nx_new+1)+i]=weight_*(4*u[2*j*(Nx+1)+2*i]
							+1*u[2*j*(Nx+1)+2*i+1]+1*u[2*j*(Nx+1)+2*i-1]
							+1*u[(2*j+1)*(Nx+1)+2*i]+1*u[(2*j-1)*(Nx+1)+2*i] 
							)/8.0;
	}
	return result;
}
}

