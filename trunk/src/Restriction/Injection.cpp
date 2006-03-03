/** \file Injection.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class Injection
 */
#include "Injection.h"
#include <stdexcept>

namespace mg
{
NumericArray Injection::restriction(
    const NumericArray& u,
    const Stencil&,
    const Prolongation&,
    const size_t nx, const size_t ny) const
{
	//if it is not possible to do standart coarsening throw an exeption
	if ((nx%2 != 0) || (ny%2 != 0))
		throw std::domain_error("u");
	const size_t nxNew = nx/2;
	const size_t nyNew = ny/2;
	NumericArray result(0.0,(nxNew+1)*(nyNew+1));
	for (size_t sy=0;sy<=nyNew;sy++)
	{
		for(size_t sx=0;sx<=nxNew;sx++)
			result[sy*(nxNew+1)+sx]=weight_*u[2*sy*(nx+1)+2*sx];
	}
	return result;
}
}
