/** \file error.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the function error
 */
#include<cmath>
#include "error.h"

namespace mg
{
	Precision error(const std::valarray<Precision>& u,const size_t Nx,
					const size_t Ny,const function2D solu)
	{
		Precision result=0;
		const Precision hx=1.0/Nx;
		const Precision hy=1.0/Ny;
		for (size_t j=1; j<Nx;j++)
			for(size_t i=1; i<Ny;i++)
			{
				Precision temp = std::fabs(solu(i*hx,j*hy)-u[j*(Nx+1)+i]);
				if ((result-temp) < 0)
					result=temp;
			}
			
		return result;
	}
}
