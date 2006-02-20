/** \file putrhs.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the function putrhs
 */
#include "putrhs.h"

namespace mg
{
	void putrhs(std::valarray<precision>& fv,const size_t Nx, const size_t Ny,
					const function2D f)
	{
		const precision hx=1.0/Nx;
		const precision hy=1.0/Ny;
		for (size_t j=0;j<=Ny;j++)
			for(size_t i=0;i<=Nx;i++)
				fv[j*(Nx+1)+i]=f(i*hx,j*hy);
	}
}
