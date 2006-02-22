/** \file putbval.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the function putbval
 */
#include "putbval.h"

namespace mg
{
	void putbval(std::valarray<Precision>& u,const size_t Nx, const size_t Ny,
					const function2D g)
	{
		const Precision hx = 1.0/Nx;
		const Precision hy = 1.0/Ny;
		for(size_t i=0;i<=Nx;i++)
		{
			u[i]=g(i*hx,0);						//first line
			u[u.size()-1-i]=g((Nx-i)*hx,1);		//last line
		}
		for (size_t i=0;i<=Ny;i++)
		{
			u[i*(Nx+1)]=g(0,i*hy);				//left boarder
			u[(i+1)*(Nx+1)-1]=g(1,i*hy);		//right boarder
		
		}
	}
}
