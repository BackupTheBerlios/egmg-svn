/** \file TransposedProlongation.cpp
 * \author Benedikt Engbroks
 * \brief Contains the implementation of the class TransposedProlongation.
 */

#include "TransposedProlongation.h"
#include <stdexcept>

namespace mg
{

NumericArray TransposedProlongation::restriction(
    const NumericArray& u,
    const Stencil& stencil,
    const Index nx,
    const Index ny) const
{
	//if it is not possible to do standard coarsening throw an exeption
	if (nx%2!=0 || ny%2!=0)
		throw std::domain_error("u");

	const Index nxNew=nx/2;
	const Index nyNew=ny/2;

	PositionArray Jx = prolongation_.getJx(C);
	PositionArray Jy = prolongation_.getJy(C);
	NumericArray result(0.0,(nxNew+1)*(nyNew+1));
	Precision scale = 0;
	//do injection on the borders
	for (Index sy=0; sy<=nyNew; sy++)
	{
		result[sy*(nxNew+1)]=u[2*sy*(nx+1)];
		result[sy*(nxNew+1)+nxNew]=u[2*sy*(nx+1)+nx];
    }
    for (Index sx=0; sx<=nxNew; sx++)
    {
		result[sx]=u[2*sx];
	    result[nyNew*(nxNew+1)+sx]=u[ny*(nx+1)+2*sx];
	}
	for (Index sy=1; sy<nyNew; sy++)
		for (Index sx=1; sx<nxNew; sx++)
		{
			Position pos = C;
			if (sy == 1)
			{
				if (sx == 1)
                    pos = SW;
				else if (sx == nxNew-1)
                    pos = SE;
				else
                    pos = S;
			}
			else if (sy == nyNew-1)
			{
				if (sx == 1)
                    pos = NW;
				else if (sx == nxNew-1)
                    pos = NE;
				else
                    pos = N;
			}
			else
			{
				if (sx == 1)
                    pos = W;
				else if (sx == nxNew-1)
                    pos = E;
				else
                    pos = C;
			}
			// Get the stencil of the adjoint prolongation
			NumericArray I = prolongation_.getI(pos, 2*sx, 2*sy, nx, ny, stencil);
			scale = 0;

			// Apply the stencil with regard to the Jx/Jy-vectors
			for (Index no=0; no<I.size(); no++)
			{
                //! \todo Why this condition?
				//if (I[no] > 0.000000001 || I[no] < -0.000000001)	
				result[sy*(nxNew+1)+sx] += I[no] * u[(2*sy+Jy[no])*(nx+1)+2*sx+Jx[no]];
				scale += I[no];
			}
			result[sy*(nxNew+1)+sx] *= weight_ / scale;
		}
	return result;
}

}

