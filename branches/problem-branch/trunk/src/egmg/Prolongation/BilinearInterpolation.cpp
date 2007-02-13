/** \file BilinearInterpolation.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class BilinearInterpolation
 */
#include "BilinearInterpolation.h"

namespace mg
{
    
NumericArray BilinearInterpolation::prolongate(
    const Problem& problem) const
{
    const Index nx = problem.getNx();
    const Index ny = problem.getNy();
    const Index nxNew = 2*nx;
    const Index nyNew = 2*ny;
    const NumericArray& u = problem.getSolution();
    NumericArray result(0.0,(nxNew+1)*(nxNew+1));

    //"interpolation" of coarse grid points
    for (Index sy=0; sy<=ny; ++sy)
        for (Index sx=0; sx<=nx; ++sx)
            result[2*sy*(nxNew+1)+2*sx]=u[sy*(nx+1)+sx];

    //interpolation of fine grid points on coarse grid lines
    for (Index sy=0; sy<=nyNew; sy+=2)
        for (Index sx=1; sx<=nxNew; sx+=2)
            result[sy*(nxNew+1)+sx]=
                1./2*(result[sy*(nxNew+1)+sx-1]
                     +result[sy*(nxNew+1)+sx+1]);
    //interpolation of fine grid points on fine grid lines
    for (Index sy=1; sy<=nyNew; sy+=2)
        for (Index sx=0; sx<=nxNew; ++sx)
            result[sy*(nxNew+1)+sx]=
                1./2*(result[(sy-1)*(nxNew+1)+sx]
                     +result[(sy+1)*(nxNew+1)+sx]);
    problem.applyBoundaryConstraint(result);
    return result;
}

}
