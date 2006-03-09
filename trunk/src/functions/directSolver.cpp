/** \file directSolver.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 *  <a href="mailto:kai.lists.berlios.de@tragetaschen.dyndns.org">Kai Ruhnau</a>
 * \brief directSolver.cpp contains the impl. of the func. directSolver.
 */

#include <algorithm>
#include "directSolver.h"

namespace mg
{

namespace
{
    
void fillBorderValues(
    NumericArray& matrix,
    const Index dimension,
    const NumericArray& u,
    NumericArray& rightSide,
    const Index nx,
    const Index ny)
{
    //border values XDIR
    for (Index sx=0; sx<=nx; ++sx)
    {
        //lower border
        matrix[sx*dimension+sx]=1;
        rightSide[sx]=u[sx];
        //upper border
        matrix[(dimension-1-sx)*dimension+(dimension-1-sx)]=1;
        rightSide[dimension-1-sx]=u[dimension-sx];
    }
    //border values YDIR
    //corners have been process in XDIR (y=1..ny-1 instead of y=0..ny)
    for (Index sy=1; sy<ny; ++sy)
    {
        //left border
        matrix[sy*(nx+1)*dimension+sy*(nx+1)]=1;
        rightSide[sy*(nx+1)]=u[sy*(nx+1)];
        //right border
        matrix[(sy*(nx+1)+nx)*dimension+(sy*(nx+1)+nx)]=1;
        rightSide[(sy*(nx+1)+nx)]=u[(sy*(nx+1)+nx)];
    }
}

void pointFillMatrix(
    NumericArray& matrix,
    const Index dimension,
    const Stencil& stencil,
    const Position postion,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny)
{
    PositionArray jX=stencil.getJx(postion);
    PositionArray jY=stencil.getJy(postion);
    NumericArray operatorL=stencil.getL(postion,sx,sy,nx,ny);
    for (Index i=0; i<operatorL.size(); ++i)
        matrix[(sy*(nx+1)+sx)*dimension+((sy+jY[i])*(nx+1)+sx+jX[i])]=operatorL[i];
}

void fillMatrix(
    NumericArray& matrix,
    const Index dimension,
    const Stencil& stencil,
    const Index nx,
    const Index ny)
{
    //corner points
    pointFillMatrix(matrix,dimension,stencil,NW,1,ny-1,nx,ny);
    pointFillMatrix(matrix,dimension,stencil,NE,nx-1,ny-1,nx,ny);
    pointFillMatrix(matrix,dimension,stencil,SE,nx-1,1,nx,ny);
    pointFillMatrix(matrix,dimension,stencil,SW,1,1,nx,ny);
    //border points
    for (Index sy=2; sy<(ny-1); ++sy)
        pointFillMatrix(matrix,dimension,stencil,W,1,sy,nx,ny);
    for (Index sy=2; sy<(ny-1); ++sy)
        pointFillMatrix(matrix,dimension,stencil,E,nx-1,sy,nx,ny);
    for (Index sx=2; sx<(nx-1); ++sx)
        pointFillMatrix(matrix,dimension,stencil,N,sx,ny-1,nx,ny);
    for (Index sx=2; sx<(nx-1); ++sx)
        pointFillMatrix(matrix,dimension,stencil,S,sx,1,nx,ny);
    //center Points
    for (Index sy=2; sy<(ny-1); ++sy)
        for (Index sx=2; sx<(nx-1); ++sx)
            pointFillMatrix(matrix,dimension,stencil,C,sx,sy,nx,ny);
}

PositionArray pivotLU(
    NumericArray& matrix,
    const Index dimension)
{
    PositionArray permutation(dimension);
    for (Index sx=0; sx<dimension; ++sx)
        permutation[sx]=sx;

    for (Index sx=0; sx<dimension; ++sx)
    {
        Index max_zeile=sx;

        for (Index sy=sx; sy<dimension; ++sy)
            if (std::abs(matrix[permutation[sy       ]*dimension+sx])
               >std::abs(matrix[permutation[max_zeile]*dimension+sx]))
                max_zeile=sy;

        std::swap(permutation[sx],permutation[max_zeile]);

        for (Index sy=sx+1; sy<dimension; ++sy)
        {
            matrix[permutation[sy]*dimension+sx]/=
                    matrix[permutation[sx]*dimension+sx];
            Precision pivot=matrix[permutation[sy]*dimension+sx];

            for (Index k=sx+1; k<dimension; k++)
                matrix[permutation[sy]*dimension+k]-=
                     pivot
                    *matrix[permutation[sx]*dimension+k];
        }
    }
    return permutation;
}

NumericArray solve(
    NumericArray& matrix,
    const NumericArray& b,
    const Index dimension)
{
    PositionArray permutation=pivotLU(matrix,dimension);

    NumericArray result(dimension);
    //solve Ly=b
    for (Index sx=0; sx<dimension; ++sx)
    {
        result[sx]=b[permutation[sx]];
        for (Index sy=0; sy<sx; ++sy)
            result[sx]-=matrix[permutation[sx]*dimension+sy]*result[sy];
    }
    //solve Rx=y
    for (Index sx=0; sx<dimension; ++sx)
    {
        for (Index sy=0; sy<sx; ++sy)
            result[dimension-1-sx]-=
                 matrix[permutation[dimension-1-sx]*dimension+dimension-1-sy]
                *result[dimension-1-sy];
        result[dimension-1-sx]/=
            matrix[(permutation[dimension-1-sx]+1)*dimension-1-sx];
    }
    return result;
}

}

void directSolver(
    NumericArray& u,
    const NumericArray& f,
    const Stencil& stencil,
    const Index nx,
    const Index ny)
{
    const Index dimension=u.size();
    NumericArray matrix(0.0,dimension*dimension);
    NumericArray rightSide=f;
    
    fillBorderValues(matrix,dimension,u,rightSide,nx,ny);

    fillMatrix(matrix,dimension,stencil,nx,ny);
    
    u=solve(matrix,rightSide,u.size());
}  
}
