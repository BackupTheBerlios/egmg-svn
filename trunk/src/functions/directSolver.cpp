/** \file directSolver.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 *  <a href="mailto:kai.lists.berlios.de@tragetaschen.dyndns.org">Kai Ruhnau</a>
 * \brief directSolver.cpp contains the impl. of the func. directSolver.
 */
#include "directSolver.h"

namespace mg
{

namespace
{
    
void fillBoarderValues(
    NumericArray& matrixA,
    const Index dimA,
    const NumericArray& u,
    NumericArray& rightSide,
    const Index nx,
    const Index ny
    )
{
    //boarder values XDIR
    for (Index sx=0; sx<=nx; ++sx)
    {
        //lower boarder
        matrixA[sx*dimA+sx]=1;
        rightSide[sx]=u[sx];
        //upper boarder
        matrixA[(dimA-1-sx)*dimA+(dimA-1-sx)]=1;
        rightSide[dimA-1-sx]=u[dimA-sx];
    }
    //boarder values YDIR
    //corners have been process in XDIR (y=1..ny-1 instead of y=0..ny)
    for (Index sy=1; sy<ny; ++sy)
    {
        //left boarder
        matrixA[sy*(nx+1)*dimA+sy*(nx+1)]=1;
        rightSide[sy*(nx+1)]=u[sy*(nx+1)];
        //right boarder
        matrixA[(sy*(nx+1)+nx)*dimA+(sy*(nx+1)+nx)]=1;
        rightSide[(sy*(nx+1)+nx)]=u[(sy*(nx+1)+nx)];
    }
}

void pointFillMatrixA(
    NumericArray& matrixA,
    const Index dimA,
    const Stencil& stencil,
    const Position postion,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny
    )
{
    PositionArray jX=stencil.getJx(postion);
    PositionArray jY=stencil.getJy(postion);
    NumericArray operatorL=stencil.getL(postion,sx,sy,nx,ny);
    for (Index i=0; i<operatorL.size(); ++i)
        matrixA[(sy*(nx+1)+sx)*dimA+((sy+jY[i])*(nx+1)+sx+jX[i])]=operatorL[i];
}

void fillMatrixA(
    NumericArray& matrixA,
    const Index dimA,
    const Stencil& stencil,
    const Index nx,
    const Index ny
)
{
    //corner points
    pointFillMatrixA(matrixA,dimA,stencil,NW,1,ny-1,nx,ny);
    pointFillMatrixA(matrixA,dimA,stencil,NE,nx-1,ny-1,nx,ny);
    pointFillMatrixA(matrixA,dimA,stencil,SE,nx-1,1,nx,ny);
    pointFillMatrixA(matrixA,dimA,stencil,SW,1,1,nx,ny);
    //boarder points
    for (Index sy=2; sy<(ny-1); ++sy)
        pointFillMatrixA(matrixA,dimA,stencil,W,1,sy,nx,ny);
    for (Index sy=2; sy<(ny-1); ++sy)
        pointFillMatrixA(matrixA,dimA,stencil,E,nx-1,sy,nx,ny);
    for (Index sx=2; sx<(nx-1); ++sx)
        pointFillMatrixA(matrixA,dimA,stencil,N,sx,ny-1,nx,ny);
    for (Index sx=2; sx<(nx-1); ++sx)
        pointFillMatrixA(matrixA,dimA,stencil,S,sx,1,nx,ny);
    //center Points
    for (Index sy=2; sy<(ny-1); ++sy)
    {
        for (Index sx=2; sx<(nx-1); ++sx)
        {
            pointFillMatrixA(matrixA,dimA,stencil,C,sx,sy,nx,ny);
        }
    }
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
    const NumericArray& matrix,
    const NumericArray& b,
    const Index dimension)
{
    NumericArray lu=matrix;
    PositionArray permutation=pivotLU(lu,dimension);

    NumericArray result(dimension);
    //solve Ly=b
    for (Index sx=0; sx<dimension; ++sx)
    {
        result[sx]=b[permutation[sx]];
        for (Index sy=0; sy<sx; ++sy)
            result[sx]-=lu[permutation[sx]*dimension+sy]*result[sy];
    }
    //solve Rx=y
    for (Index sx=0; sx<dimension; ++sx)
    {
        for (Index sy=0; sy<sx; ++sy)
            result[dimension-1-sx]-=
                 lu[permutation[dimension-1-sx]*dimension+dimension-1-sy]
                *result[dimension-1-sy];
        result[dimension-1-sx]/=
            lu[(permutation[dimension-1-sx]+1)*dimension-1-sx];
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
    const Index dimA=u.size();
    NumericArray matrixA(0.0,dimA*dimA);
    NumericArray rightSide=f;
    
    fillBoarderValues(matrixA,dimA,u,rightSide,nx,ny);

    fillMatrixA(matrixA,dimA,stencil,nx,ny);
    
    u=solve(matrixA,rightSide,u.size());
}  
}
