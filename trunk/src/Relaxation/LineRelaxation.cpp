/** \file LineRelaxation.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief LineRelaxation.cpp contains the impl. of the class LineRelaxation.
 * \see LineRelaxation.h
 */
#include <iostream>
#include "LineRelaxation.h"
#include <stdexcept>
#include "../functions/printStencil.h"

namespace mg
{

void LineRelaxation::xLRSolver(
    NumericArray& solveInto,
    const Index sy,
    const Index nx,
    const Index /*ny*/,
    NumericArray& rhs,
    NumericArray& leftDiagonal,
    NumericArray& mainDiagonal,
    const NumericArray& rightDiagonal) const
{
    for(Index sx=1; sx<nx-1; sx++)
    {
        leftDiagonal[sx-1]=leftDiagonal[sx-1]/mainDiagonal[sx-1];
        mainDiagonal[sx]-=leftDiagonal[sx-1]*rightDiagonal[sx-1];
        rhs[sx]-=leftDiagonal[sx-1]*rhs[sx-1];
    }
    // solve the linear system of equations R x = rhs and store in solveInto
    solveInto[sy*(nx+1)+(nx-1)]=rhs[nx-2]/mainDiagonal[nx-2];
    for(Index sx=nx-2; sx>0; sx--)
    {
        solveInto[sy*(nx+1)+sx]=1/mainDiagonal[sx-1]*
                (rhs[sx-1]-rightDiagonal[sx-1]*solveInto[sy*(nx+1)+sx+1]);
    }
}

void LineRelaxation::yLRSolver(
    NumericArray& solveInto,
    const Index sx,
    const Index nx,
    const Index ny,
    NumericArray& rhs,
    NumericArray& leftDiagonal,
    NumericArray& mainDiagonal,
    const NumericArray& rightDiagonal) const
{
    // LR-decomposition + transformation of the rhs vector
    for(Index sy=1; sy<ny-1; sy++)
    {
        leftDiagonal[sy-1]=leftDiagonal[sy-1]/mainDiagonal[sy-1];
        mainDiagonal[sy]-=leftDiagonal[sy-1] * rightDiagonal[sy-1];
        rhs[sy]-=leftDiagonal[sy-1]*rhs[sy-1];
    }
    // solve the linear system of equations R solveInto = rhs
    solveInto[sx+(nx+1)*(ny-1)]=rhs[ny-2]/mainDiagonal[ny-2];
    for(Index sy=ny-2; sy>0; sy--)
    {
        solveInto[sx+sy*(nx+1)]=1/mainDiagonal[sy-1]*
                (rhs[sy-1]-rightDiagonal[sy-1]*solveInto[sx+(sy+1)*(nx+1)]);
    }
}

void LineRelaxation::xLRSolver(
    NumericArray& solveInto,
    const Index sy,
    const Index nx,
    const Index /*ny*/,
    NumericArray& rhs,
    NumericArray& leftDiagonal2,
    NumericArray& leftDiagonal1,
    NumericArray& mainDiagonal,
    NumericArray& rightDiagonal1,
    const NumericArray& rightDiagonal2) const
{
    // LR-decomposition + transformation of the rhs
    for(Index sx=1; sx<nx-2; sx++)
    {
        leftDiagonal1[sx-1]=leftDiagonal1[sx-1]/mainDiagonal[sx-1];
        mainDiagonal[sx]-=leftDiagonal1[sx-1]*rightDiagonal1[sx-1];
        rightDiagonal1[sx]-=leftDiagonal1[sx-1]*rightDiagonal2[sx-1];
        rhs[sx]-=leftDiagonal1[sx-1]*rhs[sx-1];
        leftDiagonal2[sx-1]=leftDiagonal2[sx-1]/mainDiagonal[sx-1];
        leftDiagonal1[sx]-=leftDiagonal2[sx-1]*rightDiagonal1[sx-1];
        mainDiagonal[sx+1]-=leftDiagonal2[sx-1]*rightDiagonal2[sx-1];
        rhs[sx+1]-=leftDiagonal2[sx-1]*rhs[sx-1];
    }
    leftDiagonal1[nx-2-1]=leftDiagonal1[nx-2-1]/mainDiagonal[nx-2-1];
    mainDiagonal[nx-2]-=leftDiagonal1[nx-2-1]*rightDiagonal1[nx-2-1];
    rhs[nx-2]-=leftDiagonal1[nx-2-1]*rhs[nx-2-1];

    // solve the linear system of equations R u = rhs
    solveInto[nx-1+sy*(nx+1)]=rhs[nx-2]/mainDiagonal[nx-2];
    for(Index sx=nx-2; sx>1; sx--)
    {
        solveInto[sx+sy*(nx+1)]=1/mainDiagonal[sx-1]*
            (rhs[sx-1]-rightDiagonal1[sx-1]*solveInto[sx+1+sy*(nx+1)]);
        rhs[sx-2]-=rightDiagonal2[sx-2]*solveInto[sx+1+sy*(nx+1)];
    }
    solveInto[1+sy*(nx+1)]=1/mainDiagonal[0]*
        (rhs[0]-rightDiagonal1[0]*solveInto[2+sy*(nx+1)]);
}

void LineRelaxation::yLRSolver(
    NumericArray& solveInto,
    const Index sx,
    const Index nx,
    const Index ny,
    NumericArray& rhs,
    NumericArray& leftDiagonal2,
    NumericArray& leftDiagonal1,
    NumericArray& mainDiagonal,
    NumericArray& rightDiagonal1,
    const NumericArray& rightDiagonal2) const
{
    // LR-decomposition + transformation of the rhs
    for(Index sy=1; sy<ny-2; sy++)
    {
        leftDiagonal1[sy-1]=leftDiagonal1[sy-1]/mainDiagonal[sy-1];
        mainDiagonal[sy]-=leftDiagonal1[sy-1]*rightDiagonal1[sy-1];
        rightDiagonal1[sy]-=leftDiagonal1[sy-1]*rightDiagonal2[sy-1];
        rhs[sy]-=leftDiagonal1[sy-1]*rhs[sy-1];
        leftDiagonal2[sy-1]=leftDiagonal2[sy-1]/mainDiagonal[sy-1];
        leftDiagonal1[sy]-=
            leftDiagonal2[sy-1]
            *rightDiagonal1[sy-1];
        mainDiagonal[sy+1]-=leftDiagonal2[sy-1]*rightDiagonal2[sy-1];
        rhs[sy+1]-=leftDiagonal2[sy-1]*rhs[sy-1];
    }
    leftDiagonal1[ny-2-1]=leftDiagonal1[ny-2-1]/mainDiagonal[ny-2-1];
    mainDiagonal[ny-2]-=leftDiagonal1[ny-2-1]*rightDiagonal1[ny-2-1];
    rhs[ny-2]-=leftDiagonal1[ny-2-1]*rhs[ny-2-1];

    // solve the linear system of equations R u = rhs
    solveInto[sx+(ny-1)*(nx+1)]=rhs[ny-2]/mainDiagonal[ny-2];
    for(Index sy=ny-2; sy>1; sy--)
    {
        solveInto[sx+sy*(nx+1)]=1/mainDiagonal[sy-1]*
            (rhs[sy-1]-rightDiagonal1[sy-1]*solveInto[sx+(sy+1)*(nx+1)]);
        rhs[sy-2]-=rightDiagonal2[sy-2]*solveInto[sx+(sy+1)*(nx+1)];
    }
    solveInto[sx+1*(nx+1)]=1/mainDiagonal[0]*
        (rhs[0]-rightDiagonal1[0]*solveInto[sx+2*(nx+1)]);
}

void LineRelaxation::solveNinePointXLine(
    NumericArray &u,
    const NumericArray &f,
    NumericArray &solveInto,
    const Stencil &stencil,
    const Index nx,
    const Index ny,
    const Index,
	const Index sy,
	const Position east,
	const Position center,
	const Position west) const
{
    NumericArray leftDiagonal(nx-2);
    NumericArray mainDiagonal(nx-1);
    NumericArray rightDiagonal(nx-2);
    NumericArray rhs(nx-1);

    NumericArray operatorL=stencil.getL(east,1,sy,nx,ny);
    PositionArray jX=stencil.getJx(C,nx,ny);
    PositionArray jY=stencil.getJy(C,nx,ny);

    mainDiagonal[0]=operatorL[C];
    rightDiagonal[0]=operatorL[E];

	rhs[0]=f[1+sy*(nx+1)]
		-operatorL[N]*u[1+(sy+jY[N])*(nx+1)]
		-operatorL[S]*u[1+(sy+jY[S])*(nx+1)]
		-operatorL[W]*u[0+ sy       *(nx+1)];
	for(Index i=5; i<operatorL.size(); i++)
	{
		rhs[0]-=operatorL[i]*u[1+jX[i]+(sy+jY[i])*(nx+1)];
	}

	for(Index sx=2; sx<nx-1; sx++)
	{
        operatorL=stencil.getL(center,sx,sy,nx,ny);
        leftDiagonal[sx-2]=operatorL[W];
        mainDiagonal[sx-1]=operatorL[C];
        rightDiagonal[sx-1]=operatorL[E];

		rhs[sx-1]=f[sx+sy*(nx+1)]
			-operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
			-operatorL[S]*u[sx+(sy+jY[S])*(nx+1)];
		for(Index i=5; i<operatorL.size(); i++)
		{
			rhs[sx-1]-=operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
		}
	}

    operatorL=stencil.getL(west,nx-1,sy,nx,ny);
    leftDiagonal[nx-3]=operatorL[W];
    mainDiagonal[nx-2]=operatorL[C];

	rhs[nx-2]=f[nx-1+sy*(nx+1)]
		-operatorL[N]*u[nx-1+(sy+jY[N])*(nx+1)]
		-operatorL[S]*u[nx-1+(sy+jY[S])*(nx+1)]
		-operatorL[E]*u[nx  + sy       *(nx+1)];
	for(Index i=5; i<operatorL.size(); i++)
	{
		rhs[nx-2]-=operatorL[i]*u[nx-1+jX[i]+(sy+jY[i])*(nx+1)];
	}

	xLRSolver(solveInto,sy,nx,ny,rhs,leftDiagonal,mainDiagonal,rightDiagonal);
}

void LineRelaxation::solveNinePointYLine(
    NumericArray &u,
    const NumericArray &f,
    NumericArray &solveInto,
    const Stencil &stencil,
    const Index nx,
    const Index ny,
    const Index sx,
	const Index,
	const Position south,
	const Position center,
	const Position north) const
{
    NumericArray leftDiagonal(ny-2);
    NumericArray mainDiagonal(ny-1);
    NumericArray rightDiagonal(ny-2);
    NumericArray rhs(nx-1);

    NumericArray operatorL=stencil.getL(south,sx,1,nx,ny);
    PositionArray jX = stencil.getJx(C,nx,ny);
    PositionArray jY = stencil.getJy(C,nx,ny);

    mainDiagonal[0]=operatorL[C];
    rightDiagonal[0]=operatorL[N];

    rhs[0]=f[sx+(nx+1)]
        -operatorL[W]*u[sx+jX[W]+          (nx+1)]
        -operatorL[E]*u[sx+jX[E]+          (nx+1)]
        -operatorL[S]*u[sx+      (1+jY[S])*(nx+1)];
    for(Index i=5; i<operatorL.size(); i++)
    {
        rhs[0]-=operatorL[i]*u[sx+jX[i]+(1+jY[i])*(nx+1)];
    }

    for(Index sy=2; sy<ny-1; sy++)
    {
        operatorL=stencil.getL(center,sx,sy,nx,ny);
        leftDiagonal[sy-2]=operatorL[S];
        mainDiagonal[sy-1]=operatorL[C];
        rightDiagonal[sy-1]=operatorL[N];

        rhs[sy-1]=f[sx+sy*(nx+1)]
            -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
            -operatorL[E]*u[sx+jX[E]+sy*(nx+1)];
        for(Index i=5; i<operatorL.size(); i++)
        {
            rhs[sy-1]-=operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
        }
    }

    operatorL=stencil.getL(north,sx,ny-1,nx,ny);
    leftDiagonal[ny-3] = operatorL[S];
    mainDiagonal[ny-2] = operatorL[C];

    rhs[ny-2]=f[sx+(ny-1)*(nx+1)]
        -operatorL[W]*u[sx+jX[W]+(ny-1)*(nx+1)]
        -operatorL[E]*u[sx+jX[E]+(ny-1)*(nx+1)]
        -operatorL[N]*u[sx+(ny-1+jY[N])*(nx+1)];
    for(Index i=5; i<operatorL.size(); i++)
    {
        rhs[ny-2]-=operatorL[i]*u[sx+jX[i]+(ny-1+jY[i])*(nx+1)];
    }

    yLRSolver(solveInto,sx,nx,ny,rhs,leftDiagonal,mainDiagonal,rightDiagonal);
}

void LineRelaxation::solveFullXLine(
    NumericArray &u,
    const NumericArray &f,
    NumericArray &solveInto,
    const Stencil &stencil,
    const Index nx,
    const Index ny,
    const Index,
	const Index sy,
	const Position west,
	const Position center,
	const Position east) const
{
    NumericArray operatorLB=stencil.getL(west,1,sy,nx,ny);
    PositionArray jXB=stencil.getJx(west,nx,ny);
    PositionArray jYB=stencil.getJy(west,nx,ny);

    Index start=9;
    if (center!=C)
        start=8;

    NumericArray leftDiagonal2(0.0,nx-3);
    NumericArray leftDiagonal1(0.0,nx-2);
    NumericArray mainDiagonal(0.0,nx-1);
    NumericArray rightDiagonal1(0.0,nx-2);
    NumericArray rightDiagonal2(0.0,nx-3);
    NumericArray rhs(nx-1);

    mainDiagonal[0]=operatorLB[C];
    rightDiagonal1[0]=operatorLB[E];

    rhs[0]=f[1+sy*(nx+1)]
        -operatorLB[N]*u[1+(sy+jYB[N])*(nx+1)]
        -operatorLB[S]*u[1+(sy+jYB[S])*(nx+1)]
        -operatorLB[W]*u[0+ sy        *(nx+1)];

    switch (center)
    {
    case S:
        rightDiagonal2[0]=operatorLB[6];
        rhs[0]-=operatorLB[5]*u[1+(sy+jYB[5])*(nx+1)];
        break;
    case C:
        rightDiagonal2[0]=operatorLB[6];
        rhs[0]-=operatorLB[5]*u[1+(sy+jYB[5])*(nx+1)];
        rhs[0]-=operatorLB[7]*u[1+(sy+jYB[7])*(nx+1)];
        break;
    case N:
        rightDiagonal2[0]=operatorLB[5];
        rhs[0]-=operatorLB[6]*u[1+(sy+jYB[6])*(nx+1)];
        break;
    default:
        throw std::logic_error("out of range");
    }

    for(Index i=start-1; i<operatorLB.size(); i++)
    {
        rhs[0]-=operatorLB[i]*u[1+jXB[i]+(sy+jYB[i])*(nx+1)];
    }

    NumericArray operatorL=stencil.getL(center,2,sy,nx,ny);
    PositionArray jX=stencil.getJx(center,nx,ny);
    PositionArray jY=stencil.getJy(center,nx,ny);

    leftDiagonal1[0]=operatorL[W];
    mainDiagonal[1]=operatorL[C];
    rightDiagonal1[1]=operatorL[E];

    rhs[1]=f[2+sy*(nx+1)]
        -operatorL[N]*u[2+(sy+jY[N])*(nx+1)]
        -operatorL[S]*u[2+(sy+jY[S])*(nx+1)];
    switch (center)
    {
    case S:
        rightDiagonal2[1]=operatorL[7];
        rhs[1]-=operatorL[5]*u[2+jX[5]+sy       *(nx+1)];
        rhs[1]-=operatorL[6]*u[2+     (sy+jY[6])*(nx+1)];
        break;
    case C:
        rightDiagonal2[1]=operatorL[7];
        rhs[1]-=operatorL[5]*u[2+jX[5]+sy       *(nx+1)];
        rhs[1]-=operatorL[6]*u[2+     (sy+jY[6])*(nx+1)];
        rhs[1]-=operatorL[8]*u[2+     (sy+jY[8])*(nx+1)];
        break;
    case N:
        rightDiagonal2[1]=operatorL[6];
        rhs[1]-=operatorL[5]*u[2+jX[5]+sy       *(nx+1)];
        rhs[1]-=operatorL[7]*u[2+     (sy+jY[7])*(nx+1)];
        break;
    default:
        throw std::logic_error("out of range");
    }

    for(Index i=start; i<operatorL.size(); i++)
    {
        rhs[1]-=operatorL[i]*u[2+jX[i]+(sy+jY[i])*(nx+1)];
    }

    for(Index sx=3; sx<nx-2; sx++)
    {
        operatorL=stencil.getL(center,sx,sy,nx,ny);

        leftDiagonal1[sx-2]=operatorL[W];
        mainDiagonal[sx-1]=operatorL[C];
        rightDiagonal1[sx-1]=operatorL[E];

        rhs[sx-1]=f[sx+sy*(nx+1)]
            -operatorL[N]*u[sx+(sy+jY[N])*(nx+1)]
            -operatorL[S]*u[sx+(sy+jY[S])*(nx+1)];

        switch (center)
        {
        case S:
            leftDiagonal2[sx-3]=operatorL[5];
            rightDiagonal2[sx-1]=operatorL[7];
            rhs[sx-1]-=operatorL[6]*u[sx+(sy+jY[6])*(nx+1)];
            break;
        case C:
            leftDiagonal2[sx-3]=operatorL[5];
            rightDiagonal2[sx-1]=operatorL[7];
            rhs[sx-1]-=operatorL[6]*u[sx+(sy+jY[6])*(nx+1)];
            rhs[sx-1]-=operatorL[8]*u[sx+(sy+jY[8])*(nx+1)];
            break;
        case N:
            leftDiagonal2[sx-3]=operatorL[5];
            rightDiagonal2[sx-1]=operatorL[6];
            rhs[sx-1]-=operatorL[7]*u[sx+(sy+jY[7])*(nx+1)];
            break;
        default:
            throw std::logic_error("out of range");
        }

        for(Index i=start; i<operatorL.size(); i++)
        {
            rhs[sx-1]-=operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
        }
    }

    operatorL=stencil.getL(center,nx-2,sy,nx,ny);

    leftDiagonal1[nx-4]=operatorL[W];
    mainDiagonal[nx-3]=operatorL[C];
    rightDiagonal1[nx-3]=operatorL[E];

    rhs[nx-3]=f[nx-2+sy*(nx+1)]
        -operatorL[N]*u[nx-2+(sy+jY[N])*(nx+1)]
        -operatorL[S]*u[nx-2+(sy+jY[S])*(nx+1)];

    switch (center)
    {
    case S:
        leftDiagonal2[nx-5]=operatorL[5];
        rhs[nx-3]-=operatorL[6]*u[nx-2+     (sy+jY[6])*(nx+1)];
        rhs[nx-3]-=operatorL[7]*u[nx-2+jX[7]+sy       *(nx+1)];
        break;
    case C:
        leftDiagonal2[nx-5]=operatorL[5];
        rhs[nx-3]-=operatorL[6]*u[nx-2+     (sy+jY[6])*(nx+1)];
        rhs[nx-3]-=operatorL[7]*u[nx-2+jX[7]+sy       *(nx+1)];
        rhs[nx-3]-=operatorL[8]*u[nx-2+     (sy+jY[8])*(nx+1)];
        break;
    case N:
        leftDiagonal2[nx-5]=operatorL[5];
        rhs[nx-3]-=operatorL[6]*u[nx-2+jX[6]+sy      *(nx+1)];
        rhs[nx-3]-=operatorL[7]*u[nx-2+    (sy+jY[7])*(nx+1)];
        break;
    default:
        throw std::logic_error("out of range");
    }

    for(Index i=start; i<operatorL.size(); i++)
    {
        rhs[nx-3]-=operatorL[i]*u[nx-2+jX[i]+(sy+jY[i])*(nx+1)];
    }

    operatorLB=stencil.getL(east,nx-1,sy,nx,ny);
    jXB=stencil.getJx(east,nx,ny);
    jYB=stencil.getJy(east,nx,ny);

    leftDiagonal1[nx-3]=operatorLB[W];
    mainDiagonal[nx-2]=operatorLB[C];

    rhs[nx-2]=f[nx-1+sy*(nx+1)]
        -operatorLB[N]*u[nx-1+(sy+jYB[N])*(nx+1)]
        -operatorLB[S]*u[nx-1+(sy+jYB[S])*(nx+1)]
        -operatorLB[E]*u[nx+sy*(nx+1)];

    switch (center)
    {
    case S:
        leftDiagonal2[nx-4]=operatorLB[5];
        rhs[nx-2]-=operatorLB[6]*u[nx-1+(sy+jYB[6])*(nx+1)];
        break;
    case C:
        leftDiagonal2[nx-4]=operatorLB[5];
        rhs[nx-2]-=operatorLB[6]*u[nx-1+(sy+jYB[6])*(nx+1)];
        rhs[nx-2]-=operatorLB[7]*u[nx-1+(sy+jYB[7])*(nx+1)];
        break;
    case N:
        leftDiagonal2[nx-4]=operatorLB[5];
        rhs[nx-2]-=operatorLB[6]*u[nx-1+(sy+jYB[6])*(nx+1)];
        break;
    default:
        throw std::logic_error("out of range");
    }

    for(Index i=start-1; i<operatorLB.size(); i++)
    {
        rhs[nx-2]-=operatorLB[i]*u[nx-1+jXB[i]+(sy+jYB[i])*(nx+1)];
    }

    xLRSolver(solveInto,sy,nx,ny,rhs,
        leftDiagonal2,leftDiagonal1,mainDiagonal,rightDiagonal1,rightDiagonal2);
}

void LineRelaxation::solveFullYLine(
    NumericArray &u,
    const NumericArray &f,
    NumericArray &solveInto,
    const Stencil &stencil,
    const Index nx,
    const Index ny,
    const Index sx,
	const Index,
	const Position south,
	const Position center,
	const Position north) const
{
    NumericArray operatorLB=stencil.getL(south,sx,1,nx,ny);
    PositionArray jXB=stencil.getJx(south,nx,ny);
    PositionArray jYB=stencil.getJy(south,nx,ny);

    Index start=9;
    if (center!=C)
        start=8;

    NumericArray leftDiagonal2(0.0,ny-3);
    NumericArray leftDiagonal1(0.0,ny-2);
    NumericArray mainDiagonal(0.0,ny-1);
    NumericArray rightDiagonal1(0.0,ny-2);
    NumericArray rightDiagonal2(0.0,ny-3);
    NumericArray rhs(ny-1);

    mainDiagonal[0]=operatorLB[C];
    rightDiagonal1[0]=operatorLB[N];

    rhs[0]=f[sx+1*(nx+1)]
        -operatorLB[S]*u[sx+(1+jYB[S])*(nx+1)]
        -operatorLB[W]*u[sx-1+          nx+1 ]
        -operatorLB[E]*u[sx+1+          nx+1 ];

    switch (center)
    {
    case W:
        rightDiagonal2[0]=operatorLB[5];
        rhs[0]-=operatorLB[6]*u[sx+jXB[6]+1*(nx+1)];
        break;
    case C:
        rightDiagonal2[0]=operatorLB[6];
        rhs[0]-=operatorLB[5]*u[sx+jXB[6]+1*(nx+1)];
        rhs[0]-=operatorLB[7]*u[sx+jXB[7]+1*(nx+1)];
        break;
    case E:
        rightDiagonal2[0]=operatorLB[6];
        rhs[0]-=operatorLB[5]*u[sx+jXB[5]+1*(nx+1)];
        break;
    default:
        throw std::logic_error("out of range");
    }

    for(Index i=start-1; i<operatorLB.size(); i++)
    {
        rhs[0]-=operatorLB[i]*u[sx+jXB[i]+(1+jYB[i])*(nx+1)];
    }

    NumericArray operatorL=stencil.getL(center,sx,2,nx,ny);
    PositionArray jX=stencil.getJx(center,nx,ny);
    PositionArray jY=stencil.getJy(center,nx,ny);

    leftDiagonal1[0]=operatorL[S];
    mainDiagonal[1]=operatorL[C];
    rightDiagonal1[1]=operatorL[N];

    rhs[1]=f[sx+2*(nx+1)]
        -operatorL[W]*u[sx+jX[W]+2*(nx+1)]
        -operatorL[E]*u[sx+jX[E]+2*(nx+1)];
    switch (center)
    {
    case W:
        rightDiagonal2[1]=operatorL[5];
        rhs[1]-=operatorL[6]*u[sx+jX[6]+2       *(nx+1)];
        rhs[1]-=operatorL[7]*u[sx+     (2+jY[7])*(nx+1)];
        break;
    case C:
        rightDiagonal2[1]=operatorL[6];
        rhs[1]-=operatorL[5]*u[sx+jX[5]+2       *(nx+1)];
        rhs[1]-=operatorL[7]*u[sx+jX[7]+2       *(nx+1)];
        rhs[1]-=operatorL[8]*u[sx+     (2+jY[8])*(nx+1)];
        break;
    case E:
        rightDiagonal2[1]=operatorL[6];
        rhs[1]-=operatorL[5]*u[sx+jX[5]+2       *(nx+1)];
        rhs[1]-=operatorL[7]*u[sx+     (2+jY[7])*(nx+1)];
        break;
    default:
        throw std::logic_error("out of range");
    }

    for(Index i=start; i<operatorL.size(); i++)
    {
        rhs[1]-=operatorL[i]*u[sx+jX[i]+(2+jY[i])*(nx+1)];
    }

    for(Index sy=3; sy<ny-2; sy++)
    {
        operatorL=stencil.getL(center,sx,sy,nx,ny);

        leftDiagonal1[sy-2]=operatorL[S];
        mainDiagonal[sy-1]=operatorL[C];
        rightDiagonal1[sy-1]=operatorL[N];

        rhs[sy-1]=f[sx+sy*(nx+1)]
            -operatorL[W]*u[sx+jX[W]+sy*(nx+1)]
            -operatorL[E]*u[sx+jX[E]+sy*(nx+1)];

        switch (center)
        {
        case W:
            leftDiagonal2[sy-3]=operatorL[7];
            rightDiagonal2[sy-1]=operatorL[5];
            rhs[sy-1]-=operatorL[6]*u[sx+jX[6]+sy*(nx+1)];
            break;
        case C:
            leftDiagonal2[sy-3]=operatorL[8];
            rightDiagonal2[sy-1]=operatorL[6];
            rhs[sy-1]-=operatorL[5]*u[sx+jX[5]+sy*(nx+1)];
            rhs[sy-1]-=operatorL[7]*u[sx+jX[7]+sy*(nx+1)];
            break;
        case E:
            leftDiagonal2[sy-3]=operatorL[7];
            rightDiagonal2[sy-1]=operatorL[6];
            rhs[sy-1]-=operatorL[5]*u[sx+jX[5]+sy*(nx+1)];
            break;
        default:
            throw std::logic_error("out of range");
        }

        for(Index i=start; i<operatorL.size(); i++)
        {
            rhs[sy-1]-=operatorL[i]*u[sx+jX[i]+(sy+jY[i])*(nx+1)];
        }
    }

    operatorL=stencil.getL(center,sx, ny-2,nx,ny);

    leftDiagonal1[ny-4]=operatorL[S];
    mainDiagonal[ny-3]=operatorL[C];
    rightDiagonal1[ny-3]=operatorL[N];

    rhs[ny-3]=f[sx+(ny-2)*(nx+1)]
        -operatorL[W]*u[sx+jX[W]+(ny-2)*(nx+1)]
        -operatorL[E]*u[sx+jX[E]+(ny-2)*(nx+1)];

    switch (center)
    {
    case W:
        leftDiagonal2[ny-5]=operatorL[7];
        rhs[ny-3]-=operatorL[5]*u[sx+      (ny-2+jY[5])*(nx+1)];
        rhs[ny-3]-=operatorL[6]*u[sx+jX[6]+(ny-2)      *(nx+1)];
        break;
    case C:
        leftDiagonal2[ny-5]=operatorL[8];
        rhs[ny-3]-=operatorL[5]*u[sx+jX[5]+(ny-2)      *(nx+1)];
        rhs[ny-3]-=operatorL[6]*u[sx+      (ny-2+jY[6])*(nx+1)];
        rhs[ny-3]-=operatorL[7]*u[sx+jX[7]+(ny-2)      *(nx+1)];
        break;
    case E:
        leftDiagonal2[ny-5]=operatorL[7];
        rhs[ny-3]-=operatorL[5]*u[sx+jX[5]+(ny-2)      *(nx+1)];
        rhs[ny-3]-=operatorL[6]*u[sx+      (ny-2+jY[6])*(nx+1)];
        break;
    default:
        throw std::logic_error("out of range");
    }

    for(Index i=start; i<operatorL.size(); i++)
    {
        rhs[ny-3]-=operatorL[i]*u[sx+jX[i]+(ny-2+jY[i])*(nx+1)];
    }

    operatorLB=stencil.getL(north,sx,ny-1,nx,ny);
    jXB=stencil.getJx(north,nx,ny);
    jYB=stencil.getJy(north,nx,ny);

    leftDiagonal1[ny-3]=operatorLB[S];
    mainDiagonal[ny-2]=operatorLB[C];

    rhs[ny-2]=f[sx+(ny-1)*(nx+1)]
        -operatorLB[N]*u[sx+        ny   *(nx+1)]
        -operatorLB[W]*u[sx+jXB[W]+(ny-1)*(nx+1)]
        -operatorLB[E]*u[sx+jXB[E]+(ny-1)*(nx+1)];

    switch (center)
    {
    case W:
        leftDiagonal2[ny-4]=operatorLB[6];
        rhs[ny-2]-=operatorLB[5]*u[sx+jX[5]+(ny-1)*(nx+1)];
        break;
    case C:
        leftDiagonal2[ny-4]=operatorLB[7];
        rhs[ny-2]-=operatorLB[5]*u[sx+jX[5]+(ny-1)*(nx+1)];
        rhs[ny-2]-=operatorLB[6]*u[sx+jX[6]+(ny-1)*(nx+1)];
        break;
    case E:
        leftDiagonal2[ny-4]=operatorLB[6];
        rhs[ny-2]-=operatorLB[5]*u[sx+jX[5]+(ny-1)*(nx+1)];
        break;
    default:
        throw std::logic_error("out of range");
    }

    for(Index i=start-1; i<operatorLB.size(); i++)
    {
        rhs[ny-2]-=operatorLB[i]*u[sx+jXB[i]+(ny-1+jYB[i])*(nx+1)];
    }

    yLRSolver(solveInto,sx,nx,ny,rhs,
        leftDiagonal2,leftDiagonal1,mainDiagonal,rightDiagonal1,rightDiagonal2);
}

void LineRelaxation::relax(
    NumericArray &u,
    const NumericArray &f,
    const Stencil &stencil,
    const Index nx,
    const Index ny) const
{
    switch (stencil.size())
    {
    case 1:
        switch (direction_)
        {
        case ALTDIR:
            ninePointX(u,f,stencil,nx,ny);
            ninePointY(u,f,stencil,nx,ny);
            break;
        case XDIR:
            ninePointX(u,f,stencil,nx,ny);
            break;
        case YDIR:
            ninePointY(u,f,stencil,nx,ny);
            break;
        default:
            std::cerr<<"Error in direction of the line relaxation!\n";
            break;
        }
        break;
    case 2:
        switch (direction_)
        {
        case ALTDIR:
            fullX(u,f,stencil,nx,ny);
            fullY(u,f,stencil,nx,ny);
            break;
        case XDIR:
            fullX(u,f,stencil,nx,ny);
            break;
        case YDIR:
            fullY(u,f,stencil,nx,ny);
            break;
        default:
            std::cerr<<"Error in direction of the line relaxation!\n";
            break;
        }
        break;
    default:
        std::cerr << "Stencil is too big (size>2)!" << std::endl;
        break;
    }
}

}
