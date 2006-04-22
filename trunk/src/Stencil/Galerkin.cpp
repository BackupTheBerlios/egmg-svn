/** \file Galerkin.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class Galerkin.
 */

#include <cstdlib>
//#include "Galerkin.h"
#include "Stencil.h"
#include "../general/parameters.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"
#include <stack>

namespace mg
{

namespace
{
    
void GeneratePositionArrays(
    PositionArray& jX,
    PositionArray& jY,
    const Index size,
    const Position pos)
{
    //ASSERT( pos == c )
    //ASSERT( (2*size+1)*(2*size+1) == jX.size() && jX.size() == jY.size() );
    std::stack<Integer> xPositions;
    std::stack<Integer> yPositions;
    
    Integer intSize = size;
    for ( Integer k=intSize; k>=1; --k)
    {
        xPositions.push(k);
        xPositions.push(k);
        yPositions.push(k);
        yPositions.push(k);
    }
    jX[0]=0;
    jY[0]=0;
    for ( Index i=1; i<4*size+1; ++i )
    {
        if ( i%4 == 0 )
        {
            jX[i] = 0;
            jY[i] = -1*yPositions.top();
            yPositions.pop();
        }
        else if ( i%4 == 1 )
        {
            jX[i] = -1*xPositions.top();
            xPositions.pop();
            jY[i] = 0;
        }
        else if ( i%4 == 2 )
        {
            jX[i] = 0;
            jY[i] = yPositions.top();
            yPositions.pop();
        }
        else if ( i%4 == 3 )
        {
            jX[i] = xPositions.top();
            xPositions.pop();
            jY[i] = 0;
        }
    }
    //ASSERT( xPositions.empty() && yPositions.empty() );
    for ( Integer k=-intSize; k<=-1; ++k )
    {
        xPositions.push(k);
        yPositions.push(k);
    }
    for ( Integer k=1; k<=intSize; ++k )
    {
        xPositions.push(k);
        yPositions.push(k);
    }
    Index i=4*size+1;
    while( i<jX.size() )
    {
        std::stack<Integer> xP = xPositions;
        Integer yPosition = yPositions.top();
        yPositions.pop();
        while( !xP.empty() )
        {
            jX[i] = xP.top();
            xP.pop();
            jY[i] = yPosition;
            ++i;
            if ( i>=jX.size() )
                break;
        }
        //ASSERT( !yPositions.empty() );    
    }
}

Index ComputeSize(
    const Restriction& restriction,
    const Stencil& stencil)
{
    PositionArray restrictionJx = restriction.getJx();
    PositionArray restrictionJy = restriction.getJy();
    Index restSize = std::max( restrictionJx.expansion(),
                               restrictionJy.expansion() );
    return restSize + stencil.size();
}

Index ComputeSize(
    const PositionArray& jX,
    const PositionArray& jY,
    const Prolongation& prolong)
{
    PositionArray prolongJx = prolong.getJx();
    PositionArray prolongJy = prolong.getJy();
    const Index prolongSize = std::max( prolongJx.expansion(),
                                        prolongJy.expansion() );
    const Index result = std::max( jX.expansion(), jY.expansion() );
    if ( result%2 == 0 )
    {
        return result/2 + prolongSize - 1;
    }
    return ( result - 1 )/2 + prolongSize;
}

void RestritionTimesStencil(
    NumericArray& resultL,
    PositionArray& jX,
    PositionArray& jY,
    const Position pos,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny,
    const Restriction& restriction,
    const Stencil& stencil)
{
    //ASSERT( pos == C )
    const PositionArray restrictionJx = restriction.getJx();
    const PositionArray restriciionJy = restriction.getJy();
    const NumericArray restrictionI = restriction.getI(sx,sy,nx,ny,stencil);

    const Index size = ComputeSize(restriction,stencil);
    
    jX.resize((2*size+1)*(2*size+1));
    jY.resize((2*size+1)*(2*size+1));
    jX=0;
    jY=0;
    GeneratePositionArrays(jX,jY,size,C);
    resultL.resize((2*size+1)*(2*size+1));
    resultL = 0.0;

    for ( Index i=0; i<restrictionI.size(); ++i )
    {
        PositionArray stencilJx = stencil.getJx( pos );
        PositionArray stencilJy = stencil.getJy( pos );
        NumericArray stencilL = stencil.getL(
            pos,
            sx+restrictionJx[i],
            sy+restriciionJy[i],
            nx, ny);
        stencilJx+=restrictionJx[i];
        stencilJy+=restriciionJy[i];
        stencilL*=restrictionI[i];
        for ( Index j=0; j<resultL.size(); ++j )
        {
            for ( Index k=0; k<stencilL.size(); ++k)
            {
                if ( stencilJx[k] == jX[j] && stencilJy[k] == jY[j] )
                {
                    resultL[j] += stencilL[k];
                }
            }
        }
    }
}

void ProlongateStencil(
    NumericArray& resultL,
    PositionArray& resultJx,
    PositionArray& resultJy,
    const NumericArray& opL,
    const PositionArray& jX,
    const PositionArray& jY,
    const Position pos,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny,
    const Prolongation& prolongation,
    const Stencil& stencil)
{
    const PositionArray prolongJx = prolongation.getJx();
    const PositionArray prolongJy = prolongation.getJy();
    const NumericArray prolongI = prolongation.getI(sx,sy,nx,ny,stencil);
    const Index size = ComputeSize(jX,jY,prolongation);
    
    resultJx.resize((2*size+1)*(2*size+1));
    resultJy.resize((2*size+1)*(2*size+1));
    resultJx=0;
    resultJy=0;
    GeneratePositionArrays(resultJx,resultJy,size,C);
    resultL.resize((2*size+1)*(2*size+1));
    resultL = 0.0;
    
    for (Index i=0; i<opL.size(); ++i)
    {
        if (jX[i]%2 == 0 && jY[i]%2 == 0) //coars grid point
        {
            for (Index j=0; j<resultL.size(); ++j)
            {
                if (resultJx[j] == jX[i]/2 && resultJy[j] == jY[i]/2)
                {
                    resultL[j]+=opL[i];
                    break;
                }
            }
        }
        else if (jX[i]%2 == 0) //fine grid point on x-coars grid line
        {
            for ( Index j=0; j<prolongI.size(); ++j )
            {
                if ( prolongJy[j] == 0 )
                {
                    for ( Index k=0; k<resultL.size(); ++k )
                    {
                        if ( resultJx[k] == jX[i] + prolongJx[j] &&
                             resultJy[k] == jY[i] + prolongJy[j] )
                        {
                            resultL[k]+=opL[i]*prolongI[j];
                        }
                    }
                }
            }
        }
        else if (jY[i]%2 == 0) //fine grid point on y-coars grid line
        {
            for ( Index j=0; j<prolongI.size(); ++j )
            {
                if ( prolongJx[j] == 0 )
                {
                    for ( Index k=0; k<resultL.size(); ++k )
                    {
                        if ( resultJx[k] == jX[i] + prolongJx[j] &&
                             resultJy[k] == jY[i] + prolongJy[j] )
                        {
                            resultL[k]+=opL[i]*prolongI[j];
                        }
                    }
                }
            }
        }
        else
        {
            for ( Index j=0; j<prolongI.size(); ++j )
            {
                if ( prolongJx[j] != 0 && prolongJy[j] != 0 )
                {
                    for ( Index k=0; k<resultL.size(); ++k )
                    {
                        if ( resultJx[k] == jX[i] + prolongJx[j] &&
                             resultJy[k] == jY[i] + prolongJy[j] )
                        {
                            resultL[k]+=opL[i]*prolongI[j];
                        }
                    }
                }
            }
        }
    }
}

}
}



/*std::vector<PositionArray > TwoGridGalerkin::initJX_(
	const Restriction& restriction,
	const Stencil& stencil,
	const Prolongation& prolongation)
{
	std::vector<PositionArray > result;
	for (Index i=0; i<9; ++i)
		result.push_back(computeJ(i,restriction.getJx(i),stencil.getJx(i),prolongation.getJx(i)));
	return result;
}

std::vector<PositionArray > TwoGridGalerkin::initJY_(
	const Restriction& restriction,
	const Stencil& stencil,
	const Prolongation& prolongation)
{
	std::vector<PositionArray > result(9);
	for (Index i=0; i<9; ++i)
		somecall;
	return result;
}

Index TwoGridGalerkin::initSize_(
	const Restriction& restriction,
	const Stencil& stencil,
	const Prolongation& prolongation)
{
	Index resultSize=2;
	Index restrictionSize=std::max(
		std::max(
		    std::abs(restriction.getJx(C).max()),
		    std::abs(restriction.getJx(C).min())),
		std::max(
		    std::abs(restriction.getJy(C).max()),
		    std::abs(restriction.getJy(C).min())));
	Index prolongationSize=std::max(
		std::max(
		    std::abs(prolongation.getJx(C).max()),
		    std::abs(prolongation.getJx(C).min())),
		std::max(
		    std::abs(prolongation.getJy(C).max()),
		    std::abs(prolongation.getJy(C).min())));
	if (stencil.size()<=1 && restrictionSize<=1
						  && prolongationSize<=1)
		resultSize=1;
	return resultSize;
}

}

//std::vector<std::vector<PositionArray > > Galerkin::initJx_(const Stencil& stencil)
//{
//	std::vector<std::vector<PositionArray > > jx(1);
//	jx[0]->resize(9);
//	for (Index i=0; i<9; ++i)
//	{
//		PositionArray temp=stencil.getJx(static_cast<Position>(i));
//		jx[0][i]->resize(temp.size());
//		jx[0][i]=temp;
//	}
//	return jx;
//}
//
//std::vector<std::vector<PositionArray > > Galerkin::initJy_(const Stencil& stencil)
//{
//	std::vector<std::vector<PositionArray > > jy(1);
//	jy[0]->resize(9);
//	for (Index i=0; i<9; ++i)
//	{
//		PositionArray temp=stencil.getJy(static_cast<Position>(i));
//		jy[0][i].resize(temp.size());
//		jy[0][i]=temp;
//	}
//	return jy;
//}
//
//void Galerkin::updateSize_()
//{
//	if (size_==2)
//		return;
//	else
//	{
//		std::vector<const Prolongation&>::const_iterator prolongation=prolongations_.front();
//	}
//}
//void Galerkin::updateJxJy_()
//{
//	
//}
//
//NumericArray Galerkin::computeL(
//    const Position position,
//	const Index sx,
//	const Index sy,
//	const Index nx,
//	const Index ny) const
//{
//	NumericArray result;
//	return result;
//}
//
//virtual Precision Galerkin::apply(
//	const NumericArray& u,
//    const Position position,
//    const Index sx,
//    const Index sy,
//    const Index nx,
//    const Index ny) const
//{
//	Precision result=0;
//	NumericArray operatorL=getL(position,sx,sy,nx,ny);
//	PositionArray jX=getJx(position);
//	PositionArray jY=getJy(position);
//	for (Index i=0; i<operatorL.size(); ++i)
//		result+=operatorL[i]*u[(sy+jY[i])*(nx+1)+sx+jX[i]];
//	return result;
//}
//
//virtual Precision Galerkin::getCenter(
//	const Position position,
//	const Index sx,
//	const Index sy,
//	const Index nx,
//	const Index ny) const
//{
//	NumericArray operatorL=getL(position,sx,sy,nx,ny);
//	return operatorL[0];
//}
//
//virtual const NumericArray& Galerkin::getL(
//    const Position position,
//	const Index sx,
//	const Index sy,
//	const Index nx,
//	const Index ny) const
//{
//    if (data_[position].find(Quadruple(sx,sy,nx,ny))==data_[postition].end())
//		data_[postition].insert(Quadruple(sx,sy,nx,ny),computeL(postition,sx,sy,nx,ny);
//    return data_[position][Quadruple(sx,sy,nx,ny)];
//}

}*/
