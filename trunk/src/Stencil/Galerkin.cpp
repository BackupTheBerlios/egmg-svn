/** \file Galerkin.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class Galerkin.
 */

#include <cstdlib>
#include "Galerkin.h"

namespace mg
{

namespace
{

PositionArray computeJ(
	const Position pos,
	const PositionArray& restriction,
	const PositionArray& stencil,
	const PositionArray& prolongation)
{
	std::vector<PositionArray> result;
	for (Index i=0; i<restriction.size() ; ++i)
	{
		PositionArray temp(stencil);
		temp+=restriction[i];
		for (Index j=0; j<temp.size(); ++j)
		{
			if (temp[j]%2!=0)
			{
				PositionArray temp2(prolongation);
				temp2+=temp[j];
				result.push_back(temp2);
				temp[j]=0;
			}
		}
		result.push_back(temp);
	}

	return result;
}

std::vector<PositionArray > TwoGridGalerkin::initJX_(
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

}
