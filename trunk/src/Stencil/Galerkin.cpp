/** \file Galerkin.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class Galerkin.
 */

#include <cstdlib>
#include "Galerkin.h"
#include "Stencil.h"
#include "../general/parameters.h"
#include "../functions/expansion.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"
#include "../functions/generatePositionArrays.h"
#include "../functions/printStencil.h"
#include <stack>
#include <iostream>
#include <iomanip>

namespace mg
{

namespace
{

Index ComputeSize(
    const Restriction& restriction,
    const Stencil& stencil)
{
    PositionArray restrictionJx = restriction.getJx( C );
    PositionArray restrictionJy = restriction.getJy( C );
    const Index restSize = expansion(restrictionJx,restrictionJy);
    return restSize + stencil.size();
}

Index ComputeSize(
    const PositionArray& jX,
    const PositionArray& jY,
    const Prolongation& prolong)
{
    PositionArray prolongJx = prolong.getJx( C );
    PositionArray prolongJy = prolong.getJy( C );
    const Index prolongSize = expansion(prolongJx, prolongJy);
    const Index result = expansion(jX, jY);
    if ( result%2 == 0 )
    {
        return result/2 + prolongSize - 1;
    }
    return ( result - 1 )/2 + prolongSize;
}

inline Position getNewPos(
    const Position pos,
    const Integer jX,
    const Integer jY)
{
    Position newPos = C;
    switch ( pos )
    {
        case W:
            if ( jY < 0 )
            {
                newPos = W;
            }
            break;
        case N:
            if ( jX > 0 )
            {
                newPos = N;
            }
            break;
        case E:
            if ( jY > 0 )
            {
                newPos = E;
            }
            break;
        case S:
            if ( jX < 0 )
            {
                newPos = S;
            }
            break;
        case NW:
            if ( jX < 0 && jY > 0 )
            {
                newPos = NW;
            }
            break;
        case NE:
            if ( jX > 0 && jY > 0 )
            {
                newPos = NW;
            }
            break;
        case SE:
            if ( jX > 0 && jY < 0 )
            {
                newPos = NW;
            }
            break;
        case SW:
            if ( jX > 0 && jY < 0 )
            {
                newPos = NW;
            }
            break;
        default:        // case C:
            newPos = C;
    }
    return newPos;
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
    const PositionArray restrictionJx = restriction.getJx( C );
    const PositionArray restriciionJy = restriction.getJy( C );
    const NumericArray restrictionI = restriction.getI( C,sx,sy,nx,ny,stencil );

    const Index size = ComputeSize(restriction,stencil);

    generatePositionArrays(jX,jY,pos,size,2);
    resultL.resize(jX.size());
    resultL = 0.0;

    for ( Index i=0; i<restrictionI.size(); ++i )
    {
        const Position newPos = getNewPos(
            pos, restrictionJx[i], restriciionJy[i] );
        PositionArray stencilJx = stencil.getJx( newPos );
        PositionArray stencilJy = stencil.getJy( newPos );
        NumericArray stencilL = stencil.getL(
            newPos,
            2*sx+restrictionJx[i],
            2*sy+restriciionJy[i],
            2*nx, 2*ny);
    
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
    const Index size = ComputeSize(jX,jY,prolongation);
    generatePositionArrays( resultJx,resultJy,pos,size,1 );
    resultL.resize( resultJx.size() );
    resultL = 0.0;
    for (Index i=0; i<opL.size(); ++i)
    {
        const Position newPos = getNewPos(
            pos, jX[i], jY[i] );
        const PositionArray prolongJx = prolongation.getJx( newPos );
        const PositionArray prolongJy = prolongation.getJy( newPos );
        const NumericArray prolongI = prolongation.getI( 
            newPos,sx+jX[i],sy+jY[i],2*nx,2*ny,stencil );
        if (jX[i]%2 == 0 && jY[i]%2 == 0) //coars grid point
        {
            for (Index j=0; j<resultL.size(); ++j)
            {
                if (resultJx[j] == jX[i]/2 && resultJy[j] == jY[i]/2)
                {
                    resultL[j]+=prolongI[0]*opL[i];
                    break;
                }
            }
        }
        else if (jX[i]%2 == 0) //fine grid point on x-coars grid line
        {
            for ( Index j=0; j<prolongI.size(); ++j )
            {
                if ( prolongJy[j] != 0 && prolongJx[j] == 0 )
                {
                    Integer posX = (jX[i] + prolongJx[j])/2;
                    Integer posY = (jY[i] + prolongJy[j])/2;
                    for ( Index k=0; k<resultL.size(); ++k )
                    {
                        if ( resultJx[k] == posX &&
                             resultJy[k] == posY )
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
                if ( prolongJx[j] != 0 && prolongJy[j] == 0 )
                {
                    Integer posX = (jX[i] + prolongJx[j])/2;
                    Integer posY = (jY[i] + prolongJy[j])/2;
                    for ( Index k=0; k<resultL.size(); ++k )
                    {
                        if ( resultJx[k] == posX &&
                             resultJy[k] == posY )
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
                    Integer posX = (jX[i] + prolongJx[j])/2;
                    Integer posY = (jY[i] + prolongJy[j])/2;
                    for ( Index k=0; k<resultL.size(); ++k )
                    {
                        if ( resultJx[k] == posX &&
                             resultJy[k] == posY )
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

void computeGalerkin(
    NumericArray& resultL,
    PositionArray& resultJx,
    PositionArray& resultJy,
    const Position pos,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny,
    const Restriction& restriction,
    const Stencil& stencil,
    const Prolongation& prolongation)
{
    NumericArray interResultL;
    PositionArray interJx;
    PositionArray interJy;
    RestritionTimesStencil(
        interResultL,interJx,interJy,
        pos,sx,sy,nx,ny,
        restriction,
        stencil);
    printStencil(interResultL,interJx,interJy,std::cout);
    ProlongateStencil(
        resultL,resultJx,resultJy,
        interResultL,interJx,interJy,
        pos,sx,sy,nx,ny,
        prolongation,
        stencil);
    
}
}
/*
void testGenPosArray(const Index size, const Index size2)
{
    PositionArray jX;
    PositionArray jY;
    GeneratePositionArrays(jX,jY,C,size,size2);
    Integer neg = size;
    for ( mg::Integer sy=neg; sy>=-1*neg; --sy )
    {
        for ( mg::Integer sx=-1*neg; sx<=neg; ++sx )
        {
            for ( mg::Index k=0; k<jX.size(); ++k )
            {
                if ( jY[k] == sy && jX[k] == sx )
                {
                    std::cout<<std::setw(3)<<k;
                }   
            }
            std::cout<<" ";
        }
        std::cout<<std::endl<<std::endl;
    }
    std::cout<<std::endl;
    GeneratePositionArrays(jX,jY,W,size,size2);
    std::cout<<"W"<<std::endl;
    for ( mg::Integer sy=neg; sy>=-1*neg; --sy )
    {
        for ( mg::Integer sx=-1*neg; sx<=neg; ++sx )
        {
            mg::Index k=0;
            for ( k=0; k<jX.size(); ++k )
            {
                if ( jY[k] == sy && jX[k] == sx )
                {
                    std::cout<<std::setw(3)<<k;
                    break;
                }   
            }
            if ( k == jX.size() )
                std::cout<<"0x0";
            std::cout<<" ";
        }
        std::cout<<std::endl<<std::endl;
    }
    std::cout<<std::endl;
    GeneratePositionArrays(jX,jY,N,size,size2);
    std::cout<<"N"<<std::endl;
    for ( mg::Integer sy=neg; sy>=-1*neg; --sy )
    {
        for ( mg::Integer sx=-1*neg; sx<=neg; ++sx )
        {
            mg::Index k=0;
            for ( k=0; k<jX.size(); ++k )
            {
                if ( jY[k] == sy && jX[k] == sx )
                {
                    std::cout<<std::setw(3)<<k;
                    break;
                }   
            }
            if ( k == jX.size() )
                std::cout<<"0x0";
            std::cout<<" ";
        }
        std::cout<<std::endl<<std::endl;
    }
    std::cout<<std::endl;
    GeneratePositionArrays(jX,jY,E,size,size2);
    std::cout<<"E"<<std::endl;
    for ( mg::Integer sy=neg; sy>=-1*neg; --sy )
    {
        for ( mg::Integer sx=-1*neg; sx<=neg; ++sx )
        {
            mg::Index k=0;
            for ( k=0; k<jX.size(); ++k )
            {
                if ( jY[k] == sy && jX[k] == sx )
                {
                    std::cout<<std::setw(3)<<k;
                    break;
                }   
            }
            if ( k == jX.size() )
                std::cout<<"0x0";
            std::cout<<" ";
        }
        std::cout<<std::endl<<std::endl;
    }
    std::cout<<std::endl;
    GeneratePositionArrays(jX,jY,S,size,size2);
    std::cout<<"S"<<std::endl;
    for ( mg::Integer sy=neg; sy>=-1*neg; --sy )
    {
        for ( mg::Integer sx=-1*neg; sx<=neg; ++sx )
        {
            mg::Index k=0;
            for ( k=0; k<jX.size(); ++k )
            {
                if ( jY[k] == sy && jX[k] == sx )
                {
                    std::cout<<std::setw(3)<<k;
                    break;
                }   
            }
            if ( k == jX.size() )
                std::cout<<"0x0";
            std::cout<<" ";
        }
        std::cout<<std::endl<<std::endl;
    }
    std::cout<<std::endl;
    GeneratePositionArrays(jX,jY,NW,size,size2);
    std::cout<<"NW"<<std::endl;
    for ( mg::Integer sy=neg; sy>=-1*neg; --sy )
    {
        for ( mg::Integer sx=-1*neg; sx<=neg; ++sx )
        {
            mg::Index k=0;
            for ( k=0; k<jX.size(); ++k )
            {
                if ( jY[k] == sy && jX[k] == sx )
                {
                    std::cout<<std::setw(3)<<k;
                    break;
                }   
            }
            if ( k == jX.size() )
                std::cout<<"0x0";
            std::cout<<" ";
        }
        std::cout<<std::endl<<std::endl;
    }
    std::cout<<"NE"<<std::endl;
    GeneratePositionArrays(jX,jY,NE,size,size2);
    for ( mg::Integer sy=neg; sy>=-1*neg; --sy )
    {
        for ( mg::Integer sx=-1*neg; sx<=neg; ++sx )
        {
            mg::Index k=0;
            for ( k=0; k<jX.size(); ++k )
            {
                if ( jY[k] == sy && jX[k] == sx )
                {
                    std::cout<<std::setw(3)<<k;
                    break;
                }   
            }
            if ( k == jX.size() )
                std::cout<<"0x0";
            std::cout<<" ";
        }
        std::cout<<std::endl<<std::endl;
    }
    GeneratePositionArrays(jX,jY,SE,size,size2);
    std::cout<<"SE"<<std::endl;
    for ( mg::Integer sy=neg; sy>=-1*neg; --sy )
    {
        for ( mg::Integer sx=-1*neg; sx<=neg; ++sx )
        {
            mg::Index k=0;
            for ( k=0; k<jX.size(); ++k )
            {
                if ( jY[k] == sy && jX[k] == sx )
                {
                    std::cout<<std::setw(3)<<k;
                    break;
                }   
            }
            if ( k == jX.size() )
                std::cout<<"0x0";
            std::cout<<" ";
        }
        std::cout<<std::endl<<std::endl;
    }
    GeneratePositionArrays(jX,jY,SW,size,size2);
    std::cout<<"SW"<<std::endl;
    for ( mg::Integer sy=neg; sy>=-1*neg; --sy )
    {
        for ( mg::Integer sx=-1*neg; sx<=neg; ++sx )
        {
            mg::Index k=0;
            for ( k=0; k<jX.size(); ++k )
            {
                if ( jY[k] == sy && jX[k] == sx )
                {
                    std::cout<<std::setw(3)<<k;
                    break;
                }   
            }
            if ( k == jX.size() )
                std::cout<<"0x0";
            std::cout<<" ";
        }
        std::cout<<std::endl<<std::endl;
    }
    
}

}*/



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
//  std::vector<std::vector<PositionArray > > jx(1);
//  jx[0]->resize(9);
//  for (Index i=0; i<9; ++i)
//  {
//      PositionArray temp=stencil.getJx(static_cast<Position>(i));
//      jx[0][i]->resize(temp.size());
//      jx[0][i]=temp;
//  }
//  return jx;
//}
//
//std::vector<std::vector<PositionArray > > Galerkin::initJy_(const Stencil& stencil)
//{
//  std::vector<std::vector<PositionArray > > jy(1);
//  jy[0]->resize(9);
//  for (Index i=0; i<9; ++i)
//  {
//      PositionArray temp=stencil.getJy(static_cast<Position>(i));
//      jy[0][i].resize(temp.size());
//      jy[0][i]=temp;
//  }
//  return jy;
//}
//
//void Galerkin::updateSize_()
//{
//  if (size_==2)
//      return;
//  else
//  {
//      std::vector<const Prolongation&>::const_iterator prolongation=prolongations_.front();
//  }
//}
//void Galerkin::updateJxJy_()
//{
//  
//}
//
//NumericArray Galerkin::computeL(
//    const Position position,
//  const Index sx,
//  const Index sy,
//  const Index nx,
//  const Index ny) const
//{
//  NumericArray result;
//  return result;
//}
//
//virtual Precision Galerkin::apply(
//  const NumericArray& u,
//    const Position position,
//    const Index sx,
//    const Index sy,
//    const Index nx,
//    const Index ny) const
//{
//  Precision result=0;
//  NumericArray operatorL=getL(position,sx,sy,nx,ny);
//  PositionArray jX=getJx(position);
//  PositionArray jY=getJy(position);
//  for (Index i=0; i<operatorL.size(); ++i)
//      result+=operatorL[i]*u[(sy+jY[i])*(nx+1)+sx+jX[i]];
//  return result;
//}
//
//virtual Precision Galerkin::getCenter(
//  const Position position,
//  const Index sx,
//  const Index sy,
//  const Index nx,
//  const Index ny) const
//{
//  NumericArray operatorL=getL(position,sx,sy,nx,ny);
//  return operatorL[0];
//}
//
//virtual const NumericArray& Galerkin::getL(
//    const Position position,
//  const Index sx,
//  const Index sy,
//  const Index nx,
//  const Index ny) const
//{
//    if (data_[position].find(Quadruple(sx,sy,nx,ny))==data_[postition].end())
//      data_[postition].insert(Quadruple(sx,sy,nx,ny),computeL(postition,sx,sy,nx,ny);
//    return data_[position][Quadruple(sx,sy,nx,ny)];
//}

}*/
