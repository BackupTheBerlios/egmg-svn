/** \file Galerkin.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class Galerkin.
 */

#include <stack>
#include <algorithm>
#include <cmath>

#include "Galerkin.h"
#include "../functions/expansion.h"
#include "../functions/generatePositionArrays.h"

namespace mg
{
    
namespace
{

Index computeSize(
    const Restriction& restriction,
    const Stencil& stencil,
    const Prolongation& prolongation)
{
    PositionArray restrictionJx = restriction.getJx( C );
    PositionArray restrictionJy = restriction.getJy( C );
    PositionArray prolongationJx = prolongation.getJx( C );
    PositionArray prolongationJy = prolongation.getJy( C );
    const Index restSize = expansion( restrictionJx, restrictionJy );
    const Index prolongSize = expansion( prolongationJx, prolongationJy );
    const Index result = restSize+stencil.size();
    if ( result%2 == 0 )
    {
        return result-1+prolongSize;
    }
    return result+prolongSize;
}

Position getNewPos(
    const Position pos,
    const Integer jX,
    const Integer jY)
{
    Position newPos = C;
    switch ( pos )
    {
    case W:
        if ( jX < 0 )
        {
            newPos = W;
        }
        break;
    case N:
        if ( jY > 0 )
        {
            newPos = N;
        }
        break;
    case E:
        if ( jX > 0 )
        {
            newPos = E;
        }
        break;
    case S:
        if ( jY < 0 )
        {
            newPos = S;
        }
        break;
    case NW:
        if ( jX < 0 && jY > 0 )
        {
            newPos = NW;
        }
        else if ( jX < 0 )
        {
            newPos = W;
        }
        else if ( jY > 0 )
        {
            newPos = N;
        }
        break;
    case NE:
        if ( jX > 0 && jY > 0 )
        {
            newPos = NE;
        }
        else if ( jX > 0 )
        {
            newPos = E;
        }
        else if ( jY > 0 )
        {
            newPos = N;
        }
        break;
    case SE:
        if ( jX > 0 && jY < 0 )
        {
            newPos = SE;
        }
        else if ( jX > 0 )
        {
            newPos = E;
        }
        else if ( jY < 0 )
        {
            newPos = S;
        }
        break;
    case SW:
        if ( jX < 0 && jY < 0 )
        {
            newPos = SW;
        }
        else if ( jX < 0 )
        {
            newPos = W;
        }
        else if ( jY < 0 )
        {
            newPos = S;
        }
        break;
    default:        // case C:
        newPos = C;
    }
    return newPos;
}

void restritionTimesStencil(
    NumericArray& resultL,
    const PositionArray& resultJx,
    const PositionArray& resultJy,
    const Position pos,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny,
    const Restriction& restriction,
    const Stencil& stencil)
{
    const PositionArray restrictionJx = restriction.getJx( C );
    const PositionArray restriciionJy = restriction.getJy( C );
    const NumericArray restrictionI = restriction.getI( C, 2 * sx, 2 * sy, 2 * nx, 2 * ny,stencil );
    for ( Index i=0; i<restrictionI.size(); ++i )
    {
        const Position newPos = getNewPos(
            pos, restrictionJx[i], restriciionJy[i] );
        PositionArray stencilJx = stencil.getJx( newPos, 2 * nx, 2 * ny );
        PositionArray stencilJy = stencil.getJy( newPos, 2 * nx, 2 * ny );
        NumericArray stencilL = stencil.getL(
            newPos,
            2 * sx+restrictionJx[i],
            2 * sy+restriciionJy[i],
            2*nx, 2*ny);

        stencilJx+=restrictionJx[i];
        stencilJy+=restriciionJy[i];
        stencilL*=restrictionI[i];
        for ( Index j=0; j<stencilL.size(); ++j)
        {
            for ( Index k=0; k<resultL.size(); ++k )
            {
                if ( stencilJx[j] == resultJx[k] && stencilJy[j] == resultJy[k] )
                {
                    resultL[k] += stencilL[j];
                    break;
                }
            }
        }
    }
}

void prolongateStencil(
    NumericArray& resultL,
    const PositionArray& resultJx,
    const PositionArray& resultJy,
    const Position pos,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny,
    const Prolongation& prolongation,
    const Stencil& stencil)
{
    const NumericArray opL(resultL);
    for (Index i=0; i<resultL.size(); ++i)
    {
        const Position newPos = getNewPos(
            pos, resultJx[i], resultJy[i] );
        const PositionArray prolongJx = prolongation.getJx( newPos );
        const PositionArray prolongJy = prolongation.getJy( newPos );
        const NumericArray prolongI = prolongation.getI( 
            newPos,2*sx+resultJx[i],2*sy+resultJy[i],2 * nx,2 * ny,stencil );
        if ( resultJx[i]%2 == 0 ) //fine grid point on x-coars grid line
        {
            for ( Index j=0; j<prolongI.size(); ++j )
            {
                if ( prolongJy[j] != 0 && prolongJx[j] == 0 )
                {
                    Integer posX = resultJx[i];
                    Integer posY = resultJy[i] + prolongJy[j];
                    for ( Index k=0; k<resultL.size(); ++k )
                    {
                        if ( resultJx[k] == posX &&
                             resultJy[k] == posY )
                        {
                            resultL[k]+=opL[i]*prolongI[j];
                            break;
                        }
                    }
                }
            }
        }
        else if (resultJy[i]%2 == 0) //fine grid point on y-coars grid line
        {
            for ( Index j=0; j<prolongI.size(); ++j )
            {
                if ( prolongJx[j] != 0 && prolongJy[j] == 0 )
                {
                    Integer posX = resultJx[i]+prolongJx[j];
                    Integer posY = resultJy[i];
                    for ( Index k=0; k<resultL.size(); ++k )
                    {
                        if ( resultJx[k] == posX &&
                             resultJy[k] == posY )
                        {
                            resultL[k]+=opL[i]*prolongI[j];
                            break;
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
                    const Integer posX = resultJx[i] + prolongJx[j];
                    const Integer posY = resultJy[i] + prolongJy[j];
                    for ( Index k=0; k<resultL.size(); ++k )
                    {
                        if ( resultJx[k] == posX &&
                             resultJy[k] == posY )
                        {
                            resultL[k]+=opL[i]*prolongI[j];
                            break;
                        }
                    }
                }
            }
        }
    }
}

void shrinkStencil(
    NumericArray& resultL,
    PositionArray& resultJx,
    PositionArray& resultJy,
    const NumericArray& opL,
    const PositionArray& jX,
    const PositionArray& jY,
    const Position pos,
    Index newSize,
    Index newSizeToBorder)
{
    newSize = std::min( expansion( jX, jY )/2, newSize );
    newSizeToBorder = std::min( expansion( jX, jY )/2, newSizeToBorder );
    generatePositionArrays( resultJx, resultJy, pos, newSize, newSizeToBorder );
    resultL.resize(resultJx.size());
    resultL = 0.0;
    Integer sizeToBorder = static_cast<Integer>( newSizeToBorder );
    for ( Index i = 0; i < opL.size(); ++i )
    {
        if ( jX[i]%2 == 0 && jY[i]%2 == 0 )
        {
            const Integer posX = jX[i]/2;
            const Integer posY = jY[i]/2;
            const Index absPosX = std::abs( posX );
            const Index absPosY = std::abs( posY );
            if ((pos == C ||
                 pos == W && posX>=-1*sizeToBorder || 
                 pos == N && posY<=sizeToBorder ||
                 pos == E && posX<=sizeToBorder ||
                 pos == S && posY>=-1*sizeToBorder ||
                 pos == NW && posX>=-1*sizeToBorder && posY<=sizeToBorder ||
                 pos == NE && posX<=sizeToBorder && posY<=sizeToBorder ||
                 pos == SE && posX<=sizeToBorder && posY>=-1*sizeToBorder ||
                 pos == SW && posX>=-1*sizeToBorder && posY>=-1*sizeToBorder)
                 && absPosX <= newSize && absPosY <= newSize )
            {
                for ( Index j = 0; j < resultL.size(); ++j )
                {
                    if ( resultJx[j] == posX && resultJy[j] == posY )
                    {
                        resultL[j] = opL[i];
                        break;
                    }
                }
            }
            else
                resultL[0]+=opL[i];
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
    const Index resultSize = computeSize( restriction, stencil, prolongation );
    PositionArray interResultJx;
    PositionArray interResultJy;
    generatePositionArrays( interResultJx, interResultJy, pos, resultSize, 2 );
    NumericArray interResultL( 0.0, interResultJx.size() );
    restritionTimesStencil(
        interResultL,interResultJx,interResultJy,
        pos,sx,sy,nx,ny,
        restriction,
        stencil);
    prolongateStencil(
        interResultL,interResultJx,interResultJy,
        pos,sx,sy,nx,ny,
        prolongation,
        stencil);
    shrinkStencil(
        resultL,resultJx,resultJy,
        interResultL,interResultJx,interResultJy,
        pos,2,1);
}

}

Galerkin::Galerkin( const Stencil& fineGridOperator )
    : PreCalculatedStencil( fineGridOperator ),
      prolongations_(0),
      restrictions_(0) {}

Galerkin::~Galerkin() {}

void Galerkin::update()
{
    prolongations_.pop_back();
    restrictions_.pop_back();
}

void Galerkin::update(
    const Restriction& restriction,
    const Prolongation& prolongation,
    const Index nx,
    const Index ny )
{
    prolongations_.push_back(&prolongation);
    restrictions_.push_back(&restriction);
    //check if we have already calculated this level
    if ( data_.contains( currentDepth_, C, 1, 1 ) )
    {
        return;
    }
    PositionArray jX;
    PositionArray jY;
    NumericArray operatorL;
    //C
    for ( Index sx = 1; sx < nx; ++sx )
    {
        for ( Index sy = 1; sy < ny; ++sy )
        {
            computeGalerkin(
                operatorL,jX,jY,
                C,sx,sy,nx,ny,
                *restrictions_.front(),*this,*prolongations_.front());
            if ( sx == nx/2 && sy == ny/2 )
                data_.insert( currentDepth_, C, jX, jY );
            data_.insert( currentDepth_, C, sx, sy, operatorL );
        }
    }
    //W
    for ( Index sy = 1; sy < ny; ++sy )
    {
        computeGalerkin(
            operatorL,jX,jY,
            W,1,sy,nx,ny,
            *restrictions_.front(),*this,*prolongations_.front());
        if ( sy == ny/2 )
            data_.insert( currentDepth_, W, jX, jY );
        data_.insert( currentDepth_, W, 1, sy, operatorL );
    }
    //N
    for ( Index sx = 1; sx < nx; ++sx )
    {
        computeGalerkin(
            operatorL,jX,jY,
            N,sx,ny-1,nx,ny,
            *restrictions_.front(),*this,*prolongations_.front());
        if ( sx == nx/2 )
            data_.insert( currentDepth_, N, jX, jY );
        data_.insert( currentDepth_, N, sx, ny-1, operatorL );
    }
    //E
    for ( Index sy = 1; sy < ny; ++sy )
    {
        computeGalerkin(
            operatorL,jX,jY,
            E,nx-1,sy,nx,ny,
            *restrictions_.front(),*this,*prolongations_.front());
        if ( sy == ny/2 )
            data_.insert( currentDepth_, E, jX, jY );
        data_.insert( currentDepth_, E, nx-1, sy, operatorL );
    }
    //S
    for ( Index sx = 1; sx < nx; ++sx )
    {
        computeGalerkin(
            operatorL,jX,jY,
            S,sx,1,nx,ny,
            *restrictions_.front(),*this,*prolongations_.front());
        if ( sx == nx/2 )
            data_.insert( currentDepth_, S, jX, jY );
        data_.insert( currentDepth_, S, sx, 1, operatorL );
    }
    //NW
    computeGalerkin(
        operatorL,jX,jY,
        NW,1,ny-1,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
    data_.insert( currentDepth_, NW, jX, jY );
    data_.insert( currentDepth_, NW, 1, ny-1, operatorL );
    //NE
    computeGalerkin(
        operatorL,jX,jY,
        NE,nx-1,ny-1,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
    data_.insert( currentDepth_, NE, jX, jY );
    data_.insert( currentDepth_, NE, nx-1, ny-1, operatorL );
    //SE
    computeGalerkin(
        operatorL,jX,jY,
        SE,nx-1,1,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
    data_.insert( currentDepth_, SE, jX, jY );
    data_.insert( currentDepth_, SE, nx-1, 1, operatorL );
    //SW
    computeGalerkin(
        operatorL,jX,jY,
        SW,1,1,nx,ny,
        *restrictions_.front(),*this,*prolongations_.front());
    data_.insert( currentDepth_, SW, jX, jY );
    data_.insert( currentDepth_, SW, 1, 1, operatorL );
    currentDepth_ = prolongations_.size();
}

}
