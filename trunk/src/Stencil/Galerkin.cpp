/** \file Galerkin.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class Galerkin.
 */

#include <stack>

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
    const NumericArray restrictionI = restriction.getI( C,sx,sy,nx,ny,stencil );
    for ( Index i=0; i<restrictionI.size(); ++i )
    {
        const Position newPos = getNewPos(
            pos, restrictionJx[i], restriciionJy[i] );
        PositionArray stencilJx = stencil.getJx( newPos );
        PositionArray stencilJy = stencil.getJy( newPos );
        NumericArray stencilL = stencil.getL(
            newPos,
            sx+restrictionJx[i],
            sy+restriciionJy[i],
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
            newPos,sx+resultJx[i],sy+resultJy[i],2*nx,2*ny,stencil );
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
    Integer sizeToBorder = newSizeToBorder;
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

}
