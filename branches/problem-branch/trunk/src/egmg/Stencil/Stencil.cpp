/** \file Stencil.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implimentation of Stencil.
 */
#include "Stencil.h"

namespace mg
{

Stencil::Stencil() {}

Stencil::~Stencil() {}

NumericArray Stencil::getLInSize(
    const Position pos,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny,
	const Index size ) const
{
	assert( size == 1 || size == 2 );
	NumericArray result( 0.0, size == 1 ? 9 : 25 );
	
	NumericArray operatorL = getL( pos, sx, sy, nx, ny );
	PositionArray jX = getJx( pos, nx, ny );
	PositionArray jY = getJy( pos, nx, ny );

	for ( Index i = 0; i < operatorL.size(); i++ )
	{
		if( jX[ i ] == -1 )
		{
			
			if( jY[ i ] == -1 )
				result[ SW ] = operatorL[ i ];
			else if( jY[ i ] ==  0 )
				result[ W ] = operatorL[ i ];
			else if( jY[ i ] ==  1)
				result[ NW ] = operatorL[ i ];
			if ( size > 1 )
			{
				if ( jY[ i ] == -2 )
					result[ SSW ] = operatorL[ i ];
				else if( jY[ i ] ==  2 )
					result[ NNW ] = operatorL[ i ];
			}
			else if ( jY[i] == -2 || jY[i] ==  2 )
			{
				result[ C ] += operatorL[ i ];
			}
		}
		else if( jX[ i ] == 0 )
		{
			if( jY[ i ] == -1 )
				result[ S ] = operatorL[ i ];
			else if( jY[ i ] ==  0 )
				result[ C ] += operatorL[ i ];
			else if( jY[ i ] ==  1 )
				result[ N ] = operatorL[ i ];
			if ( size > 1 )
			{
				if ( jY[ i ] == -2 )
					result[ SS ] = operatorL[ i ];
				else if( jY[ i ] ==  2 )
					result[ NN ] = operatorL[ i ];
			}
			else if ( jY[ i ] == -2 || jY[ i ] ==  2 )
			{
				result[ C ] += operatorL[ i ];
			}
		}
		else if( jX[ i ] == 1 )
		{
			if( jY[ i ] == -1 )
				result[ SE ] = operatorL[ i ];
			else if( jY[ i ] ==  0 )
				result[ E ] = operatorL[ i ];
			else if( jY[ i ] ==  1 )
				result[ NE ] = operatorL[ i ];
			if ( size > 1 )
			{
				if ( jY[ i ] == -2 )
					result[ SSE ] = operatorL[ i ];
				else if( jY[ i ] ==  2 )
					result[ NNE ] = operatorL[ i ];
			}
			else if ( jY[ i ] == -2 || jY[ i ] ==  2 )
			{
				result[ C ] += operatorL[ i ];
			}
		}

		if ( size > 1 )
		{
			if ( jX[ i ] == -2 )
			{
				if ( jY[ i ] == -2 )
					result[ SSWW ] = operatorL[ i ];
				else if ( jY[ i ] == -1 )
					result[ SWW ] = operatorL[ i ];
				else if ( jY[ i ] ==  0 )
					result[ WW ] = operatorL[ i ];
				else if ( jY[ i ] ==  1 )
					result[ NWW ] = operatorL[ i ];
				else if ( jY[ i ] ==  2 )
					result[ NNWW ] = operatorL[ i ];
			}
			
			else if ( jX[ i ] == 2 )
			{
				if ( jY[ i ] == -2 )
					result[ SSEE ] = operatorL[ i ];
				else if ( jY[ i ] == -1 )
					result[ SEE ] = operatorL[ i ];
				else if ( jY[ i ] ==  0 )
					result[ EE ] = operatorL[ i ];
				else if ( jY[ i ] ==  1 )
					result[ NEE ] = operatorL[ i ];
				else if ( jY[ i ] ==  2 )
					result[ NNEE ] = operatorL[ i ];
			}
		}
		else if ( jX[ i ] == -2 || jX[ i ] == 2 )
		{
			result[ C ] += operatorL[ i ];
		}
	}
	return result;
}

void Stencil::pushTransferOperators(
    const Restriction&,
    const Prolongation&,
    const Index,
    const Index ) {}

void Stencil::popTransferOperators() {}

}
