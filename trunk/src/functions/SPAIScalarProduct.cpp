/** \file SPAIScalarProduct.cpp
 * \author André Oeckerath
 */
#ifndef SPAIScalarProduct_H_
#define SPAIScalarProduct_H_

#include <cmath>


namespace mg
{
	const Precision halfPi = asin( 1 );

	inline Precision SPAIScalarProduct
		(const int i, const int j) 
    {
		if(i != 0)
		{
			if(j != 0)
			{
				return -4.0*sin(i*halfPi)*sin(j*halfPi)/(i*j);
			}
			else
			{
				return -4.0*halfPi*sin(i*halfPi)/i;
			}
		}
		else
		{
			if(j != 0)
			{
				return -4.0*halfPi*sin(j*halfPi)/j;
			}
			else
			{
				return 12*halfPi*halfPi;
			}
		}

        
    }


}

#endif /*SPAIScalarProduct_H_*/
