/** \file WHighJacScalarProduct.cpp
 * \author Andre Oeckerath
 */
#ifndef WHighJacScalarProduct_H_
#define WHighJacScalarProduct_H_

#include <cmath>


namespace mg
{
	const Precision halfPi = asin( 1.0 );

	inline Precision WHighJacScalarProduct
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

#endif /*WHighJacScalarProduct_H_*/
