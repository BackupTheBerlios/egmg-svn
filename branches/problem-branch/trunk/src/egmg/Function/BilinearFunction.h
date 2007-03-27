#ifndef BILINEARFUNCTION_H_
#define BILINEARFUNCTION_H_

#include "Function.h"

namespace mg
{

class BilinearFunction: public Function
{
public:
    BilinearFunction();
	BilinearFunction(const NumericArray& a);

protected:
	virtual Precision apply(Precision x, Precision y) const;
private:
    NumericArray initA() const;
    const NumericArray a_;
};

}

#endif /*BILINEARFUNCTION_H_*/
