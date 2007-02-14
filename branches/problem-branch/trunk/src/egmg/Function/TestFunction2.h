#ifndef TESTFUNCTION2_H_
#define TESTFUNCTION2_H_

#include "Function.h"

namespace mg
{

class TestFunction2: public Function
{
public:
    TestFunction2(Precision c = 1.0) : c_(c) {}
protected:
	virtual Precision apply(Precision x, Precision y) const;
private:
    const Precision c_;
};

}

#endif /*TESTFUNCTION2_H_*/
