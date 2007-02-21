#ifndef TESTFUNCTION3_H_
#define TESTFUNCTION3_H_

#include "Function.h"

namespace mg
{

class TestFunction3: public Function
{
public:
    TestFunction3(Precision c = 1.0) : c_(c) {}
protected:
	virtual Precision apply(Precision x, Precision y) const;
private:
    const Precision c_;
};

}

#endif /*TESTFUNCTION3_H_*/
