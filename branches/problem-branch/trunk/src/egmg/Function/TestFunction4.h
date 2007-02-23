#ifndef TESTFUNCTION4_H_
#define TESTFUNCTION4_H_

#include "Function.h"

namespace mg
{

class TestFunction4: public Function
{
public:
    TestFunction4(Precision c = 1.0) : c_(c) {}
protected:
	virtual Precision apply(Precision x, Precision y) const;
private:
    const Precision c_;
};

}

#endif /*TESTFUNCTION4_H_*/
