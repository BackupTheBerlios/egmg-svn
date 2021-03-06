#ifndef FUNCTION_H_
#define FUNCTION_H_

#include "../general/parameters.h"

namespace mg
{

class Function
{
protected:
	virtual Precision apply(Precision x, Precision y) const =0;
public:
	Precision operator() (Precision x, Precision y) const;
    virtual ~Function() =0;
private:
	const Function& operator =(const Function& rhs);
};

}

#endif
