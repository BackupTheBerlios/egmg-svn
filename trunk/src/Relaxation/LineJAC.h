#ifndef LINEJAC_H_
#define LINEJAC_H_

#include "Relaxation.h"

namespace mg
{

class LineJAC : public mg::Relaxation
{
public:
	LineJAC();
	virtual ~LineJAC();
};

}

#endif /*LINEJAC_H_*/
